#!/usr/bin/env Rscript
# =============================================================================
# 01_pairwise_wilcoxon.R
# Pairwise Wilcoxon Signed-Rank Test Comparison Heatmaps
#
# Reads vcfeval_all_stats.tsv OR hap.py all_stats.tsv and produces:
#   - Pairwise F1 difference heatmaps (SNP + INDEL) with Wilcoxon p-value stars
#
# Adapted from: archive/caller_benchmark-main/benchmark_analysis.R (lines 306–440)
#
# Usage:
#   Rscript visualization/01_pairwise_wilcoxon.R [stats_tsv] [output_dir]
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(dplyr)
  library(cowplot)
})

# ── Command-line args ──
args <- commandArgs(trailingOnly = TRUE)

# Try to find stats file: prefer hap.py output (has SNP/INDEL breakdown)
find_stats_file <- function() {
  project_dir <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), ".."), mustWork = FALSE)
  candidates <- c(
    file.path(project_dir, "results", "eval", "all_stats.tsv"),          # hap.py output
    file.path(project_dir, "results", "eval", "vcfeval_all_stats.tsv")   # vcfeval output
  )
  for (f in candidates) {
    if (file.exists(f)) return(f)
  }
  return(NULL)
}

stats_file <- if (length(args) >= 1) args[1] else find_stats_file()
output_dir <- if (length(args) >= 2) args[2] else file.path(dirname(dirname(stats_file)), "plots")

if (is.null(stats_file) || !file.exists(stats_file)) {
  stop("Stats file not found. Usage: Rscript 01_pairwise_wilcoxon.R <stats_tsv> [output_dir]")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cat("Input:  ", stats_file, "\n")
cat("Output: ", output_dir, "\n\n")

# ── Read data ──
raw <- read.delim(stats_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cat("Columns: ", paste(colnames(raw), collapse = ", "), "\n")
cat("Rows:    ", nrow(raw), "\n\n")

# ── Detect format: hap.py vs vcfeval ──
is_happy <- "Type" %in% colnames(raw) && "Filter" %in% colnames(raw)

if (is_happy) {
  cat("Detected: hap.py format (with SNP/INDEL breakdown)\n\n")

  # Standardize column names
  if ("Caller" %in% colnames(raw))       colnames(raw)[colnames(raw) == "Caller"] <- "CallerName"
  if ("CallerAlias" %in% colnames(raw))  colnames(raw)[colnames(raw) == "CallerAlias"] <- "Caller"
  if ("METRIC.F1_Score" %in% colnames(raw)) colnames(raw)[colnames(raw) == "METRIC.F1_Score"] <- "F1"
  if ("METRIC.Precision" %in% colnames(raw)) colnames(raw)[colnames(raw) == "METRIC.Precision"] <- "Precision"
  if ("METRIC.Recall" %in% colnames(raw)) colnames(raw)[colnames(raw) == "METRIC.Recall"] <- "Recall"

  # Filter: PASS filter, overall (*) subset, main subtypes
  bdata <- raw %>%
    filter(Filter == "PASS", Subtype == "*",
           Subset == "*" | Subset == "")

  variant_types <- c("SNP", "INDEL")

} else {
  cat("Detected: vcfeval format (overall stats only)\n")
  cat("NOTE: vcfeval summary does not separate SNP/INDEL.\n")
  cat("      Heatmap will show overall F-measure only.\n\n")

  # vcfeval columns: Coverage, Caller, CallerAlias, Threshold, True-pos-baseline,
  #                  True-pos-call, False-pos, False-neg, Precision, Sensitivity, F-measure
  colnames(raw) <- gsub("[.-]", "_", colnames(raw))

  raw <- raw %>%
    rename(F1 = F_measure, Recall = Sensitivity) %>%
    mutate(F1 = as.numeric(F1),
           Precision = as.numeric(Precision),
           Recall = as.numeric(Recall))

  if ("CallerAlias" %in% colnames(raw)) {
    raw$Caller <- raw$CallerAlias
  }

  # Use "None" threshold row (unfiltered)
  bdata <- raw %>% filter(Threshold == "None" | grepl("^[0-9]", Threshold))
  variant_types <- c("Overall")
  bdata$Type <- "Overall"
}

# ── Define caller aliases and colors ──
caller_order <- c("HC", "DV", "ST", "FB", "DS")
caller_colors <- c(
  "HC" = "#E64B35",
  "DV" = "#4DBBD5",
  "ST" = "#00A087",
  "FB" = "#3C5488",
  "DS" = "#F39B7F"
)

# Filter to known callers
bdata <- bdata %>% filter(Caller %in% caller_order)
bdata$Caller <- factor(bdata$Caller, levels = caller_order)

cat("Callers found: ", paste(unique(bdata$Caller), collapse = ", "), "\n")
cat("Data points per caller:\n")
print(table(bdata$Caller, bdata$Type))

# ── Pairwise comparison function ──
draw_pairwise_heatmap <- function(df, metric_col = "F1", title = "") {
  callers <- levels(droplevels(df$Caller))
  n_callers <- length(callers)

  if (n_callers < 2) {
    cat("  SKIP: fewer than 2 callers for", title, "\n")
    return(NULL)
  }

  # Build matrices
  diff_mat <- matrix(0, nrow = n_callers, ncol = n_callers,
                     dimnames = list(callers, callers))
  pval_mat <- matrix("", nrow = n_callers, ncol = n_callers,
                     dimnames = list(callers, callers))

  for (i in callers) {
    for (j in callers) {
      if (i == j) next

      left  <- df[[metric_col]][df$Caller == i]
      right <- df[[metric_col]][df$Caller == j]

      # Ensure paired and same length
      n_obs <- min(length(left), length(right))
      if (n_obs < 2) {
        diff_mat[i, j] <- median(left, na.rm = TRUE) - median(right, na.rm = TRUE)
        pval_mat[i, j] <- ""
        next
      }

      left  <- left[1:n_obs]
      right <- right[1:n_obs]

      diff_mat[i, j] <- median(left - right, na.rm = TRUE)

      tryCatch({
        pval <- wilcox.test(left, right, paired = TRUE)$p.value
        pval_mat[i, j] <- ifelse(pval < 0.001, "***",
                           ifelse(pval < 0.01, "**",
                            ifelse(pval < 0.05, "*", "")))
      }, error = function(e) {
        pval_mat[i, j] <<- ""
      })
    }
  }

  # Melt for ggplot
  diff_df <- melt(diff_mat, varnames = c("CallerA", "CallerB"), value.name = "diff")
  pval_df <- melt(pval_mat, varnames = c("CallerA", "CallerB"), value.name = "pval")
  plot_df <- merge(diff_df, pval_df)

  # Create label with stars
  plot_df$label <- ifelse(plot_df$CallerA == plot_df$CallerB, "",
                          paste0(round(plot_df$diff, 4), "\n", plot_df$pval))

  # Symmetric color scale
  max_abs <- max(abs(plot_df$diff), na.rm = TRUE)
  if (max_abs == 0) max_abs <- 0.01

  p <- ggplot(plot_df, aes(x = CallerA, y = CallerB, fill = diff)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = label), size = 3, color = "black") +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-max_abs, max_abs),
      name = "Median\nΔ F1"
    ) +
    labs(title = title, x = "", y = "") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    coord_fixed()

  return(p)
}

# ── Generate plots ──
plots <- list()

for (vtype in variant_types) {
  cat("\nProcessing:", vtype, "\n")

  if (vtype == "Overall") {
    sub_df <- bdata
  } else {
    sub_df <- bdata %>% filter(Type == vtype)
  }

  if (nrow(sub_df) == 0) {
    cat("  SKIP: no data for", vtype, "\n")
    next
  }

  p <- draw_pairwise_heatmap(sub_df, metric_col = "F1",
                              title = paste("Pairwise F1 Comparison —", vtype))
  if (!is.null(p)) {
    plots[[vtype]] <- p

    fname <- paste0("pairwise_wilcoxon_", tolower(vtype), ".pdf")
    ggsave(file.path(output_dir, fname), p, width = 7, height = 6)
    cat("  Saved:", fname, "\n")

    fname_png <- paste0("pairwise_wilcoxon_", tolower(vtype), ".png")
    ggsave(file.path(output_dir, fname_png), p, width = 7, height = 6, dpi = 150)
  }
}

# Combined panel if we have both SNP and INDEL
if (length(plots) >= 2) {
  combined <- plot_grid(plotlist = plots, ncol = 2, labels = c("a", "b"))
  ggsave(file.path(output_dir, "pairwise_wilcoxon_combined.pdf"),
         combined, width = 14, height = 6)
  ggsave(file.path(output_dir, "pairwise_wilcoxon_combined.png"),
         combined, width = 14, height = 6, dpi = 150)
  cat("\nSaved: pairwise_wilcoxon_combined.pdf/png\n")
}

cat("\n✓ Pairwise Wilcoxon analysis complete.\n")
cat("  Output directory:", output_dir, "\n")
