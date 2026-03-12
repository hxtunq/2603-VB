#!/usr/bin/env Rscript
# =============================================================================
# 03_stratified_plots.R
# Boxplots, Line Plots, and Stratified Performance Analysis
#
# Reads hap.py extended CSV output (with stratification), generates:
#   - Line plots: F1 by GC content per caller
#   - Line plots: F1 by coverage stratum per caller
#   - Line plots: F1 by CDS distance per caller
#   - Box plots: F1 per caller, grouped by coverage
#
# Adapted from: archive/caller_benchmark-main/benchmark_analysis.R (lines 618–800)
#
# Usage:
#   Rscript visualization/03_stratified_plots.R [stats_tsv] [output_dir]
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(cowplot)
  library(stringr)
})

# ── Command-line args ──
args <- commandArgs(trailingOnly = TRUE)

project_dir <- if (file.exists("workflow.md") || dir.exists("visualization")) "." else ".."

stats_file <- if (length(args) >= 1) args[1] else file.path(project_dir, "results", "eval", "all_stats.tsv")
output_dir <- if (length(args) >= 2) args[2] else file.path(project_dir, "results", "plots")

if (!file.exists(stats_file)) {
  stop(paste("Stats file not found:", stats_file,
             "\nRun eval_happy.sh first to generate hap.py output with stratification."))
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cat("Input:  ", stats_file, "\n")
cat("Output: ", output_dir, "\n\n")

# ── Caller settings ──
caller_colors <- c(
  "HC" = "#E64B35", "DV" = "#4DBBD5", "ST" = "#00A087",
  "FB" = "#3C5488", "DS" = "#F39B7F"
)
caller_order <- c("HC", "DV", "ST", "FB", "DS")

# ── Read data ──
bdata <- read.delim(stats_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cat("Loaded", nrow(bdata), "rows\n")
cat("Columns:", paste(colnames(bdata), collapse = ", "), "\n\n")

# Standardize caller name column
if ("Caller" %in% colnames(bdata) && !"CallerAlias" %in% colnames(bdata)) {
  # Map full names to aliases
  alias_map <- c(gatk = "HC", deepvariant = "DV", strelka2 = "ST",
                 freebayes = "FB", dnascope = "DS")
  bdata$CallerAlias <- sapply(bdata$Caller, function(x) {
    if (x %in% names(alias_map)) alias_map[[x]] else x
  })
}
if ("CallerAlias" %in% colnames(bdata)) bdata$Caller <- bdata$CallerAlias

# Standardize metric columns
if ("METRIC.F1_Score" %in% colnames(bdata)) {
  bdata$F1 <- as.numeric(bdata$METRIC.F1_Score)
  bdata$Precision <- as.numeric(bdata$METRIC.Precision)
  bdata$Recall <- as.numeric(bdata$METRIC.Recall)
}

# Filter known callers
bdata <- bdata %>% filter(Caller %in% caller_order)
bdata$Caller <- factor(bdata$Caller, levels = caller_order)

# Extract coverage as numeric
bdata$CovNum <- as.numeric(gsub("[^0-9]", "", bdata$Coverage))

# ── Unique subsets for debugging ──
cat("Unique Subset values:", length(unique(bdata$Subset)), "\n")
cat(head(unique(bdata$Subset), 20), sep = "\n  ")
cat("\n\n")

# =============================================================================
# PLOT 1: Overall F1 Boxplot by Caller × Coverage
# =============================================================================

cat("--- Plot 1: Overall F1 Boxplots ---\n")

overall <- bdata %>%
  filter(Filter == "PASS", Subtype == "*",
         Subset == "*" | Subset == "")

if (nrow(overall) > 0) {
  for (vtype in c("SNP", "INDEL")) {
    sub_df <- overall %>% filter(Type == vtype)
    if (nrow(sub_df) == 0) next

    p <- ggplot(sub_df, aes(x = factor(CovNum), y = F1, fill = Caller)) +
      geom_boxplot(position = position_dodge(0.8), width = 0.7,
                   alpha = 0.85, outlier.size = 1) +
      scale_fill_manual(values = caller_colors) +
      scale_y_continuous(limits = c(max(0, min(sub_df$F1, na.rm = TRUE) - 0.05), 1)) +
      labs(title = paste("F1 Score by Coverage —", vtype),
           x = "Coverage (x)", y = "F1 Score") +
      theme_bw(base_size = 13) +
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "bottom")

    fname <- paste0("boxplot_f1_by_coverage_", tolower(vtype))
    ggsave(file.path(output_dir, paste0(fname, ".pdf")), p, width = 9, height = 6)
    ggsave(file.path(output_dir, paste0(fname, ".png")), p, width = 9, height = 6, dpi = 150)
    cat("  Saved:", fname, "\n")
  }
}

# =============================================================================
# PLOT 2: F1 by GC Content (Line Plot)
# =============================================================================

cat("\n--- Plot 2: F1 by GC Content ---\n")

# Detect GC strata from Subset names
bdata$GC_pct <- as.numeric(
  ifelse(grepl("GC_", bdata$Subset),
         str_extract(bdata$Subset, "\\d+$"),
         NA)
)

# Alternative pattern: GC_0_20, GC_20_30, etc.
gc_range <- str_match(bdata$Subset, "GC_(\\d+)_(\\d+)")
if (any(!is.na(gc_range[, 2]))) {
  bdata$GC_pct <- ifelse(!is.na(gc_range[, 2]),
                          (as.numeric(gc_range[, 2]) + as.numeric(gc_range[, 3])) / 2,
                          bdata$GC_pct)
}

by_gc <- bdata %>%
  filter(!is.na(GC_pct),
         Filter == "PASS", Subtype == "*",
         TRUTH.TOTAL > 10) %>%
  filter(!is.na(F1))

if (nrow(by_gc) > 0) {
  gc_aggr <- by_gc %>%
    group_by(Type, Caller, GC_pct) %>%
    summarise(
      F1_mean = mean(F1, na.rm = TRUE),
      F1_se = sd(F1, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )

  for (vtype in unique(gc_aggr$Type)) {
    sub_df <- gc_aggr %>% filter(Type == vtype)
    if (nrow(sub_df) == 0) next

    p <- ggplot(sub_df, aes(x = GC_pct, y = F1_mean, color = Caller)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2.5) +
      scale_color_manual(values = caller_colors) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_x_continuous(limits = c(0, 100)) +
      labs(title = paste("F1 Score by GC Content —", vtype),
           x = "GC Content (%)", y = "F1 Score") +
      theme_bw(base_size = 13) +
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "bottom")

    fname <- paste0("lineplot_f1_by_gc_", tolower(vtype))
    ggsave(file.path(output_dir, paste0(fname, ".pdf")), p, width = 8, height = 6)
    ggsave(file.path(output_dir, paste0(fname, ".png")), p, width = 8, height = 6, dpi = 150)
    cat("  Saved:", fname, "\n")
  }
} else {
  cat("  SKIP: No GC stratification data found in Subset column\n")
}

# =============================================================================
# PLOT 3: F1 by Coverage Stratum (Line Plot — normalized coverage)
# =============================================================================

cat("\n--- Plot 3: F1 by Coverage Stratum ---\n")

# Extract normalized coverage from subset names like "wgs_mean_cov_0.75_1"
cov_match <- str_match(bdata$Subset, "mean_cov_(\\d+\\.?\\d*)_(\\d+\\.?\\d*|Inf)")
bdata$cov_stratum <- ifelse(!is.na(cov_match[, 2]),
                            (as.numeric(cov_match[, 2]) +
                             ifelse(cov_match[, 3] == "Inf", 3,
                                    as.numeric(cov_match[, 3]))) / 2,
                            NA)

by_cov <- bdata %>%
  filter(!is.na(cov_stratum),
         Filter == "PASS", Subtype == "*",
         TRUTH.TOTAL > 10) %>%
  filter(!is.na(F1))

if (nrow(by_cov) > 0) {
  cov_aggr <- by_cov %>%
    group_by(Type, Caller, cov_stratum) %>%
    summarise(
      F1_mean = mean(F1, na.rm = TRUE),
      F1_se = sd(F1, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )

  for (vtype in unique(cov_aggr$Type)) {
    sub_df <- cov_aggr %>% filter(Type == vtype)
    if (nrow(sub_df) == 0) next

    p <- ggplot(sub_df, aes(x = cov_stratum, y = F1_mean, color = Caller)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2.5) +
      scale_color_manual(values = caller_colors) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(title = paste("F1 by Normalized Coverage —", vtype),
           x = "Normalized Coverage (depth / mean)", y = "F1 Score") +
      theme_bw(base_size = 13) +
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "bottom")

    fname <- paste0("lineplot_f1_by_covstratum_", tolower(vtype))
    ggsave(file.path(output_dir, paste0(fname, ".pdf")), p, width = 8, height = 6)
    ggsave(file.path(output_dir, paste0(fname, ".png")), p, width = 8, height = 6, dpi = 150)
    cat("  Saved:", fname, "\n")
  }
} else {
  cat("  SKIP: No coverage stratification data found\n")
}

# =============================================================================
# PLOT 4: F1 by CDS Vicinity Distance (Line Plot)
# =============================================================================

cat("\n--- Plot 4: F1 by CDS Distance ---\n")

cds_match <- str_match(bdata$Subset, "vicinity_(\\d+)bp")
bdata$cds_distance <- as.numeric(cds_match[, 2])

by_cds <- bdata %>%
  filter(!is.na(cds_distance),
         Filter == "PASS", Subtype == "*",
         TRUTH.TOTAL > 10) %>%
  filter(!is.na(F1))

if (nrow(by_cds) > 0) {
  cds_aggr <- by_cds %>%
    group_by(Type, Caller, cds_distance) %>%
    summarise(
      F1_mean = mean(F1, na.rm = TRUE),
      F1_se = sd(F1, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )

  for (vtype in unique(cds_aggr$Type)) {
    sub_df <- cds_aggr %>% filter(Type == vtype)
    if (nrow(sub_df) == 0) next

    p <- ggplot(sub_df, aes(x = cds_distance, y = F1_mean, color = Caller)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2.5) +
      scale_color_manual(values = caller_colors) +
      scale_y_continuous(limits = c(max(0, min(sub_df$F1_mean, na.rm = TRUE) - 0.1), 1)) +
      labs(title = paste("F1 by CDS Vicinity Distance —", vtype),
           x = "Distance from CDS boundary (bp)", y = "F1 Score") +
      theme_bw(base_size = 13) +
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "bottom")

    fname <- paste0("lineplot_f1_by_cds_distance_", tolower(vtype))
    ggsave(file.path(output_dir, paste0(fname, ".pdf")), p, width = 8, height = 6)
    ggsave(file.path(output_dir, paste0(fname, ".png")), p, width = 8, height = 6, dpi = 150)
    cat("  Saved:", fname, "\n")
  }
} else {
  cat("  SKIP: No CDS vicinity data found\n")
}

# =============================================================================
# PLOT 5: F1 Comparison in CDS Regions (Overall vs. All CDS vs. Canonical CDS)
# =============================================================================

cat("\n--- Plot 5: F1 Comparison in CDS Regions ---\n")

cds_comp <- bdata %>%
  filter(Filter == "PASS", Subtype == "*") %>%
  filter(Subset %in% c("*", "", "proteincoding_only", "CDS_canonical")) %>%
  mutate(SubsetLabel = case_when(
    Subset == "*" | Subset == "" ~ "Overall",
    Subset == "proteincoding_only" ~ "All Protein Coding CDS",
    Subset == "CDS_canonical" ~ "Canonical CDS",
    TRUE ~ "Other"
  ))

cds_comp$SubsetLabel <- factor(cds_comp$SubsetLabel, levels = c("Overall", "All Protein Coding CDS", "Canonical CDS"))

if (nrow(cds_comp) > 0) {
  for (vtype in c("SNP", "INDEL")) {
    sub_df <- cds_comp %>% filter(Type == vtype)
    if (nrow(sub_df) == 0) next

    p <- ggplot(sub_df, aes(x = SubsetLabel, y = F1, fill = Caller)) +
      geom_boxplot(position = position_dodge(0.8), width = 0.7, alpha = 0.85, outlier.size = 1) +
      scale_fill_manual(values = caller_colors) +
      scale_y_continuous(limits = c(max(0, min(sub_df$F1, na.rm = TRUE) - 0.05), 1)) +
      labs(title = paste("F1 Score in CDS Regions (Across All Coverages) —", vtype),
           x = "Genomic Region", y = "F1 Score") +
      theme_bw(base_size = 13) +
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "bottom",
            axis.text.x = element_text(angle = 15, hjust = 1, face = "bold"))

    fname <- paste0("boxplot_f1_cds_comparison_", tolower(vtype))
    ggsave(file.path(output_dir, paste0(fname, ".pdf")), p, width = 9, height = 6)
    ggsave(file.path(output_dir, paste0(fname, ".png")), p, width = 9, height = 6, dpi = 150)
    cat("  Saved:", fname, "\n")
  }
} else {
  cat("  SKIP: No CDS stratification data found for comparison\n")
}

# =============================================================================
# COMBINED FIGURE: Stratification Panel
# =============================================================================

cat("\n--- Assembling combined stratification figure ---\n")

# Collect whatever plots were generated and combine
all_saved <- list.files(output_dir, pattern = "lineplot_.*\\.pdf$")
cat("  Line plots generated:", length(all_saved), "\n")

# =============================================================================
# PLOT 6: F1 by %GC in CDS Region (Line Plot)
# =============================================================================

cat("\n--- Plot 6: F1 by %GC in CDS Region ---\n")

# Detect CDS_GC strata from Subset names: CDS_GC_0_20, CDS_GC_20_30, etc.
cds_gc_range <- str_match(bdata$Subset, "CDS_GC_(\\d+)_(\\d+)")
bdata$cds_gc_pct <- ifelse(!is.na(cds_gc_range[, 2]),
                        (as.numeric(cds_gc_range[, 2]) + as.numeric(cds_gc_range[, 3])) / 2,
                        NA)

by_cds_gc <- bdata %>%
  filter(!is.na(cds_gc_pct),
         Filter == "PASS", Subtype == "*",
         TRUTH.TOTAL > 5) %>%
  filter(!is.na(F1))

if (nrow(by_cds_gc) > 0) {
  cds_gc_aggr <- by_cds_gc %>%
    group_by(Type, Caller, cds_gc_pct) %>%
    summarise(
      F1_mean = mean(F1, na.rm = TRUE),
      F1_se = sd(F1, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )

  for (vtype in unique(cds_gc_aggr$Type)) {
    sub_df <- cds_gc_aggr %>% filter(Type == vtype)
    if (nrow(sub_df) == 0) next

    p_cds_gc <- ggplot(sub_df, aes(x = cds_gc_pct, y = F1_mean, color = Caller)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2.5) +
      scale_color_manual(values = caller_colors) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_x_continuous(limits = c(0, 100)) +
      labs(title = paste("%GC in CDS region —", vtype),
           x = "%GC in CDS region", y = "F1 Score") +
      theme_bw(base_size = 13) +
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "bottom")

    fname <- paste0("lineplot_f1_by_cds_gc_", tolower(vtype))
    ggsave(file.path(output_dir, paste0(fname, ".pdf")), p_cds_gc, width = 8, height = 6)
    ggsave(file.path(output_dir, paste0(fname, ".png")), p_cds_gc, width = 8, height = 6, dpi = 150)
    cat("  Saved:", fname, "\n")
  }
} else {
  cat("  SKIP: No CDS_GC stratification data found in Subset column\n")
}

cat("\n✓ Stratified plots complete.\n")
cat("  Output directory:", output_dir, "\n")

