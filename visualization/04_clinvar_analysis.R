#!/usr/bin/env Rscript
# =============================================================================
# 04_clinvar_analysis.R
# ClinVar Pathogenic Variant Detection Analysis
#
# Reads hap.py output stratified by ClinVar BED region, extracts performance
# on pathogenic/likely_pathogenic variants, generates bar charts showing
# % TP detection rates per caller.
#
# Usage:
#   Rscript visualization/04_clinvar_analysis.R [stats_tsv] [output_dir]
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(cowplot)
})

# ── Command-line args ──
args <- commandArgs(trailingOnly = TRUE)

project_dir <- if (file.exists("workflow.md") || dir.exists("visualization")) "." else ".."

stats_file <- if (length(args) >= 1) args[1] else file.path(project_dir, "results", "eval", "all_stats.tsv")
output_dir <- if (length(args) >= 2) args[2] else file.path(project_dir, "results", "plots")

if (!file.exists(stats_file)) {
  stop(paste("Stats file not found:", stats_file,
             "\nRun eval_happy.sh with stratification first."))
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

# Standardize caller names
if ("Caller" %in% colnames(bdata) && !"CallerAlias" %in% colnames(bdata)) {
  alias_map <- c(gatk = "HC", deepvariant = "DV", strelka2 = "ST",
                 freebayes = "FB", dnascope = "DS")
  bdata$CallerAlias <- sapply(bdata$Caller, function(x) {
    if (x %in% names(alias_map)) alias_map[[x]] else x
  })
}
if ("CallerAlias" %in% colnames(bdata)) bdata$Caller <- bdata$CallerAlias

# Standardize metric names
if ("METRIC.F1_Score" %in% colnames(bdata)) {
  bdata$F1 <- as.numeric(bdata$METRIC.F1_Score)
  bdata$Precision <- as.numeric(bdata$METRIC.Precision)
  bdata$Recall <- as.numeric(bdata$METRIC.Recall)
}
if ("TRUTH.TP" %in% colnames(bdata)) {
  bdata$TP <- as.numeric(bdata$TRUTH.TP)
  bdata$FN <- as.numeric(bdata$TRUTH.FN)
  bdata$FP <- as.numeric(bdata$QUERY.FP)
}

bdata <- bdata %>% filter(Caller %in% caller_order)
bdata$Caller <- factor(bdata$Caller, levels = caller_order)
bdata$CovNum <- as.numeric(gsub("[^0-9]", "", bdata$Coverage))

# ── Filter ClinVar stratum ──
cat("Looking for ClinVar stratum in Subset column...\n")

clinvar_data <- bdata %>%
  filter(grepl("clinvar", Subset, ignore.case = TRUE),
         Filter == "PASS",
         Subtype == "*")

if (nrow(clinvar_data) == 0) {
  cat("\n⚠ No ClinVar stratification data found in hap.py output.\n")
  cat("  Ensure you ran eval_happy.sh after prepare_clinvar.sh.\n")
  cat("  Required: ClinVar BED listed in stratification_chr22.tsv\n")
  cat("\n  Skipping ClinVar analysis.\n")
  q(status = 0)
}

cat("Found", nrow(clinvar_data), "ClinVar-stratified rows\n\n")

# ── Compute detection rate ──
clinvar_stats <- clinvar_data %>%
  mutate(
    Total = TP + FN,
    DetectionRate = ifelse(Total > 0, TP / Total * 100, NA)
  )

# ── PLOT 1: Bar chart — % ClinVar pathogenic variants detected ──
cat("--- Plot 1: ClinVar Detection Rate ---\n")

for (vtype in c("SNP", "INDEL")) {
  sub_df <- clinvar_stats %>% filter(Type == vtype)
  if (nrow(sub_df) == 0) {
    cat("  SKIP:", vtype, "— no data\n")
    next
  }

  p <- ggplot(sub_df, aes(x = Caller, y = DetectionRate, fill = Caller)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3,
             position = position_dodge(0.8), width = 0.7) +
    geom_text(aes(label = paste0(round(DetectionRate, 1), "%")),
              vjust = -0.5, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = caller_colors) +
    scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 20)) +
    facet_wrap(~CovNum, labeller = label_both, nrow = 1) +
    labs(
      title = paste("ClinVar Pathogenic Variant Detection —", vtype),
      subtitle = paste0("Variants in ClinVar (Pathogenic/Likely_pathogenic) on chr22"),
      x = "", y = "Detection Rate (%)"
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold")
    )

  fname <- paste0("clinvar_detection_rate_", tolower(vtype))
  ggsave(file.path(output_dir, paste0(fname, ".pdf")), p, width = 10, height = 6)
  ggsave(file.path(output_dir, paste0(fname, ".png")), p, width = 10, height = 6, dpi = 150)
  cat("  Saved:", fname, "\n")
}

# ── PLOT 2: Summary table of counts ──
cat("\n--- ClinVar Detection Summary ---\n")

summary_df <- clinvar_stats %>%
  group_by(Caller, Type, CovNum) %>%
  summarise(TP = sum(TP, na.rm = TRUE),
            FN = sum(FN, na.rm = TRUE),
            FP = sum(FP, na.rm = TRUE),
            Total = sum(Total, na.rm = TRUE),
            DetectionRate = ifelse(sum(Total, na.rm = TRUE) > 0,
                                   sum(TP, na.rm = TRUE) / sum(Total, na.rm = TRUE) * 100, NA),
            .groups = "drop")

cat("\n")
print(as.data.frame(summary_df))

# Save summary
write.csv(summary_df, file.path(output_dir, "clinvar_detection_summary.csv"),
          row.names = FALSE)
cat("\n  Saved: clinvar_detection_summary.csv\n")

# ── PLOT 3: Heatmap Caller × Coverage for ClinVar Recall ──
cat("\n--- Plot 3: ClinVar Recall Heatmap ---\n")

for (vtype in c("SNP", "INDEL")) {
  sub_df <- summary_df %>% filter(Type == vtype)
  if (nrow(sub_df) == 0) next

  p <- ggplot(sub_df, aes(x = factor(CovNum), y = Caller, fill = DetectionRate)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = paste0(round(DetectionRate, 1), "%\n(",
                                  TP, "/", Total, ")")),
              size = 3.5, color = "black") +
    scale_fill_gradient(low = "#FEE0D2", high = "#CB181D",
                        limits = c(0, 100), name = "Detection\nRate (%)") +
    labs(title = paste("ClinVar Detection —", vtype),
         x = "Coverage (x)", y = "") +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(face = "bold")
    )

  fname <- paste0("clinvar_heatmap_", tolower(vtype))
  ggsave(file.path(output_dir, paste0(fname, ".pdf")), p, width = 8, height = 5)
  ggsave(file.path(output_dir, paste0(fname, ".png")), p, width = 8, height = 5, dpi = 150)
  cat("  Saved:", fname, "\n")
}

cat("\n✓ ClinVar analysis complete.\n")
cat("  Output directory:", output_dir, "\n")
