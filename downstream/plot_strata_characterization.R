#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- ifelse(length(args) >= 1, args[[1]], getwd())

score_dir <- file.path(project_dir, "results", "analysis", "alphagenome")
out_dir <- file.path(project_dir, "results", "plots", "downstream")
summary_dir <- file.path(project_dir, "results", "analysis", "summary_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(score_dir, pattern = "_scored\\.tsv$", full.names = TRUE)
if (length(files) == 0) {
  stop(sprintf("No *_scored.tsv files found in %s", score_dir))
}

read_scored <- function(f) {
  dt <- fread(f)
  if (!"coverage" %in% names(dt)) {
    dt$coverage <- sub("^.*_(\\d+x)_scored\\.tsv$", "\\1", basename(f))
  }
  if (!"error_type" %in% names(dt)) {
    dt$error_type <- ifelse(grepl("_fn", dt$pattern), "FN", "FP")
  }
  dt
}

all_df <- rbindlist(lapply(files, read_scored), fill = TRUE)
if (nrow(all_df) == 0) {
  stop("Scored files are empty")
}

all_df <- all_df %>%
  mutate(
    gc_bin = ifelse(is.na(gc_bin) | gc_bin == "", "GC_UNK", gc_bin),
    in_cds_canonical = ifelse(is.na(in_cds_canonical), 0, in_cds_canonical),
    cds_label = ifelse(in_cds_canonical == 1, "CDS", "non-CDS")
  )

# 1) Heatmap-like tile: counts by pattern x gc_bin
heat_df <- all_df %>%
  group_by(coverage, error_type, pattern, gc_bin) %>%
  summarise(
    n = n(),
    med_ag = median(alphagenome_score_mean, na.rm = TRUE),
    .groups = "drop"
  )

p_heat <- ggplot(heat_df, aes(x = gc_bin, y = pattern, fill = n)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = n), size = 2.7) +
  facet_grid(error_type ~ coverage, scales = "free") +
  scale_fill_gradient(low = "#E7F2F8", high = "#0D5A8A") +
  labs(
    title = "Downstream Error Patterns by GC Bin",
    subtitle = "Cell label = variant count",
    x = "GC bin",
    y = "Pattern"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "downstream_strata_heatmap_counts.png"), p_heat, width = 14, height = 8, dpi = 150)
ggsave(file.path(out_dir, "downstream_strata_heatmap_counts.pdf"), p_heat, width = 14, height = 8)

# 2) Bar plot: composition by GC bin
bar_df <- all_df %>%
  group_by(coverage, error_type, pattern, gc_bin) %>%
  summarise(n = n(), .groups = "drop")

p_bar <- ggplot(bar_df, aes(x = gc_bin, y = n, fill = pattern)) +
  geom_col(position = "dodge") +
  facet_grid(error_type ~ coverage, scales = "free_y") +
  labs(
    title = "Pattern Composition Across GC Bins",
    x = "GC bin",
    y = "Variant count"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_dir, "downstream_strata_gc_barplot.png"), p_bar, width = 14, height = 8, dpi = 150)
ggsave(file.path(out_dir, "downstream_strata_gc_barplot.pdf"), p_bar, width = 14, height = 8)

# 3) Boxplot: AlphaGenome score by GC and CDS status
score_df <- all_df %>%
  filter(is.finite(alphagenome_score_mean)) %>%
  mutate(gc_bin = factor(gc_bin, levels = sort(unique(gc_bin))))

if (nrow(score_df) > 0) {
  p_box <- ggplot(score_df, aes(x = gc_bin, y = alphagenome_score_mean, fill = cds_label)) +
    geom_boxplot(outlier.size = 0.4, linewidth = 0.3) +
    facet_grid(error_type ~ coverage) +
    scale_fill_manual(values = c("CDS" = "#D95F0E", "non-CDS" = "#1B7837")) +
    labs(
      title = "AlphaGenome Score Distribution by GC Bin and CDS",
      x = "GC bin",
      y = "AlphaGenome score (mean)",
      fill = "Region"
    ) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(out_dir, "downstream_alphagenome_gc_cds_boxplot.png"), p_box, width = 14, height = 8, dpi = 150)
  ggsave(file.path(out_dir, "downstream_alphagenome_gc_cds_boxplot.pdf"), p_box, width = 14, height = 8)
}

# 4) Export summary table for interpretation
summary_tbl <- all_df %>%
  group_by(coverage, error_type, pattern) %>%
  summarise(
    n_variants = n(),
    n_cds = sum(in_cds_canonical == 1, na.rm = TRUE),
    pct_cds = 100 * n_cds / n_variants,
    med_gc = median(gc_pct_mid, na.rm = TRUE),
    med_ag = median(alphagenome_score_mean, na.rm = TRUE),
    p90_ag = quantile(alphagenome_score_mean, probs = 0.9, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(error_type, coverage, pattern)

fwrite(summary_tbl, file.path(summary_dir, "downstream_pattern_interpretation.csv"))

cat(sprintf("[DONE] Strata characterization plots in %s\n", out_dir))
