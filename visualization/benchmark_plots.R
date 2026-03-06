#!/usr/bin/env Rscript
# Benchmark plots per Mode + Coverage combination.

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(cowplot)
})

args <- commandArgs(trailingOnly = TRUE)
STATS_FILE <- if (length(args) >= 1) args[1] else "results/eval/all_stats.tsv"
OUTPUT_DIR <- if (length(args) >= 2) args[2] else "results/plots"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("Reading stats from:", STATS_FILE, "\n")
cat("Output directory:", OUTPUT_DIR, "\n")

bdata <- read.table(STATS_FILE, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

caller_alias <- function(x) {
  if (grepl("^deepvariant", x)) return("DV")
  if (grepl("^freebayes", x)) return("FB")
  if (grepl("^gatk", x)) return("HC")
  if (grepl("^strelka2", x)) return("ST")
  if (grepl("^dnascope", x)) return("DS")
  x
}

caller_levels <- c("HC", "DV", "ST", "FB", "DS")
caller_colors <- c(
  "HC" = "#4CAF50",
  "DV" = "#2196F3",
  "ST" = "#E91E63",
  "FB" = "#FF9800",
  "DS" = "#009688"
)

bdata$CallerFilter <- vapply(bdata$Caller, caller_alias, character(1))
bdata$CallerFilter <- factor(bdata$CallerFilter, levels = caller_levels)

pass_data <- subset(bdata, Filter == "PASS" & Subtype == "*" & Subset == "*")
summary_table <- pass_data[, c("Coverage", "Mode", "CallerFilter", "Type",
                               "METRIC.F1_Score", "METRIC.Precision", "METRIC.Recall")]
write.table(summary_table,
            file = file.path(OUTPUT_DIR, "summary_stats.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

sanitize_tag <- function(mode, coverage) {
  tag <- paste(tolower(mode), coverage, sep = "_")
  gsub("[^A-Za-z0-9_]", "_", tag)
}

plot_f1 <- function(df, mode, coverage, tag) {
  snp_df <- subset(df, Type == "SNP")
  indel_df <- subset(df, Type == "INDEL")

  if (nrow(snp_df) == 0 || nrow(indel_df) == 0) return()

  p1 <- ggplot(snp_df, aes(x = CallerFilter, y = METRIC.F1_Score, fill = CallerFilter)) +
    geom_col(color = "black", width = 0.65) +
    scale_fill_manual(values = caller_colors, drop = FALSE) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
    ylab("F1 Score") + xlab("") + ggtitle("SNP")

  p2 <- ggplot(indel_df, aes(x = CallerFilter, y = METRIC.F1_Score, fill = CallerFilter)) +
    geom_col(color = "black", width = 0.65) +
    scale_fill_manual(values = caller_colors, drop = FALSE) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_bw() + theme(panel.grid = element_blank(), legend.position = "none") +
    ylab("F1 Score") + xlab("") + ggtitle("INDEL")

  combined <- plot_grid(p1, p2, nrow = 1, labels = c("a", "b"))
  title <- ggdraw() + draw_label(sprintf("F1 Scores — %s %s", mode, coverage), fontface = "bold")
  final_plot <- plot_grid(title, combined, ncol = 1, rel_heights = c(0.1, 1))

  ggsave(file.path(OUTPUT_DIR, sprintf("fig_f1_%s.png", tag)), final_plot, width = 11, height = 5, dpi = 300)
}

plot_pr <- function(df, mode, coverage, tag) {
  snp_df <- subset(df, Type == "SNP")
  indel_df <- subset(df, Type == "INDEL")

  if (nrow(snp_df) == 0 || nrow(indel_df) == 0) return()

  p1 <- ggplot(snp_df, aes(x = METRIC.Recall, y = METRIC.Precision, color = CallerFilter)) +
    geom_point(size = 4) +
    geom_text(aes(label = CallerFilter), vjust = -1, size = 3.5) +
    scale_color_manual(values = caller_colors, drop = FALSE) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") +
    xlab("Recall") + ylab("Precision") + ggtitle("SNP")

  p2 <- ggplot(indel_df, aes(x = METRIC.Recall, y = METRIC.Precision, color = CallerFilter)) +
    geom_point(size = 4) +
    geom_text(aes(label = CallerFilter), vjust = -1, size = 3.5) +
    scale_color_manual(values = caller_colors, drop = FALSE) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none") +
    xlab("Recall") + ylab("Precision") + ggtitle("INDEL")

  combined <- plot_grid(p1, p2, nrow = 1, labels = c("a", "b"))
  title <- ggdraw() + draw_label(sprintf("Precision vs Recall — %s %s", mode, coverage), fontface = "bold")
  final_plot <- plot_grid(title, combined, ncol = 1, rel_heights = c(0.1, 1))

  ggsave(file.path(OUTPUT_DIR, sprintf("fig_precision_recall_%s.png", tag)),
         final_plot, width = 12, height = 5, dpi = 300)
}

plot_counts <- function(df, mode, coverage, tag) {
  if (!all(c("TRUTH.TP", "QUERY.FP", "TRUTH.FN") %in% colnames(df))) return()

  count_df <- df[, c("CallerFilter", "Type", "TRUTH.TP", "QUERY.FP", "TRUTH.FN")]
  count_long <- melt(count_df, id.vars = c("CallerFilter", "Type"))
  count_long$variable <- gsub("TRUTH\\.", "", gsub("QUERY\\.", "", count_long$variable))

  p <- ggplot(count_long, aes(x = CallerFilter, y = value, fill = variable)) +
    geom_col(position = "dodge", color = "black", width = 0.65) +
    facet_wrap(~Type, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = c("TP" = "#4CAF50", "FP" = "#F44336", "FN" = "#FF9800")) +
    theme_bw() + theme(panel.grid = element_blank()) +
    xlab("") + ylab("Count") +
    ggtitle(sprintf("TP / FP / FN — %s %s", mode, coverage))

  ggsave(file.path(OUTPUT_DIR, sprintf("fig_counts_%s.png", tag)), p, width = 10, height = 5, dpi = 300)
}

combos <- unique(pass_data[, c("Mode", "Coverage")])
combos <- combos[order(combos$Mode, combos$Coverage), ]

for (i in seq_len(nrow(combos))) {
  mode <- combos$Mode[i]
  coverage <- combos$Coverage[i]
  subset_df <- subset(pass_data, Mode == mode & Coverage == coverage)

  if (nrow(subset_df) == 0) next

  tag <- sanitize_tag(mode, coverage)
  plot_f1(subset_df, mode, coverage, tag)
  plot_pr(subset_df, mode, coverage, tag)
  plot_counts(subset_df, mode, coverage, tag)
}

cat("Generated plot sets for", nrow(combos), "mode/coverage combinations\n")
