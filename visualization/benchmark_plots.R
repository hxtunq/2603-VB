#!/usr/bin/env Rscript
# Benchmark Analysis - Variant Caller Comparison
# Based on: caller_benchmark-main/benchmark_analysis.R
# Adapted for simulated data with 4 callers (DV, FB, HC, ST) + BWA aligner

library(ggplot2)
library(reshape2)
library(dplyr)
library(cowplot)
library(scales)
library(stringr)

# --- Config ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  STATS_FILE <- args[1]
} else {
  STATS_FILE <- 'results/eval/all_stats.tsv'
}

if (length(args) >= 2) {
  OUTPUT_DIR <- args[2]
} else {
  OUTPUT_DIR <- 'results/plots'
}
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("Reading stats from:", STATS_FILE, "\n")
cat("Output directory:", OUTPUT_DIR, "\n")

# --- Color palette ---
mycol1 = rgb(110, 224, 255, maxColorValue = 255)
mycol2 = rgb(255, 192, 78, maxColorValue = 255)
mycol3 = rgb(255, 177, 177, maxColorValue = 255)
mycol4 = rgb(221, 221, 221, maxColorValue = 255)
caller_colors <- c('DV' = '#2196F3', 'FB' = '#FF9800', 'HC' = '#4CAF50', 'ST' = '#E91E63')

# --- Read data ---
bdata <- read.table(STATS_FILE, sep = '\t', header = TRUE)

# Caller aliases
aliases <- data.frame(
  src = c('deepvariant', 'freebayes', 'gatk', 'strelka2'),
  alias = c('DV', 'FB', 'HC', 'ST')
)
if ('Caller' %in% colnames(bdata)) {
  bdata$CallerFilter <- sapply(as.character(bdata$Caller), function(x) {
    m <- aliases[aliases$src == x, 'alias']
    if (length(m) == 0) x else m
  })
}

cat("Data loaded:", nrow(bdata), "rows\n")
cat("Callers:", paste(unique(bdata$CallerFilter), collapse = ', '), "\n")

# --- Filter data ---
pass_data <- bdata[bdata$Filter == 'PASS' & bdata$Subtype == '*' & bdata$Subset == '*', ]
snp_pass <- pass_data[pass_data$Type == 'SNP', ]
indel_pass <- pass_data[pass_data$Type == 'INDEL', ]

all_data <- bdata[bdata$Filter == 'ALL' & bdata$Subtype == '*' & bdata$Subset == '*', ]

# ======================================================================
# Figure 1: F1 Score by caller (SNP + INDEL)
# ======================================================================
cat("Generating F1 Score bar plots...\n")

f1_snp <- ggplot(snp_pass, aes(x = CallerFilter, y = METRIC.F1_Score, fill = CallerFilter)) +
  geom_bar(stat = 'identity', col = 'black', width = 0.6) +
  scale_fill_manual(values = caller_colors) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = 'none') +
  ylab('F1 Score') + xlab('') + ggtitle('SNP')

f1_indel <- ggplot(indel_pass, aes(x = CallerFilter, y = METRIC.F1_Score, fill = CallerFilter)) +
  geom_bar(stat = 'identity', col = 'black', width = 0.6) +
  scale_fill_manual(values = caller_colors) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = 'none') +
  ylab('F1 Score') + xlab('') + ggtitle('INDEL')

f1_plot <- plot_grid(f1_snp, f1_indel, nrow = 1, labels = c('a', 'b'))
ggsave(file.path(OUTPUT_DIR, 'fig1_f1_scores.png'), f1_plot, width = 10, height = 5, dpi = 300)
ggsave(file.path(OUTPUT_DIR, 'fig1_f1_scores.pdf'), f1_plot, width = 10, height = 5)


# =======================================================================
# Figure 2: Precision & Recall comparison
# =======================================================================
cat("Generating Precision & Recall plots...\n")

prec_recall_data <- pass_data[, c('CallerFilter', 'Type', 'METRIC.Precision', 'METRIC.Recall', 'METRIC.F1_Score')]

pr_snp <- ggplot(snp_pass, aes(x = METRIC.Recall, y = METRIC.Precision, col = CallerFilter)) +
  geom_point(size = 5) +
  scale_color_manual(values = caller_colors) +
  scale_x_continuous(limits = c(0.8, 1)) +
  scale_y_continuous(limits = c(0.8, 1)) +
  theme_bw() + theme(panel.grid.minor = element_blank()) +
  xlab('Recall') + ylab('Precision') + ggtitle('SNP') +
  geom_text(aes(label = CallerFilter), vjust = -1, size = 4)

pr_indel <- ggplot(indel_pass, aes(x = METRIC.Recall, y = METRIC.Precision, col = CallerFilter)) +
  geom_point(size = 5) +
  scale_color_manual(values = caller_colors) +
  scale_x_continuous(limits = c(0.5, 1)) +
  scale_y_continuous(limits = c(0.5, 1)) +
  theme_bw() + theme(panel.grid.minor = element_blank()) +
  xlab('Recall') + ylab('Precision') + ggtitle('INDEL') +
  geom_text(aes(label = CallerFilter), vjust = -1, size = 4)

pr_plot <- plot_grid(pr_snp, pr_indel, nrow = 1, labels = c('a', 'b'))
ggsave(file.path(OUTPUT_DIR, 'fig2_precision_recall.png'), pr_plot, width = 12, height = 5, dpi = 300)
ggsave(file.path(OUTPUT_DIR, 'fig2_precision_recall.pdf'), pr_plot, width = 12, height = 5)


# ======================================================================
# Figure 3: Detailed metrics (Precision, Recall, F1) grouped bar chart
# ======================================================================
cat("Generating detailed metrics bar chart...\n")

metrics_long <- melt(pass_data[, c('CallerFilter', 'Type', 'METRIC.Precision', 'METRIC.Recall', 'METRIC.F1_Score')],
                     id.vars = c('CallerFilter', 'Type'))
metrics_long$variable <- gsub('METRIC\\.', '', metrics_long$variable)

detail_snp <- ggplot(metrics_long[metrics_long$Type == 'SNP', ],
                     aes(x = CallerFilter, y = value, fill = variable)) +
  geom_bar(stat = 'identity', position = position_dodge(width = 0.7), col = 'black', width = 0.6) +
  scale_fill_manual(values = c('Precision' = mycol1, 'Recall' = mycol3, 'F1_Score' = mycol2)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = 'bottom') +
  ylab('Metric Value') + xlab('') + ggtitle('SNP')

detail_indel <- ggplot(metrics_long[metrics_long$Type == 'INDEL', ],
                       aes(x = CallerFilter, y = value, fill = variable)) +
  geom_bar(stat = 'identity', position = position_dodge(width = 0.7), col = 'black', width = 0.6) +
  scale_fill_manual(values = c('Precision' = mycol1, 'Recall' = mycol3, 'F1_Score' = mycol2)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = 'bottom') +
  ylab('Metric Value') + xlab('') + ggtitle('INDEL')

detail_plot <- plot_grid(detail_snp, detail_indel, nrow = 1, labels = c('a', 'b'))
ggsave(file.path(OUTPUT_DIR, 'fig3_detail_metrics.png'), detail_plot, width = 12, height = 5, dpi = 300)
ggsave(file.path(OUTPUT_DIR, 'fig3_detail_metrics.pdf'), detail_plot, width = 12, height = 5)


# ======================================================================
# Figure 4: TP / FP / FN counts
# ======================================================================
cat("Generating TP/FP/FN count plots...\n")

if ('TRUTH.TP' %in% colnames(pass_data) & 'QUERY.FP' %in% colnames(pass_data)) {
  count_data <- pass_data[, c('CallerFilter', 'Type', 'TRUTH.TP', 'QUERY.FP', 'TRUTH.FN')]
  count_long <- melt(count_data, id.vars = c('CallerFilter', 'Type'))
  count_long$variable <- gsub('TRUTH\\.', '', gsub('QUERY\\.', '', count_long$variable))

  count_snp <- ggplot(count_long[count_long$Type == 'SNP', ],
                      aes(x = CallerFilter, y = value, fill = variable)) +
    geom_bar(stat = 'identity', position = 'dodge', col = 'black', width = 0.6) +
    scale_fill_manual(values = c('TP' = '#4CAF50', 'FP' = '#F44336', 'FN' = '#FF9800')) +
    theme_bw() + theme(panel.grid = element_blank(), legend.position = 'bottom') +
    ylab('Count') + xlab('') + ggtitle('SNP')

  count_indel <- ggplot(count_long[count_long$Type == 'INDEL', ],
                        aes(x = CallerFilter, y = value, fill = variable)) +
    geom_bar(stat = 'identity', position = 'dodge', col = 'black', width = 0.6) +
    scale_fill_manual(values = c('TP' = '#4CAF50', 'FP' = '#F44336', 'FN' = '#FF9800')) +
    theme_bw() + theme(panel.grid = element_blank(), legend.position = 'bottom') +
    ylab('Count') + xlab('') + ggtitle('INDEL')

  count_plot <- plot_grid(count_snp, count_indel, nrow = 1, labels = c('a', 'b'))
  ggsave(file.path(OUTPUT_DIR, 'fig4_tp_fp_fn.png'), count_plot, width = 12, height = 5, dpi = 300)
  ggsave(file.path(OUTPUT_DIR, 'fig4_tp_fp_fn.pdf'), count_plot, width = 12, height = 5)
}


# ======================================================================
# Figure 5: Filtering impact (ALL vs PASS)
# ======================================================================
cat("Generating filtering impact analysis...\n")

if (nrow(all_data) > 0 & nrow(pass_data) > 0) {
  merged <- merge(all_data[, c('CallerFilter', 'Type', 'METRIC.Precision', 'METRIC.Recall')],
                  pass_data[, c('CallerFilter', 'Type', 'METRIC.Precision', 'METRIC.Recall')],
                  by = c('CallerFilter', 'Type'), suffixes = c('.ALL', '.PASS'))

  merged$PrecisionGain <- merged$METRIC.Precision.PASS - merged$METRIC.Precision.ALL
  merged$RecallLoss <- merged$METRIC.Recall.PASS - merged$METRIC.Recall.ALL

  gain_loss <- melt(merged, id.vars = c('CallerFilter', 'Type'),
                    measure.vars = c('PrecisionGain', 'RecallLoss'))

  filtr_stats <- ggplot(gain_loss, aes(x = CallerFilter, y = value, fill = variable)) +
    geom_bar(stat = 'identity', col = 'black', position = 'dodge', width = 0.6) +
    scale_fill_manual(values = c('PrecisionGain' = mycol1, 'RecallLoss' = mycol3),
                      labels = c('Precision Gain', 'Recall Loss')) +
    theme_bw() + theme(panel.grid = element_blank(), legend.position = 'top') +
    ylab('Change in metric\nvalue after filtering') +
    facet_wrap(~ Type, nrow = 1)
  ggsave(file.path(OUTPUT_DIR, 'fig5_filter_impact.png'), filtr_stats, width = 10, height = 5, dpi = 300)
  ggsave(file.path(OUTPUT_DIR, 'fig5_filter_impact.pdf'), filtr_stats, width = 10, height = 5)
}


# ======================================================================
# Figure 6: Pairwise comparison heatmap of callers
# ======================================================================
cat("Generating pairwise comparison heatmap...\n")

callers <- unique(snp_pass$CallerFilter)
n_callers <- length(callers)

if (n_callers >= 2) {
  library(lattice)
  library(colorspace)

  jet <- diverge_hsv(100)

  drawCoolHM <- function(df, p_df) {
    myPanel_a <- function(x, y, z, ...) {
      panel.levelplot(x, y, z, ...)
      panel.text(x, y, p_df[cbind(x, y)])
    }
    return(levelplot(df, col.regions = jet,
                     at = seq(min(df), max(df), length.out = 100),
                     aspect = 'fill', colorkey = list(width = 3, labels = list(cex = 1.0)),
                     scales = list(x = list(rot = 45)), xlab = list(label = ''),
                     ylab = list(label = ''), panel = myPanel_a))
  }

  # SNP F1 pairwise
  caller_pairwise <- matrix(0, nrow = n_callers, ncol = n_callers)
  rownames(caller_pairwise) <- callers
  colnames(caller_pairwise) <- callers
  pvals_callers <- matrix("", nrow = n_callers, ncol = n_callers)
  rownames(pvals_callers) <- callers
  colnames(pvals_callers) <- callers

  for (i in callers) {
    for (j in callers) {
      left <- as.numeric(snp_pass[snp_pass$CallerFilter == i, 'METRIC.F1_Score'])
      right <- as.numeric(snp_pass[snp_pass$CallerFilter == j, 'METRIC.F1_Score'])

      if (length(left) > 0 & length(right) > 0) {
        difference <- median(left - right)
        caller_pairwise[i, j] <- difference
      }
    }
  }

  png(file.path(OUTPUT_DIR, 'fig6_pairwise_snp_f1.png'), width = 800, height = 600)
  print(drawCoolHM(caller_pairwise, pvals_callers))
  dev.off()

  # INDEL F1 pairwise
  caller_pairwise_id <- matrix(0, nrow = n_callers, ncol = n_callers)
  rownames(caller_pairwise_id) <- callers
  colnames(caller_pairwise_id) <- callers
  pvals_callers_id <- matrix("", nrow = n_callers, ncol = n_callers)
  rownames(pvals_callers_id) <- callers
  colnames(pvals_callers_id) <- callers

  for (i in callers) {
    for (j in callers) {
      left <- as.numeric(indel_pass[indel_pass$CallerFilter == i, 'METRIC.F1_Score'])
      right <- as.numeric(indel_pass[indel_pass$CallerFilter == j, 'METRIC.F1_Score'])

      if (length(left) > 0 & length(right) > 0) {
        difference <- median(left - right)
        caller_pairwise_id[i, j] <- difference
      }
    }
  }

  png(file.path(OUTPUT_DIR, 'fig6_pairwise_indel_f1.png'), width = 800, height = 600)
  print(drawCoolHM(caller_pairwise_id, pvals_callers_id))
  dev.off()
}


# ======================================================================
# Summary Table
# ======================================================================
cat("Writing summary table...\n")

summary_table <- pass_data[, c('CallerFilter', 'Type', 'METRIC.F1_Score',
                                'METRIC.Precision', 'METRIC.Recall')]
write.table(summary_table, file = file.path(OUTPUT_DIR, 'summary_stats.tsv'),
            sep = '\t', row.names = FALSE, quote = FALSE)


cat("\n============================\n")
cat("All plots saved to:", OUTPUT_DIR, "\n")
cat("============================\n")
