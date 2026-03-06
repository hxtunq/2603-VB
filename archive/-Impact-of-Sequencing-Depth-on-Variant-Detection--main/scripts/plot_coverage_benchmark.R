setwd("~/giab_project/results/")
library(ggplot2)


# Figure 1 — Recall vs Coverage
# Simple performance data
df <- data.frame(
  coverage = c(2, 10, 40, 80),
  SNP_Recall = c(21.91, 74.48, 95.07, 97.15),
  Indel_Recall = c(11.76, 52.94, 85.29, 92.65)
)

# Convert to long format manually (simpler than tidyr)
df_long <- data.frame(
  coverage = rep(df$coverage, 2),
  recall = c(df$SNP_Recall, df$Indel_Recall),
  variant = rep(c("SNP", "Indel"), each = 4)
)

p1 <- ggplot(df_long, aes(x = coverage, y = recall, color = variant)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = c(2, 10, 40, 80)) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(
    title = "Variant Recall Increases with Coverage",
    subtitle = "HG002 chr22 exome benchmark",
    x = "Coverage (×)",
    y = "Recall (%)",
    color = "Variant type"
  ) +
  theme_minimal(base_size = 13)

ggsave("figure1_recall_vs_coverage.png", p1, width = 7, height = 5, dpi = 300)
print(p1)


# Figure 2 — False Negatives vs Coverage
# Error data (FN only – the most important error)
error_df <- data.frame(
  coverage = c(2, 10, 40, 80),
  SNP_FN = c(713, 233, 45, 26),
  Indel_FN = c(60, 32, 10, 5)
)

error_long <- data.frame(
  coverage = rep(error_df$coverage, 2),
  FN = c(error_df$SNP_FN, error_df$Indel_FN),
  variant = rep(c("SNP", "Indel"), each = 4)
)

p2 <- ggplot(error_long, aes(x = factor(coverage), y = FN, fill = variant)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "False Negatives Decrease with Coverage",
    subtitle = "Low coverage misses real variants",
    x = "Coverage (×)",
    y = "Number of false negatives",
    fill = "Variant type"
  ) +
  theme_minimal(base_size = 13)

ggsave("figure2_false_negatives.png", p2, width = 7, height = 5, dpi = 300)
print(p2)

# Figure 3 — Precision vs Coverage

precision_df <- data.frame(
  coverage = c(2, 10, 40, 80),
  SNP_Precision = c(60.79, 91.52, 97.64, 97.37),
  Indel_Precision = c(53.33, 72.00, 88.06, 87.67)
)

precision_long <- data.frame(
  coverage = rep(precision_df$coverage, 2),
  precision = c(precision_df$SNP_Precision, precision_df$Indel_Precision),
  variant = rep(c("SNP", "Indel"), each = 4)
)

p3 <- ggplot(precision_long, aes(x = coverage, y = precision, color = variant)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = c(2, 10, 40, 80)) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(
    title = "Variant Precision Improves with Coverage",
    subtitle = "HG002 chr22 exome benchmark",
    x = "Coverage (×)",
    y = "Precision (%)",
    color = "Variant type"
  ) +
  theme_minimal(base_size = 13)

ggsave("figure3_precision_vs_coverage.png", p3, width = 7, height = 5, dpi = 300)
print(p3)

# Figure 4 — False Positives vs Coverage

fp_df <- data.frame(
  coverage = c(2, 10, 40, 80),
  SNP_FP = c(129, 63, 21, 24),
  Indel_FP = c(7, 14, 8, 9)
)

fp_long <- data.frame(
  coverage = rep(fp_df$coverage, 2),
  FP = c(fp_df$SNP_FP, fp_df$Indel_FP),
  variant = rep(c("SNP", "Indel"), each = 4)
)

p4 <- ggplot(fp_long, aes(x = factor(coverage), y = FP, fill = variant)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "False Positives Across Coverage Levels",
    subtitle = "Higher coverage does not always reduce false positives",
    x = "Coverage (×)",
    y = "Number of false positives",
    fill = "Variant type"
  ) +
  theme_minimal(base_size = 13)

ggsave("figure4_false_positives.png", p4, width = 7, height = 5, dpi = 300)
print(p4)

