#!/usr/bin/env Rscript

library(ggplot2)

setwd("~/giab_project/results/")

# ===============================
# Read data
# ===============================
data <- read.csv("caller_comparison.csv")
tstv <- read.csv("tstv_comparison.csv")

# Convert to %
data$Recall    <- data$Recall * 100
data$Precision <- data$Precision * 100
data$F1        <- data$F1 * 100

# ===============================
# FIGURE 1 — Recall vs Precision
# ===============================
p1 <- ggplot(data, aes(x=Caller)) +
  geom_bar(aes(y=Recall, fill="Recall"), stat="identity", position="dodge") +
  geom_bar(aes(y=Precision, fill="Precision"), stat="identity", position="dodge") +
  facet_wrap(~Type) +
  labs(
    title="Recall vs Precision at 40× Coverage",
    y="Percentage (%)",
    fill="Metric"
  ) +
  ylim(0,100) +
  theme_minimal(base_size=13)

ggsave("figure1_recall_precision.png", p1, width=8, height=5, dpi=300)
print(p1)

# ===============================
# FIGURE 2 — F1 Score
# ===============================
p2 <- ggplot(data, aes(x=Caller, y=F1, fill=Caller)) +
  geom_bar(stat="identity") +
  facet_wrap(~Type) +
  labs(
    title="F1-score Comparison Across Variant Callers",
    y="F1-score (%)"
  ) +
  ylim(0,100) +
  theme_minimal(base_size=13)

ggsave("figure2_f1score.png", p2, width=8, height=5, dpi=300)
print(p2)

# ===============================
# FIGURE 3 — False Positives
# ===============================
p3 <- ggplot(data, aes(x=Caller, y=FP, fill=Caller)) +
  geom_bar(stat="identity") +
  facet_wrap(~Type, scales="free_y") +
  labs(
    title="False Positives by Variant Caller",
    y="Number of False Positives"
  ) +
  theme_minimal(base_size=13)

ggsave("figure3_false_positives.png", p3, width=8, height=5, dpi=300)
print(p3)

# ===============================
# FIGURE 4 — Ts/Tv Ratio
# ===============================
truth_tstv <- tstv$Truth_TsTv[1]

p4 <- ggplot(tstv, aes(x=Caller, y=Query_TsTv, fill=Caller)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept=truth_tstv, linetype="dashed", color="red") +
  labs(
    title="Ts/Tv Ratio Compared to Truth",
    subtitle="Dashed line = GIAB truth Ts/Tv",
    y="Ts/Tv Ratio"
  ) +
  theme_minimal(base_size=13)

ggsave("figure4_tstv.png", p4, width=8, height=5, dpi=300)
print(p4)
