#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- ifelse(length(args) >= 1, args[[1]], getwd())

vcfeval_dir <- file.path(project_dir, "results", "eval", "vcfeval")
score_dir <- file.path(project_dir, "results", "analysis", "alphagenome")
out_dir <- file.path(project_dir, "results", "plots", "downstream")
out_tab <- file.path(project_dir, "results", "analysis", "summary_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tab, recursive = TRUE, showWarnings = FALSE)

conc_files <- list.files(vcfeval_dir, pattern = "^variant_concordance_\\d+x\\.tsv$", full.names = TRUE)
if (length(conc_files) == 0) {
  stop(sprintf("No variant_concordance_<cov>x.tsv in %s", vcfeval_dir))
}

caller_cols_all <- c("HC", "DV", "ST", "FB", "DS")

# Build long caller-wise status table across coverages
long_rows <- list()
for (f in conc_files) {
  cov <- sub("^variant_concordance_(\\d+x)\\.tsv$", "\\1", basename(f))
  cov_num <- as.integer(sub("x", "", cov))
  dt <- fread(f)
  callers <- intersect(caller_cols_all, names(dt))
  if (length(callers) == 0) next

  for (cl in callers) {
    tmp <- dt[, .(variant_id, truth, called = get(cl))]
    tmp$caller <- cl
    tmp$coverage <- cov
    tmp$coverage_num <- cov_num
    tmp$status <- ifelse(tmp$truth == 1 & tmp$called == 1, "TP",
                         ifelse(tmp$truth == 1 & tmp$called == 0, "FN",
                                ifelse(tmp$truth == 0 & tmp$called == 1, "FP", "TN")))
    long_rows[[paste0(cov, "_", cl)]] <- tmp
  }
}

long_df <- rbindlist(long_rows, fill = TRUE)
if (nrow(long_df) == 0) {
  stop("No concordance data loaded")
}

# Part A: low-coverage errors resolved at high coverage
low_cov <- c(10, 20)
high_cov <- c(30, 50)

low_df <- long_df %>% filter(coverage_num %in% low_cov)
high_df <- long_df %>% filter(coverage_num %in% high_cov)

low_summary <- low_df %>%
  group_by(variant_id, caller, truth) %>%
  summarise(
    any_fn_low = any(status == "FN"),
    any_fp_low = any(status == "FP"),
    .groups = "drop"
  )

high_summary <- high_df %>%
  group_by(variant_id, caller, truth) %>%
  summarise(
    any_tp_high = any(status == "TP"),
    any_tn_high = any(status == "TN"),
    .groups = "drop"
  )

recovery <- low_summary %>%
  left_join(high_summary, by = c("variant_id", "caller", "truth")) %>%
  mutate(
    # For FN: recovered if it was FN at low cov and TP at high cov
    fn_recovered_high = as.integer(truth == 1 & any_fn_low & coalesce(any_tp_high, FALSE)),
    # For FP: recovered either if TN at high cov, OR if variant absent from
    # high-cov concordance (no caller reported it -> effectively TN)
    fp_recovered_high = as.integer(truth == 0 & any_fp_low & coalesce(any_tn_high, TRUE))
  )

recovery_summary <- recovery %>%
  group_by(caller) %>%
  summarise(
    fn_low_to_high_recovered = sum(fn_recovered_high),
    fp_low_to_high_recovered = sum(fp_recovered_high),
    .groups = "drop"
  )

fwrite(recovery_summary, file.path(out_tab, "downstream_coverage_recovery_summary.csv"))

p_recovery <- recovery_summary %>%
  pivot_longer(cols = c(fn_low_to_high_recovered, fp_low_to_high_recovered),
               names_to = "metric", values_to = "n") %>%
  ggplot(aes(x = caller, y = n, fill = metric)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("fn_low_to_high_recovered" = "#2D6DA3", "fp_low_to_high_recovered" = "#CF6A32")) +
  labs(
    title = "Coverage Effect: Low-Coverage Errors Resolved at High Coverage",
    x = "Caller",
    y = "Variant count",
    fill = "Metric"
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(out_dir, "downstream_coverage_recovery_barplot.png"), p_recovery, width = 10, height = 6, dpi = 150)
ggsave(file.path(out_dir, "downstream_coverage_recovery_barplot.pdf"), p_recovery, width = 10, height = 6)

# Part B: per-variant trajectory across all coverages
traj <- long_df %>%
  arrange(variant_id, caller, coverage_num) %>%
  group_by(variant_id, caller) %>%
  summarise(
    n_status_changes = sum(status != lag(status), na.rm = TRUE),
    statuses = paste(status, collapse = "->"),
    .groups = "drop"
  )

# Optional impact ranking from downstream scored tables
score_files <- list.files(score_dir, pattern = "_scored\\.tsv$", full.names = TRUE)
if (length(score_files) > 0) {
  scores <- rbindlist(lapply(score_files, fread), fill = TRUE)
  impact <- scores %>%
    group_by(variant_id) %>%
    summarise(max_ag = max(alphagenome_score_max, na.rm = TRUE), .groups = "drop")
  traj <- traj %>% left_join(impact, by = "variant_id")
} else {
  traj$max_ag <- NA_real_
}

# pick top unstable variants with preference to high impact
traj_rank <- traj %>%
  mutate(max_ag = ifelse(is.finite(max_ag), max_ag, -1)) %>%
  arrange(desc(n_status_changes), desc(max_ag))

selected_variants <- unique(head(traj_rank$variant_id, 30))

plot_df <- long_df %>%
  filter(variant_id %in% selected_variants) %>%
  mutate(
    status_num = case_when(
      status == "TN" ~ 0,
      status == "TP" ~ 1,
      status == "FN" ~ 2,
      status == "FP" ~ 3,
      TRUE ~ NA_real_
    )
  )

fwrite(plot_df, file.path(out_tab, "downstream_variant_trajectory_table.csv"))

p_traj <- ggplot(plot_df, aes(x = coverage_num, y = variant_id, fill = status)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_wrap(~ caller, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("TN" = "#F0F0F0", "TP" = "#2E8B57", "FN" = "#2D6DA3", "FP" = "#CF6A32")) +
  labs(
    title = "Coverage Trajectory of Error Status (Top Unstable Variants)",
    x = "Coverage (x)",
    y = "Variant",
    fill = "Status"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 6)
  )

ggsave(file.path(out_dir, "downstream_coverage_trajectory_heatmap.png"), p_traj, width = 11, height = 14, dpi = 160)
ggsave(file.path(out_dir, "downstream_coverage_trajectory_heatmap.pdf"), p_traj, width = 11, height = 14)

cat(sprintf("[DONE] Coverage sensitivity outputs in %s and %s\n", out_dir, out_tab))
