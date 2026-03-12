#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(UpSetR)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- ifelse(length(args) >= 1, args[[1]], getwd())
concordance_dir <- ifelse(length(args) >= 2, args[[2]], file.path(project_dir, "results", "eval", "vcfeval"))

vcfeval_dir <- concordance_dir
scored_dir <- file.path(project_dir, "results", "analysis", "alphagenome")
out_dir <- file.path(project_dir, "results", "plots", "downstream")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(vcfeval_dir)) {
  stop(sprintf("Missing directory: %s", vcfeval_dir))
}

caller_order <- c("HC", "DV", "ST", "FB", "DS")

to_matrix <- function(dt) {
  dt <- unique(dt)
  mat <- dcast(dt, id ~ set, fun.aggregate = length, value.var = "set")
  sets <- setdiff(names(mat), "id")
  mat[, (sets) := lapply(.SD, function(x) as.integer(x > 0)), .SDcols = sets]
  list(mat = mat, sets = sets)
}

special_counts <- function(dt, sets) {
  m <- to_matrix(dt)$mat
  if (nrow(m) == 0) {
    return(list(all5 = 0L, single = data.frame(caller = sets, n = 0L)))
  }
  m[, n_called := rowSums(.SD), .SDcols = sets]
  all5 <- nrow(m[n_called == length(sets)])
  one <- m[n_called == 1]
  single <- data.frame(caller = sets, n = sapply(sets, function(s) sum(one[[s]] == 1L)))
  list(all5 = all5, single = single)
}

median_score_for <- function(pattern, cov) {
  f <- file.path(scored_dir, sprintf("%s_%s_scored.tsv", pattern, cov))
  if (!file.exists(f)) return(NA_real_)
  df <- fread(f)
  if (!"alphagenome_score_mean" %in% names(df) || nrow(df) == 0) return(NA_real_)
  median(df$alphagenome_score_mean, na.rm = TRUE)
}

coverage_files <- list.files(vcfeval_dir, pattern = "^upset_long_callset_\\d+x\\.csv$", full.names = FALSE)
coverages <- sort(unique(gsub("^upset_long_callset_(\\d+x)\\.csv$", "\\1", coverage_files)))

if (length(coverages) == 0) {
  stop("No upset_long_callset_<cov>x.csv files found")
}

for (cov in coverages) {
  fp_file <- file.path(vcfeval_dir, sprintf("upset_long_callset_%s.csv", cov))
  fn_file <- file.path(vcfeval_dir, sprintf("upset_long_baseline_%s.csv", cov))
  if (!file.exists(fp_file) || !file.exists(fn_file)) {
    next
  }

  fp <- fread(fp_file, header = FALSE, col.names = c("set", "id"))
  fn <- fread(fn_file, header = FALSE, col.names = c("set", "id"))

  sets <- intersect(caller_order, sort(unique(c(fp$set, fn$set))))
  if (length(sets) < 2) next

  if (nrow(fp) > 0) {
    x <- to_matrix(fp)
    sp <- special_counts(fp, sets)
    med_all5 <- median_score_for("all_5_fp", cov)

    sub_txt <- sprintf(
      "All-5-FP=%d | Single-caller FP: HC=%d DV=%d ST=%d FB=%d DS=%d | median AG(all_5_fp)=%.4f",
      sp$all5,
      sp$single$n[sp$single$caller == "HC"],
      sp$single$n[sp$single$caller == "DV"],
      sp$single$n[sp$single$caller == "ST"],
      sp$single$n[sp$single$caller == "FB"],
      sp$single$n[sp$single$caller == "DS"],
      ifelse(is.finite(med_all5), med_all5, NA_real_)
    )

    pdf(file.path(out_dir, sprintf("downstream_upset_fp_%s.pdf", cov)), width = 11, height = 8)
    upset(
      as.data.frame(x$mat),
      sets = sets,
      order.by = "freq",
      nintersects = 60,
      keep.order = TRUE,
      main.bar.color = "#CC5A32",
      sets.bar.color = "#3870A0",
      mainbar.y.label = sprintf("FP intersection size (%s)", cov),
      sets.x.label = "FP per caller"
    )
    title(sub = sub_txt, cex.sub = 0.8)
    dev.off()

    png(file.path(out_dir, sprintf("downstream_upset_fp_%s.png", cov)), width = 1400, height = 1000, res = 140)
    upset(
      as.data.frame(x$mat),
      sets = sets,
      order.by = "freq",
      nintersects = 60,
      keep.order = TRUE,
      main.bar.color = "#CC5A32",
      sets.bar.color = "#3870A0",
      mainbar.y.label = sprintf("FP intersection size (%s)", cov),
      sets.x.label = "FP per caller"
    )
    title(sub = sub_txt, cex.sub = 0.8)
    dev.off()
  }

  if (nrow(fn) > 0) {
    x <- to_matrix(fn)
    sp <- special_counts(fn, sets)
    med_all5 <- median_score_for("all_5_fn", cov)

    sub_txt <- sprintf(
      "All-5-FN=%d | Single-caller FN: HC=%d DV=%d ST=%d FB=%d DS=%d | median AG(all_5_fn)=%.4f",
      sp$all5,
      sp$single$n[sp$single$caller == "HC"],
      sp$single$n[sp$single$caller == "DV"],
      sp$single$n[sp$single$caller == "ST"],
      sp$single$n[sp$single$caller == "FB"],
      sp$single$n[sp$single$caller == "DS"],
      ifelse(is.finite(med_all5), med_all5, NA_real_)
    )

    pdf(file.path(out_dir, sprintf("downstream_upset_fn_%s.pdf", cov)), width = 11, height = 8)
    upset(
      as.data.frame(x$mat),
      sets = sets,
      order.by = "freq",
      nintersects = 60,
      keep.order = TRUE,
      main.bar.color = "#325A9A",
      sets.bar.color = "#2D8A6E",
      mainbar.y.label = sprintf("FN intersection size (%s)", cov),
      sets.x.label = "FN per caller"
    )
    title(sub = sub_txt, cex.sub = 0.8)
    dev.off()

    png(file.path(out_dir, sprintf("downstream_upset_fn_%s.png", cov)), width = 1400, height = 1000, res = 140)
    upset(
      as.data.frame(x$mat),
      sets = sets,
      order.by = "freq",
      nintersects = 60,
      keep.order = TRUE,
      main.bar.color = "#325A9A",
      sets.bar.color = "#2D8A6E",
      mainbar.y.label = sprintf("FN intersection size (%s)", cov),
      sets.x.label = "FN per caller"
    )
    title(sub = sub_txt, cex.sub = 0.8)
    dev.off()
  }
}

cat(sprintf("[DONE] Downstream UpSet plots in %s\n", out_dir))
