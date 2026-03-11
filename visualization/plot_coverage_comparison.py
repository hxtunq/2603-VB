#!/usr/bin/env python3
"""
Coverage-impact visualisation for variant calling benchmarks.

Reads the aggregated vcfeval stats TSV and generates plots that highlight
how each caller's performance changes across coverages (10x → 50x).

Plots generated:
  1. Line plot   – F1 / Precision / Recall vs Coverage per caller
  2. Delta bars  – Metric improvement from lowest to highest coverage
  3. Radar chart – Multi-metric comparison per coverage
  4. Stacked area – FN + FP composition across coverages

Usage:
  python visualization/plot_coverage_comparison.py [STATS_TSV] [OUTPUT_DIR]
"""

import os
import sys
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd

# ── Defaults ────────────────────────────────────────────────────────────────
STATS_FILE = sys.argv[1] if len(sys.argv) > 1 else "results/eval/vcfeval_all_stats.tsv"
OUTPUT_DIR = sys.argv[2] if len(sys.argv) > 2 else "results/plots"
os.makedirs(OUTPUT_DIR, exist_ok=True)

COLORS = {
    "HC": "#4CAF50",
    "DV": "#2196F3",
    "ST": "#E91E63",
    "FB": "#FF9800",
    "DS": "#009688",
    "DSF": "#7B1FA2",
}
MARKERS = {"HC": "o", "DV": "s", "ST": "D", "FB": "^", "DS": "v", "DSF": "P"}
CALLER_ORDER = ["HC", "DV", "ST", "FB", "DS", "DSF"]


def alias_caller(name: str) -> str:
    name = name.lower().strip()
    if name.startswith("deepvariant"):
        return "DV"
    if name.startswith("freebayes"):
        return "FB"
    if name.startswith("gatk"):
        return "HC"
    if name.startswith("strelka2"):
        return "ST"
    if name.startswith("dnascope_fastq"):
        return "DSF"
    if name.startswith("dnascope"):
        return "DS"
    return name.upper()


def parse_coverage(cov_str: str) -> int:
    """Extract numeric coverage from strings like '30x' or '30'."""
    m = re.search(r"(\d+)", str(cov_str))
    return int(m.group(1)) if m else 0


# ── Load data ───────────────────────────────────────────────────────────────
print(f"Reading: {STATS_FILE}")
df = pd.read_csv(STATS_FILE, sep="\t")

# Normalise column names — vcfeval summary.txt uses different header names
col_map = {
    "Sensitivity": "Recall",
    "F-measure": "F1",
}
df.rename(columns=col_map, inplace=True)

# Ensure we have a CallerAlias column
if "CallerAlias" not in df.columns and "Caller" in df.columns:
    df["CallerAlias"] = df["Caller"].map(alias_caller)

# Parse coverage as integer
df["cov_int"] = df["Coverage"].apply(parse_coverage)

# Keep only the "None" threshold row (overall metrics per caller)
# vcfeval summary.txt has rows for different score thresholds;
# the "None" row is the unfiltered overall metric.
if "Threshold" in df.columns:
    df_overall = df[df["Threshold"].astype(str).str.strip().str.lower() == "none"].copy()
    if df_overall.empty:
        # Fallback: keep last threshold row per group
        df_overall = df.groupby(["Coverage", "Caller"]).last().reset_index()
else:
    df_overall = df.copy()

# Cast numeric columns
for col in ["Precision", "Recall", "F1", "True-pos-baseline", "True-pos-call",
            "False-pos", "False-neg"]:
    if col in df_overall.columns:
        df_overall[col] = pd.to_numeric(df_overall[col], errors="coerce")

coverages = sorted(df_overall["cov_int"].unique())
callers = [c for c in CALLER_ORDER if c in df_overall["CallerAlias"].unique()]

print(f"Callers:   {callers}")
print(f"Coverages: {coverages}")


# =============================================================================
# PLOT 1: Line plot — F1 / Precision / Recall vs Coverage
# =============================================================================
def plot_metric_lines() -> None:
    """One figure with 3 panels (F1, Precision, Recall), one line per caller."""
    metrics = [("F1", "F1 Score"), ("Precision", "Precision"), ("Recall", "Recall")]
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)

    for ax, (col, label) in zip(axes, metrics):
        if col not in df_overall.columns:
            ax.axis("off")
            continue
        for caller in callers:
            sub = df_overall[df_overall["CallerAlias"] == caller].sort_values("cov_int")
            ax.plot(sub["cov_int"], sub[col],
                    color=COLORS.get(caller, "#999"),
                    marker=MARKERS.get(caller, "o"),
                    linewidth=2.2, markersize=8,
                    label=caller)
        ax.set_xlabel("Coverage (×)", fontsize=12)
        ax.set_ylabel(label, fontsize=12)
        ax.set_title(label, fontsize=13, fontweight="bold")
        ax.set_xticks(coverages)
        ax.set_ylim(0, 1.05)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", alpha=0.3)
        if ax == axes[-1]:
            ax.legend(fontsize=10, title="Caller", title_fontsize=11)

    plt.suptitle("Performance vs Coverage", fontsize=15, fontweight="bold", y=1.02)
    plt.tight_layout()
    out = os.path.join(OUTPUT_DIR, "coverage_lines.png")
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# =============================================================================
# PLOT 2: Delta bar chart — Improvement from min-cov to max-cov
# =============================================================================
def plot_delta_bars() -> None:
    """Grouped bars showing Δ(metric) from lowest to highest coverage."""
    min_cov, max_cov = min(coverages), max(coverages)
    metrics = [("F1", "ΔF1"), ("Precision", "ΔPrecision"), ("Recall", "ΔRecall")]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    bar_width = 0.6

    for ax, (col, label) in zip(axes, metrics):
        if col not in df_overall.columns:
            ax.axis("off")
            continue
        deltas = []
        colors_list = []
        labels_list = []
        for caller in callers:
            low = df_overall[(df_overall["CallerAlias"] == caller) &
                             (df_overall["cov_int"] == min_cov)][col]
            high = df_overall[(df_overall["CallerAlias"] == caller) &
                              (df_overall["cov_int"] == max_cov)][col]
            delta = (high.values[0] - low.values[0]) if (len(low) > 0 and len(high) > 0) else 0
            deltas.append(delta)
            colors_list.append(COLORS.get(caller, "#999"))
            labels_list.append(caller)

        x = range(len(callers))
        bars = ax.bar(x, deltas, bar_width, color=colors_list, edgecolor="black", linewidth=0.8)
        ax.set_xticks(list(x))
        ax.set_xticklabels(labels_list, fontsize=11)
        ax.set_ylabel(label, fontsize=12)
        ax.set_title(f"{label} ({min_cov}x → {max_cov}x)", fontsize=12, fontweight="bold")
        ax.axhline(0, color="grey", linewidth=0.8, linestyle="--")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # Annotate values
        for bar, val in zip(bars, deltas):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.002,
                    f"{val:+.4f}", ha="center", va="bottom", fontsize=9, fontweight="bold")

    plt.suptitle(f"Metric Improvement: {min_cov}x → {max_cov}x",
                 fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    out = os.path.join(OUTPUT_DIR, "coverage_delta_bars.png")
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# =============================================================================
# PLOT 3: Radar chart — Multi-metric per coverage
# =============================================================================
def plot_radar() -> None:
    """One radar per coverage; each caller is a polygon."""
    radar_metrics = [c for c in ["F1", "Precision", "Recall"] if c in df_overall.columns]
    if len(radar_metrics) < 2:
        return

    n_metrics = len(radar_metrics)
    angles = np.linspace(0, 2 * np.pi, n_metrics, endpoint=False).tolist()
    angles += angles[:1]  # close polygon

    n_covs = len(coverages)
    fig, axes = plt.subplots(1, n_covs, figsize=(5 * n_covs, 5),
                             subplot_kw=dict(projection="polar"))
    if n_covs == 1:
        axes = [axes]

    for ax, cov in zip(axes, coverages):
        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)
        ax.set_rlabel_position(30)
        ax.set_ylim(0, 1.05)
        ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_yticklabels(["0.2", "0.4", "0.6", "0.8", "1.0"], fontsize=8)
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(radar_metrics, fontsize=11, fontweight="bold")
        ax.set_title(f"{cov}x", fontsize=13, fontweight="bold", pad=20)

        for caller in callers:
            sub = df_overall[(df_overall["CallerAlias"] == caller) &
                             (df_overall["cov_int"] == cov)]
            if sub.empty:
                continue
            vals = [sub[m].values[0] if m in sub.columns else 0 for m in radar_metrics]
            vals += vals[:1]
            ax.plot(angles, vals, color=COLORS.get(caller, "#999"),
                    linewidth=2, marker="o", markersize=5, label=caller)
            ax.fill(angles, vals, color=COLORS.get(caller, "#999"), alpha=0.1)

        if ax == axes[-1]:
            ax.legend(loc="upper right", bbox_to_anchor=(1.35, 1.1), fontsize=9)

    plt.suptitle("Caller Profiles by Coverage", fontsize=15, fontweight="bold", y=1.05)
    plt.tight_layout()
    out = os.path.join(OUTPUT_DIR, "coverage_radar.png")
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# =============================================================================
# PLOT 4: Stacked area — FN + FP across coverages
# =============================================================================
def plot_fn_fp_area() -> None:
    """One panel per caller; stacked areas for FN and FP across coverages."""
    if "False-neg" not in df_overall.columns or "False-pos" not in df_overall.columns:
        print("  SKIP stacked area — no FN/FP columns")
        return

    n_callers = len(callers)
    fig, axes = plt.subplots(1, n_callers, figsize=(5 * n_callers, 5), sharey=True)
    if n_callers == 1:
        axes = [axes]

    for ax, caller in zip(axes, callers):
        sub = df_overall[df_overall["CallerAlias"] == caller].sort_values("cov_int")
        if sub.empty:
            ax.axis("off")
            continue
        covs = sub["cov_int"].values
        fn = sub["False-neg"].values.astype(float)
        fp = sub["False-pos"].values.astype(float)

        ax.fill_between(covs, 0, fn, color="#FF9800", alpha=0.7, label="FN")
        ax.fill_between(covs, fn, fn + fp, color="#F44336", alpha=0.7, label="FP")
        ax.plot(covs, fn + fp, color="black", linewidth=1.5)

        ax.set_xlabel("Coverage (×)", fontsize=11)
        ax.set_title(caller, fontsize=13, fontweight="bold",
                     color=COLORS.get(caller, "#333"))
        ax.set_xticks(coverages)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        if ax == axes[0]:
            ax.set_ylabel("Variant Count", fontsize=11)
            ax.legend(fontsize=10)

    plt.suptitle("Error Composition (FN + FP) vs Coverage",
                 fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    out = os.path.join(OUTPUT_DIR, "coverage_fn_fp_area.png")
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# =============================================================================
# PLOT 5: Heatmap — F1 callers × coverage (quick overview)
# =============================================================================
def plot_f1_heatmap() -> None:
    """Single heatmap: callers (y) × coverages (x) → F1 score."""
    if "F1" not in df_overall.columns:
        return

    matrix = []
    for caller in callers:
        row = []
        for cov in coverages:
            cell = df_overall[(df_overall["CallerAlias"] == caller) &
                              (df_overall["cov_int"] == cov)]
            row.append(cell["F1"].values[0] if len(cell) > 0 else 0)
        matrix.append(row)

    data = np.array(matrix)
    fig, ax = plt.subplots(figsize=(max(8, 2 * len(coverages)), max(4, len(callers))))
    im = ax.imshow(data, cmap="RdYlGn", vmin=0.0, vmax=1.0, aspect="auto")

    ax.set_xticks(range(len(coverages)))
    ax.set_yticks(range(len(callers)))
    ax.set_xticklabels([f"{c}x" for c in coverages], fontsize=12)
    ax.set_yticklabels(callers, fontsize=12)
    ax.set_xlabel("Coverage", fontsize=12)

    for i in range(len(callers)):
        for j in range(len(coverages)):
            val = data[i, j]
            color = "white" if val < 0.5 else "black"
            ax.text(j, i, f"{val:.4f}", ha="center", va="center",
                    fontsize=11, color=color, fontweight="bold")

    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("F1 Score", fontsize=11)

    ax.set_title("F1 Score — Callers × Coverage", fontsize=14, fontweight="bold")
    plt.tight_layout()
    out = os.path.join(OUTPUT_DIR, "coverage_f1_heatmap.png")
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


# =============================================================================
# Main
# =============================================================================
print("\nGenerating coverage comparison plots...")
plot_metric_lines()
plot_delta_bars()
plot_radar()
plot_fn_fp_area()
plot_f1_heatmap()
print(f"\nAll coverage comparison plots saved to: {OUTPUT_DIR}/")
