#!/usr/bin/env python3
"""
Quick summary plots for the variant calling benchmark.
Generates one plot set per Mode + Coverage combination.
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


STATS_FILE = sys.argv[1] if len(sys.argv) > 1 else "results/eval/all_stats.tsv"
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


def alias_caller(name: str) -> str:
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
    return name


def sanitize(mode: str, coverage: str) -> str:
    return "".join(ch if ch.isalnum() or ch == "_" else "_" for ch in f"{mode.lower()}_{coverage}")


def ordered_subset(df: pd.DataFrame) -> pd.DataFrame:
    order = {"HC": 0, "DV": 1, "ST": 2, "FB": 3, "DS": 4, "DSF": 5}
    return df.assign(_order=df["CallerAlias"].map(order).fillna(99)).sort_values("_order")


def plot_f1(df: pd.DataFrame, mode: str, coverage: str, tag: str) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, vtype in zip(axes, ["SNP", "INDEL"]):
        subset = ordered_subset(df[df["Type"] == vtype])
        if subset.empty:
            ax.axis("off")
            continue
        colors = [COLORS.get(c, "#999999") for c in subset["CallerAlias"]]
        ax.bar(subset["CallerAlias"], subset["METRIC.F1_Score"], color=colors, edgecolor="black", width=0.6)
        ax.set_ylim(0, 1.05)
        ax.set_title(vtype)
        ax.set_ylabel("F1 Score")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    plt.suptitle(f"F1 Scores — {mode} {coverage}", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"summary_f1_{tag}.png"), dpi=300, bbox_inches="tight")
    plt.close()


def plot_precision_recall(df: pd.DataFrame, mode: str, coverage: str, tag: str) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, vtype in zip(axes, ["SNP", "INDEL"]):
        subset = ordered_subset(df[df["Type"] == vtype])
        if subset.empty:
            ax.axis("off")
            continue
        x = range(len(subset))
        ax.bar([i - 0.18 for i in x], subset["METRIC.Precision"], 0.36, label="Precision",
               color="#64B5F6", edgecolor="black")
        ax.bar([i + 0.18 for i in x], subset["METRIC.Recall"], 0.36, label="Recall",
               color="#FFB74D", edgecolor="black")
        ax.set_xticks(list(x))
        ax.set_xticklabels(subset["CallerAlias"])
        ax.set_ylim(0, 1.05)
        ax.set_title(vtype)
        ax.set_ylabel("Value")
        ax.legend()
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    plt.suptitle(f"Precision & Recall — {mode} {coverage}", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"summary_precision_recall_{tag}.png"), dpi=300, bbox_inches="tight")
    plt.close()


def plot_counts(df: pd.DataFrame, mode: str, coverage: str, tag: str) -> None:
    required = {"TRUTH.TP", "QUERY.FP", "TRUTH.FN"}
    if not required.issubset(df.columns):
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, vtype in zip(axes, ["SNP", "INDEL"]):
        subset = ordered_subset(df[df["Type"] == vtype])
        if subset.empty:
            ax.axis("off")
            continue
        x = range(len(subset))
        ax.bar([i - 0.25 for i in x], subset["TRUTH.TP"], 0.25, label="TP", color="#4CAF50", edgecolor="black")
        ax.bar([i for i in x], subset["QUERY.FP"], 0.25, label="FP", color="#F44336", edgecolor="black")
        ax.bar([i + 0.25 for i in x], subset["TRUTH.FN"], 0.25, label="FN", color="#FF9800", edgecolor="black")
        ax.set_xticks(list(x))
        ax.set_xticklabels(subset["CallerAlias"])
        ax.set_title(vtype)
        ax.set_ylabel("Count")
        ax.legend()
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    plt.suptitle(f"TP / FP / FN — {mode} {coverage}", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"summary_counts_{tag}.png"), dpi=300, bbox_inches="tight")
    plt.close()


# ===========================================================================
# Cross-coverage plots (all coverages in single figure)
# ===========================================================================

def plot_grouped_bars(df: pd.DataFrame, mode: str) -> None:
    """Grouped bar charts: F1/Recall/Precision per caller, grouped by coverage."""
    coverages = sorted(df["Coverage"].unique(), key=lambda x: int(str(x).replace("x", "")))
    callers = [c for c in ["HC", "DV", "ST", "FB", "DS", "DSF"] if c in df["CallerAlias"].unique()]

    for vtype in ["SNP", "INDEL"]:
        subset = df[df["Type"] == vtype]
        if subset.empty:
            continue

        metrics = ["METRIC.F1_Score", "METRIC.Precision", "METRIC.Recall"]
        labels = ["F1 Score", "Precision", "Recall"]
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))

        for ax, metric, label in zip(axes, metrics, labels):
            n_callers = len(callers)
            n_covs = len(coverages)
            bar_width = 0.8 / n_covs
            cov_cmap = plt.cm.viridis([(i + 0.5) / n_covs for i in range(n_covs)])

            for ci, cov in enumerate(coverages):
                cov_data = subset[subset["Coverage"] == cov]
                vals = []
                for caller in callers:
                    row = cov_data[cov_data["CallerAlias"] == caller]
                    vals.append(row[metric].values[0] if len(row) > 0 else 0)

                positions = [i + ci * bar_width for i in range(n_callers)]
                ax.bar(positions, vals, bar_width * 0.9, label=str(cov),
                       color=cov_cmap[ci], edgecolor="black", linewidth=0.5)

            ax.set_xticks([i + (n_covs - 1) * bar_width / 2 for i in range(n_callers)])
            ax.set_xticklabels(callers, fontsize=11)
            ax.set_ylim(0, 1.05)
            ax.set_ylabel(label, fontsize=12)
            ax.set_title(label, fontsize=13, fontweight="bold")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            if ax == axes[-1]:
                ax.legend(title="Coverage", fontsize=9, title_fontsize=10)

        plt.suptitle(f"{vtype} — Metrics by Caller × Coverage ({mode})",
                     fontsize=15, fontweight="bold", y=1.02)
        plt.tight_layout()
        plt.savefig(os.path.join(OUTPUT_DIR, f"grouped_bars_{vtype.lower()}_{mode.lower()}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()


def plot_f1_heatmap(df: pd.DataFrame, mode: str) -> None:
    """Heatmap: callers × coverage → F1 score, separate for SNP and INDEL."""
    callers = [c for c in ["HC", "DV", "ST", "FB", "DS", "DSF"] if c in df["CallerAlias"].unique()]
    coverages = sorted(df["Coverage"].unique(), key=lambda x: int(str(x).replace("x", "")))
    cov_labels = [str(c) for c in coverages]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax, vtype in zip(axes, ["SNP", "INDEL"]):
        subset = df[df["Type"] == vtype]
        matrix = []
        for caller in callers:
            row = []
            for cov in coverages:
                cell = subset[(subset["CallerAlias"] == caller) & (subset["Coverage"] == cov)]
                row.append(cell["METRIC.F1_Score"].values[0] if len(cell) > 0 else 0)
            matrix.append(row)

        data = np.array(matrix)
        im = ax.imshow(data, cmap="RdYlGn", vmin=0.0, vmax=1.0, aspect="auto")

        ax.set_xticks(range(len(cov_labels)))
        ax.set_yticks(range(len(callers)))
        ax.set_xticklabels(cov_labels, fontsize=11)
        ax.set_yticklabels(callers, fontsize=11)
        ax.set_xlabel("Coverage", fontsize=12)

        for i in range(len(callers)):
            for j in range(len(cov_labels)):
                val = data[i, j]
                color = "white" if val < 0.5 else "black"
                ax.text(j, i, f"{val:.3f}", ha="center", va="center",
                        fontsize=10, color=color, fontweight="bold")

        ax.set_title(vtype, fontsize=13, fontweight="bold")

    cbar = fig.colorbar(im, ax=axes, shrink=0.8, pad=0.02)
    cbar.set_label("F1 Score", fontsize=11)
    plt.suptitle(f"F1 Score Heatmap — Callers × Coverage ({mode})",
                 fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"f1_heatmap_{mode.lower()}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()


def plot_snp_vs_indel_comparison(df: pd.DataFrame, mode: str) -> None:
    """Side-by-side SNP vs INDEL F1 per caller, across coverages."""
    callers = [c for c in ["HC", "DV", "ST", "FB", "DS", "DSF"] if c in df["CallerAlias"].unique()]
    coverages = sorted(df["Coverage"].unique(), key=lambda x: int(str(x).replace("x", "")))

    n_covs = len(coverages)
    fig, axes = plt.subplots(1, n_covs, figsize=(5 * n_covs, 5), sharey=True)
    if n_covs == 1:
        axes = [axes]

    for ax, cov in zip(axes, coverages):
        cov_data = df[df["Coverage"] == cov]
        snp_vals = []
        indel_vals = []
        for caller in callers:
            snp = cov_data[(cov_data["CallerAlias"] == caller) & (cov_data["Type"] == "SNP")]
            indel = cov_data[(cov_data["CallerAlias"] == caller) & (cov_data["Type"] == "INDEL")]
            snp_vals.append(snp["METRIC.F1_Score"].values[0] if len(snp) > 0 else 0)
            indel_vals.append(indel["METRIC.F1_Score"].values[0] if len(indel) > 0 else 0)

        x = range(len(callers))
        ax.bar([i - 0.2 for i in x], snp_vals, 0.38, label="SNP",
               color="#4CAF50", edgecolor="black", linewidth=0.5)
        ax.bar([i + 0.2 for i in x], indel_vals, 0.38, label="INDEL",
               color="#FF9800", edgecolor="black", linewidth=0.5)

        ax.set_xticks(list(x))
        ax.set_xticklabels(callers, fontsize=11)
        ax.set_ylim(0, 1.1)
        ax.set_title(f"{cov}", fontsize=12, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if ax == axes[0]:
            ax.set_ylabel("F1 Score", fontsize=12)
            ax.legend(fontsize=10)

    plt.suptitle(f"SNP vs INDEL F1 Comparison ({mode})", fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"snp_vs_indel_{mode.lower()}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()


def plot_variant_type_breakdown(df: pd.DataFrame, mode: str) -> None:
    """Stacked breakdown: TP/FP/FN by variant type per caller across coverages."""
    required = {"TRUTH.TP", "QUERY.FP", "TRUTH.FN"}
    if not required.issubset(df.columns):
        return

    callers = [c for c in ["HC", "DV", "ST", "FB", "DS", "DSF"] if c in df["CallerAlias"].unique()]
    coverages = sorted(df["Coverage"].unique(), key=lambda x: int(str(x).replace("x", "")))

    fig, axes = plt.subplots(len(coverages), 2, figsize=(14, 4 * len(coverages)), sharey="row")
    if len(coverages) == 1:
        axes = axes.reshape(1, -1)

    for row_idx, cov in enumerate(coverages):
        cov_data = df[df["Coverage"] == cov]
        for col_idx, vtype in enumerate(["SNP", "INDEL"]):
            ax = axes[row_idx, col_idx]
            vdata = cov_data[cov_data["Type"] == vtype]

            tp_vals, fp_vals, fn_vals = [], [], []
            for caller in callers:
                cdata = vdata[vdata["CallerAlias"] == caller]
                tp_vals.append(cdata["TRUTH.TP"].values[0] if len(cdata) > 0 else 0)
                fp_vals.append(cdata["QUERY.FP"].values[0] if len(cdata) > 0 else 0)
                fn_vals.append(cdata["TRUTH.FN"].values[0] if len(cdata) > 0 else 0)

            x = range(len(callers))
            ax.bar([i - 0.25 for i in x], tp_vals, 0.25, label="TP",
                   color="#4CAF50", edgecolor="black", linewidth=0.5)
            ax.bar(list(x), fp_vals, 0.25, label="FP",
                   color="#F44336", edgecolor="black", linewidth=0.5)
            ax.bar([i + 0.25 for i in x], fn_vals, 0.25, label="FN",
                   color="#FF9800", edgecolor="black", linewidth=0.5)

            ax.set_xticks(list(x))
            ax.set_xticklabels(callers, fontsize=10)
            ax.set_title(f"{vtype} — {cov}", fontsize=11, fontweight="bold")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            if row_idx == 0 and col_idx == 0:
                ax.legend(fontsize=9)
            if col_idx == 0:
                ax.set_ylabel("Count", fontsize=11)

    plt.suptitle(f"Variant Type Breakdown — TP / FP / FN ({mode})",
                 fontsize=14, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"variant_breakdown_{mode.lower()}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

    # Also output a summary table
    summary_rows = []
    for cov in coverages:
        cov_data = df[df["Coverage"] == cov]
        for vtype in ["SNP", "INDEL"]:
            vdata = cov_data[cov_data["Type"] == vtype]
            for caller in callers:
                cdata = vdata[vdata["CallerAlias"] == caller]
                if len(cdata) > 0:
                    row = cdata.iloc[0]
                    summary_rows.append({
                        "Coverage": cov, "Type": vtype, "Caller": caller,
                        "TP": int(row.get("TRUTH.TP", 0)),
                        "FP": int(row.get("QUERY.FP", 0)),
                        "FN": int(row.get("TRUTH.FN", 0)),
                        "F1": round(row.get("METRIC.F1_Score", 0), 4),
                        "Precision": round(row.get("METRIC.Precision", 0), 4),
                        "Recall": round(row.get("METRIC.Recall", 0), 4),
                    })

    if summary_rows:
        breakdown_df = pd.DataFrame(summary_rows)
        breakdown_path = os.path.join(OUTPUT_DIR, f"variant_breakdown_{mode.lower()}.tsv")
        breakdown_df.to_csv(breakdown_path, sep="\t", index=False)
        print(f"  Variant breakdown table: {breakdown_path}")


# ===========================================================================
# Main execution
# ===========================================================================

print(f"Reading: {STATS_FILE}")
df = pd.read_csv(STATS_FILE, sep="\t")
df["CallerAlias"] = df["Caller"].map(alias_caller)

mask = (df["Filter"] == "PASS") & (df["Subtype"] == "*") & (df["Subset"] == "*")
pass_df = df.loc[mask].copy()

summary_cols = ["Coverage", "Mode", "CallerAlias", "Type", "METRIC.F1_Score", "METRIC.Precision", "METRIC.Recall"]
pass_df[summary_cols].to_csv(os.path.join(OUTPUT_DIR, "summary_stats_python.tsv"), sep="\t", index=False)

combinations = pass_df[["Mode", "Coverage"]].drop_duplicates().sort_values(["Mode", "Coverage"])
print(f"Combinations: {len(combinations)}")

# --- Per-coverage plots (existing) ---
for row in combinations.itertuples(index=False):
    subset = pass_df[(pass_df["Mode"] == row.Mode) & (pass_df["Coverage"] == row.Coverage)].copy()
    tag = sanitize(row.Mode, row.Coverage)
    plot_f1(subset, row.Mode, row.Coverage, tag)
    plot_precision_recall(subset, row.Mode, row.Coverage, tag)
    plot_counts(subset, row.Mode, row.Coverage, tag)

# --- Cross-coverage plots (new) ---
modes = pass_df["Mode"].unique()
for mode in modes:
    mode_df = pass_df[pass_df["Mode"] == mode].copy()
    print(f"\nGenerating cross-coverage plots for mode: {mode}")
    plot_grouped_bars(mode_df, mode)
    plot_f1_heatmap(mode_df, mode)
    plot_snp_vs_indel_comparison(mode_df, mode)
    plot_variant_type_breakdown(mode_df, mode)

print(f"\nAll plots saved to: {OUTPUT_DIR}/")
