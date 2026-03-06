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
    if name.startswith("dnascope"):
        return "DS"
    return name


def sanitize(mode: str, coverage: str) -> str:
    return "".join(ch if ch.isalnum() or ch == "_" else "_" for ch in f"{mode.lower()}_{coverage}")


def ordered_subset(df: pd.DataFrame) -> pd.DataFrame:
    order = {"HC": 0, "DV": 1, "ST": 2, "FB": 3, "DS": 4}
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


print(f"Reading: {STATS_FILE}")
df = pd.read_csv(STATS_FILE, sep="\t")
df["CallerAlias"] = df["Caller"].map(alias_caller)

mask = (df["Filter"] == "PASS") & (df["Subtype"] == "*") & (df["Subset"] == "*")
pass_df = df.loc[mask].copy()

summary_cols = ["Coverage", "Mode", "CallerAlias", "Type", "METRIC.F1_Score", "METRIC.Precision", "METRIC.Recall"]
pass_df[summary_cols].to_csv(os.path.join(OUTPUT_DIR, "summary_stats_python.tsv"), sep="\t", index=False)

combinations = pass_df[["Mode", "Coverage"]].drop_duplicates().sort_values(["Mode", "Coverage"])
print(f"Combinations: {len(combinations)}")

for row in combinations.itertuples(index=False):
    subset = pass_df[(pass_df["Mode"] == row.Mode) & (pass_df["Coverage"] == row.Coverage)].copy()
    tag = sanitize(row.Mode, row.Coverage)
    plot_f1(subset, row.Mode, row.Coverage, tag)
    plot_precision_recall(subset, row.Mode, row.Coverage, tag)
    plot_counts(subset, row.Mode, row.Coverage, tag)

print(f"Plots saved to: {OUTPUT_DIR}/")
