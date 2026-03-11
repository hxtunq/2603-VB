#!/usr/bin/env python3
"""
Concordance heatmap visualization.

Usage:
    python plot_concordance.py <concordance_matrix.tsv> <output_dir>

Example:
    python visualization/plot_concordance.py \
        results/eval/concordance/concordance_matrix.tsv \
        results/plots
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    import seaborn as sns
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False


CALLER_ORDER = ["HC", "DV", "ST", "FB", "DS"]
CALLER_LABELS = {
    "HC": "GATK HC",
    "DV": "DeepVariant",
    "ST": "Strelka2",
    "FB": "FreeBayes",
    "DS": "DNAscope",
}


def build_matrix(df: pd.DataFrame, callers: list) -> pd.DataFrame:
    """Build a square concordance matrix from long-form data."""
    matrix = pd.DataFrame(1.0, index=callers, columns=callers)
    for _, row in df.iterrows():
        a, b = row["CallerA"], row["CallerB"]
        if a in callers and b in callers:
            matrix.loc[a, b] = row["Concordance"]
            matrix.loc[b, a] = row["Concordance"]
    return matrix


def plot_heatmap_seaborn(matrix: pd.DataFrame, title: str, output_path: str):
    """Plot heatmap using seaborn."""
    labels = [CALLER_LABELS.get(c, c) for c in matrix.index]
    fig, ax = plt.subplots(figsize=(8, 6.5))

    annot = matrix.values.copy()
    annot_text = np.vectorize(lambda x: f"{x:.3f}")(annot)

    sns.heatmap(
        matrix.values,
        annot=annot_text,
        fmt="",
        xticklabels=labels,
        yticklabels=labels,
        cmap="YlOrRd",
        vmin=0.0,
        vmax=1.0,
        linewidths=0.5,
        linecolor="white",
        square=True,
        cbar_kws={"label": "Concordance Rate", "shrink": 0.8},
        ax=ax,
    )
    ax.set_title(title, fontsize=14, fontweight="bold", pad=15)
    ax.tick_params(axis="both", labelsize=11)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_heatmap_matplotlib(matrix: pd.DataFrame, title: str, output_path: str):
    """Fallback heatmap using matplotlib only."""
    labels = [CALLER_LABELS.get(c, c) for c in matrix.index]
    data = matrix.values
    n = len(labels)

    fig, ax = plt.subplots(figsize=(8, 6.5))
    im = ax.imshow(data, cmap="YlOrRd", vmin=0.0, vmax=1.0, aspect="equal")

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(labels, fontsize=11, rotation=45, ha="right")
    ax.set_yticklabels(labels, fontsize=11)

    for i in range(n):
        for j in range(n):
            color = "white" if data[i, j] > 0.7 else "black"
            ax.text(j, i, f"{data[i, j]:.3f}", ha="center", va="center",
                    fontsize=10, color=color, fontweight="bold")

    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("Concordance Rate", fontsize=11)

    ax.set_title(title, fontsize=14, fontweight="bold", pad=15)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_heatmap(matrix, title, output_path):
    if HAS_SEABORN:
        plot_heatmap_seaborn(matrix, title, output_path)
    else:
        plot_heatmap_matplotlib(matrix, title, output_path)


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <concordance_matrix.tsv> <output_dir>")
        sys.exit(1)

    tsv_path = sys.argv[1]
    output_dir = sys.argv[2]
    os.makedirs(output_dir, exist_ok=True)

    print(f"Reading: {tsv_path}")
    df = pd.read_csv(tsv_path, sep="\t")

    # Determine available callers
    available = [c for c in CALLER_ORDER if c in df["CallerA"].unique() or c in df["CallerB"].unique()]
    if not available:
        print("ERROR: No recognized callers found in the concordance matrix.")
        sys.exit(1)

    coverages = sorted(df["Coverage"].unique(), key=lambda x: int(str(x).replace("x", "")))
    print(f"Callers: {available}")
    print(f"Coverages: {coverages}")

    # Per-coverage heatmaps
    for cov in coverages:
        cov_df = df[df["Coverage"] == cov]
        matrix = build_matrix(cov_df, available)
        title = f"Caller Concordance — {cov}"
        output_path = os.path.join(output_dir, f"concordance_heatmap_{cov}.png")
        plot_heatmap(matrix, title, output_path)
        print(f"  Saved: {output_path}")

    # Average across coverages
    if len(coverages) > 1:
        avg_df = df.groupby(["CallerA", "CallerB"])["Concordance"].mean().reset_index()
        matrix = build_matrix(avg_df, available)
        title = "Caller Concordance — Average Across Coverages"
        output_path = os.path.join(output_dir, "concordance_heatmap_average.png")
        plot_heatmap(matrix, title, output_path)
        print(f"  Saved: {output_path}")

    print(f"\nAll concordance heatmaps saved to: {output_dir}/")


if __name__ == "__main__":
    main()
