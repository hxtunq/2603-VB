#!/usr/bin/env python3
"""
Quick summary plots for variant calling benchmark.
Reads the all_stats.tsv from gather_stats.sh and generates summary bar charts.

Usage:
    python plot_summary.py [stats_file] [output_dir]
    python plot_summary.py results/eval/all_stats.tsv results/plots
"""

import sys
import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# --- Config ---
STATS_FILE = sys.argv[1] if len(sys.argv) > 1 else 'results/eval/all_stats.tsv'
OUTPUT_DIR = sys.argv[2] if len(sys.argv) > 2 else 'results/plots'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Caller colors
COLORS = {
    'DV': '#2196F3', 'deepvariant': '#2196F3',
    'FB': '#FF9800', 'freebayes': '#FF9800',
    'HC': '#4CAF50', 'gatk': '#4CAF50',
    'ST': '#E91E63', 'strelka2': '#E91E63',
}
ALIASES = {'deepvariant': 'DV', 'freebayes': 'FB', 'gatk': 'HC', 'strelka2': 'ST'}

print(f"Reading: {STATS_FILE}")
df = pd.read_csv(STATS_FILE, sep='\t')

# Normalize caller names
if 'Caller' in df.columns:
    df['Caller'] = df['Caller'].map(ALIASES).fillna(df['Caller'])
elif 'CallerFilter' in df.columns:
    df['Caller'] = df['CallerFilter'].map(ALIASES).fillna(df['CallerFilter'])

# Filter PASS + overall
mask = (df['Filter'] == 'PASS') & (df['Subtype'] == '*') & (df['Subset'] == '*')
pass_df = df[mask].copy()

print(f"Callers: {pass_df['Caller'].unique()}")
print(f"Types: {pass_df['Type'].unique()}")


# ======================================================================
# Plot 1: F1 Score bar chart
# ======================================================================
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for ax, vtype in zip(axes, ['SNP', 'INDEL']):
    subset = pass_df[pass_df['Type'] == vtype].sort_values('Caller')
    colors = [COLORS.get(c, '#999') for c in subset['Caller']]
    ax.bar(subset['Caller'], subset['METRIC.F1_Score'], color=colors, edgecolor='black', width=0.5)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel('F1 Score')
    ax.set_title(f'{vtype}')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # Add value labels
    for i, (c, v) in enumerate(zip(subset['Caller'], subset['METRIC.F1_Score'])):
        ax.text(i, v + 0.02, f'{v:.4f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

plt.suptitle('F1 Score by Variant Caller', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'summary_f1_scores.png'), dpi=300, bbox_inches='tight')
plt.close()


# ======================================================================
# Plot 2: Precision & Recall grouped bar chart
# ======================================================================
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for ax, vtype in zip(axes, ['SNP', 'INDEL']):
    subset = pass_df[pass_df['Type'] == vtype].sort_values('Caller')
    callers = list(subset['Caller'])
    x = range(len(callers))
    w = 0.35

    prec = list(subset['METRIC.Precision'])
    rec = list(subset['METRIC.Recall'])

    bars1 = ax.bar([i - w/2 for i in x], prec, w, label='Precision', color='#64B5F6', edgecolor='black')
    bars2 = ax.bar([i + w/2 for i in x], rec, w, label='Recall', color='#FFB74D', edgecolor='black')

    ax.set_xticks(list(x))
    ax.set_xticklabels(callers)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel('Value')
    ax.set_title(f'{vtype}')
    ax.legend()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

plt.suptitle('Precision & Recall by Variant Caller', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'summary_precision_recall.png'), dpi=300, bbox_inches='tight')
plt.close()


# ======================================================================
# Plot 3: TP / FP / FN counts
# ======================================================================
if 'TRUTH.TP' in pass_df.columns:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for ax, vtype in zip(axes, ['SNP', 'INDEL']):
        subset = pass_df[pass_df['Type'] == vtype].sort_values('Caller')
        callers = list(subset['Caller'])
        x = range(len(callers))
        w = 0.25

        tp = list(subset['TRUTH.TP'])
        fp = list(subset['QUERY.FP'])
        fn = list(subset['TRUTH.FN'])

        ax.bar([i - w for i in x], tp, w, label='TP', color='#4CAF50', edgecolor='black')
        ax.bar([i for i in x], fp, w, label='FP', color='#F44336', edgecolor='black')
        ax.bar([i + w for i in x], fn, w, label='FN', color='#FF9800', edgecolor='black')

        ax.set_xticks(list(x))
        ax.set_xticklabels(callers)
        ax.set_ylabel('Count')
        ax.set_title(f'{vtype}')
        ax.legend()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.suptitle('TP / FP / FN Counts', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'summary_tp_fp_fn.png'), dpi=300, bbox_inches='tight')
    plt.close()


# ======================================================================
# Summary text
# ======================================================================
print("\n" + "=" * 60)
print("BENCHMARK SUMMARY")
print("=" * 60)
for vtype in ['SNP', 'INDEL']:
    print(f"\n--- {vtype} ---")
    subset = pass_df[pass_df['Type'] == vtype].sort_values('METRIC.F1_Score', ascending=False)
    for _, row in subset.iterrows():
        print(f"  {row['Caller']:4s}  F1={row['METRIC.F1_Score']:.4f}  "
              f"Prec={row['METRIC.Precision']:.4f}  Rec={row['METRIC.Recall']:.4f}")

print(f"\nPlots saved to: {OUTPUT_DIR}/")
