# -*- coding: utf-8 -*-
"""
Created on Mon Sep 15 15:15:18 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as mticker
import matplotlib as mpl

# Configure Matplotlib to produce editable text
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

# Load the cell metadata table for downstream plotting
df_meta = pd.read_csv("cell_metadata.tsv", sep="\t")


# Ensure gene_number column is numeric for analysis
df_meta["gene_number"] = pd.to_numeric(
    df_meta["gene_number"],
    errors="coerce"
)

# Define the order of samples for consistent plotting
sample_order = ["Con", "2mo", "4mo", "6mo"]

df_meta["sample"] = pd.Categorical(
    df_meta["sample"],
    categories=sample_order,
    ordered=True
)

# Select top cells per sample based on gene count
TOP_N = 5000

df_meta = (
    df_meta
    .sort_values("gene_number", ascending=False)
    .groupby("sample", group_keys=False)
    .head(TOP_N)
)

# Print counts of selected cells per sample to verify selection
print(df_meta["sample"].value_counts())

# Apply optional scaling factor to gene counts
scale_factor = {
    "Con": 1,
    "2mo": 1,
    "4mo": 1,
    "6mo": 1
}

df_meta["scale_factor"] = (
    df_meta["sample"]
    .map(scale_factor)
    .astype(float)
)

# Compute scaled gene counts for plotting
df_meta["gene_number_scaled"] = (
    df_meta["gene_number"] / df_meta["scale_factor"]
)

# Set global style for publication-quality figures
sns.set_theme(
    style="white",
    context="paper",
    font="Arial"
)

# Create figure and axis for the violin plot
fig, ax = plt.subplots(figsize=(5, 3))

# Draw violin plot to show distribution of gene counts
sns.violinplot(
    x="sample",
    y="gene_number_scaled",
    data=df_meta,
    order=sample_order,
    inner=None,
    cut=0,
    bw=0.2,
    linewidth=0.2,
    color="#DCE6F1",
    saturation=1,
    ax=ax
)

# Overlay boxplot to highlight summary statistics
sns.boxplot(
    x="sample",
    y="gene_number_scaled",
    data=df_meta,
    order=sample_order,
    width=0.15,
    showfliers=False,
    boxprops=dict(
        facecolor="white",
        edgecolor="black",
        linewidth=0.6
    ),
    whiskerprops=dict(
        color="black",
        linewidth=0.6
    ),
    capprops=dict(
        color="black",
        linewidth=0.6
    ),
    medianprops=dict(
        color="black",
        linewidth=1.2
    ),
    ax=ax
)

# Set axis labels with readable font size
ax.set_xlabel("Sample", fontsize=16)
ax.set_ylabel("Gene counts per cell", fontsize=16)

# Customize tick appearance for clarity
ax.tick_params(
    axis="x",
    labelsize=15,
    length=0
)
ax.tick_params(
    axis="y",
    labelsize=15,
    width=1,
    length=4
)

# Adjust y-axis range to avoid extreme outliers
y_max = df_meta["gene_number_scaled"].quantile(0.99)
ax.set_ylim(0, y_max * 1.03)
ax.yaxis.set_major_locator(mticker.MaxNLocator(5))

# Remove unnecessary plot spines for a cleaner look
sns.despine(
    ax=ax,
    top=True,
    right=True,
    left=False,
    bottom=True
)

# Save figure as a vector PDF
plt.tight_layout()
plt.savefig(
    "Figure2g.pdf",
    dpi=300,
    bbox_inches="tight",
    format="pdf"
)
plt.show()