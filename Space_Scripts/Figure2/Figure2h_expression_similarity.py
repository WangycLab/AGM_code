# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 14:19:46 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
import numpy as np
import itertools
import seaborn as sns
import matplotlib.pyplot as plt

# Ensure PDFs are AI-editable by setting font type
plt.rcParams["pdf.fonttype"] = 42

# Load the gene-level logCPM expression matrix
df = pd.read_csv("gene_logCPM_matrix.csv", index_col=0)

# Harmonize sample names and arrange columns in desired order
sample_map = {
    "Con": "Con",
    "2mo": "2mo",
    "6mo": "6mo",
    "4mo": "4mo"
}
df = df.rename(columns=sample_map)

desired_order = ["Con", "2mo", "4mo", "6mo"]
df = df[desired_order]
samples = df.columns.tolist()

# Initialize an empty pairwise correlation matrix
corr_matrix = pd.DataFrame(index=samples, columns=samples, dtype=float)

# Define a function to select top 10 percent expressed genes in each sample
def top10_genes(series):
    threshold = np.percentile(series, 90)
    return series[series >= threshold].index

# Compute correlations based on shared highly expressed genes between sample pairs
for s1, s2 in itertools.combinations(samples, 2):
    top_s1 = top10_genes(df[s1])
    top_s2 = top10_genes(df[s2])
    inter = top_s1.intersection(top_s2)

    if len(inter) < 3:
        corr = np.nan
    else:
        corr = df.loc[inter, s1].corr(df.loc[inter, s2])

    corr_matrix.loc[s1, s2] = corr
    corr_matrix.loc[s2, s1] = corr

# Fill the diagonal with 1 to represent self-correlation
np.fill_diagonal(corr_matrix.values, 1)

# Plot a heatmap showing sample similarity based on gene expression
plt.figure(figsize=(5, 4))
ax = sns.heatmap(
    corr_matrix.astype(float),
    cmap="mako",
    vmin=0,
    vmax=1,
    annot=True,
    fmt=".2f",
    square=True,
    linewidths=0.6,
    cbar_kws={"label": "Pearson correlation"},
    annot_kws={"size": 15}
)

# Customize colorbar font size
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=15)
cbar.ax.yaxis.label.set_size(16)

# Set axis labels and tick formatting
plt.xticks(rotation=0, ha="center", fontsize=15)
plt.yticks(rotation=0, ha="right", fontsize=15)

# Add title to the heatmap
plt.title("Gene expression-based similarity", fontsize=16, pad=15)

# Save the figure as PDF with AI-editable text
plt.tight_layout()
plt.savefig("Figure2h.pdf", dpi=300)
plt.show()