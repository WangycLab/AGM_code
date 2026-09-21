# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 18:53:18 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency
import warnings
warnings.filterwarnings("ignore")

plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})

# Cluster colors and order
cluster_color_map = {
    "Cofactor metabolism|0": "#A6CEE3",
    "Redox homeostasis|1": "#98D277",
    "Uncharacterized function|2": "#F16667",
    "Nutrient transport|3": "#FE982C",
    "Mobile genetic elements|4": "#7D54A5",
    "Protein stress response|5": "#B15928"
}
cluster_order = list(cluster_color_map)

# Load and filter metadata
df = pd.read_csv("cell_metadata_annotated.tsv", sep="\t").dropna(subset=["timepoint", "cluster"])
df["cluster"] = pd.Categorical(df["cluster"].astype(str), categories=cluster_order, ordered=True)
astro_order = ["L-60", "L-30", "FD30", "FD90", "FD150", "R+1", "R+7", "R+14"]
df = df[df["timepoint"].isin(astro_order)]

# Build contingency table
ct_table = pd.crosstab(df["timepoint"], df["cluster"]).reindex(index=astro_order, columns=cluster_order, fill_value=0)

# Chi-square test and standardized Pearson residuals
chi2, p_chi2, dof, expected = chi2_contingency(ct_table)
residual = (ct_table - expected) / np.sqrt(expected)
print(f"Chi-square = {chi2:.3f}, df = {dof}, p = {p_chi2:.2e}")

# Plot heatmap
residual_T = residual.T
annot_mat_T = np.where(np.abs(residual_T.values) > 1.96, "*", "")

fig, ax = plt.subplots(figsize=(7, 3.7), dpi=300)
sns.heatmap(
    residual_T, ax=ax, cmap="RdBu_r", vmin=-3, vmax=3, center=0,
    linewidths=0.3, cbar_kws={"label": "Std. Pearson Residual", "shrink": 0.5, "aspect": 15, "pad": 0.05}
)

# Add significance markers
for i, j in zip(*np.where(annot_mat_T == "*")):
    ax.text(j + 0.5, i + 0.8, "*", fontsize=30, ha="center", va="center", color="white", fontweight="bold")

cbar = ax.collections[0].colorbar
cbar.set_label("Std. Pearson Residual", fontsize=18)
cbar.ax.tick_params(labelsize=15)

ax.set_xlabel("Timepoint", fontsize=20)
ax.set_ylabel("Cluster", fontsize=20)

ax.set_xticks(np.arange(len(astro_order)) + 0.5)
ax.set_xticklabels(astro_order, fontsize=20, rotation=45)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=20, rotation=0)

ax.tick_params(axis="x", labelsize=15)
ax.tick_params(axis="y", labelsize=15)

plt.tight_layout()
plt.savefig("Figure6b.pdf", bbox_inches="tight")
plt.show()