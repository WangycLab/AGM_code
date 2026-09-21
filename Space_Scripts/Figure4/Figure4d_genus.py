# -*- coding: utf-8 -*-
"""
Created on Sun Aug 24 18:32:42 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency
import warnings
warnings.filterwarnings("ignore")

# Plot settings
plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42

# Cluster colors and order
cluster_color_map = {
    "Oxidative stress response|0": "#A6CEE3", "Propionate production|1": "#3B8ABE",
    "Cofactor biosynthesis|2": "#72B29C", "Central carbon metabolism|3": "#84C868",
    "SCFA metabolism|4": "#4F9F3B", "Oxidative stress adaptation|5": "#EC9A91",
    "Cell envelope remodeling|6": "#E93E3F", "Host interface adaptation|7": "#F06C45",
    "Environmental sensing|8": "#FDAC4F", "Anaerobic metabolism|9": "#FB820F",
    "Envelope and motility|10": "#D1AAB7", "Metabolic flexibility|11": "#8C66AF",
    "Propionate utilization|12": "#A99099", "Anaerobic energy metabolism|13": "#EEDB80",
    "Organic acid metabolism|14": "#B15928"
}
cluster_order = list(cluster_color_map.keys())

# Load and prepare metadata
df = pd.read_csv("cell_metadata_annotated.tsv", sep="\t").dropna(subset=["genus", "cluster"])
df["cluster"] = pd.Categorical(df["cluster"].astype(str), categories=cluster_order, ordered=True)
genus_order = df["genus"].value_counts().index.tolist()
df = df[df["genus"].isin(genus_order)]

# Build genus × cluster table
ct_table = pd.crosstab(df["genus"], df["cluster"]).reindex(index=genus_order, columns=cluster_order, fill_value=0)

# Chi-square test and standardized Pearson residuals
chi2, p_chi2, dof, expected = chi2_contingency(ct_table)
residual = (ct_table - expected) / np.sqrt(expected)
print(f"χ² = {chi2:.3f}, df = {dof}, p-value = {p_chi2:.2e}")

# Prepare heatmap
residual_T = residual.T
annot_mat_T = np.where(np.abs(residual.values) > 1.96, "*", "").T

fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
sns.heatmap(
    residual_T, cmap="RdBu_r", vmin=-3, vmax=3, center=0, linewidths=0.3, ax=ax,
    cbar_kws={"label": "Std. Pearson Residual", "shrink": 0.5, "aspect": 15, "pad": 0.03}
)

# Add significance marks
for i, j in zip(*np.where(annot_mat_T == "*")):
    ax.text(j + 0.5, i + 0.85, "*", fontsize=25, ha="center", va="center", color="white", fontweight="bold")

cbar = ax.collections[0].colorbar
cbar.set_label("Std. Pearson Residual", fontsize=18, labelpad=15)
cbar.ax.tick_params(labelsize=15)

ax.set_xlabel("Genus", fontsize=25, labelpad=10)
ax.set_ylabel("Cluster", fontsize=25, labelpad=10)
ax.tick_params(axis="x", rotation=90, labelsize=16)
ax.tick_params(axis="y", rotation=0, labelsize=14)

plt.tight_layout()
plt.savefig("Figure4d.pdf", dpi=300, bbox_inches="tight")
plt.show()
