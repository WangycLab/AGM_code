# -*- coding: utf-8 -*-
"""
Created on Mon Dec 29 23:44:07 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Plot settings
plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "Arial"})

# Load data
df_scaled = pd.read_csv("species_abundance_profile_scaled.csv", index_col=0)
groups = {"Astronaut 1": ["L_1", "FD_1", "R_1"], "Astronaut 2": ["L_2", "FD_2", "R_2"], "Astronaut 3": ["L_3", "FD_3", "R_3"]}

res = pd.read_csv("FD_vs_L_results.csv", index_col=0)
sig_species = res.index[res["logFC"].abs() > 0.5]
heatmap_df = df_scaled.loc[df_scaled.index.intersection(sig_species)].copy()

# Build column order with spacers
cols_order, group_boundaries, spacer_cols = [], [], []
for i, (grp_name, grp_cols) in enumerate(groups.items()):
    cols_order.extend(grp_cols)
    if i < len(groups) - 1:
        spacer_name = f"_spacer_{i}"
        heatmap_df[spacer_name] = np.nan
        cols_order.append(spacer_name)
        spacer_cols.append(spacer_name)
        group_boundaries.append(len(cols_order) - 1)

heatmap_df = heatmap_df[cols_order]
sns.set(style="white")
mask = heatmap_df.isna()

# Generate heatmap
g = sns.clustermap(heatmap_df.fillna(0), row_cluster=True, col_cluster=False, cmap="vlag", linewidths=1, linecolor="white", figsize=(10, max(6, 0.3 * heatmap_df.shape[0])), dendrogram_ratio=(0.15, 0.0), cbar_kws={"label": "Scaled Abundance"}, mask=mask)
ax = g.ax_heatmap

# Format species labels
for label in ax.get_yticklabels():
    label.set_fontstyle("italic")
    label.set_fontsize(20)

g.cax.set_position([1.2, 0.3, 0.01, 0.5])
g.cax.tick_params(labelsize=20)
g.cax.set_ylabel("Scaled abundance", fontsize=20)
ax.tick_params(axis="y", labelsize=18)

# Add stage labels
group_labels = ["L", "FD", "R"] * len(groups)
tick_positions = []
for i, (grp_name, grp_cols) in enumerate(groups.items()):
    start = 0 if i == 0 else group_boundaries[i - 1] + 1
    tick_positions.extend([start + j + 0.5 for j in range(3)])

ax.set_xticks(tick_positions)
ax.set_xticklabels(group_labels, rotation=0, ha="center", fontsize=23)

# Add astronaut labels
for i, (grp_name, grp_cols) in enumerate(groups.items()):
    start = 0 if i == 0 else group_boundaries[i - 1] + 1
    mid = np.mean([start + j for j in range(len(grp_cols))])
    ax.text(mid + 0.5, -0.5, grp_name, ha="center", va="center", fontsize=23)

ax.set_xlabel("Stages", fontsize=30, labelpad=5)
ax.set_ylabel("Species", fontsize=30, labelpad=5)

plt.savefig("Figure3f.pdf", dpi=300, bbox_inches="tight")
plt.show()