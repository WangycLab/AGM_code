# -*- coding: utf-8 -*-
"""
Created on Wed Jan 7 19:01:32 2026
@author: ZHENG XINGHAI
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist
from matplotlib.patches import Patch
import matplotlib as mpl

# Configure matplotlib to embed editable fonts
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.sans-serif'] = ['Arial']

# Use the current working directory
os.chdir(os.getcwd())

# Load expression matrix and prepare pseudotime index
df = pd.read_csv("species_pseudotime_deg_expression.csv")
df.rename(columns={df.columns[0]: "barcode"}, inplace=True)
df = df.set_index("pseudotime")

if "barcode" in df.columns:
    df = df.drop(columns="barcode")

# Filter genes based on expression prevalence and total abundance
expr_cells = (df > 0).sum(axis=0)
expr_sum = df.sum(axis=0)

valid_genes = df.columns[
    (expr_cells >= 100) & (expr_sum >= 500)
]

df = df[valid_genes]
print(f"Genes kept after filtering: {df.shape[1]}")

# Convert continuous pseudotime into discrete pseudotimepoints
pseudotime_bins = np.arange(0, 120, 1)
pseudotime_labels = pseudotime_bins[:-1] + 1

df_indexed = df.copy()
df_indexed["pseudotimepoint"] = pd.cut(
    df_indexed.index,
    bins=pseudotime_bins,
    labels=pseudotime_labels,
    right=False
)

new_df = df_indexed.groupby("pseudotimepoint", observed=False).mean()

# Apply log transformation and z-score normalization
log_df = np.log1p(new_df)
scaled_df = log_df.sub(log_df.mean(axis=0), axis=1).div(log_df.std(axis=0), axis=1)
scaled_df = scaled_df.fillna(0)

# Perform hierarchical clustering of genes using correlation distance
gene_matrix = scaled_df.T.values
distance_matrix = pdist(gene_matrix, metric="correlation")
linkage_matrix = linkage(distance_matrix, method="average")

n_clusters = 10
gene_groups = fcluster(linkage_matrix, t=n_clusters, criterion="maxclust")

# Reindex cluster labels to start from zero
unique_groups = np.unique(gene_groups)
group_mapping = {old: new for new, old in enumerate(unique_groups)}
gene_clusters = [group_mapping[g] for g in gene_groups]

# Save gene cluster assignments
cluster_df = pd.DataFrame({
    "gene": scaled_df.columns,
    "cluster": gene_clusters
})

cluster_df.to_csv("gene_clusters.csv", index=False)
print("Saved gene_clusters.csv")

# Generate cluster color annotations
n_clusters = len(set(gene_clusters))
palette = sns.color_palette("Set2", n_clusters)
group_color_map = {grp: palette[i] for i, grp in enumerate(sorted(set(gene_clusters)))}
row_colors = [group_color_map[g] for g in gene_clusters]

# Draw clustered heatmap with fixed pseudotime ordering
g = sns.clustermap(
    scaled_df.T,
    row_linkage=linkage_matrix,
    col_cluster=False,
    row_colors=row_colors,
    cmap="rocket_r",
    xticklabels=False,
    yticklabels=False,
    figsize=(15, 20),
    cbar_kws={"label": "Scaled Expression Level"},
    dendrogram_ratio=(0.12, 0.02),
    cbar_pos=(1.01, 0.7, 0.01, 0.2)
)

# Redraw dendrogram to control line properties
dendrogram(
    linkage_matrix,
    ax=g.ax_row_dendrogram,
    orientation="left",
    no_labels=True
)

g.ax_row_dendrogram.invert_yaxis()

for line in g.ax_row_dendrogram.lines:
    line.set_linewidth(0.1)

# Adjust the width and position of cluster annotation bars
pos = g.ax_row_colors.get_position()
g.ax_row_colors.set_position([
    pos.x0 - 0.003,
    pos.y0,
    pos.width * 0.8,
    pos.height
])

for line in g.ax_row_dendrogram.lines:
    line.set_linewidth(0.1)

for coll in g.ax_row_dendrogram.collections:
    coll.set_linewidth(0.1)

# Configure pseudotime axis labels
ax = g.ax_heatmap
xtick_positions = np.linspace(0, scaled_df.shape[0]-1, 13)
xtick_labels = [str(i) for i in range(0, 121, 10)]

ax.set_xticks(xtick_positions)
ax.set_xticklabels(xtick_labels, fontsize=25)

ax.set_xlabel("Pseudo-timepoint", fontsize=30, labelpad=15)
ax.set_ylabel("Pseudotime-associated DEGs", fontsize=30, labelpad=20)

# Adjust colorbar label and tick formatting
g.cax.tick_params(labelsize=25)
g.cax.yaxis.label.set_size(30)

# Add cluster legend to the heatmap
legend_handles = [
    Patch(color=group_color_map[g], label=f"Cluster {g}")
    for g in sorted(group_color_map.keys())
]

g.ax_heatmap.legend(
    handles=legend_handles,
    bbox_to_anchor=(1.02, 0.03),
    loc="lower left",
    borderaxespad=0,
    fontsize=20,
    title="Clusters",
    title_fontsize=25
)

# Save the final figure as a vector PDF
plt.savefig("Figure5e.pdf", format="pdf", dpi=300, bbox_inches="tight")
plt.show()