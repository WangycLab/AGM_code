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
import matplotlib as mpl
from sklearn.cluster import KMeans
from matplotlib.patches import Patch

# Plot settings
mpl.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "Arial", "font.sans-serif": ["Arial"]})
os.chdir(os.getcwd())

# Load expression matrix
df = pd.read_csv("species_pseudotime_deg_expression.csv")
df.rename(columns={df.columns[0]: "barcode"}, inplace=True)
df = df.set_index("pseudotime")
if "barcode" in df.columns: df = df.drop(columns="barcode")

# Filter genes
expr_cells = (df > 0).sum(axis=0); expr_sum = df.sum(axis=0)
valid_genes = df.columns[(expr_cells >= 100) & (expr_sum >= 500)]
df = df[valid_genes]

# Convert pseudotime to pseudotimepoint
pseudotime_bins = np.arange(0, 57, 1); pseudotime_labels = pseudotime_bins[:-1] + 1
df["pseudotimepoint"] = pd.cut(df.index, bins=pseudotime_bins, labels=pseudotime_labels, right=False)
new_df = df.groupby("pseudotimepoint", observed=False).mean()

# Log transform and Z-score
log_df = np.log1p(new_df)
scaled_df = log_df.sub(log_df.mean(axis=0), axis=1).div(log_df.std(axis=0), axis=1).fillna(0)
gene_matrix = scaled_df.T.values

# Determine optimal k
K = list(range(2, 16)); inertias = []
for k in K: inertias.append(KMeans(n_clusters=k, random_state=0, n_init=20).fit(gene_matrix).inertia_)

points = np.column_stack((K, inertias)); first, last = points[0], points[-1]
line_vec = last - first; line_vec /= np.linalg.norm(line_vec)
vec_from_first = points - first
projection = np.outer(np.dot(vec_from_first, line_vec), line_vec)
distances = np.linalg.norm(vec_from_first - projection, axis=1)
elbow_index = np.argmax(distances); optimal_k = K[elbow_index]

# Plot elbow curve
plt.figure(figsize=(4, 4))
plt.plot(K, inertias, "-o", color="black", linewidth=2, markersize=6)
plt.scatter(optimal_k, inertias[elbow_index], color="red", s=90, zorder=10)
plt.axvline(optimal_k, color="red", linestyle="--", linewidth=1.8)
plt.axhline(inertias[elbow_index], color="red", linestyle="--", linewidth=1.2, alpha=0.6)
plt.annotate(f"Optimal k = {optimal_k}", xy=(optimal_k, inertias[elbow_index]), xytext=(optimal_k + 0.8, inertias[elbow_index] * 1.05), fontsize=12, color="red", arrowprops=dict(arrowstyle="->", color="red", lw=1.5))
plt.xlabel("Number of clusters (k)", fontsize=14); plt.ylabel("Within-cluster SSE", fontsize=14)
plt.xticks(K); plt.tight_layout()
plt.show()

# K-means clustering
n_clusters = optimal_k
kmeans = KMeans(n_clusters=n_clusters, random_state=0, n_init=20)
gene_clusters = kmeans.fit_predict(gene_matrix)

# Save cluster assignments
cluster_df = pd.DataFrame({"gene": scaled_df.columns, "cluster": gene_clusters})
cluster_df.to_csv("gene_clusters.csv", index=False)

# Sort genes by cluster
cluster_df = cluster_df.sort_values("cluster")
scaled_df = scaled_df[cluster_df["gene"]]
gene_clusters = cluster_df["cluster"].tolist()

# Cluster colors
palette = sns.color_palette("Set2", n_clusters)
group_color_map = {g: palette[i] for i, g in enumerate(sorted(set(gene_clusters)))}
row_colors = [group_color_map[g] for g in gene_clusters]

# Draw heatmap
g = sns.clustermap(scaled_df.T, row_cluster=False, col_cluster=False, row_colors=row_colors, cmap="rocket_r", xticklabels=False, yticklabels=False, figsize=(4, 18), cbar_kws={"label": "Scaled Expression Level"}, dendrogram_ratio=(0.001, 0.02), cbar_pos=(1.01, 0.70, 0.015, 0.20))
g.ax_row_dendrogram.set_visible(False)

# Move cluster color bar
pos = g.ax_row_colors.get_position()
g.ax_row_colors.set_position([pos.x0 - 0.01, pos.y0, pos.width, pos.height])

# Configure x-axis
ax = g.ax_heatmap
xtick_positions = np.linspace(0, scaled_df.shape[0] - 1, 3)
xtick_labels = [str(i) for i in range(0, 57, 20)]
ax.set_xticks(xtick_positions); ax.set_xticklabels(xtick_labels, fontsize=25)
ax.set_xlabel("Pseudo-timepoint", fontsize=30, labelpad=15)
ax.set_ylabel("Pseudotime-associated DEGs", fontsize=30, labelpad=20)

# Configure colorbar
g.cax.tick_params(labelsize=25); g.cax.yaxis.label.set_size(30)

# Add cluster legend
legend_handles = [Patch(color=group_color_map[g], label=f"Cluster {g}") for g in sorted(group_color_map)]
g.ax_heatmap.legend(handles=legend_handles, bbox_to_anchor=(1.02, 0.03), loc="lower left", borderaxespad=0, fontsize=20, title="Clusters", title_fontsize=25)

# Save heatmap
plt.savefig("FigureS5d.pdf", format="pdf", dpi=300, bbox_inches="tight")
plt.show()