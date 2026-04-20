# -*- coding: utf-8 -*-
"""
Created on Sun Aug 24 18:32:42 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist

# Configure matplotlib to export editable vector text
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams["font.family"] = "Arial"

# Load input datasets for bubble plot and gene annotation
fraction_df = pd.read_csv("gene_fraction_by_cluster.csv", index_col=0)
expr_df = pd.read_csv("gene_mean_expr_by_cluster.csv", index_col=0)
annot_df = pd.read_csv("all_clusters_DEGs_annotated.csv")

# Filter valid COG annotations and generate merge key
annot_df = annot_df[~annot_df['COG_ID'].isin(["-"])]
annot_df = annot_df[~annot_df['COG_ID'].str.contains(",", na=False)]
annot_df['gene_key'] = annot_df['gene'].str.replace("_", "-")

# Merge gene fraction and expression matrices with annotations
fraction_df = fraction_df.merge(annot_df[['gene_key', 'COG_ID']], left_index=True, right_on='gene_key', how='inner')
expr_df = expr_df.merge(annot_df[['gene_key', 'COG_ID']], left_index=True, right_on='gene_key', how='inner')

# Convert matrices to long format for visualization
fraction_melt = fraction_df.melt(id_vars=['gene_key', 'COG_ID'], var_name='cluster', value_name='fraction')
expr_melt = expr_df.melt(id_vars=['gene_key', 'COG_ID'], var_name='cluster', value_name='expression')
merged = pd.merge(fraction_melt, expr_melt, on=['gene_key', 'COG_ID', 'cluster'])

# Identify top marker genes per cluster based on fraction values
top_genes = []
for clust in merged['cluster'].unique():
    sub = merged[merged['cluster'] == clust]
    top = sub.groupby('gene_key')['fraction'].max().sort_values(ascending=False).head(5).index
    top_genes.extend(top)

# Prepare plotting dataframe and clip values to valid range
plot_df = merged[merged['gene_key'].isin(top_genes)].copy()
plot_df['gene_key_mod'] = plot_df['gene_key']
plot_df['fraction'] = plot_df['fraction'].clip(upper=1)
plot_df['expression'] = plot_df['expression'].clip(upper=1)

# Perform hierarchical clustering of genes based on fraction profiles
pivot_fraction = plot_df.pivot_table(index='gene_key_mod', columns='cluster', values='fraction', fill_value=0)
dist_mat = pdist(pivot_fraction.values, metric='euclidean')
Z = linkage(dist_mat, method='ward')
dendro = dendrogram(Z, no_plot=True)
ordered_genes = pivot_fraction.index[dendro['leaves']]

# Apply hierarchical ordering to gene axis
plot_df['gene_key_mod'] = pd.Categorical(plot_df['gene_key_mod'], categories=ordered_genes, ordered=True)
plot_df = plot_df.sort_values(['cluster', 'gene_key_mod'])

# Load cell metadata for species composition analysis
df_meta = pd.read_csv("cell_metadata.tsv", sep="\t")
df_meta['cluster'] = df_meta['cluster'].apply(lambda x: f"Cluster_{x}")

# Compute cluster by species proportion table
count_table = pd.crosstab(df_meta['cluster'], df_meta['species'])
prop_table = count_table.div(count_table.sum(axis=1), axis=0)

# Retain dominant species per cluster
def keep_top_n(row, n=3):
    top_species = row.sort_values(ascending=False).head(n).index
    new_row = pd.Series(0, index=row.index)
    new_row[top_species] = row[top_species]
    return new_row

prop_table = prop_table.apply(keep_top_n, axis=1)
prop_table = prop_table.loc[:, (prop_table != 0).any(axis=0)]

# Define cluster ordering for consistent visualization
cluster_order = [f"Cluster_{i}" for i in range(23)]
prop_table = prop_table.reindex(cluster_order)

# Generate color palette for stacked bar visualization
num_colors = len(prop_table.columns)
palette_bar = sns.color_palette("cubehelix", num_colors)

# Apply predefined cluster ordering
sorted_clusters = [np.int64(10), np.int64(9), np.int64(7), np.int64(11), np.int64(19), np.int64(20), np.int64(5), np.int64(13), np.int64(6), np.int64(1), np.int64(17), np.int64(12), np.int64(15), np.int64(21), np.int64(0), np.int64(2), np.int64(4), np.int64(14), np.int64(22), np.int64(8), np.int64(3), np.int64(16)]
sorted_cluster_names = [f"Cluster_{i}" for i in sorted_clusters]
prop_table = prop_table.reindex(sorted_cluster_names)

# Ensure cluster column uses ordered categorical type
plot_df['cluster'] = pd.Categorical(plot_df['cluster'], categories=sorted_cluster_names, ordered=True)
plot_df = plot_df.sort_values('cluster')

# Normalize fraction and expression values within each cluster
cols_to_scale = ['fraction', 'expression']
plot_df[cols_to_scale] = plot_df.groupby('cluster')[cols_to_scale].transform(lambda x: (x - x.min()) / (x.max() - x.min()))

# Create figure layout for bubble plot and stacked bar plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(25, 16), sharey=True, gridspec_kw={'width_ratios':[3, 1]})

# Generate bubble plot showing marker gene patterns
scatter = sns.scatterplot(
    data=plot_df,
    x='gene_key_mod',
    y='cluster',
    size='fraction',
    hue='expression',
    sizes=(10, 300),
    palette='crest',
    edgecolor='black',
    ax=ax1
)

# Customize bubble plot axis appearance
ax1.set_xlabel("Marker genes", fontsize=30)
ax1.set_ylabel("Clusters", fontsize=30)
ax1.tick_params(axis='x', rotation=90, labelsize=18)
ax1.tick_params(axis='y', labelsize=18)

# Apply italic style to gene labels
for label in ax1.get_xticklabels():
    label.set_fontstyle('italic')

# Modify legend labels and formatting
handles, labels = ax1.get_legend_handles_labels()
new_labels = []
new_handles = []

for h, l in zip(handles, labels):
    if l == 'fraction':
        new_labels.append("S. Exp. R")
        new_handles.append(h)
    elif l == 'expression':
        new_labels.append("S. Exp. L")
        new_handles.append(h)
    elif l not in ['size', 'hue']:
        new_labels.append(l)
        new_handles.append(h)

leg = ax1.legend(new_handles, new_labels, title="", ncol=2, fontsize=18,
                 title_fontsize=25, bbox_to_anchor=(0.05, -0.4), loc='upper left')

# Adjust legend text sizes
for text in leg.get_texts():
    if text.get_text() in ["S. Exp. L", "S. Exp. R"]:
        text.set_fontsize(25)
    else:
        text.set_fontsize(18)

# Draw stacked horizontal bar plot of species composition
bottom = np.zeros(len(prop_table))
handles = []
labels = []

for i, species in enumerate(prop_table.columns):
    values = prop_table[species].values
    bar = ax2.barh(prop_table.index, values, left=bottom,
                   color=palette_bar[i % len(palette_bar)],
                   height=0.6,
                   edgecolor="none")
    bottom += values
    handles.append(bar[0])
    labels.append(species)

# Customize stacked bar plot appearance
ax2.set_yticks(range(len(prop_table)))
ax2.set_yticklabels(prop_table.index, fontsize=15)
ax2.set_xlabel("Species ratio", fontsize=30)
ax2.set_ylabel("")
ax2.set_xticks([])
ax2.set_xticklabels([])

ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.spines['left'].set_visible(True)

# Add species legend
legend2 = ax2.legend(handles, labels, title="Core species",
                     fontsize=15, title_fontsize=25,
                     labelspacing=0.3,
                     ncol=1,
                     bbox_to_anchor=(-0.1, -0.05),
                     loc='upper left')

# Apply italic formatting to species names
for text in legend2.get_texts():
    text.set_fontstyle('italic')

for patch in legend2.get_patches():
    patch.set_edgecolor("none")

# Finalize layout and export figure with editable text
plt.tight_layout()
plt.savefig("Figure4e.pdf", dpi=300, bbox_inches='tight')
plt.show()