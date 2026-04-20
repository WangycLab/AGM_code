# -*- coding: utf-8 -*-
"""
Created on Mon Jan 5 19:13:28 2026
@author: ZHENG XINGHAI
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
import matplotlib as mpl

# Configure matplotlib to embed editable fonts
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.sans-serif'] = ['Arial']

# Load hub gene fraction, expression, and annotation tables
fraction_df = pd.read_csv("hub_gene_fraction_by_timepoint.csv", index_col=0)
expr_df = pd.read_csv("hub_gene_mean_expr_by_timepoint.csv", index_col=0)
annot_df = pd.read_csv("hub_genes_top50_anno.csv")

# Select hub genes with valid COG annotation and highest module connectivity
annot_df = annot_df[annot_df['COG_ID'] != "-"]
top_genes_df = (
    annot_df
    .sort_values(['module', 'kME'], ascending=[True, False])
    .groupby('module')
    .head(6)
)
hub_genes = top_genes_df['gene_name'].tolist()

print("Number of selected genes:", len(hub_genes))

# Filter fraction and expression matrices to selected hub genes
fraction_df = fraction_df.loc[fraction_df.index.isin(hub_genes)]
expr_df = expr_df.loc[expr_df.index.isin(hub_genes)]

if fraction_df.shape[0] == 0:
    raise ValueError("No genes passed the COG_ID and kME filtering criteria.")

# Apply row-wise min–max normalization to scale values between 0 and 1
scaled_fraction_df = fraction_df.apply(
    lambda x: (x - x.min()) / (x.max() - x.min()) if x.max() != x.min() else 0,
    axis=1
)
scaled_expr_df = expr_df.apply(
    lambda x: (x - x.min()) / (x.max() - x.min()) if x.max() != x.min() else 0,
    axis=1
)

# Perform hierarchical clustering to determine gene display order
clust_df = scaled_fraction_df.dropna()
dist = pdist(clust_df.values, metric="euclidean")
Z = linkage(dist, method="average")
gene_order = clust_df.index[leaves_list(Z)].tolist()

# Reshape fraction and expression matrices into long format and merge
fraction_melt = scaled_fraction_df.reset_index().melt(
    id_vars='index', var_name='cluster', value_name='fraction'
).rename(columns={'index': 'gene'})

expr_melt = scaled_expr_df.reset_index().melt(
    id_vars='index', var_name='cluster', value_name='expression'
).rename(columns={'index': 'gene'})

merged = pd.merge(fraction_melt, expr_melt, on=['gene', 'cluster'])

# Define ordered timepoints and gene ordering for plotting
timepoint_order = ["L-60", "L-30", "FD30", "FD90", "FD150", "R+1", "R+7", "R+14"]
merged['cluster'] = pd.Categorical(merged['cluster'], categories=timepoint_order, ordered=True)
merged['gene'] = pd.Categorical(merged['gene'], categories=gene_order, ordered=True)
merged = merged.sort_values(['gene', 'cluster'])

# Limit scaled values within the 0–1 range
merged['fraction'] = merged['fraction'].clip(0, 1)
merged['expression'] = merged['expression'].clip(0, 1)

# Generate bubble plot showing expression level and expressing cell fraction
fig, ax = plt.subplots(figsize=(7.5, 16))

sns.scatterplot(
    data=merged,
    x='cluster',
    y='gene',
    size='fraction',
    hue='expression',
    sizes=(10, 200),
    palette='flare',
    edgecolor='none',
    linewidth=0,
    ax=ax
)

ax.set_xlabel("Timepoint", fontsize=30)
ax.set_ylabel("Hub gene", fontsize=30)

ax.tick_params(axis='x', rotation=90, labelsize=18)
ax.tick_params(axis='y', labelsize=15)

for label in ax.get_yticklabels():
    label.set_fontstyle('italic')

# Adjust legend labels and layout for publication-quality figure
leg = ax.get_legend()

for text in leg.get_texts():
    if text.get_text() == "expression":
        text.set_text("S. Exp. R")
        text.set_fontsize(18)
    elif text.get_text() == "fraction":
        text.set_text("S. Exp. L")
        text.set_fontsize(18)
    else:
        text.set_fontsize(15)

leg.set_frame_on(True)
leg.set_bbox_to_anchor((1.05, 0.5))
leg._loc = 6

# Save the figure as editable vector PDF
plt.tight_layout()
plt.savefig("Figure6c.pdf", dpi=300, bbox_inches='tight')
plt.show()