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

mpl.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "Arial", "font.sans-serif": ["Arial"]})

# Load data
fraction_df = pd.read_csv("hub_gene_fraction_by_timepoint.csv", index_col=0)
expr_df = pd.read_csv("hub_gene_mean_expr_by_timepoint.csv", index_col=0)
annot_df = pd.read_csv("hub_genes_top50_anno.csv")

# Select top 3 hub genes per module
top_genes_df = annot_df[annot_df["COG_ID"] != "-"].sort_values(["module", "kME"], ascending=[True, False]).groupby("module").head(3)
hub_genes = top_genes_df["gene_name"].tolist()
print(f"Selected genes: {len(hub_genes)}")

# Create COG labels
label_df = top_genes_df[["gene_name", "COG_ID", "module"]].copy()
label_df["COG_ID_label"] = label_df.groupby("COG_ID").cumcount().add(1).astype(str)
cog_counts = label_df["COG_ID"].value_counts()
label_df["COG_ID_label"] = label_df.apply(lambda x: x["COG_ID"] if cog_counts[x["COG_ID"]] == 1 else f'{x["COG_ID"]}({x["COG_ID_label"]})', axis=1)
gene_to_cog = dict(zip(label_df["gene_name"], label_df["COG_ID_label"]))
cog_to_module = dict(zip(label_df["COG_ID_label"], label_df["module"]))

# Define module colors
module_colors = {"blue": "#377EB8", "brown": "#A65628", "green": "#4DAF4A", "red": "#E41A1C", "turquoise": "#00A6A6", "yellow": "#FFD92F"}

# Filter and normalize matrices
fraction_df = fraction_df.loc[fraction_df.index.isin(hub_genes)]
expr_df = expr_df.loc[expr_df.index.isin(hub_genes)]
if fraction_df.empty:
    raise ValueError("No genes passed filtering.")

scaled_fraction_df = fraction_df.apply(lambda x: (x - x.min()) / (x.max() - x.min()) if x.max() != x.min() else 0, axis=1)
scaled_expr_df = expr_df.apply(lambda x: (x - x.min()) / (x.max() - x.min()) if x.max() != x.min() else 0, axis=1)

# Determine gene order by hierarchical clustering
clust_df = scaled_fraction_df.dropna()
gene_order = clust_df.index[leaves_list(linkage(pdist(clust_df.values, metric="euclidean"), method="average"))].tolist()

# Prepare plotting data
fraction_melt = scaled_fraction_df.reset_index().melt(id_vars="index", var_name="cluster", value_name="fraction").rename(columns={"index": "gene"})
expr_melt = scaled_expr_df.reset_index().melt(id_vars="index", var_name="cluster", value_name="expression").rename(columns={"index": "gene"})
merged = pd.merge(fraction_melt, expr_melt, on=["gene", "cluster"])
merged["COG_ID"] = merged["gene"].map(gene_to_cog)
merged["gene_order"] = merged["gene"].map({g: i for i, g in enumerate(gene_order)})
merged = merged.sort_values(["gene_order", "cluster"])

timepoint_order = ["L-60", "L-30", "FD30", "FD90", "FD150", "R+1", "R+7", "R+14"]
merged["cluster"] = pd.Categorical(merged["cluster"], categories=timepoint_order, ordered=True)
merged["COG_ID"] = pd.Categorical(merged["COG_ID"], categories=[gene_to_cog[g] for g in gene_order], ordered=True)
merged[["fraction", "expression"]] = merged[["fraction", "expression"]].clip(0, 1)

# Plot
fig, ax = plt.subplots(figsize=(7, 5))
sns.scatterplot(data=merged, x="cluster", y="COG_ID", size="fraction", hue="expression", sizes=(10, 100), palette="flare", edgecolor="none", linewidth=0, ax=ax)
ax.set_xlabel("Timepoint", fontsize=25)
ax.set_ylabel("COG items", fontsize=25)
ax.tick_params(axis="x", rotation=90, labelsize=18)
ax.tick_params(axis="y", labelsize=15)

# Color COG labels by module
for label in ax.get_yticklabels():
    label.set_fontstyle("italic")
    label.set_fontweight("bold")
    label.set_color(module_colors.get(cog_to_module.get(label.get_text()), "black"))

# Customize legend
leg = ax.get_legend()
for text in leg.get_texts():
    text.set_text({"expression": "S. Exp. R", "fraction": "S. Exp. L"}.get(text.get_text(), text.get_text()))
    text.set_fontsize(18 if text.get_text() in ["S. Exp. R", "S. Exp. L"] else 15)
leg.set_frame_on(True)
leg.set_bbox_to_anchor((1.05, 0.5))
leg._loc = 6

# Save figure
plt.tight_layout()
plt.savefig("Figure6f.pdf", dpi=300, bbox_inches="tight")
plt.show()