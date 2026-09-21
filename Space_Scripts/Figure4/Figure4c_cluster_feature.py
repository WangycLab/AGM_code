# -*- coding: utf-8 -*-
"""
Created on Sun Aug 24 18:32:42 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Plot settings
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.family"] = "Arial"

# Load data
fraction_df = pd.read_csv("gene_fraction_by_cluster.csv", index_col=0)
expr_df = pd.read_csv("gene_mean_expr_by_cluster.csv", index_col=0)
annot_df = pd.read_csv("all_clusters_DEGs_annotated.csv")

# Filter and merge COG annotations
annot_df = annot_df[~annot_df["COG_ID"].isin(["-"]) & ~annot_df["COG_ID"].str.contains(",", na=False)].copy()
annot_df["gene_key"] = annot_df["gene"].str.replace("_", "-", regex=False)
fraction_df = fraction_df.merge(annot_df[["gene_key", "COG_ID"]], left_index=True, right_on="gene_key", how="inner")
expr_df = expr_df.merge(annot_df[["gene_key", "COG_ID"]], left_index=True, right_on="gene_key", how="inner")

# Convert to long format
fraction_melt = fraction_df.melt(id_vars=["gene_key", "COG_ID"], var_name="cluster", value_name="fraction")
expr_melt = expr_df.melt(id_vars=["gene_key", "COG_ID"], var_name="cluster", value_name="expression")
merged = pd.merge(fraction_melt, expr_melt, on=["gene_key", "COG_ID", "cluster"])

# Selected genes
selected_genes = [
    "MGYG000002478-03950", "BVU-RS19770", "MGYG000004797-03948", "MGYG000002478-00488", "BVU-RS13285",
    "MGYG000002284-01625", "MGYG000002284-01626", "MGYG000002284-00923", "MGYG000002284-00002", "MGYG000002284-01630",
    "MGYG000004775-02249", "MGYG000003210-00467", "MGYG000001556-01599", "MGYG000004491-00230", "MGYG000003115-01051",
    "MGYG000002438-03104", "MGYG000002438-00285", "MGYG000002438-02238", "MGYG000002438-03535", "MGYG000002438-00512",
    "MGYG000002284-00549", "MGYG000004755-01461", "MGYG000002221-00833",
    "MGYG000001763-00916", "MGYG000000215-02273", "MGYG000001763-02187", "MGYG000002967-00452", "MGYG000000215-00564",
    "MGYG000001263-04486", "MGYG000000243-00572", "MGYG000002527-00970", "MGYG000000198-03147",
    "MGYG000001763-01556", "MGYG000000215-01894", "MGYG000001763-00774", "MGYG000002478-01478", "MGYG000000243-00574",
    "MGYG000001045-03032", "MGYG000000092-01437"
]

# Species colors
species_color_map = {
    "Agathobacter faecis": "#8A97A8", "Agathobacter rectalis": "#B4C2D3",
    "Bacteroides eggerthii": "#B4D1AC", "Bacteroides fragilis": "#D8A27E", "Bacteroides intestinalis": "#E8C8B0",
    "Bacteroides ovatus": "#B8A878", "Bacteroides stercoris": "#CBB987", "Bacteroides thetaiotaomicron": "#D8C89A",
    "Bacteroides uniformis": "#A89B68", "Bacteroides xylanisolvens": "#E2D2B0",
    "CAG-81 sp900066785": "#C8C2BC", "Clostridium_Q sp003024715": "#C898A8",
    "Enterocloster bolteae": "#BCA8C8", "Enterocloster sp000431375": "#CDB6D6", "Enterocloster sp001517625": "#A888B8",
    "Escherichia coli_D": "#9674A8", "Faecalibacterium prausnitzii_C": "#7FA6D8",
    "Fusicatenibacter saccharivorans": "#D0B8D8", "Fusobacterium_A mortiferum": "#8F8F78", "Fusobacterium_A varium": "#B1B194",
    "Megamonas funiformis": "#80A8A8", "Parabacteroides distasonis": "#A8C8C8", "Parabacteroides merdae": "#D88888",
    "Parasutterella excrementihominis": "#F0B8B8", "Phascolarctobacterium faecium": "#C9A6D8",
    "Phocaeicola dorei": "#E8D090", "Phocaeicola massiliensis": "#D8C4A0", "Phocaeicola vulgatus": "#A88C78",
    "Prevotella stercorea": "#E0C8B8", "Roseburia intestinalis": "#8298B0", "Roseburia sp900552665": "#B0C8D8",
    "Sutterella wadsworthensis": "#98BC98", "UBA7182 sp003480725": "#C48898", "Others": "#BDBDBD"
}

# Select genes and assign unique COG labels
plot_df = merged[merged["gene_key"].isin(selected_genes)].copy()
gene_cog_map = plot_df[["gene_key", "COG_ID"]].drop_duplicates()
cog_counts = gene_cog_map["COG_ID"].value_counts()
cog_index, unique_map = {}, {}

for _, row in gene_cog_map.iterrows():
    gene, cog = row["gene_key"], row["COG_ID"]
    if cog_counts[cog] > 1:
        cog_index[cog] = cog_index.get(cog, 0) + 1
        unique_map[gene] = f"{cog} ({cog_index[cog]})"
    else:
        unique_map[gene] = cog

plot_df["COG_label"] = plot_df["gene_key"].map(unique_map)
cog_order = [unique_map[g] for g in selected_genes if g in unique_map]
plot_df["COG_label"] = pd.Categorical(plot_df["COG_label"], categories=cog_order, ordered=True)

# Clip and sort values
plot_df["fraction"] = plot_df["fraction"].clip(upper=1)
plot_df["expression"] = plot_df["expression"].clip(upper=1)
plot_df = plot_df.sort_values(["cluster", "COG_label"])

# Species composition
df_meta = pd.read_csv("cell_metadata.tsv", sep="\t")
df_meta["cluster"] = df_meta["cluster"].apply(lambda x: f"Cluster_{x}")
count_table = pd.crosstab(df_meta["cluster"], df_meta["species"])
prop_table = count_table.div(count_table.sum(axis=1), axis=0)

def keep_top_n(row, n=3):
    top_species = row.nlargest(n).index
    new_row = pd.Series(0, index=row.index)
    new_row[top_species] = row[top_species]
    return new_row

prop_table = prop_table.apply(keep_top_n, axis=1)
prop_table = prop_table.loc[:, (prop_table != 0).any(axis=0)]
cluster_order = [f"Cluster_{i}" for i in range(15)]
prop_table = prop_table.reindex(cluster_order)

# Normalize within each cluster
plot_df[["fraction", "expression"]] = plot_df.groupby("cluster")[["fraction", "expression"]].transform(lambda x: (x - x.min()) / (x.max() - x.min()))

# Cluster annotation
cluster_annotation = {str(i): f"C{i}" for i in range(15)}
plot_df["cluster"] = plot_df["cluster"].astype(str).str.replace("Cluster_", "", regex=False)
cluster_order = [str(i) for i in range(15)]
cluster_label_order = [cluster_annotation[i] for i in cluster_order]
plot_df["cluster_label"] = pd.Categorical(plot_df["cluster"].map(cluster_annotation), categories=cluster_label_order, ordered=True)

# Prepare species panel
prop_table.index = prop_table.index.astype(str).str.replace("Cluster_", "", regex=False)
prop_table = prop_table.reindex(cluster_order)
prop_table.index = [cluster_annotation[i] for i in cluster_order]

# Create figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 12), sharey=True, gridspec_kw={"width_ratios": [3, 1]})

# Bubble plot
scatter = sns.scatterplot(data=plot_df, x="COG_label", y="cluster_label", size="fraction", hue="expression", sizes=(10, 200), palette="crest", edgecolor="black", ax=ax1)
ax1.set_xlabel("COG IDs", fontsize=30)
ax1.set_ylabel("Functional clusters", fontsize=30)
ax1.tick_params(axis="x", rotation=90, labelsize=20)
ax1.tick_params(axis="y", labelsize=20)

# Bubble legend
handles, labels = ax1.get_legend_handles_labels()
new_labels, new_handles = [], []

for h, l in zip(handles, labels):
    if l == "fraction":
        new_labels.append("S. Exp. R"); new_handles.append(h)
    elif l == "expression":
        new_labels.append("S. Exp. L"); new_handles.append(h)
    elif l not in ["size", "hue"]:
        new_labels.append(l); new_handles.append(h)

leg = ax1.legend(new_handles, new_labels, title="", ncol=2, fontsize=18, bbox_to_anchor=(0.05, -0.5), loc="upper left")
for text in leg.get_texts(): text.set_fontsize(25 if text.get_text() in ["S. Exp. L", "S. Exp. R"] else 18)

# Species composition
bottom = np.zeros(len(prop_table))
handles, labels = [], []

for species in prop_table.columns:
    values = prop_table[species].values
    bar = ax2.barh(prop_table.index, values, left=bottom, color=species_color_map.get(species, "#CCCCCC"), height=0.6, edgecolor="none")
    bottom += values
    handles.append(bar[0]); labels.append(species)

ax2.set_yticks(range(len(prop_table)))
ax2.set_yticklabels(prop_table.index, fontsize=20)
ax2.set_xlabel("Species ratio", fontsize=30)
ax2.set_xticks([0, 0.25, 0.5, 0.75, 1])
ax2.set_xticklabels(["0", "0.25", "0.5", "0.75", "1.0"], fontsize=20)
ax2.set_ylabel("")
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.spines["bottom"].set_visible(False)

# Species legend
legend2 = ax2.legend(handles, labels, title="Core species", fontsize=18, title_fontsize=25, labelspacing=0.3, bbox_to_anchor=(-0.1, -0.15), loc="upper left")
for text in legend2.get_texts(): text.set_fontstyle("italic")
for patch in legend2.get_patches(): patch.set_edgecolor("none")

# Save figure
plt.tight_layout()
plt.savefig("Figure4c.pdf", dpi=300, bbox_inches="tight")
plt.show()