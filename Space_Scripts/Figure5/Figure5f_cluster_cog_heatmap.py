# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 18:53:18 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
import os
import glob
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

# Set global plotting parameters for consistent visual style
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "Arial"
AXIS_LABEL_SIZE = 25
TICK_LABEL_SIZE = 18
CBAR_WIDTH = 0.1
CBAR_LENGTH = 1
CBAR_TICK_SIZE = 15
CBAR_LABEL_SIZE = 20

# Define fixed order of clusters for visualization
sample_order = ["1", "4", "6"]

# Locate input CSV files containing COG enrichment results
input_dir = "COG_enrichment_by_cluster"
files = glob.glob(os.path.join(input_dir, "COG_enrichment_cluster_*.csv"))

q_value_dict = {}
rich_factor_dict = {}

# Read and filter each cluster enrichment file
for file in files:
    cluster_name = (
        os.path.basename(file)
        .replace("COG_enrichment_results_", "")
        .replace(".csv", "")
    )

    df = pd.read_csv(file)

    # Keep only significant entries with positive rich factor
    df = df[(df["overlap"] > 0) & (df["significant"]) & (df["rich_factor"] >= 0.0005)]

    q_value_dict[cluster_name] = df.set_index("COG_ID")["q_value"]
    rich_factor_dict[cluster_name] = df.set_index("COG_ID")["rich_factor"]

# Combine per-cluster results into DataFrames
q_value_df = pd.DataFrame(q_value_dict).fillna(1)
q_value_df.columns = [col.replace('COG_enrichment_cluster_', '') for col in q_value_df.columns]
rich_factor_df = pd.DataFrame(rich_factor_dict).fillna(0)
rich_factor_df.columns = [col.replace('COG_enrichment_cluster_', '') for col in rich_factor_df.columns]

# Reorder columns to ensure consistent cluster order
rich_factor_df = rich_factor_df.reindex(columns=sample_order)
q_value_df = q_value_df.reindex(columns=sample_order)

# Keep only COGs that are significant in at least one cluster
significant_mask = (q_value_df < 0.05).any(axis=1)
q_value_df = q_value_df[significant_mask]
rich_factor_df = rich_factor_df.loc[significant_mask]

# Load COG functional annotation file
cog_anno = pd.read_csv(
    "cog-24.def.tab",
    sep="\t",
    header=None,
    usecols=[0, 2],
    on_bad_lines="skip",
    engine="python"
)
cog_anno.columns = ["COG_ID", "Annotation"]
cog_anno_dict = cog_anno.set_index("COG_ID")["Annotation"].to_dict()

# Map COG_ID to descriptive annotation names
rich_factor_df.index = rich_factor_df.index.map(lambda x: cog_anno_dict.get(x, x))
q_value_df.index = q_value_df.index.map(lambda x: cog_anno_dict.get(x, x))

# Preview filtered q_value DataFrame
print("Filtered q_value_df index (only significant COGs):")
print(q_value_df.index)

# Define abbreviation dictionary for shorter y-axis labels
abbreviation_dict = {
    'Glyceraldehyde-3-phosphate dehydrogenase/erythrose-4-phosphate dehydrogenase':'GAPDH/E4P dehydrogenase',
    'Small heat shock protein IbpA, HSP20 family':'Small heat shock protein IbpA',
    'DNA-directed RNA polymerase, alpha subunit/40 kD subunit':'RNA polymerase alpha subunit',
    '6-phosphogluconate dehydrogenase':'6-phosphogluconate dehydrogenase',
    '6-phosphogluconolactonase/Glucosamine-6-phosphate isomerase/deaminase':'6-phosphogluconolactonase',
    'Glucose-6-phosphate 1-dehydrogenase':'Glucose-6-phosphate dehydrogenase',
    'DNA polymerase B elongation subunit':'DNA polymerase B',
    'Proline dehydrogenase':'Proline dehydrogenase',
    'Phosphoglycerate mutase (BPG-dependent)':'Phosphoglycerate mutase',
    'Pyruvate:ferredoxin oxidoreductase or related 2-oxoacid:ferredoxin oxidoreductase, alpha subunit':'PFOR alpha subunit',
    'RNA recognition motif (RRM) domain':'RNA recognition motif',
    'Multidrug efflux pump subunit AcrB':'AcrB efflux pump',
    'Pyruvate:ferredoxin oxidoreductase or related 2-oxoacid:ferredoxin oxidoreductase, beta subunit':'PFOR beta subunit',
    'Pyruvate:ferredoxin oxidoreductase or related 2-oxoacid:ferredoxin oxidoreductase, gamma subunit':'PFOR gamma subunit',
    'Superfamily I DNA and/or RNA helicase':'SF1 helicase',
    'Outer membrane receptor protein, Fe transport':'Outer membrane Fe receptor',
    'Fructose-bisphosphate aldolase class Ia, DhnA family':'Fructose-bisphosphate aldolase',
    'Phosphoenolpyruvate carboxykinase, ATP-dependent':'PEP carboxykinase',
    'Uncharacterized conserved protein YfaS, alpha-2-macroglobulin family':'YfaS protein',
    'Uncharacterized lipoprotein NlpE involved in copper resistance':'Lipoprotein NlpE',
    'Cell division protein ZapB, interacts with FtsZ':'Cell division protein ZapB',
    'Cell division protein ZipA, interacts with FtsZ':'Cell division protein ZipA',
    'Periplasmic ligand-binding sensor domain':'Periplasmic ligand sensor',
    'Retron-type reverse transcriptase':'Retron reverse transcriptase',
    'Accessory protein of IS66 transposable element':'IS66 accessory protein',
    'Transposase':'Transposase',
    'Putative alpha-1,2-mannosidase, GH92 family':'Alpha-1,2-mannosidase',
    'Outer membrane cobalamin receptor protein BtuB':'Cobalamin receptor BtuB',
    'Delta 1-pyrroline-5-carboxylate dehydrogenase':'P5C dehydrogenase',
    'Predicted exporter':'Predicted exporter',
    'TPR-like repeat domain':'TPR repeat protein',
    'Uncharacterized conserved protein':'Conserved hypothetical protein',
    'Outer membrane receptor for ferrienterochelin and colicins':'Ferrienterochelin receptor',
    'Outer membrane receptor for ferric coprogen and ferric-rhodotorulic acid':'Coprogen receptor',
    'Outer membrane receptor for monomeric catechols':'Catechol receptor',
    'Uncharacterized N-terminal domain of tricorn protease, contains WD40 repeats':'Tricorn protease WD40',
    'Phage-related tail fiber protein':'Phage tail fiber protein',
    'Protein involved in initiation of plasmid replication':'Plasmid replication initiator',
    'Self-loading helicase or an inactivated derivative':'Self-loading helicase'
}

# Apply abbreviations to y-axis labels
rich_factor_df.index = rich_factor_df.index.map(lambda x: abbreviation_dict.get(x, x))
q_value_df.index = q_value_df.index.map(lambda x: abbreviation_dict.get(x, x))

# Normalize rich factor values between 0 and 1 per cluster
rich_factor_df_norm = rich_factor_df.copy()
for col in rich_factor_df_norm.columns:
    col_min = rich_factor_df_norm[col].min()
    col_max = rich_factor_df_norm[col].max()
    if col_max > col_min:
        rich_factor_df_norm[col] = (rich_factor_df_norm[col] - col_min) / (col_max - col_min)
    else:
        rich_factor_df_norm[col] = 0

# Determine hierarchical row order for heatmap clustering
cg = sns.clustermap(
    rich_factor_df_norm,
    method="average",
    metric="euclidean",
    row_cluster=True,
    col_cluster=False,
    cmap="Reds",
    figsize=(1, 1)
)
row_order = cg.dendrogram_row.reordered_ind
plt.close()

# Reorder rows in normalized matrix and q-value matrix
rich_factor_df_norm = rich_factor_df_norm.iloc[row_order, :]
q_value_df = q_value_df.iloc[row_order, :]
q_value_df = q_value_df[rich_factor_df_norm.columns]

# Plot heatmap with significance marks
plt.figure(figsize=(1.5, 15))
ax = sns.heatmap(
    rich_factor_df_norm,
    cmap="Reds",
    linewidths=1,
    linecolor="grey",
    annot=False,
    cbar_kws={
        "label": "Scaled rich factor",
        "fraction": CBAR_WIDTH,
        "shrink": CBAR_LENGTH
    }
)

# Overlay stars for statistically significant enrichments
for i in range(rich_factor_df_norm.shape[0]):
    for j in range(rich_factor_df_norm.shape[1]):
        if q_value_df.iloc[i, j] < 0.05:
            ax.text(
                j + 0.5,
                i + 0.85,
                "*",
                ha="center",
                va="center",
                color="black",
                fontsize=30,
                fontweight="bold"
            )

# Customize axis labels, ticks, and colorbar
ax.set_xlabel("Cluster", fontsize=AXIS_LABEL_SIZE)
ax.set_ylabel("COG function", fontsize=AXIS_LABEL_SIZE)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=TICK_LABEL_SIZE)
ax.tick_params(axis="x", labelsize=TICK_LABEL_SIZE)
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=CBAR_TICK_SIZE)
cbar.set_label("Scaled rich factor", fontsize=CBAR_LABEL_SIZE)

# Save the figure as PDF
plt.tight_layout()
plt.savefig("Figure5f.pdf", dpi=300, bbox_inches="tight")
plt.show()