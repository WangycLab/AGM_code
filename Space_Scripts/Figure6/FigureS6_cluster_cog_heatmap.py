# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 18:53:18 2025
@author: ZHENG XINGHAI
"""

import os
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "Arial"})
AXIS_LABEL_SIZE, TICK_LABEL_SIZE = 25, 20
CBAR_WIDTH, CBAR_LENGTH, CBAR_TICK_SIZE, CBAR_LABEL_SIZE = 0.1, 0.3, 20, 25
sample_order = ["0", "1", "2", "3", "4", "5"]

# Load enrichment results
input_dir = "COG_enrichment_by_cluster"
files = glob.glob(os.path.join(input_dir, "COG_enrichment_cluster_*.csv"))
q_value_dict, rich_factor_dict = {}, {}

for file in files:
    cluster_name = os.path.basename(file).replace("COG_enrichment_cluster_", "").replace(".csv", "")
    df = pd.read_csv(file)
    df = df[(df["overlap"] > 0) & df["significant"] & (df["rich_factor"] > 0.0002)]
    q_value_dict[cluster_name] = df.set_index("COG_ID")["q_value"]
    rich_factor_dict[cluster_name] = df.set_index("COG_ID")["rich_factor"]

q_value_df = pd.DataFrame(q_value_dict).fillna(1)
rich_factor_df = pd.DataFrame(rich_factor_dict).fillna(0)
q_value_df.columns = q_value_df.columns.str.replace("COG_enrichment_cluster_", "", regex=False)
rich_factor_df.columns = rich_factor_df.columns.str.replace("COG_enrichment_cluster_", "", regex=False)
q_value_df = q_value_df.reindex(columns=sample_order)
rich_factor_df = rich_factor_df.reindex(columns=sample_order)

# Keep COGs significant in at least one cluster
significant_mask = (q_value_df < 0.05).any(axis=1)
q_value_df = q_value_df.loc[significant_mask]
rich_factor_df = rich_factor_df.loc[significant_mask]

# Load COG annotations
cog_anno = pd.read_csv("cog-24.def.tab", sep="\t", header=None, usecols=[0, 2], on_bad_lines="skip", engine="python")
cog_anno.columns = ["COG_ID", "Annotation"]
cog_anno_dict = cog_anno.set_index("COG_ID")["Annotation"].to_dict()
q_value_df.index = q_value_df.index.map(lambda x: cog_anno_dict.get(x, x))
rich_factor_df.index = rich_factor_df.index.map(lambda x: cog_anno_dict.get(x, x))

# Define abbreviations
abbreviation_dict = {
    "Malate/lactate dehydrogenase": "Malate/lactate dehydrogenase",
    "Translation elongation factor EF-Tu, a GTPase": "EF-Tu elongation factor",
    "Small heat shock protein IbpA, HSP20 family": "Small heat shock protein IbpA",
    "DNA-directed RNA polymerase, beta' subunit/160 kD subunit": "RNA polymerase beta subunit",
    "Ribosomal protein S9": "Ribosomal protein S9",
    "Anthranilate/para-aminobenzoate synthases component I": "Anthranilate/PABA synthase",
    "Enolase": "Enolase",
    "Inorganic pyrophosphatase": "Inorganic pyrophosphatase",
    "Molecular chaperone, HSP90 family": "HSP90 chaperone",
    "Glutamyl-tRNA reductase": "Glutamyl-tRNA reductase",
    "Chaperonin GroEL (HSP60 family)": "GroEL chaperonin",
    "Leucyl-tRNA synthetase": "Leu-tRNA synthetase",
    "Translation initiation factor IF-2, a GTPase": "IF-2 initiation factor",
    "Superoxide dismutase": "Superoxide dismutase",
    "Triacylglycerol esterase/lipase EstA, alpha/beta hydrolase fold": "EstA lipase",
    "ABC-type sulfate/molybdate transport systems, ATPase component": "ABC sulfate/molybdate transporter",
    "ABC-type dipeptide/oligopeptide/nickel transport system, ATPase component": "ABC peptide/Ni transporter",
    "Neutral trehalase": "Neutral trehalase",
    "Epoxyqueuosine reductase QueH (tRNA modification)": "QueH tRNA modification enzyme",
    "Siroheme synthase (precorrin-2 oxidase/ferrochelatase domain)": "Siroheme synthase",
    "D-alanyl-lipoteichoic acid acyltransferase DltB, MBOAT superfamily": "DltB lipoteichoic acid acyltransferase",
    "Lipoprotein Med, regulator of KinD/Spo0A, PBP1-ABC superfamily, includes NupN": "Med Spo0A regulator",
    "Peroxiredoxin": "Peroxiredoxin",
    "Uncharacterized conserved protein YfaS, alpha-2-macroglobulin family": "YfaS protein",
    "Cell division protein ZipA, interacts with FtsZ": "ZipA cell division protein",
    "Iron-binding extracellular (periplasmic) sensor domain CHASE4": "CHASE4 iron sensor",
    "Retron-type reverse transcriptase": "Retron reverse transcriptase",
    "Accessory protein of IS66 transposable element": "IS66 accessory protein",
    "Extracellular (periplasmic) sensor domain CHASE (specificity unknown)": "CHASE extracellular sensor",
    "Transposase": "Transposase",
    "Predicted ATP-dependent endonuclease of the OLD family, contains P-loop ATPase and TOPRIM domains": "OLD family ATP-dependent nuclease",
    "Predicted ATPase": "Predicted ATPase",
    "Predicted phospholipase, patatin/cPLA2 family": "Patatin-like phospholipase",
    "ABC-type transport system involved in cytochrome bd biosynthesis, ATPase and permease components": "ABC cytochrome bd transporter",
    "Extracytoplasmic sensor domain CHASE3 (specificity unknown)": "CHASE3 extracellular sensor",
    "Uncharacterized conserved protein, heparinase superfamily": "Heparinase-like protein",
    "Protein involved in initiation of plasmid replication": "Plasmid replication protein"
}

q_value_df.index = q_value_df.index.map(lambda x: abbreviation_dict.get(x, x))
rich_factor_df.index = rich_factor_df.index.map(lambda x: abbreviation_dict.get(x, x))

# Normalize rich factors by cluster
rich_factor_df_norm = rich_factor_df.apply(lambda col: (col - col.min()) / (col.max() - col.min()) if col.max() > col.min() else 0)

# Determine hierarchical row order
cg = sns.clustermap(rich_factor_df_norm, method="average", metric="euclidean", row_cluster=True, col_cluster=False, cmap="Reds", figsize=(1, 1))
row_order = cg.dendrogram_row.reordered_ind
plt.close()

rich_factor_df_norm = rich_factor_df_norm.iloc[row_order]
q_value_df = q_value_df.iloc[row_order][rich_factor_df_norm.columns]

# Plot heatmap
fig, ax = plt.subplots(figsize=(3, 16))
sns.heatmap(
    rich_factor_df_norm, ax=ax, cmap="Reds", linewidths=1, linecolor="grey",
    annot=False, cbar_kws={"label": "Scaled rich factor", "fraction": CBAR_WIDTH, "shrink": CBAR_LENGTH}
)

# Add significance markers
for i, j in zip(*((q_value_df < 0.05).to_numpy().nonzero())):
    ax.text(j + 0.5, i + 0.85, "*", ha="center", va="center", color="black", fontsize=30, fontweight="bold")

ax.set_xlabel("Cluster", fontsize=AXIS_LABEL_SIZE)
ax.set_ylabel("COG function", fontsize=AXIS_LABEL_SIZE)
ax.tick_params(axis="both", labelsize=TICK_LABEL_SIZE)
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=CBAR_TICK_SIZE)
cbar.set_label("Scaled rich factor", fontsize=CBAR_LABEL_SIZE)

plt.tight_layout()
plt.savefig("FigureS6.pdf", dpi=300, bbox_inches="tight")
plt.show()

# Print a concise summary
sig_counts = (q_value_df < 0.05).sum()
print(f"COG enrichment completed: {len(q_value_df)} terms across {len(q_value_df.columns)} clusters.")
print("Significant terms per cluster:", sig_counts.to_dict())