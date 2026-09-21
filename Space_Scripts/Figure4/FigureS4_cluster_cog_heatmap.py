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

# Plot settings
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "Arial"
AXIS_LABEL_SIZE, TICK_LABEL_SIZE, CBAR_WIDTH, CBAR_LENGTH, CBAR_TICK_SIZE, CBAR_LABEL_SIZE = 25, 20, 0.1, 0.3, 20, 25
sample_order = [str(i) for i in range(15)]

# Load COG enrichment results
input_dir = "COG_enrichment_by_cluster"
files = glob.glob(os.path.join(input_dir, "COG_enrichment_cluster_*.csv"))
q_value_dict, rich_factor_dict = {}, {}

for file in files:
    cluster_name = os.path.basename(file).replace("COG_enrichment_results_", "").replace(".csv", "")
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
q_value_df, rich_factor_df = q_value_df.loc[significant_mask], rich_factor_df.loc[significant_mask]

# Load COG annotations
cog_anno = pd.read_csv("cog-24.def.tab", sep="\t", header=None, usecols=[0, 2], on_bad_lines="skip", engine="python")
cog_anno.columns = ["COG_ID", "Annotation"]
cog_anno_dict = cog_anno.set_index("COG_ID")["Annotation"].to_dict()
rich_factor_df.index = rich_factor_df.index.map(lambda x: cog_anno_dict.get(x, x))
q_value_df.index = q_value_df.index.map(lambda x: cog_anno_dict.get(x, x))

# Shorten COG labels
abbreviation_dict = {
    "Malate/lactate dehydrogenase": "Malate/lactate dehydrogenase",
    "Translation elongation factor EF-Tu, a GTPase": "EF-Tu elongation factor",
    "Glyceraldehyde-3-phosphate dehydrogenase/erythrose-4-phosphate dehydrogenase": "GAPDH/E4P dehydrogenase",
    "Ribosomal protein S9": "Ribosomal protein S9",
    "Anthranilate/para-aminobenzoate synthases component I": "Anthranilate/PABA synthase",
    "Enolase": "Enolase",
    "Preprotein translocase subunit SecY": "SecY translocase",
    "Inorganic pyrophosphatase": "Inorganic pyrophosphatase",
    "6-phosphogluconate dehydrogenase": "6-phosphogluconate dehydrogenase",
    "Citrate synthase": "Citrate synthase",
    "Glutamyl-tRNA reductase": "Glutamyl-tRNA reductase",
    "Propionyl CoA:succinate CoA transferase": "Propionyl-CoA transferase",
    "Phosphoenolpyruvate synthase/pyruvate phosphate dikinase": "PEP synthase/Pyruvate PDK",
    "mRNA degradation ribonuclease J1/J2": "RNase J1/J2",
    "Superoxide dismutase": "Superoxide dismutase",
    "Alanine dehydrogenase (includes sporulation protein SpoVN)": "Alanine dehydrogenase",
    "Aspartate ammonia-lyase": "Aspartate ammonia-lyase",
    "Aconitase A": "Aconitase",
    "Triacylglycerol esterase/lipase EstA, alpha/beta hydrolase fold": "EstA lipase",
    "Superfamily I DNA and/or RNA helicase": "SF1 helicase",
    "CO dehydrogenase/acetyl-CoA synthase alpha subunit": "CODH/ACS alpha subunit",
    "Flagellar motor protein MotB": "Flagellar protein MotB",
    "FAD:protein FMN transferase ApbE": "FMN transferase ApbE",
    "Endolytic transglycosylase MltG, terminates peptidoglycan polymerization": "MltG transglycosylase",
    "TRAP-type C4-dicarboxylate transport system, large permease component": "TRAP dicarboxylate transporter",
    "Siroheme synthase (precorrin-2 oxidase/ferrochelatase domain)": "Siroheme synthase",
    "Fructose-bisphosphate aldolase class Ia, DhnA family": "Fructose-bisphosphate aldolase",
    "Phosphoenolpyruvate carboxykinase, ATP-dependent": "PEP carboxykinase",
    "LD-carboxypeptidase LdcB, LAS superfamily": "LD-carboxypeptidase LdcB",
    "Pyruvate-formate lyase": "Pyruvate-formate lyase",
    "Methylmalonyl-CoA mutase, N-terminal domain/subunit": "Methylmalonyl-CoA mutase N-term",
    "Peroxiredoxin": "Peroxiredoxin",
    "Sirohydrochlorin ferrochelatase": "Ferrochelatase",
    "Methylmalonyl-CoA mutase, C-terminal domain/subunit (cobalamin-binding)": "Methylmalonyl-CoA mutase C-term",
    "Mg/Co/Ni transporter MgtE (contains CBS domain)": "MgtE metal transporter",
    "Uncharacterized conserved protein YfaS, alpha-2-macroglobulin family": "YfaS protein",
    "Aldehyde:ferredoxin oxidoreductase": "Aldehyde ferredoxin oxidoreductase",
    "N-acetyl-anhydromuramyl-L-alanine amidase AmpD": "AmpD amidase",
    "PepSY domain containing protein, regulator of zincin peptidase activity": "PepSY regulator",
    "Predicted periplasmic protein, DUF2271 domain": "DUF2271 periplasmic protein",
    "Transposase, IS1182 family": "IS1182 transposase",
    "Signal transduction histidine kinase regulating C4-dicarboxylate transport system": "C4-dicarboxylate sensor kinase",
    "Mg2+ and Co2+ transporter CorB, contains DUF21, CBS pair, and CorC-HlyC domains": "CorB metal transporter",
    "Predicted phospholipase, patatin/cPLA2 family": "Patatin-like phospholipase",
    "Endo-beta-N-acetylglucosaminidase D": "Endo-beta-N-acetylglucosaminidase",
    "Outer membrane receptor for monomeric catechols": "Catechol receptor",
    "Acetyl-CoA carboxylase, carboxyltransferase component": "Acetyl-CoA carboxylase",
    "Pyruvate/oxaloacetate carboxyltransferase": "Pyruvate carboxyltransferase",
    "Signaling protein combining a Ser/Thr protein kinase domain and the GUN4/Ycf53 porphyrin-binding domain": "Ser/Thr kinase-GUN4 signaling protein"
}

rich_factor_df.index = rich_factor_df.index.map(lambda x: abbreviation_dict.get(x, x))
q_value_df.index = q_value_df.index.map(lambda x: abbreviation_dict.get(x, x))

# Normalize rich factors
rich_factor_df_norm = rich_factor_df.copy()
for col in rich_factor_df_norm.columns:
    col_min, col_max = rich_factor_df_norm[col].min(), rich_factor_df_norm[col].max()
    rich_factor_df_norm[col] = (rich_factor_df_norm[col] - col_min) / (col_max - col_min) if col_max > col_min else 0

# Cluster rows
cg = sns.clustermap(rich_factor_df_norm, method="average", metric="euclidean", row_cluster=True, col_cluster=False, cmap="Reds", figsize=(1, 1))
row_order = cg.dendrogram_row.reordered_ind
plt.close()
rich_factor_df_norm = rich_factor_df_norm.iloc[row_order, :]
q_value_df = q_value_df.iloc[row_order, :][rich_factor_df_norm.columns]

# Plot heatmap
plt.figure(figsize=(12, 20))
ax = sns.heatmap(rich_factor_df_norm, cmap="Reds", linewidths=1, linecolor="grey", annot=False, cbar_kws={"label": "Scaled rich factor", "fraction": CBAR_WIDTH, "shrink": CBAR_LENGTH})

for i in range(rich_factor_df_norm.shape[0]):
    for j in range(rich_factor_df_norm.shape[1]):
        if q_value_df.iloc[i, j] < 0.05:
            ax.text(j + 0.5, i + 0.85, "*", ha="center", va="center", color="black", fontsize=30, fontweight="bold")

ax.set_xlabel("Cluster", fontsize=AXIS_LABEL_SIZE)
ax.set_ylabel("COG function", fontsize=AXIS_LABEL_SIZE)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=TICK_LABEL_SIZE)
ax.tick_params(axis="x", labelsize=TICK_LABEL_SIZE)
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=CBAR_TICK_SIZE)
cbar.set_label("Scaled rich factor", fontsize=CBAR_LABEL_SIZE)

# Save figure
plt.tight_layout()
plt.savefig("FigureS4.pdf", dpi=300, bbox_inches="tight")
plt.show()