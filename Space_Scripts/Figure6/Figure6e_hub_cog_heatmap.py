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

# Plot settings
matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "Arial"})
AXIS_LABEL_SIZE, TICK_LABEL_SIZE = 25, 18
CBAR_WIDTH, CBAR_LENGTH, CBAR_TICK_SIZE, CBAR_LABEL_SIZE = 0.1, 1, 15, 20
module_order = ["blue", "brown", "green", "turquoise", "yellow"]

# Load enrichment results
input_dir = "COG_enrichment_by_module"
files = glob.glob(os.path.join(input_dir, "COG_enrichment_module_*.csv"))
q_value_dict, rich_factor_dict = {}, {}

for file in files:
    module = os.path.basename(file).replace("COG_enrichment_module_", "").replace(".csv", "")
    df = pd.read_csv(file)
    df = df[(df["overlap"] > 0) & df["significant"] & (df["rich_factor"] > 0.0005)]
    q_value_dict[module] = df.set_index("COG_ID")["q_value"]
    rich_factor_dict[module] = df.set_index("COG_ID")["rich_factor"]

q_value_df = pd.DataFrame(q_value_dict).reindex(columns=module_order).fillna(1)
rich_factor_df = pd.DataFrame(rich_factor_dict).reindex(columns=module_order).fillna(0)

# Keep significantly enriched COGs
sig_mask = (q_value_df < 0.05).any(axis=1)
q_value_df, rich_factor_df = q_value_df.loc[sig_mask], rich_factor_df.loc[sig_mask]

# Load COG annotations
cog_anno = pd.read_csv("cog-24.def.tab", sep="\t", header=None, usecols=[0, 2], on_bad_lines="skip", engine="python")
cog_anno.columns = ["COG_ID", "Annotation"]
cog_anno_dict = cog_anno.set_index("COG_ID")["Annotation"].to_dict()

rich_factor_df.index = rich_factor_df.index.map(lambda x: cog_anno_dict.get(x, x))
q_value_df.index = q_value_df.index.map(lambda x: cog_anno_dict.get(x, x))

# Abbreviate functional descriptions
abbreviation_dict = {
    "Glyceraldehyde-3-phosphate dehydrogenase/erythrose-4-phosphate dehydrogenase": "GAPDH/E4PDH",
    "Small heat shock protein IbpA, HSP20 family": "Small HSP IbpA",
    "6-phosphogluconate dehydrogenase": "6PG dehydrogenase",
    "Glucose-6-phosphate 1-dehydrogenase": "G6P dehydrogenase",
    "DNA polymerase B elongation subunit": "DNA polymerase B",
    "Chaperonin GroEL (HSP60 family)": "Chaperonin GroEL",
    "Isocitrate dehydrogenase": "Isocitrate dehydrogenase",
    "Phosphoglycerate mutase (BPG-dependent)": "Phosphoglycerate mutase",
    "Superoxide dismutase": "Superoxide dismutase",
    "Superfamily I DNA and/or RNA helicase": "SF1 helicase",
    "ABC-type dipeptide/oligopeptide/nickel transport system, ATPase component": "ABC peptide/Ni ATPase",
    "Cobalamin biosynthesis protein CobN, Mg-chelatase": "CobN cobalamin synthase",
    "Redox-sensitive bicupin YhaK, pirin superfamily": "Redox bicupin YhaK",
    "Fructose-bisphosphate aldolase class Ia, DhnA family": "FBP aldolase",
    "Methylmalonyl-CoA mutase, N-terminal domain/subunit": "MCM N-terminal",
    "Glucuronate isomerase": "Glucuronate isomerase",
    "Threonine aldolase": "Threonine aldolase",
    "Peroxiredoxin": "Peroxiredoxin",
    "N-acetylglucosaminyl deacetylase, LmbE family": "GlcNAc deacetylase",
    "Methylmalonyl-CoA mutase, C-terminal domain/subunit (cobalamin-binding)": "MCM C-terminal",
    "Uncharacterized membrane protein YphA, DoxX/SURF4 family": "YphA membrane protein",
    "Ferritin-like DNA-binding protein, DPS (DNA Protection under Starvation) family": "Ferritin-like DPS",
    "Antitoxin component YwqK of the YwqJK toxin-antitoxin module": "YwqK antitoxin",
    "Uncharacterized lipoprotein NlpE involved in copper resistance": "Lipoprotein NlpE",
    "Putative hemolysin": "Hemolysin",
    "Outer membrane protein OmpW": "Outer membrane OmpW",
    "Cell division protein ZapB, interacts with FtsZ": "Cell division ZapB",
    "Heat shock protein HslJ": "HslJ heat shock protein",
    "Periplasmic ligand-binding sensor domain": "Periplasmic sensor",
    "Iron-binding extracellular (periplasmic) sensor domain CHASE4": "CHASE4 iron sensor",
    "Accessory protein of IS66 transposable element": "IS66 accessory protein",
    "Extracellular (periplasmic) sensor domain CHASE (specificity unknown)": "CHASE sensor domain",
    "Transposase": "Transposase",
    "FRM2-like oxidoreductase, nitroreductase family": "FRM2 oxidoreductase",
    "Opacity protein LomR and related surface antigens": "LomR opacity protein",
    "Transposase, IS1182 family": "IS1182 transposase",
    "Toxin component of the Txe-Axe toxin-antitoxin module, Txe/YoeB family": "Txe/YoeB toxin",
    "Mu-like prophage I protein": "Mu-like prophage protein",
    "Adenine-specific DNA methylase, N12 class": "DNA adenine methylase",
    "Outer membrane receptor for monomeric catechols": "Catechol OM receptor",
    "Phage-related tail fiber protein": "Phage tail fiber",
    "Preprotein translocase subunit Sec63": "Translocase Sec63",
    "Protein involved in initiation of plasmid replication": "Plasmid replication initiator",
    "Self-loading helicase or an inactivated derivative": "Self-loading helicase"
}

rich_factor_df.index = rich_factor_df.index.map(lambda x: abbreviation_dict.get(x, x))
q_value_df.index = q_value_df.index.map(lambda x: abbreviation_dict.get(x, x))

# Normalize each module to 0–1
rich_factor_df_norm = rich_factor_df.apply(lambda col: (col - col.min()) / (col.max() - col.min()) if col.max() > col.min() else 0)

# Cluster rows
cg = sns.clustermap(rich_factor_df_norm, method="average", metric="euclidean", row_cluster=True, col_cluster=False, cmap="Reds", figsize=(1, 1))
row_order = cg.dendrogram_row.reordered_ind
plt.close(cg.fig)
rich_factor_df_norm = rich_factor_df_norm.iloc[row_order]
q_value_df = q_value_df.iloc[row_order]

# Plot heatmap
fig, ax = plt.subplots(figsize=(2, 15))
sns.heatmap(rich_factor_df_norm, cmap="Reds", linewidths=1, linecolor="grey", annot=False, cbar_kws={"label": "Rich factor", "fraction": CBAR_WIDTH, "shrink": CBAR_LENGTH}, ax=ax)

# Add significance stars
for i, j in zip(*((q_value_df < 0.05).to_numpy()).nonzero()):
    ax.text(j + 0.5, i + 0.95, "*", ha="center", va="center", color="black", fontsize=30, fontweight="bold")

# Customize axes
ax.set(xlabel="Module", ylabel="COG function")
ax.set_xlabel("Module", fontsize=AXIS_LABEL_SIZE)
ax.set_ylabel("COG function", fontsize=AXIS_LABEL_SIZE)
ax.tick_params(axis="x", labelsize=TICK_LABEL_SIZE, rotation=90)
ax.tick_params(axis="y", labelsize=TICK_LABEL_SIZE)

# Customize colorbar
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=CBAR_TICK_SIZE)
cbar.set_label("Scaled rich factor", fontsize=CBAR_LABEL_SIZE)

# Save figure
plt.tight_layout()
plt.savefig("Figure6e.pdf", dpi=300, bbox_inches="tight")
plt.show()