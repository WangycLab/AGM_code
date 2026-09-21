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
sample_order = ["Phocaeicola", "Bacteroides", "All_species"]

# Load COG enrichment results
input_dir = "COG_enrichment"
files = glob.glob(os.path.join(input_dir, "*.csv"))
q_value_dict, rich_factor_dict = {}, {}

for file in files:
    sample = os.path.basename(file).replace("COG_enrichment_", "").replace(".csv", "")
    df = pd.read_csv(file)
    df = df[(df["overlap"] > 0) & df["significant"] & (df["rich_factor"] > 0.0002)]
    q_value_dict[sample] = df.set_index("COG_ID")["q_value"]
    rich_factor_dict[sample] = df.set_index("COG_ID")["rich_factor"]

q_value_df = pd.DataFrame(q_value_dict).fillna(1).reindex(columns=sample_order)
rich_factor_df = pd.DataFrame(rich_factor_dict).fillna(0).reindex(columns=sample_order)
significant_mask = (q_value_df < 0.05).any(axis=1)
q_value_df, rich_factor_df = q_value_df.loc[significant_mask], rich_factor_df.loc[significant_mask]

# Load COG annotations
cog_anno = pd.read_csv("cog-24.def.tab", sep="\t", header=None, usecols=[0, 2], on_bad_lines="skip", engine="python")
cog_anno.columns = ["COG_ID", "Annotation"]
cog_anno_dict = cog_anno.set_index("COG_ID")["Annotation"].to_dict()
q_value_df.index = q_value_df.index.map(lambda x: cog_anno_dict.get(x, x))
rich_factor_df.index = rich_factor_df.index.map(lambda x: cog_anno_dict.get(x, x))

abbreviation_dict = {
    "Translation elongation factor EF-Tu, a GTPase": "EF-Tu elongation factor",
    "Glyceraldehyde-3-phosphate dehydrogenase/erythrose-4-phosphate dehydrogenase": "GAPDH/E4PDH",
    "Imidazoleglycerol phosphate synthase glutamine amidotransferase subunit HisH": "HisH glutamine amidase",
    "Inorganic pyrophosphatase": "Inorganic PPase",
    "Ribosomal protein L10": "Ribosomal protein L10",
    "Glutamyl-tRNA reductase": "Glutamyl-tRNA reductase",
    "Na+/proline symporter": "Na+/proline symporter",
    "Superoxide dismutase": "Superoxide dismutase",
    "Signal transduction histidine kinase": "Histidine kinase",
    "RNA recognition motif (RRM) domain": "RRM RNA-binding domain",
    "DNA-binding response regulator, OmpR family, contains REC and winged-helix (wHTH) domain": "OmpR response regulator",
    "C-terminal processing protease CtpA/Prc, contains a PDZ domain": "CtpA/Prc protease",
    "Triacylglycerol esterase/lipase EstA, alpha/beta hydrolase fold": "EstA lipase",
    "Phosphomannomutase": "Phosphomannomutase",
    "Pseudouridylate synthase RsuA/RluF, specific for 16S rRNA U516, 23S rRNA U2604/U2605/U2457, and tRNA(Tyr)-35": "RsuA/RluF pseudouridylase",
    "Outer membrane receptor protein, Fe transport": "Fe OM receptor",
    "Epoxyqueuosine reductase QueH (tRNA modification)": "QueH tRNA modification",
    "Siroheme synthase (precorrin-2 oxidase/ferrochelatase domain)": "Siroheme synthase",
    "LD-carboxypeptidase LdcB, LAS superfamily": "LdcB carboxypeptidase",
    "Stereoselective (R,S)-S-adenosylmethionine hydrolase (adenosine-forming)": "SAM hydrolase",
    "K+-sensing histidine kinase KdpD": "KdpD K+ sensor kinase",
    "AraC-type DNA-binding domain and AraC-containing proteins": "AraC transcription regulator",
    "Mg/Co/Ni transporter MgtE (contains CBS domain)": "MgtE metal transporter",
    "Outer membrane protein OmpA and related peptidoglycan-associated (lipo)proteins": "OmpA outer membrane protein",
    "Cell division protein ZapB, interacts with FtsZ": "Cell division ZapB",
    "Beta-galactosidase/beta-glucuronidase": "Beta-gal/glucuronidase",
    "Periplasmic ligand-binding sensor domain": "Periplasmic sensor",
    "Putative alpha-1,2-mannosidase, GH92 family": "GH92 alpha-mannosidase",
    "Alanine-alpha-ketoisovalerate (or valine-pyruvate) aminotransferase": "Alanine aminotransferase",
    "Outer membrane cobalamin receptor protein BtuB": "BtuB cobalamin receptor",
    "TPR-like repeat domain": "TPR repeat domain",
    "Predicted phospholipase, patatin/cPLA2 family": "Patatin phospholipase",
    "Outer membrane receptor for ferrienterochelin and colicins": "Ferrienterochelin OM receptor",
    "Outer membrane receptor for ferric coprogen and ferric-rhodotorulic acid": "Ferric siderophore receptor",
    "Outer membrane receptor for monomeric catechols": "Catechol OM receptor",
    "Uncharacterized N-terminal domain of tricorn protease, contains WD40 repeats": "Tricorn protease WD40",
    "Transcriptional regulator GlxA, contains an amidase domain and an AraC-type DNA-binding HTH domain": "GlxA AraC regulator",
    "Predicted metal-binding protein, DUF2284 domain": "DUF2284 metal-binding protein",
    "Stage VI sporulation protein SpoVID (LysM domain), required for assembly of the spore coat": "SpoVID spore coat protein"
}

q_value_df.index = q_value_df.index.map(lambda x: abbreviation_dict.get(x, x))
rich_factor_df.index = rich_factor_df.index.map(lambda x: abbreviation_dict.get(x, x))

# Normalize rich factors
rich_factor_df_norm = rich_factor_df.apply(lambda x: (x - x.min()) / (x.max() - x.min()) if x.max() > x.min() else 0)

# Cluster rows
cg = sns.clustermap(rich_factor_df_norm, method="average", metric="euclidean", row_cluster=True, col_cluster=False, cmap="Reds", figsize=(1, 1))
row_order = cg.dendrogram_row.reordered_ind
plt.close()
rich_factor_df_norm, q_value_df = rich_factor_df_norm.iloc[row_order], q_value_df.iloc[row_order]
q_value_df = q_value_df[rich_factor_df_norm.columns]

# Plot heatmap
plt.figure(figsize=(1.2, 12))
ax = sns.heatmap(rich_factor_df_norm, cmap="Reds", linewidths=1, linecolor="grey", annot=False, cbar_kws={"label": "Scaled rich factor", "fraction": CBAR_WIDTH, "shrink": CBAR_LENGTH})

for i in range(rich_factor_df_norm.shape[0]):
    for j in range(rich_factor_df_norm.shape[1]):
        if q_value_df.iloc[i, j] < 0.05:
            ax.text(j + 0.5, i + 0.9, "*", ha="center", va="center", color="black", fontsize=30, fontweight="bold")

ax.set_xlabel("Species", fontsize=AXIS_LABEL_SIZE)
ax.set_ylabel("COG function", fontsize=AXIS_LABEL_SIZE)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=TICK_LABEL_SIZE)
ax.tick_params(axis="x", labelsize=TICK_LABEL_SIZE)
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=CBAR_TICK_SIZE)
cbar.set_label("Scaled rich factor", fontsize=CBAR_LABEL_SIZE)

plt.tight_layout()
plt.savefig("Figure5f.pdf", dpi=300, bbox_inches="tight")
plt.show()