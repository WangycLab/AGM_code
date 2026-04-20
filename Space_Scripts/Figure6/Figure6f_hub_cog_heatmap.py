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

# Define plotting parameters and module order
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "Arial"
AXIS_LABEL_SIZE = 25
TICK_LABEL_SIZE = 18
CBAR_WIDTH = 0.1
CBAR_LENGTH = 1
CBAR_TICK_SIZE = 15
CBAR_LABEL_SIZE = 20
module_order = ["blue", "brown", "green", "turquoise", "yellow"]

# Load module enrichment result files
input_dir = "COG_enrichment_by_module"
files = glob.glob(os.path.join(input_dir, "COG_enrichment_module_*.csv"))

q_value_dict = {}
rich_factor_dict = {}

# Extract q-value and rich factor information from each module result
for file in files:
    module_name = os.path.basename(file).replace("COG_enrichment_module_", "").replace(".csv", "")
    df = pd.read_csv(file)

    df = df[(df["overlap"] > 0) & df["significant"] & (df["rich_factor"] >= 0.0006)]

    q_value_dict[module_name] = df.set_index("COG_ID")["q_value"]
    rich_factor_dict[module_name] = df.set_index("COG_ID")["rich_factor"]

# Merge module matrices into unified dataframes
q_value_df = pd.DataFrame(q_value_dict).fillna(1)
rich_factor_df = pd.DataFrame(rich_factor_dict).fillna(0)

rich_factor_df = rich_factor_df.reindex(columns=module_order)
q_value_df = q_value_df.reindex(columns=module_order)

# Keep COG rows with at least one significant enrichment result
sig_mask = (q_value_df < 0.05).any(axis=1)
q_value_df = q_value_df.loc[sig_mask]
rich_factor_df = rich_factor_df.loc[sig_mask]

# Load COG functional annotation table
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

# Replace COG identifiers with functional descriptions
rich_factor_df.index = rich_factor_df.index.map(lambda x: cog_anno_dict.get(x, x))
q_value_df.index = q_value_df.index.map(lambda x: cog_anno_dict.get(x, x))

# Abbreviation dictionary for y-axis labels
abbreviation_dict = {
    'Small heat shock protein IbpA, HSP20 family': 'Small HSP IbpA',
    'DNA-directed RNA polymerase, alpha subunit/40 kD subunit': 'RNA polymerase alpha',
    '6-phosphogluconate dehydrogenase': '6PG dehydrogenase',
    'Glucose-6-phosphate 1-dehydrogenase': 'G6P dehydrogenase',
    'Glycine cleavage system protein P (pyridoxal-binding), N-terminal domain': 'Glycine cleavage P N-term',
    'DNA polymerase B elongation subunit': 'DNA polymerase B',
    'Chaperonin GroEL (HSP60 family)': 'Chaperonin GroEL',
    'Isocitrate dehydrogenase': 'Isocitrate dehydrogenase',
    'Phosphoglycerate mutase (BPG-dependent)': 'Phosphoglycerate mutase',
    'Superoxide dismutase': 'Superoxide dismutase',
    'RNA recognition motif (RRM) domain': 'RRM domain protein',
    'Biopolymer transport protein ExbD': 'ExbD transporter',
    'Glycine cleavage system protein P (pyridoxal-binding), C-terminal domain': 'Glycine cleavage P C-term',
    'Superfamily I DNA and/or RNA helicase': 'SF1 helicase',
    'Cobalamin biosynthesis protein CobN, Mg-chelatase': 'CobN cobalamin synthase',
    'Redox-sensitive bicupin YhaK, pirin superfamily': 'Redox bicupin YhaK',
    'Fructose-bisphosphate aldolase class Ia, DhnA family': 'FBP aldolase',
    'Methylmalonyl-CoA mutase, N-terminal domain/subunit': 'MCM N-terminal',
    'Glucuronate isomerase': 'Glucuronate isomerase',
    'Peroxiredoxin': 'Peroxiredoxin',
    'L-arabinose isomerase': 'L-arabinose isomerase',
    'Methylmalonyl-CoA mutase, C-terminal domain/subunit (cobalamin-binding)': 'MCM C-terminal',
    'Peroxiredoxin family protein': 'Peroxiredoxin family',
    'Ferritin-like DNA-binding protein, DPS (DNA Protection under Starvation) family': 'Ferritin-like DPS',
    'Predicted kinase related to galactokinase and mevalonate kinase': 'Galactokinase-like kinase',
    'Archaeal/vacuolar-type H+-ATPase subunit H': 'V-type ATPase subunit H',
    'Uncharacterized lipoprotein NlpE involved in copper resistance': 'Lipoprotein NlpE',
    'Putative hemolysin': 'Hemolysin',
    'Outer membrane protein OmpW': 'Outer membrane OmpW',
    'Periplasmic ligand-binding sensor domain': 'Periplasmic sensor',
    'Transposase': 'Transposase',
    'Alpha-L-arabinofuranosidase': 'Alpha-L-arabinofuranosidase',
    'FRM2-like oxidoreductase, nitroreductase family': 'FRM2 oxidoreductase',
    'Transposase, IS1182 family': 'IS1182 transposase',
    'Phage tail tape-measure protein, controls tail length': 'Phage tape-measure protein',
    'ATP-responsive regulator of NAD biosynthesis': 'ATP-regulated NAD regulator',
    'Toxin component of the Txe-Axe toxin-antitoxin module, Txe/YoeB family': 'Txe/YoeB toxin',
    'Uncharacterized membrane permease YidK, sodium:solute symporter family': 'Membrane permease YidK',
    'Mu-like prophage I protein': 'Mu-like prophage protein',
    'Adenine-specific DNA methylase, N12 class': 'DNA adenine methylase',
    'Outer membrane receptor for ferric coprogen and ferric-rhodotorulic acid': 'Ferric siderophore receptor',
    'L-rhamnose isomerase': 'L-rhamnose isomerase',
    'Phage-related tail fiber protein': 'Phage tail fiber',
    'Preprotein translocase subunit Sec63': 'Translocase Sec63',
    'Predicted small metal-binding protein': 'Metal-binding protein',
    'Protein involved in initiation of plasmid replication': 'Plasmid replication initiator',
    'Self-loading helicase or an inactivated derivative': 'Self-loading helicase'
}

# Apply abbreviation dictionary to y-axis labels
rich_factor_df.index = rich_factor_df.index.map(lambda x: abbreviation_dict.get(x, x))
q_value_df.index = q_value_df.index.map(lambda x: abbreviation_dict.get(x, x))

# Normalize each column 0-1
rich_factor_df_norm = rich_factor_df.copy()
for col in rich_factor_df_norm.columns:
    col_min = rich_factor_df_norm[col].min()
    col_max = rich_factor_df_norm[col].max()
    if col_max > col_min:
        rich_factor_df_norm[col] = (rich_factor_df_norm[col] - col_min) / (col_max - col_min)
    else:
        rich_factor_df_norm[col] = 0

# Cluster rows to determine order
cg = sns.clustermap(rich_factor_df_norm, method="average", metric="euclidean", row_cluster=True, col_cluster=False, cmap="Reds", figsize=(1,1))
row_order = cg.dendrogram_row.reordered_ind
plt.close()
rich_factor_df_norm = rich_factor_df_norm.iloc[row_order, :]
q_value_df = q_value_df.iloc[row_order, :]

# Plot heatmap
plt.figure(figsize=(2, 15))
ax = sns.heatmap(
    rich_factor_df_norm,
    cmap="Reds",
    linewidths=1,
    linecolor="grey",
    annot=False,
    cbar_kws={"label": "Rich factor", "fraction": CBAR_WIDTH, "shrink": CBAR_LENGTH}
)

# Add significance stars
for i in range(rich_factor_df_norm.shape[0]):
    for j in range(rich_factor_df_norm.shape[1]):
        if q_value_df.iloc[i, j] < 0.05:
            ax.text(
                j + 0.5,
                i + 0.95,
                "*",
                ha="center",
                va="center",
                color="black",
                fontsize=30,
                fontweight="bold"
            )

# Set axis labels and ticks
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
plt.savefig("Figure6f.pdf", dpi=300, bbox_inches="tight")
plt.show()