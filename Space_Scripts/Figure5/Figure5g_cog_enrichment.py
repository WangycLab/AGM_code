# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 16:50:33 2025
@author: ZHENG XINGHAI
"""

import os
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

# Set paths
current_dir = os.path.dirname(os.path.abspath(__file__))
desktop_dir = os.path.dirname(os.path.dirname(current_dir))
anno_file = os.path.join(desktop_dir, "reference_genome_gene_annotation.tsv")
gene_file = os.path.join(current_dir, "Figure5g_annotated.csv")

# Create output folder
output_dir = os.path.join(current_dir, "COG_enrichment_FD30")
os.makedirs(output_dir, exist_ok=True)

# Load and standardize gene identifiers
df_anno = pd.read_csv(anno_file, sep="\t", index_col="ID")
df_gene = pd.read_csv(gene_file, index_col=0)
df_anno.index = df_anno.index.astype(str).str.replace("-", "_", regex=False)
df_gene.index = df_gene.index.astype(str).str.replace("-", "_", regex=False)

species_cols = [
    "Bacteroides_ovatus", "Bacteroides_thetaiotaomicron",
    "Parabacteroides_distasonis", "Phocaeicola_dorei",
    "Phocaeicola_vulgatus", "Prevotella_stercorea"
]
df_gene[species_cols] = df_gene[species_cols].map(lambda x: str(x).upper() == "TRUE")

# Define gene sets
type1_genes = df_gene.loc[
    df_gene["Phocaeicola_vulgatus"] & df_gene["Phocaeicola_dorei"] &
    ~df_gene[["Bacteroides_ovatus", "Bacteroides_thetaiotaomicron",
              "Parabacteroides_distasonis", "Prevotella_stercorea"]].any(axis=1)
].index

type2_genes = df_gene.loc[
    df_gene["Bacteroides_ovatus"] & df_gene["Bacteroides_thetaiotaomicron"] &
    ~df_gene[["Parabacteroides_distasonis", "Phocaeicola_dorei",
              "Phocaeicola_vulgatus", "Prevotella_stercorea"]].any(axis=1)
].index

gene_sets = {
    "type1_Phocaeicola": set(type1_genes),
    "type2_Bacteroides": set(type2_genes)
}
print("Gene numbers:", {k: len(v) for k, v in gene_sets.items()})

# Build COG dictionary
cog_dict = {}
for gene_id, row in df_anno.iterrows():
    cog_ids = str(row["COG_ID"]).strip()
    if cog_ids == "-" or cog_ids.lower() == "nan":
        continue
    for cog in cog_ids.split(","):
        cog = cog.strip()
        if cog:
            cog_dict.setdefault(cog, []).append(gene_id)

N = len(df_anno)

# COG enrichment
def COG_enrichment(gene_selected, output_file):
    K, results = len(gene_selected), []

    for cog, genes in cog_dict.items():
        genes_set = set(genes)
        n, k = len(genes_set), len(genes_set & gene_selected)
        if k == 0:
            continue
        results.append([cog, n, k, hypergeom(N, n, K).sf(k - 1)])

    if not results:
        return

    df_res = pd.DataFrame(results, columns=["COG_ID", "COG_size", "overlap", "p_value"])
    df_res["q_value"] = multipletests(df_res["p_value"], method="fdr_bh")[1]
    df_res["significant"] = df_res["q_value"] <= 0.05
    df_res["rich_factor"] = df_res["overlap"] / df_res["COG_size"]
    df_res.sort_values("q_value").to_csv(output_file, index=False)

# Run enrichment and save results
output_dict = {
    "type1_Phocaeicola": os.path.join(output_dir, "COG_enrichment_Phocaeicola.csv"),
    "type2_Bacteroides": os.path.join(output_dir, "COG_enrichment_Bacteroides.csv")
}

for name, genes in gene_sets.items():
    COG_enrichment(genes, output_dict[name])

print(f"All COG enrichment completed. Results saved to: {output_dir}")
