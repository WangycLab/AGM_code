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
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
project_root = os.path.dirname(os.path.dirname(script_dir))
anno_file = os.path.join(project_root, "reference_genome_gene_annotation.tsv")
hub_file = os.path.join(script_dir, "hub_genes_top50.csv")
output_dir = os.path.join(script_dir, "COG_enrichment_by_module")
os.makedirs(output_dir, exist_ok=True)

# Load and standardize gene identifiers
df_anno = pd.read_csv(anno_file, sep="\t")
df_hub = pd.read_csv(hub_file)
df_hub["gene"] = df_hub["gene_name"].str.replace("-", "_", regex=False)

# Build COG-to-gene mapping
cog_dict = {}
for _, row in df_anno.iterrows():
    gene_id, cog_ids = row["ID"], str(row["COG_ID"]).strip()
    if cog_ids == "-" or cog_ids.lower() == "nan":
        continue
    for cog in cog_ids.split(","):
        cog = cog.strip()
        if cog:
            cog_dict.setdefault(cog, []).append(gene_id)

N = df_anno["ID"].nunique()

# Perform enrichment analysis for each module
for mod in sorted(df_hub["module"].unique()):
    gene_selected = set(df_hub.loc[df_hub["module"] == mod, "gene"].astype(str))
    K, results = len(gene_selected), []

    for cog, genes in cog_dict.items():
        genes_set = set(genes)
        n, k = len(genes_set), len(genes_set & gene_selected)
        if k == 0:
            continue
        pval = hypergeom.sf(k - 1, N, n, K)
        results.append([cog, n, k, pval])

    if not results:
        continue

    df_res = pd.DataFrame(results, columns=["COG_ID", "COG_size", "overlap", "p_value"])
    df_res["q_value"] = multipletests(df_res["p_value"], method="fdr_bh")[1]
    df_res["significant"] = df_res["q_value"] <= 0.05
    df_res["rich_factor"] = df_res["overlap"] / df_res["COG_size"]
    df_res = df_res.sort_values("q_value")

    out_file = os.path.join(output_dir, f"COG_enrichment_module_{mod}.csv")
    df_res.to_csv(out_file, index=False)
    print(f"Module {mod} results saved.")

print(f"All module enrichment results saved in: {output_dir}")