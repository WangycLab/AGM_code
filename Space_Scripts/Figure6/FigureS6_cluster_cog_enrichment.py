# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 16:50:33 2025
@author: ZHENG XINGHAI
"""

import os
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

# Define paths
current_dir = os.getcwd()
project_root = current_dir.replace(os.path.join("Space_Scripts", "Figure6"), "")
anno_file = os.path.join(project_root, "reference_genome_gene_annotation.tsv")
deg_file = os.path.join(current_dir, "all_clusters_DEGs_annotated.csv")
output_dir = os.path.join(current_dir, "COG_enrichment_by_cluster")
filtered_deg_file = os.path.join(current_dir, "all_clusters_DEGs_filted.csv")
os.makedirs(output_dir, exist_ok=True)

# Load and filter annotations
df_anno = pd.read_csv(anno_file, sep="\t", index_col="ID")
df_anno.index = df_anno.index.str.replace("-", "_", regex=False)

df_deg = pd.read_csv(deg_file)
df_deg["gene"] = df_deg["gene"].str.replace("-", "_", regex=False)
if "p_val_adj" in df_deg.columns:
    df_deg = df_deg[df_deg["p_val_adj"] <= 0.05].copy()

df_deg = df_deg[df_deg["COG_ID"].notna() & (df_deg["COG_ID"] != "-")].copy()

# Retain top 10 genes per cluster
df_deg = (
    df_deg.sort_values(["cluster", "cells_expressing", "avg_log2FC"], ascending=[True, False, False])
    .groupby("cluster", group_keys=False)
    .head(10)
    .set_index("gene")
)

# Build COG-to-gene mapping
cog_dict = {}
for gene_id, cog_ids in df_anno["COG_ID"].items():
    for cog in str(cog_ids).split(","):
        cog = cog.strip()
        if cog and cog != "-" and cog.lower() != "nan":
            cog_dict.setdefault(cog, []).append(gene_id)

# Hypergeometric enrichment analysis
N = df_anno.shape[0]
clusters = sorted(df_deg["cluster"].unique())

for cl in clusters:
    gene_selected = set(df_deg.loc[df_deg["cluster"] == cl].index)
    K = len(gene_selected)
    results = []

    for cog, genes in cog_dict.items():
        genes_set = set(genes)
        n, k = len(genes_set), len(genes_set & gene_selected)
        if k > 0:
            results.append([cog, n, k, hypergeom.sf(k - 1, N, n, K)])

    if results:
        df_res = pd.DataFrame(results, columns=["COG_ID", "COG_size", "overlap", "p_value"])
        df_res["q_value"] = multipletests(df_res["p_value"], method="fdr_bh")[1]
        df_res["significant"] = df_res["q_value"] <= 0.05
        df_res["rich_factor"] = df_res["overlap"] / df_res["COG_size"]
        df_res.to_csv(os.path.join(output_dir, f"COG_enrichment_cluster_{cl}.csv"), index=False)

# Save filtered DEG table
df_deg.reset_index().to_csv(filtered_deg_file, index=False)

print(f"COG enrichment completed: {len(clusters)} clusters, {len(df_deg)} genes retained.")