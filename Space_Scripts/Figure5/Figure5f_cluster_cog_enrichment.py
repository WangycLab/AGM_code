# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 16:50:33 2025
@author: ZHENG XINGHAI
"""

import os
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

# Determine the project root directory
current_dir = os.getcwd()
project_root = current_dir.replace(os.path.join("Space_Scripts", "Figure5"), "")

# Define input files and output directory
anno_file = os.path.join(project_root, "reference_genome_gene_annotation.tsv")
deg_file = os.path.join(current_dir, "gene_clusters_annotated.csv")
output_dir = os.path.join(current_dir, "COG_enrichment_by_cluster")
os.makedirs(output_dir, exist_ok=True)

# Load annotation table and clustered gene table
df_anno = pd.read_csv(anno_file, sep="\t", index_col="ID")
df_deg = pd.read_csv(deg_file, sep=",", index_col="gene")

# Standardize gene identifiers by replacing hyphens with underscores
df_deg.index = df_deg.index.str.replace("-", "_", regex=False)
df_anno.index = df_anno.index.str.replace("-", "_", regex=False)

# Filter genes by adjusted p value threshold when available
if "p_val_adj" in df_deg.columns:
    df_deg = df_deg[df_deg["p_val_adj"] <= 0.05].copy()

# Build a mapping dictionary from COG identifiers to gene lists
cog_dict = {}
for gene_id, row in df_anno.iterrows():
    cog_ids = str(row["COG_ID"]).strip()
    if cog_ids == "-" or cog_ids.lower() == "nan":
        continue

    for cog in cog_ids.split(","):
        cog = cog.strip()
        if cog == "":
            continue
        cog_dict.setdefault(cog, []).append(gene_id)

# Define the background gene universe size
N = df_anno.shape[0]

# Perform hypergeometric enrichment analysis for each gene cluster
clusters = sorted(df_deg["cluster"].unique())

for cl in clusters:

    gene_selected = set(df_deg[df_deg["cluster"] == cl].index)
    K = len(gene_selected)

    results = []

    for cog, genes in cog_dict.items():
        genes_set = set(genes)

        n = len(genes_set)
        k = len(genes_set & gene_selected)

        if k == 0:
            continue

        rv = hypergeom(N, n, K)
        pval = rv.sf(k - 1)

        results.append([cog, n, k, pval])

    if results:

        df_res = pd.DataFrame(
            results,
            columns=["COG_ID", "COG_size", "overlap", "p_value"]
        )

        reject, pvals_corrected, _, _ = multipletests(
            df_res["p_value"],
            method="fdr_bh"
        )

        df_res["q_value"] = pvals_corrected
        df_res["significant"] = reject
        df_res["rich_factor"] = df_res["overlap"] / df_res["COG_size"]

        out_file = os.path.join(
            output_dir,
            f"COG_enrichment_cluster_{cl}.csv"
        )

        df_res.to_csv(out_file, index=False)

        print(f"Cluster {cl} results saved to {out_file}")

# Print completion message after all enrichment results are generated
print(f"All results have been saved to directory {output_dir}")