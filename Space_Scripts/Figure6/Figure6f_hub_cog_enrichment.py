# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 16:50:33 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import os

# Set working directory to the current script location
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# Define project root
project_root = os.path.dirname(os.path.dirname(script_dir))

# Define input files and output directory
anno_file = os.path.join(project_root, "reference_genome_gene_annotation.tsv")
hub_file = os.path.join(script_dir, "hub_genes_top50.csv")
output_dir = os.path.join(script_dir, "COG_enrichment_by_module")
os.makedirs(output_dir, exist_ok=True)

# Load genome annotation and hub gene tables
df_anno = pd.read_csv(anno_file, sep="\t")
df_hub = pd.read_csv(hub_file)

# Standardize gene identifiers to match annotation format
df_hub["gene"] = df_hub["gene_name"].str.replace("-", "_", regex=False)

# Construct COG to gene mapping dictionary
cog_dict = {}
for _, row in df_anno.iterrows():
    gene_id = row["ID"]
    cog_ids = str(row["COG_ID"]).strip()

    if cog_ids == "-" or cog_ids.lower() == "nan":
        continue

    for cog in cog_ids.split(","):
        cog = cog.strip()
        if cog == "":
            continue
        cog_dict.setdefault(cog, []).append(gene_id)

# Determine total number of background genes
N = df_anno["ID"].nunique()

# Perform hypergeometric enrichment analysis for each module
modules = sorted(df_hub["module"].unique())

for mod in modules:

    gene_selected = set(df_hub[df_hub["module"] == mod]["gene"].astype(str).tolist())
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

        df_res = pd.DataFrame(results, columns=["COG_ID", "COG_size", "overlap", "p_value"])

        reject, pvals_corrected, _, _ = multipletests(df_res["p_value"], method="fdr_bh")
        df_res["q_value"] = pvals_corrected
        df_res["significant"] = reject

        df_res["rich_factor"] = df_res["overlap"] / df_res["COG_size"]

        df_res = df_res.sort_values("q_value")

        out_file = os.path.join(output_dir, f"COG_enrichment_module_{mod}.csv")
        df_res.to_csv(out_file, index=False)

        print(f"Module {mod} results saved: {out_file}")

# Print completion message after all modules are processed
print(f"All module enrichment results saved in folder: {output_dir}")