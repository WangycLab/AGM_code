# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 14:43:00 2025
@author: ZHENG XINGHAI
"""

import os
import pandas as pd

# Set working directory and project root
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
project_root = os.path.dirname(os.path.dirname(script_dir))

# Load and standardize gene identifiers
hub_df = pd.read_csv(os.path.join(script_dir, "hub_genes_top50.csv"))
hub_df["gene"] = hub_df["gene_name"].str.replace("-", "_", regex=False)
anno_df = pd.read_csv(os.path.join(project_root, "reference_genome_gene_annotation.tsv"), sep="\t").rename(columns={"ID": "gene"})

# Merge and save
merged_df = pd.merge(hub_df, anno_df, on="gene", how="left")
merged_df.to_csv(os.path.join(script_dir, "hub_genes_top50_anno.csv"), index=False)