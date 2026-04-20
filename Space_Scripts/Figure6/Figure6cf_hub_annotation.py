# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 14:43:00 2025
@author: ZHENG XINGHAI
"""

import os
import pandas as pd

# Set working directory to the directory where the script is located
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# Determine project root by moving two levels up
project_root = os.path.dirname(os.path.dirname(script_dir))

# Load hub gene table and normalize gene identifiers
hub_file = os.path.join(script_dir, "hub_genes_top50.csv")
hub_df = pd.read_csv(hub_file)
hub_df["gene"] = hub_df["gene_name"].str.replace("-", "_")

# Load reference genome gene annotation from the project root directory
anno_file = os.path.join(project_root, "reference_genome_gene_annotation.tsv")
anno_df = pd.read_csv(anno_file, sep="\t")
anno_df = anno_df.rename(columns={"ID": "gene"})

# Merge hub gene information with genome annotation
merged_df = pd.merge(hub_df, anno_df, on="gene", how="left")

# Save annotated hub gene table to the current directory
output_file = os.path.join(script_dir, "hub_genes_top50_anno.csv")
merged_df.to_csv(output_file, index=False)