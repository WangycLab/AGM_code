# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 14:43:00 2025
@author: ZHENG XINGHAI
"""

import os
import pandas as pd

# Set paths
current_dir = os.getcwd()
project_root = current_dir.replace(os.path.join("Space_Scripts", "Figure4"), "")
deg_file = os.path.join(current_dir, "all_clusters_DEGs.tsv")
anno_file = os.path.join(project_root, "reference_genome_gene_annotation.tsv")
output_file = os.path.join(current_dir, "all_clusters_DEGs_annotated.csv")

# Load and standardize data
deg_df = pd.read_csv(deg_file, sep="\t")
deg_df["gene"] = deg_df["gene"].astype(str).str.replace("-", "_", regex=False)
anno_df = pd.read_csv(anno_file, sep="\t").rename(columns={"ID": "gene"})

# Merge and save
merged_df = pd.merge(deg_df, anno_df, on="gene", how="left")
merged_df.to_csv(output_file, index=False)
print(f"Annotation completed. Output saved to {output_file}")