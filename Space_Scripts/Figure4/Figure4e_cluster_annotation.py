# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 14:43:00 2025
@author: ZHENG XINGHAI
"""

import os
import pandas as pd

# Determine the current working directory where the script is executed
current_dir = os.getcwd()

# Derive the project root by removing the script folder path
project_root = current_dir.replace(os.path.join("Space_Scripts", "Figure4"), "")

# Define input and output file paths
deg_file = os.path.join(current_dir, "all_clusters_DEGs.tsv")
anno_file = os.path.join(project_root, "reference_genome_gene_annotation.tsv")

# Load the differential expression gene table
deg_df = pd.read_csv(deg_file, sep="\t")

# Standardize gene identifiers by replacing hyphens with underscores
deg_df["gene"] = deg_df["gene"].astype(str).str.replace("-", "_", regex=False)

# Load the reference genome annotation table
anno_df = pd.read_csv(anno_file, sep="\t")

# Rename the annotation identifier column to match the DEG table
anno_df = anno_df.rename(columns={"ID": "gene"})

# Merge DEG results with genome annotations based on gene identifiers
merged_df = pd.merge(deg_df, anno_df, on="gene", how="left")

# Save the annotated DEG table to the working directory
output_file = os.path.join(current_dir, "all_clusters_DEGs_annotated.csv")
merged_df.to_csv(output_file, index=False)

# Print a completion message with the output file location
print(f"Annotation completed. Output saved to {output_file}")