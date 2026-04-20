# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 14:43:00 2025
@author: ZHENG XINGHAI
"""

import os
import pandas as pd

# Determine the project root directory
current_dir = os.getcwd()
project_root = current_dir.replace(os.path.join("Space_Scripts", "Figure5"), "")

# Load gene cluster assignment table and set gene as index
deg_file = os.path.join(current_dir, "gene_clusters.csv")
deg_df = pd.read_csv(deg_file, sep=",")
deg_df.set_index("gene", inplace=True)

# Standardize gene identifiers by replacing hyphens with underscores
deg_df.index = deg_df.index.str.replace("-", "_")

# Load reference genome annotation table from the project root directory
anno_file = os.path.join(project_root, "reference_genome_gene_annotation.tsv")
anno_df = pd.read_csv(anno_file, sep="\t")
anno_df.set_index("ID", inplace=True)
anno_df.index = anno_df.index.str.replace("-", "_")

# Merge cluster assignments with gene annotation using gene identifiers
merged_df = deg_df.join(anno_df, how="left")

# Save the annotated gene cluster table
output_file = os.path.join(current_dir, "gene_clusters_annotated.csv")
merged_df.to_csv(output_file)

# Print completion message after annotation
print(f"Annotation completed. Output saved to {output_file}")