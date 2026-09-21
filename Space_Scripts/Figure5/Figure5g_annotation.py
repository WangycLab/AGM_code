# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 14:43:00 2025
@author: ZHENG XINGHAI
"""

import os
import pandas as pd

# Set paths
current_dir = os.path.dirname(os.path.abspath(__file__))
desktop_dir = os.path.dirname(os.path.dirname(current_dir))

# Load and standardize gene identifiers
deg_file = os.path.join(current_dir, "Figure5g.csv")
deg_df = pd.read_csv(deg_file, index_col=0)
deg_df.index = deg_df.index.astype(str).str.replace("-", "_", regex=False)

anno_file = os.path.join(desktop_dir, "reference_genome_gene_annotation.tsv")
anno_df = pd.read_csv(anno_file, sep="\t", index_col="ID")
anno_df.index = anno_df.index.astype(str).str.replace("-", "_", regex=False)

# Merge gene matrix with annotations
merged_df = deg_df.join(anno_df, how="left")

# Save annotated table
output_file = os.path.join(current_dir, "Figure5g_annotated.csv")
merged_df.to_csv(output_file)

print(f"Annotation completed: {output_file}")