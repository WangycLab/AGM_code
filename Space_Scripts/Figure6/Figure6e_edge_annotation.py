# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 18:58:08 2026
@author: ZHENG XINGHAI
"""

import pandas as pd

# Load input files
edges = pd.read_csv("combined_WGCNA_edges.csv", sep=",")
nodes = pd.read_csv("combined_WGCNA_nodes.csv", sep=",")
hub_df = pd.read_csv("hub_genes_top10.csv", sep=",")

# Build gene-to-module mapping
gene_to_module = dict(zip(nodes["gene_name"], nodes["module"]))

# Assign module info to edges
edges["Gene1_module"] = edges["Gene1"].map(gene_to_module)
edges["Gene2_module"] = edges["Gene2"].map(gene_to_module)

# Remove edges involving grey module
edges = edges[(edges["Gene1_module"] != "grey") & (edges["Gene2_module"] != "grey")].copy()

# Define hub gene set
hub_genes = set(hub_df["gene_name"].astype(str))

# Annotate hub edges (Gene1 only)
edges["isHub"] = edges["Gene1"].astype(str).isin(hub_genes)
edges["isHub"] = edges["isHub"].map({True: "TRUE", False: "FALSE"})

# Normalize edge weight (Z-score + Min–Max)
w_mean = edges["Weight"].mean()
w_std = edges["Weight"].std(ddof=0)
edges["Weight_z"] = (edges["Weight"] - w_mean) / w_std
z_min = edges["Weight_z"].min()
z_max = edges["Weight_z"].max()
edges["Weight_norm"] = (edges["Weight_z"] - z_min) / (z_max - z_min)

# Export annotated edges
edges.to_csv("combined_WGCNA_edges_anno.csv", index=False)