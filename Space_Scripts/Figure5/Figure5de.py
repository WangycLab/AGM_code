# -*- coding: utf-8 -*-
"""
Created on Wed Jan 7 19:01:32 2026
@author: ZHENG XINGHAI
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from upsetplot import UpSet, from_indicators

# Plot settings
matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "Arial"})

input_dir = "for_upset"

# Load gene sets for each species
species_gene_dict = {}
for file in os.listdir(input_dir):
    if not file.endswith(".csv"): continue
    species_name = os.path.splitext(file)[0]
    file_path = os.path.join(input_dir, file)
    df = pd.read_csv(file_path, dtype=str)
    if "gene" not in df.columns:
        print(f"Warning: {file} does not contain gene column")
        continue
    species_gene_dict[species_name] = set(df["gene"].dropna().astype(str).unique())

all_genes = sorted(set().union(*species_gene_dict.values()))
gene_matrix = pd.DataFrame({species: pd.Index(all_genes).isin(genes) for species, genes in species_gene_dict.items()}, index=all_genes)
gene_matrix.to_csv(f"Figure5f.csv")

# Plot UpSet for genes shared by at least two species
gene_matrix_filtered = gene_matrix[gene_matrix.sum(axis=1) >= 2]
upset_data = from_indicators(gene_matrix_filtered.columns, gene_matrix_filtered)
plt.figure(figsize=(14, 10))
UpSet(upset_data, subset_size="count", show_counts=True, sort_by="cardinality", facecolor="steelblue", min_subset_size=20).plot()
plt.savefig(f"Figure5d.pdf", bbox_inches="tight", dpi=300)
plt.close()

# Export gene intersections
intersection_df = pd.DataFrame([{"gene": gene, "species_number": row.sum(), "species": ";".join(row.index[row.astype(bool)])} for gene, row in gene_matrix.iterrows()])
intersection_df.to_csv(f"gene_intersections.csv", index=False)

# Define clusters for each species
species_cluster_dict = {
    "Bacteroides_ovatus": ["0", "2", "3"],
    "Bacteroides_thetaiotaomicron": ["0", "1", "2"],
    "Phocaeicola_dorei": ["2", "3", "4"],
    "Prevotella_stercorea": ["1", "2", "3"],
    "Parabacteroides_distasonis": ["5"],
    "Phocaeicola_vulgatus": ["1", "3"]
}

# Extract genes from selected clusters
selected_cluster_gene_dict = {}
for species_name, cluster_list in species_cluster_dict.items():
    file_path = os.path.join(input_dir, f"{species_name}.csv")
    if not os.path.exists(file_path):
        print(f"Warning: {file_path} not found")
        continue
    df = pd.read_csv(file_path, dtype=str)
    if "gene" not in df.columns or "cluster" not in df.columns:
        print(f"Warning: {species_name} missing gene or cluster column")
        continue
    selected_genes = set(df.loc[df["cluster"].astype(str).isin(cluster_list), "gene"].dropna().astype(str).unique())
    selected_cluster_gene_dict[species_name] = selected_genes

all_selected_genes = sorted(set().union(*selected_cluster_gene_dict.values()))
selected_matrix = pd.DataFrame({species: pd.Index(all_selected_genes).isin(genes) for species, genes in selected_cluster_gene_dict.items()}, index=all_selected_genes)
selected_matrix.to_csv(f"Figure5g.csv")

# Plot UpSet for genes shared by at least two species
selected_matrix_filtered = selected_matrix[selected_matrix.sum(axis=1) >= 2]
selected_upset_data = from_indicators(selected_matrix_filtered.columns, selected_matrix_filtered)
plt.figure(figsize=(14, 10))
UpSet(selected_upset_data, subset_size="count", show_counts=True, sort_by="cardinality", facecolor="darkorange", min_subset_size=10).plot()
plt.savefig(f"Figure5e.pdf", bbox_inches="tight", dpi=300)
plt.close()

print("Finished!")