# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 16:22:45 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import jensenshannon
from kneed import KneeLocator
import glob
import os
import seaborn as sns
from upsetplot import UpSet, from_indicators

# Ensure PDFs are AI-editable by setting font type
plt.rcParams["pdf.fonttype"] = 42

# Define a function to remove initial decreasing part of a similarity curve
def remove_decreasing_part(data, associated_data):
    for i in range(1, len(data)):
        if data[i] > data[i - 1]:
            return data[i:], associated_data[i:]
    return data, associated_data

# Initialize lists to store species compositions and sample names
all_compositions = []
sample_names = []

# Find all filtered taxonomy CSV files for processing
file_list = glob.glob('store_species_filted/*.csv')
print(file_list)
os.makedirs("jsd_curve", exist_ok=True)

# Set default plotting style
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 15

# Load summary table of cell counts
summary_path = "summary_cell_counts.csv"
summary_df = pd.read_csv(summary_path)

# Dictionary to store final selected cell numbers based on elbow points
final_selected_dict = {}

# Calculate Jensen-Shannon similarity curves and identify elbow points
for file in file_list:
    df = pd.read_csv(file)
    df = df.sort_values(by="new_est_reads", ascending=False).reset_index(drop=True)

    similarities = []
    cutoff0 = int(np.ceil(len(df) * 0.01))
    df0 = df.iloc[:cutoff0]
    dataframe0 = df0.groupby("name")["new_est_reads"].sum() / df0["new_est_reads"].sum()
    prev_df = dataframe0

    for pct in range(5, 101, 5):
        cutoff = int(np.ceil(len(df) * pct / 100))
        dfN = df.iloc[:cutoff]
        dataframeN = dfN.groupby("name")["new_est_reads"].sum()
        dataframeN = dataframeN / dataframeN.sum()
        aligned = pd.concat([prev_df, dataframeN], axis=1).fillna(0)
        p, q = aligned.iloc[:, 0].values, aligned.iloc[:, 1].values
        sim = 1 - jensenshannon(p, q)
        similarities.append((pct, sim))
        prev_df = dataframeN

    similarity_df = pd.DataFrame(similarities, columns=["percent", "similarity"])
    x = similarity_df["percent"].tolist()
    y = similarity_df["similarity"].tolist()
    y_filted, x_filted = remove_decreasing_part(y, x)
    knee = KneeLocator(x_filted, y_filted, curve="concave", direction="increasing")
    elbow_pct = knee.knee

    cutoff_elbow = int(np.ceil(len(df) * elbow_pct / 100))
    sample_name = os.path.basename(file).replace("_sc_taxonomy_filted.csv", "")
    final_selected_dict[sample_name] = cutoff_elbow

    # Extract species composition at elbow cutoff
    df_elbow = df.iloc[:cutoff_elbow]
    composition_df = df_elbow.groupby("name")["new_est_reads"].sum()
    composition_df = composition_df / composition_df.sum()
    composition_df.name = sample_name
    all_compositions.append(composition_df)
    sample_names.append(sample_name)

    # Plot and save the JSD similarity curve for this sample
    plt.figure(figsize=(8, 5))
    plt.plot(x, y, marker='o', markersize=5, label="Similarity")
    plt.axvline(elbow_pct, color='r', linestyle='--', label=f"Elbow: {elbow_pct}%")
    plt.xlabel("Top N% of Cells", fontsize=20)
    plt.ylabel("Similarity with Previous Step", fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.title(f"{sample_name}", fontsize=25)
    plt.legend(fontsize=15)
    plt.tight_layout()
    plt.savefig(f"jsd_curve/{sample_name}_jsd_curve.pdf")
    plt.show()
    plt.close()

# Update summary table with final selected cell counts
summary_df["Final_selected"] = summary_df["Sample"].map(final_selected_dict)
output_summary_path = "summary_cell_counts_completed.csv"
summary_df.to_csv(output_summary_path, index=False)
print(f"Updated summary saved to {output_summary_path}")

# Define desired order of samples
desired_order = ["Con", "2mo", "4mo", "6mo"]

# Merge species compositions into a single DataFrame
super_composition_df = pd.concat(all_compositions, axis=1).fillna(0)
existing_order = [col for col in desired_order if col in super_composition_df.columns]
super_composition_df = super_composition_df[existing_order]
super_composition_df.to_csv("species_abundance_profile.csv")

# Compute Pearson correlation between samples
data = super_composition_df.T
corr_df = data.T.corr(method="pearson")

# Plot a heatmap of sample similarity based on species abundance
plt.figure(figsize=(5, 4))
ax = sns.heatmap(
    corr_df.astype(float),
    cmap="viridis_r",
    vmin=0,
    vmax=1,
    annot=True,
    fmt=".2f",
    square=True,
    linewidths=0.6,
    cbar_kws={"label": "Pearson correlation"},
    annot_kws={"size": 15}
)
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=15)
cbar.ax.yaxis.label.set_size(16)
plt.xticks(rotation=0, ha="center", fontsize=15)
plt.yticks(rotation=0, ha="right", fontsize=15)
plt.title("Species abundance-based similarity", fontsize=16, pad=15)
plt.tight_layout()
plt.savefig("Figure2f.pdf", dpi=300)
plt.show()

# Prepare binary presence/absence matrix for UpSet plot
data_filtered = super_composition_df.copy()
data_filtered = data_filtered.applymap(lambda x: x if x > 0 else 0)
binary_df = (data_filtered > 0)

# Create and plot UpSet diagram of shared species
reversed_columns = list(reversed(binary_df.columns.tolist()))
upset_data = from_indicators(reversed_columns, binary_df)
plt.figure(figsize=(12, 8))
upset = UpSet(
    upset_data,
    show_counts=True,
    min_subset_size=1,
    sort_by='cardinality',
    sort_categories_by=None,
    element_size=30,
    facecolor="#3B5B92"
)
upset.plot()

# Set font size for all text elements in the plot
for text in plt.gcf().findobj(match=plt.Text):
    text.set_fontsize(18)

plt.tight_layout()
plt.savefig("Figure2e.pdf", dpi=300)
plt.show()