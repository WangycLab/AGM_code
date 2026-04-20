# -*- coding: utf-8 -*-
"""
Created on Wed Sep 17 16:35:46 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler

# Configure matplotlib so that exported PDF text remains editable
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'

# Define a cosine similarity function between two vectors
def cosine_sim(a, b):
    if np.allclose(a, 0) or np.allclose(b, 0):
        return 0.0
    return float(np.dot(a, b) / (norm(a) * norm(b)))

# Calculate similarity and stability of species distributions across groups
def calc_similarity(df, stats, group_col, cluster_totals):
    similarities = []
    stabilities = []
    groups = sorted(df[group_col].unique())

    species_counts = (
        df.groupby(["species", group_col, "cluster"])
        .size()
        .reset_index(name="count")
    )

    species_counts = species_counts.merge(cluster_totals, on="cluster", how="left")
    species_counts["percent"] = (
        species_counts["count"] / species_counts["cluster_total"]
    )

    for species in stats["species"]:
        sub = species_counts[species_counts["species"] == species]

        mat = sub.pivot_table(
            index="cluster",
            columns=group_col,
            values="percent",
            fill_value=0
        )

        for g in groups:
            if g not in mat.columns:
                mat[g] = 0.0
        mat = mat[groups]

        vectors = [mat[g].values for g in groups]

        sims = []
        for i in range(len(vectors)):
            for j in range(i + 1, len(vectors)):
                sims.append(cosine_sim(vectors[i], vectors[j]))

        if sims:
            similarities.append(float(np.mean(sims)))
            stabilities.append(float(np.std(sims)))
        else:
            similarities.append(0.0)
            stabilities.append(0.0)

    return similarities, stabilities


# Load the cell-level metadata table
df = pd.read_csv("cell_metadata.tsv", sep="\t")

# Compute cell and cluster statistics for each species
stats = (
    df.groupby("species")
    .agg(
        cell_count=("barcode", "count"),
        cluster_count=("cluster", "nunique")
    )
    .reset_index()
)

# Calculate total cell counts for each cluster
cluster_totals = (
    df.groupby("cluster")
    .size()
    .reset_index(name="cluster_total")
)

# Compute similarity and stability metrics for different grouping variables
stats["astronaut_similarity"], stats["astronaut_stability"] = calc_similarity(
    df, stats, "astronaut", cluster_totals
)
stats["stage_similarity"], stats["stage_stability"] = calc_similarity(
    df, stats, "stage", cluster_totals
)
stats["timepoint_similarity"], stats["timepoint_stability"] = calc_similarity(
    df, stats, "timepoint", cluster_totals
)

# Save similarity statistics to disk
output_file = "species_similarity_stats.csv"
stats.to_csv(output_file, index=False)
print(f"Results saved to {output_file}")

# Filter species with sufficient cell numbers and sort by abundance
stats_filtered = stats[stats["cell_count"] >= 3000].copy()
stats_filtered = stats_filtered.sort_values("cell_count", ascending=False)

# Normalize similarity scores using Min-Max scaling
similarity_cols = [
    "astronaut_similarity",
    "stage_similarity",
    "timepoint_similarity"
]
scaler = MinMaxScaler()
stats_filtered[similarity_cols] = scaler.fit_transform(
    stats_filtered[similarity_cols]
)

# Convert similarity metrics into long format for bubble plotting
bubble_df = pd.melt(
    stats_filtered,
    id_vars=["species"],
    value_vars=similarity_cols,
    var_name="group",
    value_name="Similarity"
)

# Attach stability values corresponding to each similarity measurement
stability_map = {
    "astronaut_similarity": "astronaut_stability",
    "stage_similarity": "stage_stability",
    "timepoint_similarity": "timepoint_stability"
}
bubble_df["Stability"] = bubble_df.apply(
    lambda row: stats_filtered.loc[
        stats_filtered["species"] == row["species"],
        stability_map[row["group"]]
    ].values[0],
    axis=1
)

# Clean group names for visualization
bubble_df["group"] = bubble_df["group"].str.replace("_similarity", "")

# Configure seaborn plotting style and context
sns.set_style("dark")
sns.set_context("talk")

# Define font settings for axes, ticks, and legend
font_axis = {'fontname': 'Arial', 'fontsize': 28}
font_tick = {'labelsize': 20}
font_legend = {'size': 18}

# Create a multi-panel figure showing abundance and similarity metrics
fig, axes = plt.subplots(
    1, 3,
    figsize=(9, 12),
    sharey=True,
    gridspec_kw={'width_ratios': [1, 1, 2], 'wspace': 0.1}
)

# Plot cell counts per species
sns.barplot(
    data=stats_filtered,
    x="cell_count",
    y="species",
    ax=axes[0],
    color="#4C72B0",
    edgecolor=None
)
axes[0].set_xlabel("nCells", **font_axis)
axes[0].set_ylabel("Species", **font_axis)
axes[0].tick_params(axis="x", **font_tick)
axes[0].tick_params(axis="y", **font_tick)

# Plot cluster counts per species
sns.barplot(
    data=stats_filtered,
    x="cluster_count",
    y="species",
    ax=axes[1],
    color="#008080",
    edgecolor=None
)
axes[1].set_xlabel("nClusters", **font_axis)
axes[1].set_ylabel("")
axes[1].tick_params(axis="x", **font_tick)
axes[1].tick_params(axis="y", **font_tick)

# Plot similarity and stability as a bubble chart
sns.scatterplot(
    data=bubble_df,
    x="group",
    y="species",
    size="Stability",
    hue="Similarity",
    sizes=(50, 400),
    palette="mako_r",
    ax=axes[2],
    legend="brief",
    edgecolor=None
)
axes[2].set_xlabel("Group", **font_axis)
axes[2].set_ylabel("")
axes[2].tick_params(axis="x", **font_tick)
axes[2].tick_params(axis="y", **font_tick)

# Customize x-axis labels for grouping variables
axes[2].set_xticks([0, 1, 2])
axes[2].set_xticklabels(
    ["Astronaut", "Stage", "Timepoint"],
    fontname="Arial",
    fontsize=20
)
axes[2].set_xlim(-0.5, 2.5)

# Apply italic formatting to species names
species_labels = stats_filtered["species"].unique()
for ax in axes:
    ax.set_yticks(range(len(species_labels)))
    ax.set_yticklabels(
        species_labels,
        fontstyle="italic",
        fontname="Arial",
        fontsize=20
    )

# Remove axis spines for a cleaner appearance
for ax in axes:
    for spine in ax.spines.values():
        spine.set_visible(False)

# Adjust legend position and font properties
legend = axes[2].legend(
    bbox_to_anchor=(1.03, 1),
    loc="upper left",
    title="",
    labelspacing=0.5,
    handleheight=0.5,
    handlelength=0.5
)
for text in legend.get_texts():
    text.set_fontsize(font_legend["size"])
    text.set_family("Arial")

# Save the figure as an Illustrator-editable PDF and display it
plt.tight_layout()
plt.savefig(
    "Figure4d.pdf",
    dpi=300,
    bbox_inches="tight"
)
plt.show()