# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 16:08:05 2025
@author: ZHENG XINGHAI
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
from matplotlib.lines import Line2D

# Configure matplotlib to export vector graphics with editable text
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['pdf.use14corefonts'] = False
plt.rcParams["font.family"] = "Arial"

# Define the number of clusters and species displayed in the plot
N_CLUSTERS = 23
TOP_SPECIES = 20

# Define radii of concentric rings from outer to inner
R_OUTER = 1.00
R1 = 0.97
R2 = 0.94
R3 = 0.91
R4 = 0.88
R_INNER = 0.00

# Define chronological order of timepoints
timepoint_order = ["L-60", "L-30", "FD30", "FD90", "FD150", "R+1", "R+7", "R+14"]

# Load the cell metadata table
df = pd.read_csv("cell_metadata.tsv", sep="\t")

# Select the most abundant clusters
cluster_counts = df["cluster"].value_counts()
top_clusters = cluster_counts.nlargest(N_CLUSTERS).index.astype(str).tolist()

# Filter dataset to retain only selected clusters
df = df[df["cluster"].astype(str).isin(top_clusters)].copy()
df["cluster"] = df["cluster"].astype(str)

# Identify globally abundant species
top_species_list = df["species"].value_counts().nlargest(TOP_SPECIES).index.tolist()

# Collapse low abundance species into an "Others" category
df["species_filtered"] = np.where(
    df["species"].isin(top_species_list),
    df["species"],
    "Others"
)

# Clean astronaut labels and retain observed order
df["astronaut"] = df["astronaut"].astype(str).str.strip()
astronaut_levels = df["astronaut"].dropna().unique().tolist()

# Retain timepoints present in the dataset
timepoint_levels = [tp for tp in timepoint_order if tp in df["timepoint"].unique()]

# Define color palette for clusters
cluster_color_map = {
    "0": "#A6CEE3","1": "#6CA9CF","2": "#3385BB","3": "#4693A8",
    "4": "#84BF96","5": "#A3D77F","6": "#6DBD57","7": "#37A22F",
    "8": "#7F9D55","9": "#D49B84","10": "#F57C7C","11": "#EB4647",
    "12": "#E42622","13": "#F06C45","14": "#FBB268","15": "#FDA848",
    "16": "#FE8D19","17": "#F48829","18": "#DE9E83","19": "#C6ADD3",
    "20": "#9D7BBA","21": "#754AA0","22": "#977899","23": "#D6CA99",
    "24": "#F3E587","25": "#D29F57","26": "#B15928"
}

# Define color palette for species categories
species_color_map = {
    "Phocaeicola vulgatus": "#8DD3C7",
    "Phocaeicola dorei": "#CFECBB",
    "Prevotella stercorea": "#F4F4B9",
    "Megamonas funiformis": "#CFCCCF",
    "Parabacteroides distasonis": "#D1A7B9",
    "Bacteroides xylanisolvens": "#F4867C",
    "UBA3766 sp900769175": "#86B1CD",
    "Bacteroides ovatus": "#C0979F",
    "Phocaeicola massiliensis": "#CEB28B",
    "CAG-81 sp900066535": "#EDBC63",
    "Agathobacter rectalis": "#C2D567",
    "Bacteroides fragilis": "#CDD796",
    "Roseburia sp900552665": "#F8CDDE",
    "Fusicatenibacter saccharivorans": "#E9D3DE",
    "Dorea_A longicatena_B": "#D5CFD6",
    "Agathobacter faecis": "#C59CC5",
    "Escherichia coli_D": "#C09CBF",
    "Bacteroides cellulosilyticus": "#C9DAC3",
    "Parasutterella excrementihominis": "#E1EBA0",
    "Sutterella wadsworthensis": "#FFED6F",
    "Others": "#E5E5E5"
}

# Define color palette for astronauts
astronaut_color_map = {
    "Astronaut 1": "#FF9A8B",
    "Astronaut 2": "#FFE0AC",
    "Astronaut 3": "#BBDED6"
}

# Define color palette for timepoints
timepoint_color_map = {
    "L-60": "#F4F1DE",
    "L-30": "#E6D8AD",
    "FD30": "#5F6F75",
    "FD90": "#81B295",
    "FD150": "#A3BFA8",
    "R+1": "#9A031E",
    "R+7": "#BA3A26",
    "R+14": "#E07A5F"
}

# Ensure species legend includes the Others category
species_keys = top_species_list + (["Others"] if "Others" not in top_species_list else [])

# Compute proportional composition within each cluster
cluster_to_species = {}
cluster_to_astronaut = {}
cluster_to_timepoint = {}

for c in top_clusters:
    sub = df[df["cluster"] == c]

    sp = (
        sub["species_filtered"]
        .value_counts(normalize=True)
        .reindex(species_keys, fill_value=0.0)
        .to_dict()
    )
    cluster_to_species[c] = sp

    if len(astronaut_levels) > 0:
        a = (
            sub["astronaut"]
            .value_counts(normalize=True)
            .reindex(astronaut_levels, fill_value=0.0)
            .to_dict()
        )
    else:
        a = {}
    cluster_to_astronaut[c] = a

    if len(timepoint_levels) > 0:
        tp = (
            sub["timepoint"]
            .value_counts(normalize=True)
            .reindex(timepoint_levels, fill_value=0.0)
            .to_dict()
        )
    else:
        tp = {}
    cluster_to_timepoint[c] = tp

# Draw proportional segments inside a circular ring
def draw_ring(ax, start_deg, end_deg, r_inner, r_outer, parts_dict, color_map):
    angle_span = end_deg - start_deg
    current = start_deg
    for k, frac in parts_dict.items():
        if frac <= 0:
            continue
        theta1 = current
        theta2 = current + frac * angle_span
        wedge = Wedge(
            (0,0),
            r_outer,
            theta1,
            theta2,
            width=(r_outer - r_inner),
            facecolor=color_map.get(k,(0.8,0.8,0.8)),
            edgecolor="white",
            linewidth=0.1
        )
        ax.add_patch(wedge)
        current = theta2

# Initialize circular plotting canvas
fig, ax = plt.subplots(figsize=(12,12), subplot_kw=dict(aspect="equal"))
ax.set_xlim(-1.05,1.05)
ax.set_ylim(-1.05,1.05)
ax.axis("off")

# Determine angular width per cluster
n = len(top_clusters)
angle_per_cluster = 360.0 / n

# Draw outer ring representing cluster identity
for i,c in enumerate(top_clusters):
    theta1 = i * angle_per_cluster
    theta2 = (i + 1) * angle_per_cluster
    wedge = Wedge(
        (0,0),R_OUTER,
        theta1,theta2,
        width=(R_OUTER-R1),
        facecolor=cluster_color_map.get(c,"grey"),
        edgecolor="white",
        linewidth=0.8
    )
    ax.add_patch(wedge)

# Draw species composition ring
for i,c in enumerate(top_clusters):
    draw_ring(ax,
              i*angle_per_cluster,
              (i+1)*angle_per_cluster,
              R1,R2,
              cluster_to_species[c],
              species_color_map)

# Draw astronaut distribution ring
for i,c in enumerate(top_clusters):
    draw_ring(ax,
              i*angle_per_cluster,
              (i+1)*angle_per_cluster,
              R2,R3,
              cluster_to_astronaut[c],
              astronaut_color_map)

# Draw timepoint distribution ring
for i,c in enumerate(top_clusters):
    draw_ring(ax,
              i*angle_per_cluster,
              (i+1)*angle_per_cluster,
              R3,R4,
              cluster_to_timepoint[c],
              timepoint_color_map)

# Optionally draw central white core
if R_INNER > 0:
    core = Wedge((0,0),R4,0,360,width=R4-R_INNER,
                 facecolor="white",edgecolor="white")
    ax.add_patch(core)

# Save circular annotation figure
plt.savefig("Figure4c_circular.pdf",dpi=300,bbox_inches="tight")
plt.close(fig)

# Create legend figure
legend_fig, legend_ax = plt.subplots(figsize=(15,10))
legend_ax.axis("off")

legends = []

# Prepare cluster legend
cluster_handles = [
    Line2D([0],[0],
           marker='o',
           linestyle='',
           markersize=8,
           markerfacecolor=cluster_color_map.get(c,"grey"),
           markeredgecolor='none',
           label=f"Cluster {c}")
    for c in top_clusters
]
legends.append(("Cluster",cluster_handles,5))


# Prepare species legend
species_handles = [
    Line2D([0],[0],
           marker='o',
           linestyle='',
           markersize=8,
           markerfacecolor=species_color_map.get(s,"grey"),
           markeredgecolor='none',
           label=s)
    for s in species_keys
]
legends.append(("Species",species_handles,3))


# Prepare astronaut legend
astronaut_handles = [
    Line2D([0],[0],
           marker='o',
           linestyle='',
           markersize=8,
           markerfacecolor=astronaut_color_map.get(a,"grey"),
           markeredgecolor='none',
           label=a)
    for a in astronaut_levels
]

if astronaut_handles:
    legends.append(("Astronaut",astronaut_handles,3))


# Prepare timepoint legend
timepoint_handles = [
    Line2D([0],[0],
           marker='o',
           linestyle='',
           markersize=8,
           markerfacecolor=timepoint_color_map.get(t,"grey"),
           markeredgecolor='none',
           label=t)
    for t in timepoint_order if t in timepoint_levels
]

if timepoint_handles:
    legends.append(("Timepoint",timepoint_handles,4))


# Arrange legends vertically
y0 = 1.0

for title,handles,ncol in legends:

    leg = legend_ax.legend(
        handles=handles,
        title=title,
        loc="upper center",
        bbox_to_anchor=(0.5,y0),
        ncol=ncol,
        frameon=False,
        fontsize=15,
        title_fontsize=25
    )

    legend_ax.add_artist(leg)

    # Set species names in italic style while keeping text editable
    if title == "Species":
        for text in leg.get_texts():
            text.set_fontstyle('italic')

    if title == "Cluster":
        y0 -= 0.28
    elif title == "Species":
        y0 -= 0.36
    elif title == "Astronaut":
        y0 -= 0.12
    else:
        y0 -= 0.1


# Save legend figure
legend_fig.savefig(
    "Figure4c_legend.pdf",
    dpi=300,
    bbox_inches="tight"
)

plt.close(legend_fig)