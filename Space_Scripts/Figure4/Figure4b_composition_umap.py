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

# Plot settings
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.use14corefonts"] = False
plt.rcParams["font.family"] = "Arial"

N_CLUSTERS, TOP_SPECIES = 15, 20
R_OUTER, R1, R2, R3, R4, R_INNER = 1.00, 0.97, 0.94, 0.91, 0.88, 0.00
timepoint_order = ["L-60", "L-30", "FD30", "FD90", "FD150", "R+1", "R+7", "R+14"]

# Load and filter metadata
df = pd.read_csv("cell_metadata_annotated.tsv", sep="\t")
top_clusters = df["cluster"].value_counts().nlargest(N_CLUSTERS).index.astype(str).tolist()
df = df[df["cluster"].astype(str).isin(top_clusters)].copy()
df["cluster"] = df["cluster"].astype(str)
top_species_list = df["species"].value_counts().nlargest(TOP_SPECIES).index.tolist()
df["species_filtered"] = np.where(df["species"].isin(top_species_list), df["species"], "Others")
df["astronaut"] = df["astronaut"].astype(str).str.strip()
astronaut_levels = df["astronaut"].dropna().unique().tolist()
timepoint_levels = [tp for tp in timepoint_order if tp in df["timepoint"].unique()]

# Define color maps
cluster_color_map = {
    "Oxidative stress response|0": "#A6CEE3", "Propionate production|1": "#3B8ABE", "Cofactor biosynthesis|2": "#72B29C",
    "Central carbon metabolism|3": "#84C868", "SCFA metabolism|4": "#4F9F3B", "Oxidative stress adaptation|5": "#EC9A91",
    "Cell envelope remodeling|6": "#E93E3F", "Host interface adaptation|7": "#F06C45", "Environmental sensing|8": "#FDAC4F",
    "Anaerobic metabolism|9": "#FB820F", "Envelope and motility|10": "#D1AAB7", "Metabolic flexibility|11": "#8C66AF",
    "Propionate utilization|12": "#A99099", "Anaerobic energy metabolism|13": "#EEDB80", "Organic acid metabolism|14": "#B15928"
}

species_color_map = {
    "Agathobacter faecis": "#8A97A8", "Agathobacter rectalis": "#B4C2D3", "Bacteroides eggerthii": "#B4D1AC",
    "Bacteroides fragilis": "#D8A27E", "Bacteroides intestinalis": "#E8C8B0", "Bacteroides ovatus": "#B8A878",
    "Bacteroides stercoris": "#CBB987", "Bacteroides thetaiotaomicron": "#D8C89A", "Bacteroides uniformis": "#A89B68",
    "Bacteroides xylanisolvens": "#E2D2B0", "CAG-81 sp900066785": "#C8C2BC", "Clostridium_Q sp003024715": "#C898A8",
    "Enterocloster bolteae": "#BCA8C8", "Enterocloster sp000431375": "#CDB6D6", "Enterocloster sp001517625": "#A888B8",
    "Escherichia coli_D": "#9674A8", "Faecalibacterium prausnitzii_C": "#7FA6D8", "Fusicatenibacter saccharivorans": "#D0B8D8",
    "Fusobacterium_A mortiferum": "#8F8F78", "Fusobacterium_A varium": "#B1B194", "Megamonas funiformis": "#80A8A8",
    "Parabacteroides distasonis": "#A8C8C8", "Parabacteroides merdae": "#D88888", "Parasutterella excrementihominis": "#F0B8B8",
    "Phascolarctobacterium faecium": "#C9A6D8", "Phocaeicola dorei": "#E8D090", "Phocaeicola massiliensis": "#D8C4A0",
    "Phocaeicola vulgatus": "#A88C78", "Prevotella stercorea": "#E0C8B8", "Roseburia intestinalis": "#8298B0",
    "Roseburia sp900552665": "#B0C8D8", "Sutterella wadsworthensis": "#98BC98", "UBA7182 sp003480725": "#C48898",
    "Others": "#BDBDBD"
}

astronaut_color_map = {"Astronaut 1": "#FF9A8B", "Astronaut 2": "#FFE0AC", "Astronaut 3": "#BBDED6"}
timepoint_color_map = {"L-60": "#F4F1DE", "L-30": "#E6D8AD", "FD30": "#5F6F75", "FD90": "#81B295", "FD150": "#A3BFA8", "R+1": "#9A031E", "R+7": "#BA3A26", "R+14": "#E07A5F"}
species_keys = top_species_list + (["Others"] if "Others" not in top_species_list else [])

# Calculate within-cluster compositions
cluster_to_species, cluster_to_astronaut, cluster_to_timepoint = {}, {}, {}
for c in top_clusters:
    sub = df[df["cluster"] == c]
    cluster_to_species[c] = sub["species_filtered"].value_counts(normalize=True).reindex(species_keys, fill_value=0.0).to_dict()
    cluster_to_astronaut[c] = sub["astronaut"].value_counts(normalize=True).reindex(astronaut_levels, fill_value=0.0).to_dict() if astronaut_levels else {}
    cluster_to_timepoint[c] = sub["timepoint"].value_counts(normalize=True).reindex(timepoint_levels, fill_value=0.0).to_dict() if timepoint_levels else {}

# Draw proportional ring segments
def draw_ring(ax, start_deg, end_deg, r_inner, r_outer, parts_dict, color_map):
    current = start_deg
    for k, frac in parts_dict.items():
        if frac <= 0: continue
        theta1, theta2 = current, current + frac * (end_deg - start_deg)
        ax.add_patch(Wedge((0, 0), r_outer, theta1, theta2, width=r_outer - r_inner, facecolor=color_map.get(k, (0.8, 0.8, 0.8)), edgecolor="white", linewidth=0.1))
        current = theta2

# Create circular plot
fig, ax = plt.subplots(figsize=(12, 12), subplot_kw=dict(aspect="equal"))
ax.set_xlim(-1.05, 1.05); ax.set_ylim(-1.05, 1.05); ax.axis("off")
n, angle_per_cluster = len(top_clusters), 360.0 / len(top_clusters)

for i, c in enumerate(top_clusters):
    theta1, theta2 = i * angle_per_cluster, (i + 1) * angle_per_cluster
    ax.add_patch(Wedge((0, 0), R_OUTER, theta1, theta2, width=R_OUTER - R1, facecolor=cluster_color_map.get(c, "grey"), edgecolor="white", linewidth=0.8))
    draw_ring(ax, theta1, theta2, R1, R2, cluster_to_species[c], species_color_map)
    draw_ring(ax, theta1, theta2, R2, R3, cluster_to_astronaut[c], astronaut_color_map)
    draw_ring(ax, theta1, theta2, R3, R4, cluster_to_timepoint[c], timepoint_color_map)

if R_INNER > 0:
    ax.add_patch(Wedge((0, 0), R4, 0, 360, width=R4 - R_INNER, facecolor="white", edgecolor="white"))

plt.savefig("Figure4b_circular.pdf", dpi=300, bbox_inches="tight")
plt.close(fig)

# Create legends
legend_fig, legend_ax = plt.subplots(figsize=(15, 10))
legend_ax.axis("off")
legends = [
    ("Cluster", [Line2D([0], [0], marker="o", linestyle="", markersize=8, markerfacecolor=cluster_color_map.get(c, "grey"), markeredgecolor="none", label=c) for c in top_clusters], 3),
    ("Species", [Line2D([0], [0], marker="o", linestyle="", markersize=8, markerfacecolor=species_color_map.get(s, "grey"), markeredgecolor="none", label=s) for s in species_keys], 3),
    ("Astronaut", [Line2D([0], [0], marker="o", linestyle="", markersize=8, markerfacecolor=astronaut_color_map.get(a, "grey"), markeredgecolor="none", label=a) for a in astronaut_levels], 3),
    ("Timepoint", [Line2D([0], [0], marker="o", linestyle="", markersize=8, markerfacecolor=timepoint_color_map.get(t, "grey"), markeredgecolor="none", label=t) for t in timepoint_order if t in timepoint_levels], 4)
]
y0 = 1.0
for title, handles, ncol in legends:
    if not handles: continue
    leg = legend_ax.legend(handles=handles, title=title, loc="upper center", bbox_to_anchor=(0.5, y0), ncol=ncol, frameon=False, fontsize=15, title_fontsize=25)
    legend_ax.add_artist(leg)
    if title == "Species":
        for text in leg.get_texts(): text.set_fontstyle("italic")
    y0 -= {"Cluster": 0.28, "Species": 0.36, "Astronaut": 0.12, "Timepoint": 0.1}[title]

legend_fig.savefig("Figure4b_legend.pdf", dpi=300, bbox_inches="tight")
plt.close(legend_fig)