# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 16:22:45 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import jensenshannon
import seaborn as sns
from sklearn.manifold import MDS
from matplotlib.lines import Line2D
import matplotlib as mpl

# Use vector-friendly fonts
mpl.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "Arial",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "axes.edgecolor": "black",
    "axes.linewidth": 0.8,
    "font.size": 15,
    "axes.titlesize": 15,
    "axes.labelsize": 15,
    "xtick.labelsize": 15,
    "ytick.labelsize": 15,
    "legend.fontsize": 15,
    "lines.linewidth": 2.0,
    "lines.markersize": 5,
    "xtick.direction": "out",
    "ytick.direction": "out",
})

# Load merged species abundance across samples
super_composition_df = pd.read_csv("species_abundance_profile.csv", index_col=0)
data = super_composition_df.T 

# Compute similarity matrix using Jensen–Shannon similarity
samples = data.index
n = len(samples)
similarity_matrix = np.zeros((n, n))

for i in range(n):
    for j in range(n):
        d = jensenshannon(data.iloc[i], data.iloc[j])
        similarity_matrix[i, j] = 1 - d 

similarity_df = pd.DataFrame(similarity_matrix, index=samples, columns=samples)

# Perform MDS embedding for visualization
mds = MDS(n_components=2, dissimilarity="precomputed", random_state=1)
coords = mds.fit_transform(1 - similarity_df.values)

# Define color mapping for timepoints
time_colors = {
    "L-60": "#F4F1DE", "L-30": "#E6D8AD", "FD30": "#5F6F75",
    "FD90": "#81B295", "FD150": "#A3BFA8", "R+1": "#9A031E",
    "R+7": "#BA3A26", "R+14": "#E07A5F"
}

# Define marker mapping for astronaut replicates
astronaut_markers = {"_1": "o", "_2": "s", "_3": "^"}

# Create figure
plt.figure(figsize=(12, 8))

# Plot MDS scatter points
for i, label in enumerate(similarity_df.index):
    x, y = coords[i, 0], coords[i, 1]
    marker = next((m for key, m in astronaut_markers.items() if label.endswith(key)), "o")
    time_key = label.split("_")[0]
    color = time_colors.get(time_key, "silver")
    plt.scatter(x, y, c=color, s=500, marker=marker, edgecolor="k", linewidth=0.1)

plt.xlabel("MDS1", fontsize=30)
plt.ylabel("MDS2", fontsize=30)
plt.xticks(fontsize=28)
plt.yticks(fontsize=28)
plt.grid(False)

# Build astronaut legend
shape_legend = [Line2D([0], [0], marker=m, color="w", label=f"{a}",
                       markerfacecolor="gray", markersize=15)
                for a, m in zip(["1","2","3"], ["o", "s", "^"])]

# Build timepoint legend
color_legend = [Line2D([0], [0], marker="o", color="w", label=tp,
                       markerfacecolor=col, markersize=15)
                for tp, col in time_colors.items()]

# Add both legends manually to the same axes
ax = plt.gca()
legend1 = ax.legend(handles=shape_legend, title="Astronaut", title_fontsize=26,
                    fontsize=23, loc="upper right", bbox_to_anchor=(1.32, 0.99))
ax.add_artist(legend1)  # Keep this legend visible
ax.legend(handles=color_legend, title="Timepoint", title_fontsize=26,
          fontsize=23, loc="upper right", bbox_to_anchor=(1.35, 0.66))

plt.tight_layout()
plt.savefig("Figure3d.pdf", dpi=300, bbox_inches="tight")
plt.show()

# Define heatmap colormap
cmap_choice = "inferno"
astronaut_colors = ["#FF9A8B", "#FFE0AC", "#BBDED6"]
samples_per_astronaut = 8
astronauts = ["Astronaut 1", "Astronaut 2", "Astronaut 3"]

# Fixed order of timepoints for row annotation
y_time_labels = ["L-60","L-30","FD30","FD90","FD150","R+1","R+7","R+14"]*3

plt.figure(figsize=(12, 12))
ax = sns.heatmap(similarity_df, cmap=cmap_choice, annot=False, square=True, linewidths=0,
                 cbar=True, cbar_kws={"shrink":0.78, "orientation":"vertical", "pad":0.03, "aspect":30})

# Customize colorbar
cbar = ax.collections[0].colorbar
cbar.set_label("Jensen-Shannon similarity", fontsize=28)
cbar.ax.tick_params(labelsize=25, width=3, length=8)
cbar.ax.set_aspect(30)
cbar.ax.set_position([0.93, 0.25, 0.015, 0.5])

ax.set_xticks([])
ax.set_xticklabels([])

# Add astronaut annotation bar above heatmap
bar_ax = ax.inset_axes([0, 1.02, 1, 0.03])
bar_ax.axis("off")
num_samples = similarity_df.shape[1]
for i, color in enumerate(astronaut_colors):
    start = i * samples_per_astronaut / num_samples
    width = samples_per_astronaut / num_samples
    bar_ax.add_patch(plt.Rectangle((start, 0), width, 1, color=color))
    bar_ax.text(start + width/2, 1.1, astronauts[i], ha="center", va="bottom", fontsize=25)

ax.set_yticks([])
ax.set_yticklabels([])

# Add timepoint annotation bar on the left
bar_ax_y = ax.inset_axes([-0.05, 0, 0.03, 1])
bar_ax_y.axis("off")
num_rows = similarity_df.shape[0]
for idx, time_label in enumerate(y_time_labels):
    color = time_colors.get(time_label, "#CCCCCC")
    bar_ax_y.add_patch(plt.Rectangle((0, 1-(idx+1)/num_rows), 1, 1/num_rows, color=color))
    bar_ax_y.text(-0.1, 1-(idx+0.5)/num_rows, time_label, va="center", ha="right", fontsize=18)

# Draw separation lines between timepoints
split_chars = ["L", "F", "R"]
positions_red = [i+1 for i in range(num_rows-1) if similarity_df.index[i][0] != similarity_df.index[i+1][0]]
for pos in positions_red:
    ax.hlines(pos, *ax.get_xlim(), color="white", linewidth=1)
    ax.vlines(pos, *ax.get_ylim(), color="white", linewidth=1)

# Draw separation lines between astronaut groups
sample_groups = ["_1","_2","_3"]
current_pos = 0
group_positions = []
for sample in sample_groups:
    current_pos += (similarity_df.index.str.endswith(sample)).sum()
    group_positions.append(current_pos)

for pos in group_positions[:-1]:
    ax.hlines(pos, *ax.get_xlim(), color="gold", linewidth=1.5)
    ax.vlines(pos, *ax.get_ylim(), color="gold", linewidth=1.5)

plt.tight_layout()
plt.savefig("Figure3c.pdf", dpi=300, bbox_inches="tight")
plt.show()