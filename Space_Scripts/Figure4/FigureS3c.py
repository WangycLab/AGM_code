# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 14:43:00 2025
@author: ZHENG XINGHAI
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load data
df = pd.read_csv("cell_metadata0.tsv", sep="\t")
species_counts = df["species"].value_counts()
thresholds = np.array([1, 2, 3, 5, 10, 20, 30, 50, 100, 200, 300, 500, 1000, 2000, 3000, 5000])
species_numbers = np.array([(species_counts >= th).sum() for th in thresholds])

# Detect elbow point
x, y = thresholds.astype(float), species_numbers.astype(float)
x_norm = (x - x.min()) / (x.max() - x.min())
y_norm = (y - y.min()) / (y.max() - y.min())
line_start, line_end = np.array([x_norm[0], y_norm[0]]), np.array([x_norm[-1], y_norm[-1]])
line_vec = line_end - line_start
line_vec /= np.linalg.norm(line_vec)
vecs = np.column_stack((x_norm, y_norm)) - line_start
proj_vec = np.outer(vecs @ line_vec, line_vec)
dist = np.linalg.norm(vecs - proj_vec, axis=1)
elbow_idx = np.argmax(dist)
elbow_x, elbow_y = thresholds[elbow_idx], species_numbers[elbow_idx]

# Plot
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.family"] = "Arial"

fig, ax = plt.subplots(figsize=(5, 4.5))
ax.plot(thresholds, species_numbers, color="#2C3E50", linewidth=2, marker="o", markersize=5, markerfacecolor="white", markeredgewidth=1.5)
ax.scatter(elbow_x, elbow_y, color="#8B0000", s=60, zorder=5)
ax.axvline(elbow_x, linestyle="--", color="#8B0000", linewidth=2, alpha=0.6)
ax.text(elbow_x * 1.5, elbow_y, f"Elbow point\nnCell={elbow_x}\nnSpecies={elbow_y}", fontsize=13, color="#8B0000", va="bottom", ha="left")
ax.set_xlabel("Cell number threshold", fontsize=15)
ax.set_ylabel("Number of species", fontsize=15)
major_ticks = [1, 300, 1000, 2000, 3000, 5000]
ax.set_xticks(major_ticks)
ax.set_xticklabels(major_ticks)
ax.tick_params(axis="both", labelsize=13, width=1)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_linewidth(1)
ax.spines["bottom"].set_linewidth(1)
ax.grid(False)

plt.tight_layout()
plt.savefig("FigureS3c.pdf", dpi=300, bbox_inches="tight")
plt.show()
print(f"Elbow point: nCell={elbow_x}, nSpecies={elbow_y}")