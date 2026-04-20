# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 21:19:56 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

# Configure matplotlib to export vector figures with editable text
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['font.family'] = 'Arial'

# Load pseudotime values for all cells
df_pseudo = pd.read_csv("species_pseudotime.csv", sep=",")

# Load cell metadata containing cluster, timepoint and stage information
df_meta = pd.read_csv("Phocaeicola_vulgatus_metadata.tsv", sep="\t")

# Merge pseudotime values with metadata using cell barcodes
df_merged = df_pseudo.merge(
    df_meta[["barcode", "cluster", "timepoint", "stage"]],
    left_on="cell",
    right_on="barcode",
    how="left"
)

# Remove redundant barcode column after merging
df_merged = df_merged.drop(columns=["barcode"])

# Save merged pseudotime and metadata table for downstream analysis
df_merged.to_csv("species_meta_pseudotime.csv", index=False)

# Initialize visualization style for pseudotime plots
sns.set_style("dark")
sns.set_context("talk")

# Create violin plot showing pseudotime distribution across timepoints
plt.figure(figsize=(10, 6))

order = ["L-60","L-30","FD30","FD90","FD150","R+1","R+7","R+14"]
uniform_color = "tomato"

sns.violinplot(
    data=df_merged,
    x="timepoint",
    y="pseudotime",
    inner="box",
    order=order,
    color=uniform_color
)

sns.stripplot(
    data=df_merged,
    x="timepoint",
    y="pseudotime",
    color="black",
    alpha=0.2,
    jitter=0.2,
    size=1,
    order=order
)

plt.xlabel("Timepoint", fontsize=25)
plt.ylabel("Pseudotime", fontsize=25)

plt.xticks(fontsize=20, rotation=0)
plt.yticks(fontsize=20)

plt.tight_layout()
plt.savefig("Figure5d_timepoint.pdf", format="pdf")
plt.show()

# Create violin plot showing pseudotime distribution across biological stages
sns.set_style("dark")
sns.set_context("talk")

plt.figure(figsize=(5, 6))

order = ["L","FD","R"]
uniform_color = "tomato"

sns.violinplot(
    data=df_merged,
    x="stage",
    y="pseudotime",
    inner="box",
    order=order,
    color=uniform_color
)

sns.stripplot(
    data=df_merged,
    x="stage",
    y="pseudotime",
    color="black",
    alpha=0.2,
    jitter=0.2,
    size=1,
    order=order
)

plt.xlabel("Stage", fontsize=25)
plt.ylabel("Pseudotime", fontsize=25)

plt.xticks(fontsize=20, rotation=0)
plt.yticks(fontsize=20)

plt.tight_layout()
plt.savefig("Figure5d_stage.pdf", format="pdf")
plt.show()