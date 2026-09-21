# -*- coding: utf-8 -*-
"""
Created on Wed Jan 7 19:01:32 2026
@author: ZHENG XINGHAI
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Plot settings
mpl.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "Arial"})

input_dir = "species_meta_pseudotime"
species_order = ["Phocaeicola_vulgatus", "Prevotella_stercorea", "Phocaeicola_dorei", "Parabacteroides_distasonis", "Bacteroides_thetaiotaomicron", "Bacteroides_ovatus"]

species_color_map = {
    "Agathobacter_faecis": "#8A97A8", "Agathobacter_rectalis": "#B4C2D3",
    "Bacteroides_eggerthii": "#B4D1AC", "Bacteroides_fragilis": "#D8A27E", "Bacteroides_intestinalis": "#E8C8B0",
    "Bacteroides_ovatus": "#B8A878", "Bacteroides_stercoris": "#CBB987", "Bacteroides_thetaiotaomicron": "#D8C89A",
    "Bacteroides_uniformis": "#A89B68", "Bacteroides_xylanisolvens": "#E2D2B0", "CAG-81_sp900066785": "#C8C2BC",
    "Clostridium_Q_sp003024715": "#C898A8", "Enterocloster_bolteae": "#BCA8C8", "Enterocloster_sp000431375": "#CDB6D6",
    "Enterocloster_sp001517625": "#A888B8", "Escherichia_coli_D": "#9674A8", "Faecalibacterium_prausnitzii_C": "#7FA6D8",
    "Fusicatenibacter_saccharivorans": "#D0B8D8", "Fusobacterium_A_mortiferum": "#8F8F78", "Fusobacterium_A_varium": "#B1B194",
    "Megamonas_funiformis": "#80A8A8", "Parabacteroides_distasonis": "#A8C8C8", "Parabacteroides_merdae": "#D88888",
    "Parasutterella_excrementihominis": "#F0B8B8", "Phascolarctobacterium_faecium": "#C9A6D8", "Phocaeicola_dorei": "#E8D090",
    "Phocaeicola_massiliensis": "#D8C4A0", "Phocaeicola_vulgatus": "#A88C78", "Prevotella_stercorea": "#E0C8B8",
    "Roseburia_intestinalis": "#8298B0", "Roseburia_sp900552665": "#B0C8D8", "Sutterella_wadsworthensis": "#98BC98",
    "UBA7182_sp003480725": "#C48898", "Others": "#BDBDBD"
}

timepoint_color_map = {
    "L-60": "#F4F1DE", "L-30": "#E6D8AD", "FD30": "#5F6F75", "FD90": "#81B295",
    "FD150": "#A3BFA8", "R+1": "#9A031E", "R+7": "#BA3A26", "R+14": "#E07A5F"
}
timepoint_order = ["FD30", "FD90", "FD150"]

# Collect available files
tsv_files = [os.path.join(input_dir, f"{species}.tsv") for species in species_order if os.path.exists(os.path.join(input_dir, f"{species}.tsv"))]
if not tsv_files: raise FileNotFoundError("No TSV files found!")

# Create figure
fig, axes = plt.subplots(1, len(tsv_files), figsize=(3 * len(tsv_files), 3), sharey=True)
if len(tsv_files) == 1: axes = [axes]

# Plot each species
for i, file in enumerate(tsv_files):
    ax = axes[i]
    species = os.path.splitext(os.path.basename(file))[0]
    df = pd.read_csv(file, sep="\t")
    df["Timepoint"] = df["timepoint"].astype(str)
    df = df[df["Timepoint"].isin(timepoint_order)].copy()
    if len(df) == 0: continue

    pt = df["pseudotime"].astype(float).values
    df["Pseudo_norm"] = 50 if pt.max() == pt.min() else ((pt - pt.min()) / (pt.max() - pt.min())) * 100

    violin_data = []
    for tp in timepoint_order:
        vals = df.loc[df["Timepoint"] == tp, "Pseudo_norm"].values
        violin_data.append(np.repeat(vals, 2) if len(vals) == 1 else vals)

    vp = ax.violinplot(violin_data, positions=[1, 2, 3], widths=0.6, showmeans=False, showmedians=True, showextrema=False)
    fill_color = species_color_map.get(species, "#BDBDBD")

    for body, tp in zip(vp["bodies"], timepoint_order):
        body.set_facecolor(fill_color)
        body.set_edgecolor(timepoint_color_map[tp])
        body.set_linewidth(1)
        body.set_alpha(0.8)

    vp["cmedians"].set_colors([timepoint_color_map[tp] for tp in timepoint_order])
    vp["cmedians"].set_linewidth(1)

    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(timepoint_order, fontsize=15)
    ax.set_ylim(0, 100)
    ax.set_xlim(0.5, 3.5)
    ax.set_title(species.replace("_", " "), fontsize=13, fontstyle="italic")
    ax.tick_params(axis="x", labelsize=15)

    if i == 0:
        ax.set_ylabel("Normalized pseudotime", fontsize=15)
        ax.tick_params(axis="y", labelsize=15)
    else:
        ax.tick_params(axis="y", left=False, labelleft=False)
        ax.spines["left"].set_visible(False)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

# Save figure
plt.subplots_adjust(wspace=0.08)
plt.savefig("FigureS5c.pdf", dpi=300, bbox_inches="tight")
plt.show()