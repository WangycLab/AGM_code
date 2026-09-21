# -*- coding: utf-8 -*-
"""
Created on Thu Jan 1 20:49:41 2026
@author: ZHENG XINGHAI
"""

import os
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.spatial.distance import jensenshannon
from kneed import KneeLocator

# Elbow detection based on read counts
def jsd_elbow_from_sorted_df(df_sorted):
    full_df = df_sorted.groupby("name")["new_est_reads"].sum()
    full_df /= full_df.sum()
    pct_values = list(range(95, 4, -5))
    similarities = []
    for pct in pct_values:
        cutoff = int(np.ceil(len(df_sorted) * pct / 100))
        sub_df = df_sorted.iloc[:cutoff].groupby("name")["new_est_reads"].sum()
        sub_df /= sub_df.sum()
        aligned = pd.concat([full_df, sub_df], axis=1).fillna(0)
        similarities.append(1 - jensenshannon(aligned.iloc[:, 0].values, aligned.iloc[:, 1].values))
    similarity_df = pd.DataFrame({"percent": pct_values, "similarity": similarities}).sort_values("percent")
    knee = KneeLocator(similarity_df["percent"].values, similarity_df["similarity"].values, curve="concave", direction="increasing")
    return knee.knee, similarity_df

# Elbow detection based on purity threshold
def jsd_elbow_purity_threshold(df, step=0.05):
    full_df = df.groupby("name")["new_est_reads"].sum()
    full_df /= full_df.sum()
    thresholds = np.round(np.arange(0.2, 0.8 + step, step), 2).tolist()
    similarities = []
    for t in thresholds:
        df_remaining = df[df["fraction_total_reads"] >= t]
        if df_remaining.empty:
            break
        sub_df = df_remaining.groupby("name")["new_est_reads"].sum()
        sub_df /= sub_df.sum()
        aligned = pd.concat([full_df, sub_df], axis=1).fillna(0)
        similarities.append(1 - jensenshannon(aligned.iloc[:, 0].values, aligned.iloc[:, 1].values))
    purity_jsd_df = pd.DataFrame({"purity_threshold": thresholds[:len(similarities)], "similarity": similarities}).sort_values("purity_threshold")
    knee = KneeLocator(purity_jsd_df["purity_threshold"].values, purity_jsd_df["similarity"].values, curve="concave", direction="decreasing")
    return knee.knee, purity_jsd_df

# Set paths and sample order
input_dir, output_curve_dir, summary_path = "space_species_filted", "jsd_curve", "summary_cell_counts.csv"
os.makedirs(output_curve_dir, exist_ok=True)
timepoints, replicates = ["L-60", "L-30", "FD30", "FD90", "FD150", "R+1", "R+7", "R+14"], ["_1", "_2", "_3"]
sample_order = [tp + rep for rep in replicates for tp in timepoints]

# Determine optimal filtering thresholds
summary_df = pd.read_csv(summary_path)
all_compositions, final_selected_dict, optimal_read_cutoff, optimal_purity_cutoff = [], {}, {}, {}

for sample_name in sample_order:
    sample_file = os.path.join(input_dir, f"{sample_name}_sc_taxonomy_filted.csv")
    if not os.path.exists(sample_file):
        continue
    df = pd.read_csv(sample_file)
    df["total_reads"] = np.ceil(df["new_est_reads"] / df["fraction_total_reads"].replace(0, np.nan)).astype("Int64")
    df_read_sorted = df.sort_values("total_reads", ascending=False).reset_index(drop=True)
    read_pct, read_jsd_df = jsd_elbow_from_sorted_df(df_read_sorted)
    read_cutoff_idx = int(np.ceil(len(df_read_sorted) * read_pct / 100))
    read_threshold = df_read_sorted.iloc[read_cutoff_idx - 1]["total_reads"]
    optimal_read_cutoff[sample_name] = {"read_pct": read_pct, "read_threshold": read_threshold}
    purity_elbow, purity_jsd_df = jsd_elbow_purity_threshold(df, step=0.05)
    optimal_purity_cutoff[sample_name] = {"purity_threshold": purity_elbow}
    df_final = df[(df["total_reads"] >= read_threshold) & (df["fraction_total_reads"] >= purity_elbow)]
    final_selected_dict[sample_name] = len(df_final)
    composition_df = df_final.groupby("name")["new_est_reads"].sum() / df_final["new_est_reads"].sum()
    composition_df.name = sample_name
    all_compositions.append(composition_df)

# Update summary and composition profiles
summary_df["Final_selected"] = summary_df["Sample"].map(final_selected_dict)
summary_df["Best_read_pct"] = summary_df["Sample"].map(lambda x: optimal_read_cutoff.get(x, {}).get("read_pct"))
summary_df["Best_read_threshold"] = summary_df["Sample"].map(lambda x: optimal_read_cutoff.get(x, {}).get("read_threshold"))
summary_df["Best_purity_threshold"] = summary_df["Sample"].map(lambda x: optimal_purity_cutoff.get(x, {}).get("purity_threshold"))
summary_df.to_csv("summary_cell_counts_completed.csv", index=False)

super_composition_df = pd.concat(all_compositions, axis=1).fillna(0)
super_composition_df = super_composition_df[[x for x in sample_order if x in super_composition_df.columns]]
super_composition_df.to_csv("species_abundance_profile.csv")

# Plot settings
mpl.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white", "axes.edgecolor": "black", "axes.linewidth": 0.8,
    "font.family": "Arial", "font.size": 15, "axes.titlesize": 15, "axes.labelsize": 15,
    "xtick.labelsize": 15, "ytick.labelsize": 15, "legend.fontsize": 15,
    "lines.linewidth": 2.0, "lines.markersize": 5, "xtick.direction": "out", "ytick.direction": "out",
    "pdf.fonttype": 42, "ps.fonttype": 42
})

# Plot read-complexity curves
fig, axes = plt.subplots(8, 3, figsize=(18, 24), sharex=True, sharey=True)
for i, tp in enumerate(timepoints):
    for j, rep in enumerate(replicates):
        ax, sample_name = axes[i, j], f"{tp}{rep}"
        sample_file = os.path.join(input_dir, f"{sample_name}_sc_taxonomy_filted.csv")
        if not os.path.exists(sample_file):
            ax.axis("off")
            continue
        df_sample = pd.read_csv(sample_file)
        df_sample["total_reads"] = df_sample["new_est_reads"] / df_sample["fraction_total_reads"].replace(0, np.nan)
        df_read_sorted = df_sample.sort_values("total_reads", ascending=False).reset_index(drop=True)
        read_pct, read_jsd_df = jsd_elbow_from_sorted_df(df_read_sorted)
        ax.plot(read_jsd_df["percent"], read_jsd_df["similarity"], linewidth=3)
        ax.axvline(read_pct, linestyle="--", linewidth=2)
        ax.text(0.95, 0.05, f"Elbow point = {read_pct:.0f}%", transform=ax.transAxes, ha="right", va="bottom", fontsize=15)
        ax.set_title(sample_name, fontsize=25, pad=10)
        if i == 7: ax.set_xlabel("Top N% Cells", fontsize=15)
        if j == 0: ax.set_ylabel("Jensen–Shannon similarity", fontsize=15)

plt.tight_layout()
plt.savefig(os.path.join(output_curve_dir, "species_similarity_elbow_gene_counts.pdf"), dpi=300)
plt.close()

# Plot purity-threshold curves
fig, axes = plt.subplots(8, 3, figsize=(18, 24), sharex=True, sharey=True)
for i, tp in enumerate(timepoints):
    for j, rep in enumerate(replicates):
        ax, sample_name = axes[i, j], f"{tp}{rep}"
        sample_file = os.path.join(input_dir, f"{sample_name}_sc_taxonomy_filted.csv")
        if not os.path.exists(sample_file):
            ax.axis("off")
            continue
        df_sample = pd.read_csv(sample_file)
        purity_elbow, purity_jsd_df = jsd_elbow_purity_threshold(df_sample, step=0.05)
        ax.plot(purity_jsd_df["purity_threshold"], purity_jsd_df["similarity"], linewidth=3)
        ax.axvline(purity_elbow, linestyle="--", linewidth=2)
        ax.text(0.95, 0.05, f"Elbow point = {purity_elbow:.2f}", transform=ax.transAxes, ha="right", va="bottom", fontsize=15)
        ax.set_title(sample_name, fontsize=25, pad=10)
        if i == 7: ax.set_xlabel("Purity threshold", fontsize=15)
        if j == 0: ax.set_ylabel("Jensen–Shannon similarity", fontsize=15)

plt.tight_layout()
plt.savefig(os.path.join(output_curve_dir, "species_similarity_elbow_cell_purity.pdf"), dpi=300)
plt.close()

print("Analysis finished.")