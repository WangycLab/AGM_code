# -*- coding: utf-8 -*-
"""
Created on Wed Jan 7 19:01:32 2026
@author: ZHENG XINGHAI
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

# Load and merge metadata
metadata = pd.read_csv("cell_metadata.tsv", sep="\t")
pseudotime = pd.read_csv("species_pseudotime.csv").rename(columns={"cell": "barcode"})
metadata = metadata.merge(pseudotime[["barcode", "pseudotime"]], on="barcode", how="left")
metadata.to_csv("species_metadata_pseudotime.tsv", sep="\t", index=False)

# Clean data
df = metadata.copy()
df["pseudotime"] = pd.to_numeric(df["pseudotime"], errors="coerce")
df = df.dropna(subset=["pseudotime", "timepoint"]).copy()
df["timepoint"] = df["timepoint"].astype(str)

timepoint_order = ["L-60", "L-30", "FD30", "FD90", "FD150", "R+1", "R+7", "R+14"]
df = df[df["timepoint"].isin(timepoint_order)].copy()

time_colors = {
    "L-60": "#F4F1DE", "L-30": "#E6D8AD", "FD30": "#5F6F75", "FD90": "#81B295",
    "FD150": "#A3BFA8", "R+1": "#9A031E", "R+7": "#BA3A26", "R+14": "#E07A5F"
}

# Divide pseudotime into 10 equal-width bins
n_bins = 10
min_pseudotime, max_pseudotime = df["pseudotime"].min(), df["pseudotime"].max()
bins = np.linspace(min_pseudotime, max_pseudotime, n_bins + 1)
df["pseudotime_bin"] = pd.cut(df["pseudotime"], bins=bins, labels=False, include_lowest=True) + 1

# Calculate cell proportions
cell_counts = df.groupby(["pseudotime_bin", "timepoint"], observed=False).size().reset_index(name="cell_count")
all_combinations = pd.MultiIndex.from_product([range(1, n_bins + 1), timepoint_order], names=["pseudotime_bin", "timepoint"]).to_frame(index=False)
cell_counts = all_combinations.merge(cell_counts, on=["pseudotime_bin", "timepoint"], how="left")
cell_counts["cell_count"] = cell_counts["cell_count"].fillna(0)
cell_counts["proportion"] = cell_counts["cell_count"] / cell_counts.groupby("pseudotime_bin")["cell_count"].transform("sum")
cell_counts["bin_start"] = cell_counts["pseudotime_bin"].map(lambda x: bins[x - 1])
cell_counts["bin_end"] = cell_counts["pseudotime_bin"].map(lambda x: bins[x])
cell_counts["bin_midpoint"] = (cell_counts["bin_start"] + cell_counts["bin_end"]) / 2

# Calculate Spearman correlations
correlation_df = pd.DataFrame([
    {"timepoint": tp, "Spearman_rho": spearmanr(sub["bin_midpoint"], sub["proportion"]).statistic, "p_value": spearmanr(sub["bin_midpoint"], sub["proportion"]).pvalue}
    for tp in timepoint_order
    for sub in [cell_counts[cell_counts["timepoint"] == tp].sort_values("pseudotime_bin")]
])
correlation_df.to_csv("pseudotime_timepoint_spearman.csv", index=False)

# Prepare plotting data
plot_df = cell_counts.pivot(index="pseudotime_bin", columns="timepoint", values="proportion").fillna(0)[timepoint_order]

# Plot settings
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "axes.unicode_minus": False})

fig, ax = plt.subplots(figsize=(10, 6))
bottom = np.zeros(n_bins)

for tp in timepoint_order:
    values = plot_df[tp].values
    ax.bar(np.arange(1, n_bins + 1), values, bottom=bottom, width=0.95, color=time_colors[tp], edgecolor="none", label=tp)
    bottom += values

# Format axes
ax.set_xlabel("Pseudotime", fontsize=20)
ax.set_xticks(np.arange(1, n_bins + 1))
ax.set_xticklabels([f"{bins[i]:.1f}–{bins[i + 1]:.1f}" for i in range(n_bins)], rotation=45, ha="right", fontsize=15)
ax.set_ylabel("Proportion of cells", fontsize=20)
ax.set_ylim(0, 1)
ax.set_yticks(np.arange(0, 1.01, 0.2))
ax.set_yticklabels(["0%", "20%", "40%", "60%", "80%", "100%"])
ax.tick_params(axis="both", labelsize=15)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.yaxis.grid(True, linestyle="--", linewidth=0.6, alpha=0.3)
ax.set_axisbelow(True)
ax.legend(title="Timepoint", fontsize=13, title_fontsize=15, frameon=False, bbox_to_anchor=(1, 1), loc="upper left")

plt.tight_layout()
plt.savefig("Figure5c.pdf", bbox_inches="tight")
plt.show()

print(f"Finished: {len(metadata)} cells, {metadata['pseudotime'].notna().sum()} matched pseudotime values; output saved.")