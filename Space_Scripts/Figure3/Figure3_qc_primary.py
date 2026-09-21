# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 15:46:16 2025
@author: ZHENG XINGHAI
"""

import os
import glob
import pandas as pd
from scipy.stats import ttest_1samp

# Set directories
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(BASE_DIR, "space_species_filted")
INPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(BASE_DIR)), "Space_Matrix")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Process taxonomy reports
summary_list = []
for filepath in glob.glob(os.path.join(INPUT_DIR, "*_sc_taxonomy.report")):
    filename = os.path.basename(filepath)
    parts = os.path.splitext(filename)[0].split("_")
    prefix = "_".join(parts[:2])
    df = pd.read_csv(filepath, sep="\t")
    n_initial = len(df)

    def test_significant(row):
        _, p = ttest_1samp([row["fraction_total_reads2"], row["fraction_total_reads3"]], row["fraction_total_reads"], alternative="less")
        return p <= 0.05

    df["significant"] = df.apply(test_significant, axis=1)
    df_significant = df[df["significant"]].sort_values("new_est_reads", ascending=False)
    n_after_primary = len(df_significant)
    output_path = os.path.join(OUTPUT_DIR, f"{prefix}_sc_taxonomy_filted.csv")
    df_significant.to_csv(output_path, index=False)
    summary_list.append({"Sample": prefix, "Initial": n_initial, "After_primary": n_after_primary})
    print(f"{prefix}: {n_initial} -> {n_after_primary}")

# Build ordered summary
summary_df = pd.DataFrame(summary_list)
samples = [f"{t}{batch}" for batch in ["_1", "_2", "_3"] for t in ["L-60", "L-30", "FD30", "FD90", "FD150", "R+1", "R+7", "R+14"]]
summary_df["Sample"] = pd.Categorical(summary_df["Sample"], categories=samples, ordered=True)
summary_df = summary_df.sort_values("Sample")

# Save summary
summary_csv_path = os.path.join(BASE_DIR, "summary_cell_counts.csv")
summary_df.to_csv(summary_csv_path, index=False)
print(f"Processed {len(summary_df)} samples; summary saved to {summary_csv_path}")