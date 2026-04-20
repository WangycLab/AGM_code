# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 15:46:16 2025
@author: ZHENG XINGHAI
"""

import os
import glob
import pandas as pd
from scipy.stats import ttest_1samp

# Set up project directories for input and output
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(BASE_DIR, "store_species_filted")

# Dynamically set INPUT_DIR
INPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(BASE_DIR)), "Store_Matrix")
print(INPUT_DIR)

# Create output directory if it does not exist
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Initialize list to store summary statistics of cell counts
summary_list = []

# Iterate through all taxonomy report files in the input folder
for filepath in glob.glob(os.path.join(INPUT_DIR, "*_sc_taxonomy.report")):

    filename = os.path.basename(filepath)
    base = os.path.splitext(filename)[0]
    parts = base.split("_")
    prefix = "_".join(parts[:1]) if len(parts) >= 2 else parts[0]
    print(f"Processing {prefix}...")

    # Load taxonomy data and record the initial number of cells
    df = pd.read_csv(filepath, sep="\t")
    n_initial = df.shape[0]

    # Filter cells with fraction_total_reads >= 0.5
    df_filtered = df[df["fraction_total_reads"] >= 0.5]
    n_after_fraction = df_filtered.shape[0]

    # Define a one-sample t-test function to check enrichment significance
    def test_significant(row):
        other_values = [row["fraction_total_reads2"], row["fraction_total_reads3"]]
        _, p_value = ttest_1samp(other_values, row["fraction_total_reads"], alternative="less")
        return p_value < 0.05

    # Apply the significance test and keep only significant cells
    df_filtered["significant"] = df_filtered.apply(test_significant, axis=1)
    df_significant = df_filtered[df_filtered["significant"]].sort_values(
        by="new_est_reads",
        ascending=False
    )
    n_after_significant = df_significant.shape[0]

    # Save the filtered significant cells for each sample
    output_path = os.path.join(OUTPUT_DIR, f"{prefix}_sc_taxonomy_filted.csv")
    df_significant.to_csv(output_path, index=False)

    print(
        f"{filename} -> Initial: {n_initial}, "
        f"After fraction>=0.5: {n_after_fraction}, "
        f"After significance: {n_after_significant}"
    )
    print(f"Saved results to {output_path}")

    # Append statistics for summary table
    summary_list.append({
        "Sample": prefix,
        "Initial": n_initial,
        "After_primary": n_after_significant
    })

# Generate summary table of cell counts and save as CSV
summary_df = pd.DataFrame(summary_list)
summary_csv_path = os.path.join(BASE_DIR, "summary_cell_counts.csv")
summary_df.to_csv(summary_csv_path, index=False)

print(f"Summary saved to {summary_csv_path}")