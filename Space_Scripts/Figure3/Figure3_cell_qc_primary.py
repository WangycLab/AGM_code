# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 15:46:16 2025
@author: ZHENG XINGHAI
"""

import os
import glob
import pandas as pd
from scipy.stats import ttest_1samp

# Set project directories
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(BASE_DIR, "space_species_filted")
INPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(BASE_DIR)), "Space_Matrix")

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Store cell count statistics
summary_list = []

# Process all taxonomy report files
for filepath in glob.glob(os.path.join(INPUT_DIR, "*_sc_taxonomy.report")):
    filename = os.path.basename(filepath)
    base = os.path.splitext(filename)[0]
    parts = base.split("_")
    prefix = parts[0] + "_" + parts[1]

    print(f"Processing {prefix}...")

    # Load taxonomy table
    df = pd.read_csv(filepath, sep="\t")
    n_initial = df.shape[0]

    # Identify statistically supported assignments
    def test_significant(row):
        other_values = [
            row["fraction_total_reads2"],
            row["fraction_total_reads3"]
        ]
        _, p_value = ttest_1samp(
            other_values,
            row["fraction_total_reads"],
            alternative="less"
        )
        return p_value < 0.05

    # Filter significant rows
    df["significant"] = df.apply(test_significant, axis=1)
    df_significant = (
        df[df["significant"]]
        .sort_values(by="new_est_reads", ascending=False)
    )

    n_after_primary = df_significant.shape[0]

    # Save filtered results
    output_path = os.path.join(
        OUTPUT_DIR,
        f"{prefix}_sc_taxonomy_filted.csv"
    )
    df_significant.to_csv(output_path, index=False)

    print(
        f"{filename} -> Initial: {n_initial}, "
        f"After primary: {n_after_primary}"
    )
    print(f"Saved results to {output_path}")

    # Record sample statistics
    summary_list.append({
        "Sample": prefix,
        "Initial": n_initial,
        "After_primary": n_after_primary
    })

# Build summary table
summary_df = pd.DataFrame(summary_list)

# Define sample order
samples = []
for batch in ["_1", "_2", "_3"]:
    for t in ["L-60", "L-30", "FD30", "FD90", "FD150", "R+1", "R+7", "R+14"]:
        samples.append(f"{t}{batch}")

# Apply ordered sample labels
summary_df["Sample"] = pd.Categorical(
    summary_df["Sample"],
    categories=samples,
    ordered=True
)

# Sort samples
summary_df = summary_df.sort_values("Sample")

# Export summary
summary_csv_path = os.path.join(BASE_DIR, "summary_cell_counts.csv")
summary_df.to_csv(summary_csv_path, index=False)

print(f"Summary saved to {summary_csv_path}")