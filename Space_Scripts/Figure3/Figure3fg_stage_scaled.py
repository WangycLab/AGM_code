# -*- coding: utf-8 -*-
"""
Created on Mon Dec 29 23:44:07 2025
@author: ZHENG XINGHAI
"""

import pandas as pd

# Load the original species abundance profile table
df = pd.read_csv("species_abundance_profile.csv", index_col=0)

# Define stage-level sample groups corresponding to each astronaut replicate
new_columns = {
    "L_1": ["L-60_1", "L-30_1"],
    "FD_1": ["FD30_1", "FD90_1", "FD150_1"],
    "R_1": ["R+1_1", "R+7_1", "R+14_1"],
    "L_2": ["L-60_2", "L-30_2"],
    "FD_2": ["FD30_2", "FD90_2", "FD150_2"],
    "R_2": ["R+1_2", "R+7_2", "R+14_2"],
    "L_3": ["L-60_3", "L-30_3"],
    "FD_3": ["FD30_3", "FD90_3", "FD150_3"],
    "R_3": ["R+1_3", "R+7_3", "R+14_3"]
}

# Compute stage-level average abundance values for each species
df_new = pd.DataFrame(index=df.index)
for new_col, old_cols in new_columns.items():
    df_new[new_col] = df[old_cols].mean(axis=1)

# Save the stage-averaged species abundance table
df_new.to_csv("species_abundance_profile_stage.csv")

# Define astronaut-specific stage groups for scaling
groups = [
    ["L_1", "FD_1", "R_1"],
    ["L_2", "FD_2", "R_2"],
    ["L_3", "FD_3", "R_3"]
]

# Initialize a dataframe to store the scaled abundance values
df_scaled = pd.DataFrame(index=df_new.index)

# Apply min-max scaling for each species within each astronaut stage group
for grp in groups:
    sub_df = df_new[grp].copy()
    row_min = sub_df.min(axis=1)
    row_max = sub_df.max(axis=1)
    row_range = row_max - row_min
    row_range_replaced = row_range.replace(0, 1)
    sub_scaled = sub_df.sub(row_min, axis=0).div(row_range_replaced, axis=0)
    sub_scaled[row_range == 0] = 0
    df_scaled[grp] = sub_scaled

# Save the scaled species abundance table
df_scaled.to_csv("species_abundance_profile_scaled.csv")