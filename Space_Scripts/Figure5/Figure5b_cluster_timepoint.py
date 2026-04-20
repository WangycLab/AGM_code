# -*- coding: utf-8 -*-
"""
Created on Mon Jan 5 23:25:54 2026
@author: ZHENG XINGHAI
"""

import pandas as pd

# Load the metadata table containing cell annotations
data = pd.read_csv("Phocaeicola_vulgatus_metadata.tsv", sep="\t")

# Define the biological order of timepoints and convert the column to an ordered categorical variable
timepoint_order = ["L-60", "L-30", "FD30", "FD90", "FD150", "R+1", "R+7", "R+14"]
data['timepoint'] = pd.Categorical(
    data['timepoint'],
    categories=timepoint_order,
    ordered=True
)

# Count the number of cells for each combination of timepoint and cluster
grouped = (
    data
    .groupby(['timepoint', 'cluster'])
    .size()
    .reset_index(name='count')
)

# Calculate the total number of cells within each timepoint
total_per_timepoint = grouped.groupby('timepoint')['count'].transform('sum')

# Compute the relative abundance of each cluster within each timepoint
grouped['percentage'] = grouped['count'] / total_per_timepoint

# Round percentage values for consistent numerical formatting
grouped['percentage'] = grouped['percentage'].round(5)

# Reshape the table into a matrix with clusters as rows and timepoints as columns
percentage_matrix = (
    grouped
    .pivot(index='cluster', columns='timepoint', values='percentage')
    .fillna(0)
)

# Print the resulting percentage matrix for quick inspection
print(percentage_matrix)

# Save the cluster-by-timepoint percentage matrix to a CSV file
percentage_matrix.to_csv("timepoint_cluster_percentage.csv", float_format='%.5f')