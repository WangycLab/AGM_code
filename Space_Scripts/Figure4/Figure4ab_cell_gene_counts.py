# -*- coding: utf-8 -*-
"""
Created on Mon Sep 15 15:15:18 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as mticker

# Configure matplotlib so that PDF text remains editable
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Load the cell metadata table
df_meta = pd.read_csv("cell_metadata.tsv", sep="\t")

# Define the desired order of sampling timepoints
x_labels = ["L-60","L-30","FD30","FD90","FD150","R+1","R+7","R+14"]

# Extract timepoint labels and sort samples according to the defined order
def sort_samples(df):
    order_dict = {label: i for i, label in enumerate(x_labels)}
    df['timepoint_simple'] = df['sample'].str.extract(r'([A-Z0-9\+\-]+)_')[0]
    df = df[df['timepoint_simple'].isin(x_labels)]
    df = df.sort_values(by='timepoint_simple', key=lambda x: x.map(order_dict))
    return df

# Split metadata into three astronauts and apply sorting
astro01 = sort_samples(df_meta[df_meta['astronaut']=="Astronaut 1"])
astro02 = sort_samples(df_meta[df_meta['astronaut']=="Astronaut 2"])
astro03 = sort_samples(df_meta[df_meta['astronaut']=="Astronaut 3"])

# Define consistent colors for each astronaut
astronaut_color_map = {
    "Astronaut 1": "#FF9A8B",
    "Astronaut 2": "#FFE0AC",
    "Astronaut 3": "#BBDED6"
}

# Set y-axis ranges for gene count visualization
y_limits = {
    "Astronaut 1": (-10, 200),
    "Astronaut 2": (-10, 200),
    "Astronaut 3": (-10, 200)
}

# Apply seaborn visualization style and global font
sns.set_theme(style="whitegrid", context="talk")
plt.rcParams['font.family'] = 'Arial'

# Create stacked subplots for gene count distributions
fig, axes = plt.subplots(3, 1, figsize=(12, 18), sharex=True)

# Draw violin plots with overlaid boxplots to show gene count distributions
def draw_violin(ax, df, astronaut):
    violin_color = astronaut_color_map[astronaut]
    y_min, y_max = y_limits[astronaut]

    sns.violinplot(
        x='timepoint_simple', y='gene_number',
        data=df,
        order=x_labels,
        inner=None,
        linewidth=0.1,
        width=1,
        color=violin_color,
        ax=ax
    )

    sns.boxplot(
        x='timepoint_simple', y='gene_number',
        data=df,
        order=x_labels,
        width=0.08,
        showcaps=True,
        showfliers=False,
        linewidth=0.1,
        boxprops={
            'facecolor': 'white',
            'edgecolor': 'black',
            'linewidth': 0.1,
            'alpha': 0.8
        },
        whiskerprops={'color': 'black', 'linewidth': 0.1},
        capprops={'color': 'black', 'linewidth': 0.1},
        medianprops={'color': 'red', 'linewidth': 0.1},
        ax=ax
    )

    ax.set_ylabel("Gene counts", fontsize=40)
    ax.set_ylim(y_min, y_max)
    ax.yaxis.set_major_locator(mticker.MultipleLocator(30))
    ax.tick_params(axis="y", labelsize=30)

    ax.set_xlabel("")
    ax.set_xticklabels([])

    ax.text(0.01, 0.95, astronaut, transform=ax.transAxes,
            fontsize=40, va='top', ha='left',
            bbox=dict(facecolor='white', alpha=0.5, edgecolor='none', pad=2))

    ax.set_facecolor('#f9f9f9')

    ax.yaxis.grid(True, color='gray', linestyle='--', linewidth=0.1, alpha=0.5)
    ax.xaxis.grid(False)

    sns.despine(ax=ax, top=True, right=True)

# Generate violin plots for the three astronauts
draw_violin(axes[0], astro01, "Astronaut 1")
draw_violin(axes[1], astro02, "Astronaut 2")
draw_violin(axes[2], astro03, "Astronaut 3")

axes[2].set_xlabel("Sample", fontsize=40)
axes[2].set_xticklabels(x_labels, fontsize=30)

fig.patch.set_facecolor('#f9f9f9')
plt.tight_layout()
plt.savefig("Figure4b.pdf", dpi=300)
plt.show()

# Calculate the number of cells detected at each timepoint
def count_cells(df):
    count_df = (
        df.groupby('timepoint_simple')
          .size()
          .reindex(x_labels)
          .reset_index()
    )
    count_df.columns = ['timepoint_simple', 'cell_number']
    return count_df

astro01_count = count_cells(astro01)
astro02_count = count_cells(astro02)
astro03_count = count_cells(astro03)

# Create stacked subplots for cell count comparisons
fig_bar, axes_bar = plt.subplots(3, 1, figsize=(12, 18), sharex=True)

# Draw bar plots showing the number of cells per timepoint
def draw_bar(ax, count_df, astronaut):
    bar_color = astronaut_color_map[astronaut]

    sns.barplot(
        x='timepoint_simple',
        y='cell_number',
        data=count_df,
        order=x_labels,
        color=bar_color,
        width=0.3,
        edgecolor='black',
        linewidth=0.1,
        ax=ax
    )

    ax.set_ylim(0, 25000)
    ax.yaxis.set_major_locator(mticker.MultipleLocator(5000))

    ax.set_ylabel("Cell counts", fontsize=40)
    ax.tick_params(axis="y", labelsize=30)
    ax.set_xlabel("")
    ax.set_xticklabels([])

    ax.text(0.01, 0.95, astronaut, transform=ax.transAxes,
            fontsize=40, va='top', ha='left',
            bbox=dict(facecolor='white', alpha=0.5, edgecolor='none', pad=2))

    ax.set_facecolor('#f9f9f9')
    ax.yaxis.grid(True, linestyle='--', linewidth=0.1, alpha=0.5)
    ax.xaxis.grid(False)

    sns.despine(ax=ax, top=True, right=True)

# Generate bar plots for the three astronauts
draw_bar(axes_bar[0], astro01_count, "Astronaut 1")
draw_bar(axes_bar[1], astro02_count, "Astronaut 2")
draw_bar(axes_bar[2], astro03_count, "Astronaut 3")

axes_bar[2].set_xlabel("Sample", fontsize=40)
axes_bar[2].set_xticklabels(x_labels, fontsize=30)

fig_bar.patch.set_facecolor('#f9f9f9')
plt.tight_layout()
plt.savefig("Figure4a.pdf", dpi=300)
plt.show()