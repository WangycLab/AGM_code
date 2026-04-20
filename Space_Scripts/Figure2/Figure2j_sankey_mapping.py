# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 14:19:46 2025
@author: ZHENG XINGHAI
"""

import pandas as pd
import plotly.graph_objects as go
from matplotlib.colors import to_rgba
import matplotlib as mpl

# Configure Matplotlib to generate PDFs with editable text
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

# Set input and output file paths
input_file = "cell_metadata.tsv"
output_pdf = "Figure2j.pdf"

# Load the cell metadata table for plotting
df = pd.read_csv(input_file, sep="\t")

# Extract unique clusters and samples
clusters = df['cluster'].astype(str).unique().tolist()
samples = df['sample'].astype(str).unique().tolist()
nodes = clusters + samples  # Combine clusters and samples into nodes

# Map node names to indices for Sankey plotting
node_indices = {name: i for i, name in enumerate(nodes)}

# Compute the number of cells per cluster → sample combination
link_df = df.groupby(['cluster', 'sample']).size().reset_index(name='count')
link_df['cluster'] = link_df['cluster'].astype(str)
link_df['sample'] = link_df['sample'].astype(str)
source_indices = link_df['cluster'].map(node_indices)
target_indices = link_df['sample'].map(node_indices)
values = link_df['count']

# Define custom colors for clusters
cluster_color_map = {
    "0": "#A6CEE3", "1": "#438EC0", "2": "#63A8A0", "3": "#98D277",
    "4": "#3BA432", "5": "#B89B74", "6": "#F16667", "7": "#E62F27",
    "8": "#F9A963", "9": "#FE982C", "10": "#ED8F47", "11": "#C3AAD2",
    "12": "#7D54A5", "13": "#B9A499", "14": "#EAD27A", "15": "#B15928"
}

# Define custom colors for samples
sample_color_map = {
    "Con": "#8DD3C7", "2mo": "#E69394", "4mo": "#E4C264", "6mo": "#D9D9D9"
}

# Assign colors to nodes: clusters first, then samples
node_colors = [cluster_color_map.get(c, "#CCCCCC") for c in clusters] + \
              [sample_color_map.get(s, "lightgray") for s in samples]

# Assign semi-transparent colors to links based on cluster colors
link_colors = [
    f'rgba({int(r*255)},{int(g*255)},{int(b*255)},0.6)'
    for r, g, b, _ in [to_rgba(cluster_color_map.get(str(c), "#CCCCCC")) for c in link_df['cluster']]
]

# Create Sankey diagram
fig = go.Figure(data=[go.Sankey(
    arrangement="snap",
    node=dict(
        pad=25,
        thickness=20,
        line=dict(color="gray", width=0.1),
        label=nodes,  # Node labels are editable in Illustrator
        color=node_colors
    ),
    link=dict(
        source=source_indices,
        target=target_indices,
        value=values,
        color=link_colors
    )
)])

# Update layout for publication-quality figure
fig.update_layout(
    title_text="Cluster → Sample",
    title_font_size=80,
    title_x=0.5,
    font_size=80,
    font_family="Arial",
    margin=dict(t=150)
)

# Save figure as PDF
try:
    fig.write_image(output_pdf, width=800, height=2400)
    print(f"PDF successfully saved to: {output_pdf}")
except Exception as e:
    print("Error saving PDF. Make sure 'kaleido' is installed.")
    print(e)

# Display interactive figure in browser
fig.show()