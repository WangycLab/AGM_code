# ===============================================
# Author : ZHENG XINGHAI
# Date   : 2025-11-30
# Project: Astronaut mscRNA-seq Analysis
# ===============================================

library(tidyverse)
library(ggalluvial)
library(RColorBrewer)

# Use the current working directory and load the percentage matrix
setwd(getwd())
df_raw <- read.csv("timepoint_cluster_percentage.csv", row.names = 1, check.names = FALSE)

# Define a color palette for clusters
palette_set3 <- brewer.pal(8, "Set3")

# Create a reusable function to generate an alluvial plot
make_alluvial <- function(
    species_df,
    time_cols,
    time_labels,
    species_colors,
    show_x_axis = TRUE,
    show_legend = TRUE
) {
  
  # Convert the wide matrix into long format for ggplot
  df_long <- species_df %>%
    rownames_to_column(var = "Cluster") %>%
    pivot_longer(
      cols = -Cluster,
      names_to = "Time",
      values_to = "Percentage"
    ) %>%
    mutate(
      Time = factor(Time, levels = time_cols),
      Cluster = factor(Cluster, levels = 0:7)
    )
  
  # Build the alluvial visualization
  p <- ggplot(
    df_long,
    aes(
      x = Time,
      stratum = Cluster,
      alluvium = Cluster,
      y = Percentage,
      fill = Cluster
    )
  ) +
    geom_flow(stat = "alluvium", lode.guidance = "forward", color = "gray40", alpha = 0.3) +
    geom_stratum(color = "black", size = 0.3, alpha = 0.8) +
    scale_fill_manual(values = species_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.001, 0.02))) +
    scale_x_discrete(expand = expansion(add = c(0.01, 0.01)), labels = time_labels) +
    theme_minimal(base_size = 15) +
    theme(
      legend.position = if (show_legend) "right" else "none",
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 20),
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(size = 20),
      panel.grid = element_blank()
    ) +
    guides(fill = guide_legend(title = "Cluster", ncol = 1, byrow = TRUE))
  
  # Control the display of x-axis labels
  if (show_x_axis) {
    p <- p +
      xlab("Timepoints") +
      theme(
        axis.text.x = element_text(size = 15),
        axis.title.x = element_text(size = 20)
      )
  } else {
    p <- p +
      theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank()
      )
  }
  
  return(p)
}

# Define the labels corresponding to each timepoint
time_labels <- c("L-60", "L-30", "FD30", "FD90", "FD150", "R+1", "R+7", "R+14")

# Generate the alluvial plot using the full matrix
p1 <- make_alluvial(
  df_raw,
  time_cols = colnames(df_raw),
  time_labels = time_labels,
  species_colors = palette_set3,
  show_x_axis = TRUE,
  show_legend = TRUE
)

# Save the resulting plot as a PDF file
ggsave("Figure5b.pdf", plot = p1, width = 8, height = 5)