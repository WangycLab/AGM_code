# ===============================================
# Author : ZHENG XINGHAI
# Date   : 2025-11-30
# Project: Astronaut mscRNA-seq Analysis
# ===============================================

library(tidyverse)  
library(ggalluvial)  
library(RColorBrewer)  
library(patchwork)  

# Set working directory to the script location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  

# Load species abundance data
df_raw <- read.csv("species_abundance_profile.csv", row.names = 1, check.names = FALSE)  
df_raw$Species <- rownames(df_raw)  

# Define timepoint columns and compute mean abundance
all_time_cols <- colnames(df_raw)[1:24]  
df_raw$mean_abundance <- rowMeans(df_raw[, all_time_cols], na.rm = TRUE)  

# Select top N species based on mean abundance
topN <- 25  
top_species_ordered <- df_raw %>%
  arrange(desc(mean_abundance)) %>%
  slice(1:topN) %>%
  pull(Species)  

# Aggregate top species and merge remaining as Others
species_df <- df_raw %>%
  mutate(Species_group = ifelse(Species %in% top_species_ordered, Species, "Others")) %>%
  group_by(Species_group) %>%
  summarise(across(all_of(all_time_cols), sum, na.rm = TRUE), .groups = "drop")  

# Set factor levels for plotting
species_df$Species_group <- factor(species_df$Species_group, levels = c(top_species_ordered, "Others"))  

# Define color palette for species groups
palette_raw <- c(brewer.pal(12, "Paired"), brewer.pal(12, "Set3"), brewer.pal(8, "Dark2"), brewer.pal(9, "Set1"))  
species_colors <- c(setNames(palette_raw[1:topN], top_species_ordered), Others = "gray70")  

# Function to create an alluvial plot
make_alluvial <- function(species_df, time_cols, time_labels, species_colors, show_x_axis = TRUE, show_legend = TRUE, subtitle = NULL) {
  df_long <- species_df %>%
    select(Species_group, all_of(time_cols)) %>%
    pivot_longer(cols = all_of(time_cols), names_to = "Time", values_to = "Abundance") %>%
    mutate(alluvium = Species_group, Time = factor(Time, levels = time_cols),
           Species_group = factor(Species_group, levels = names(species_colors)))  
  
  p <- ggplot(df_long, aes(x = Time, stratum = Species_group, alluvium = alluvium, y = Abundance, fill = Species_group)) +
    geom_flow(stat = "alluvium", lode.guidance = "forward", color = "gray40", alpha = 0.5) +
    geom_stratum(color = "black", size = 0.3, alpha = 0.9) +
    scale_fill_manual(values = species_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.001, 0.02))) +
    scale_x_discrete(expand = expansion(add = c(0.01, 0.01)), labels = time_labels) +
    theme_minimal(base_size = 15) +
    theme(legend.position = if (show_legend) "right" else "none", legend.text = element_text(face = "italic", size = 12),
          legend.title = element_text(size = 16), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 16),
          panel.grid = element_blank()) +
    guides(fill = guide_legend(title = "Species", ncol = 1, byrow = TRUE))  
  
  if (!is.null(subtitle)) p <- p + labs(subtitle = subtitle)  
  if (show_x_axis) p <- p + xlab("Timepoints") + theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 16)) else
    p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank())  
  
  return(p)
}

# Define timepoint labels
time_labels <- c("L-60", "L-30", "FD30", "FD90", "FD150", "R+1", "R+7", "R+14")  

# Generate alluvial plots for each astronaut
p1 <- make_alluvial(species_df, time_cols = colnames(df_raw)[1:8], time_labels = time_labels, species_colors = species_colors,
                    show_x_axis = FALSE, show_legend = FALSE, subtitle = "Astronaut 1")  
p2 <- make_alluvial(species_df, time_cols = colnames(df_raw)[9:16], time_labels = time_labels, species_colors = species_colors,
                    show_x_axis = FALSE, show_legend = FALSE, subtitle = "Astronaut 2")  
p3 <- make_alluvial(species_df, time_cols = colnames(df_raw)[17:24], time_labels = time_labels, species_colors = species_colors,
                    show_x_axis = TRUE, show_legend = TRUE, subtitle = "Astronaut 3")  

# Combine plots into a single figure
combined_plot <- p1 / p2 / p3 + plot_layout(ncol = 1, guides = "collect")  

# Display the combined figure
combined_plot  

# Save the figure as an editable PDF
ggsave(filename = "Figure3e.pdf", plot = combined_plot, device = cairo_pdf, width = 7.5, height = 8, units = "in")  