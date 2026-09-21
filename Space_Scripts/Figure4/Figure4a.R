# ===============================================
# Author : ZHENG XINGHAI
# Date   : 2025-11-30
# Project: Astronaut mscRNA-seq Analysis
# ===============================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)
library(tidyverse)
library(tidyr)
library(ggrepel)
library(glmGamPoi)
library(RColorBrewer)
library(ggsci)
library(viridis)
library(Matrix)
library(monocle3)
library(SeuratWrappers)
library(patchwork)
library(irlba)

# Set working and project directories
current_dir <- trimws(getwd())
setwd(current_dir)
cat("Working directory:", current_dir, "\n")

parts <- strsplit(current_dir, .Platform$file.sep)[[1]]
project_root <- if (length(parts) >= 2 && all(tail(parts, 2) == c("Space_Scripts", "Figure4"))) {
  paste(parts[1:(length(parts) - 2)], collapse = .Platform$file.sep)
} else current_dir

data_path <- file.path(project_root, "Space_Matrix")
cat("Data path:", data_path, "\n")

# Sample order and astronaut colors
x_labels <- c("L-60", "L-30", "FD30", "FD90", "FD150", "R+1", "R+7", "R+14")

astronaut_color_map <- c(
  "Astronaut 1" = "#FF9A8B",
  "Astronaut 2" = "#FFE0AC",
  "Astronaut 3" = "#BBDED6"
)

# Read all samples and calculate the number of detected genes per cell
sample_list <- list.dirs(data_path, recursive = FALSE, full.names = FALSE)
plot_df <- data.frame()

for (sample_name in sample_list) {
  cat("Processing:", sample_name, "\n")
  sample_dir <- file.path(data_path, sample_name)
  expr_mat <- tryCatch(Read10X(sample_dir), error = function(e) NULL)
  
  if (is.null(expr_mat)) {
    cat("Failed:", sample_name, "\n")
    next
  }
  
  detected_genes <- Matrix::colSums(expr_mat > 0)
  
  # Extract astronaut identity from the sample suffix
  astronaut <- case_when(
    grepl("_1$", sample_name) ~ "Astronaut 1",
    grepl("_2$", sample_name) ~ "Astronaut 2",
    grepl("_3$", sample_name) ~ "Astronaut 3",
    TRUE ~ "Unknown"
  )
  
  sample_clean <- gsub("_[123]$", "", sample_name)
  plot_df <- rbind(plot_df, data.frame(sample = sample_clean, astronaut = astronaut, detected_genes = detected_genes))
}

# Filter cells and set plotting order
plot_df <- plot_df %>%
  filter(detected_genes > 0) %>%
  mutate(
    sample = factor(sample, levels = x_labels),
    astronaut = factor(astronaut, levels = c("Astronaut 1", "Astronaut 2", "Astronaut 3"))
  )

# Calculate median gene counts for labels
med_df <- plot_df %>%
  group_by(sample, astronaut) %>%
  summarise(med = median(detected_genes), .groups = "drop")

# Plot gene counts across samples
p1 <- ggplot(plot_df, aes(x = sample, y = detected_genes, fill = astronaut)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.2, color = NA, alpha = 0.95) +
  geom_boxplot(width = 0.10, outlier.shape = NA, linewidth = 0.25, fill = "white", color = "black") +
  stat_summary(fun = median, geom = "point", size = 1.2, color = "black") +
  geom_text(
    data = med_df,
    aes(x = sample, y = med + 20, label = round(med, 1)),
    inherit.aes = FALSE, size = 5, color = "black"
  ) +
  facet_wrap(~astronaut, ncol = 1) +
  scale_fill_manual(values = astronaut_color_map) +
  scale_y_continuous(breaks = seq(0, 300, 50), limits = c(0, 300), expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = "Sample", y = "Gene counts") +
  theme_classic(base_size = 16) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 16),
    panel.spacing = unit(0.3, "cm"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(face = "bold", size = 16),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    legend.position = "none"
  )

print(p1)

# Save figure
ggsave("Figure4a.pdf", p1, width = 8, height = 8)