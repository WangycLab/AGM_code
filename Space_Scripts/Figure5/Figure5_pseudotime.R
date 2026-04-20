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
library(hdWGCNA)
library(cowplot)
library(patchwork)
library(WGCNA)

# Set the working directory to the current folder
current_dir <- trimws(getwd())
setwd(current_dir)
cat("Working directory set to:", current_dir, "\n")

# Determine the project root directory based on folder structure
parts <- strsplit(current_dir, .Platform$file.sep)[[1]]
if (length(parts) >= 2 && all(tail(parts, 2) == c("Space_Scripts", "Figure5"))) {
  project_root <- paste(parts[1:(length(parts)-2)], collapse = .Platform$file.sep)
} else {
  project_root <- current_dir
}
cat("Project root:", project_root, "\n")

# Define the path to the data directory relative to the project root
data_path <- file.path(project_root, "Space_Matrix")
cat("Data path set to:", data_path, "\n")

# Set a global cutoff for cell purity
GLOBAL_PURITY_CUTOFF <- 0.5

# Retrieve the list of sample directories
name_list <- list.dirs(data_path, full.names = FALSE, recursive = FALSE)
seurat_list <- list()

# Initialize containers for summary statistics and expression counts
initial_stats <- data.frame(Sample = character(), Genes = integer(), Cells = integer())
after_filter_stats <- data.frame(Sample = character(), Genes = integer(), Cells = integer())
after_seurat_stats <- data.frame(Sample = character(), Genes = integer(), Cells = integer())
expr_count_raw <- list()
expr_count_filtered <- list()
gene_number_list <- list()

# Define a t-test-based filtering function for cell enrichment
pass_ttest <- function(fraction_total_reads, fraction_total_reads2, fraction_total_reads3) {
  values <- as.numeric(c(fraction_total_reads2, fraction_total_reads3))
  fraction_total_reads <- as.numeric(fraction_total_reads)
  values <- values[!is.na(values)]
  if (length(values) < 2 || var(values) == 0 || is.na(fraction_total_reads)) {
    return(FALSE)
  }
  test <- t.test(values, mu = fraction_total_reads)
  test$p.value < 0.05 && fraction_total_reads > max(values)
}

# Main loop to process each sample
for (sample_name in name_list) {
  cat("Processing sample:", sample_name, "\n")
  
  # Define file paths
  report_S_file <- file.path(data_path, paste0(sample_name, "_sc_taxonomy_greedy_S.report"))
  report_G_file <- file.path(data_path, paste0(sample_name, "_sc_taxonomy_greedy_G.report"))
  taxonomy_file <- file.path(data_path, paste0(sample_name, "_sc_taxonomy.report"))
  if (!file.exists(report_S_file) || !file.exists(report_G_file) || !file.exists(taxonomy_file)) next
  
  # Read taxonomy files
  bacteria_info <- read.delim(report_S_file, header = TRUE)
  colnames(bacteria_info) <- c("BC","name","taxonomy_id","taxonomy_lvl","reads","all_reads")
  bacteria_info_G <- read.delim(report_G_file, header = TRUE)
  colnames(bacteria_info_G) <- c("BC","name","taxonomy_id","taxonomy_lvl","reads","all_reads")
  tax_df <- read.delim(taxonomy_file, header = TRUE, stringsAsFactors = FALSE)
  
  # Read expression matrix
  sample_dir <- file.path(data_path, sample_name)
  expr_mat <- tryCatch(Read10X(sample_dir), error = function(e) NULL)
  if (is.null(expr_mat) || ncol(expr_mat) == 0) next
  
  # Record initial stats
  initial_stats <- rbind(initial_stats,
                         data.frame(Sample = sample_name,
                                    Genes = nrow(expr_mat),
                                    Cells = ncol(expr_mat)))
  expr_count_raw[[sample_name]] <- data.frame(cell_id = seq_len(ncol(expr_mat)),
                                              count = Matrix::colSums(expr_mat != 0))
  
  # Filter cells by global purity and t-test
  filtered_tax_df <- tax_df[tax_df$fraction_total_reads >= GLOBAL_PURITY_CUTOFF, ]
  filtered_tax_df <- filtered_tax_df[
    apply(filtered_tax_df[, c("fraction_total_reads","fraction_total_reads2","fraction_total_reads3")],
          1, function(x) pass_ttest(x[1], x[2], x[3])), ]
  keep_cells <- intersect(colnames(expr_mat), filtered_tax_df$barcode)
  expr_mat <- expr_mat[, keep_cells, drop = FALSE]
  
  # Record filtered stats
  after_filter_stats <- rbind(after_filter_stats,
                              data.frame(Sample = sample_name,
                                         Genes = nrow(expr_mat),
                                         Cells = ncol(expr_mat)))
  if (ncol(expr_mat) == 0) next
  
  # Create Seurat object
  sample_obj <- CreateSeuratObject(expr_mat, min.cells = 10, min.features = 10)
  BC_info   <- bacteria_info[match(colnames(sample_obj), bacteria_info$BC), ]
  BC_info_G <- bacteria_info_G[match(colnames(sample_obj), bacteria_info_G$BC), ]
  sample_obj$species_info <- BC_info$name
  sample_obj$genus_info <- BC_info_G$name
  sample_obj$sample <- sample_name
  parts <- strsplit(sample_name, "_")[[1]]
  timepoint <- parts[1]
  astronaut <- paste("Astronaut", parts[2])
  stage <- gsub("[0-9+-]", "", timepoint)
  sample_obj$timepoint <- factor(timepoint)
  sample_obj$astronaut <- factor(astronaut)
  sample_obj$stage <- factor(stage)
  sample_obj$orig.ident <- sample_name
  sample_obj$gene_number <- Matrix::colSums(GetAssayData(sample_obj, slot = "counts") != 0)
  
  # Store per-cell gene numbers
  gene_number_list[[sample_name]] <- data.frame(cell_id = colnames(sample_obj),
                                                gene_number = sample_obj$gene_number,
                                                sample = sample_name)
  
  # Record Seurat stats
  after_seurat_stats <- rbind(after_seurat_stats,
                              data.frame(Sample = sample_name,
                                         Genes = nrow(sample_obj),
                                         Cells = ncol(sample_obj)))
  expr_count_filtered[[sample_name]] <- data.frame(cell_id = seq_len(ncol(sample_obj)),
                                                   count = sample_obj$gene_number)
  seurat_list[[sample_name]] <- sample_obj
  cat("Finished:", sample_name,
      "| global purity cutoff = 0.5",
      "| cells =", ncol(sample_obj), "\n")
}

# Save summary statistics to CSV files
write.csv(initial_stats, "initial_stats.csv", row.names = FALSE)
write.csv(after_filter_stats, "after_filter_stats.csv", row.names = FALSE)
write.csv(after_seurat_stats, "after_seurat_stats.csv", row.names = FALSE)

# Add sample prefixes to barcodes to avoid conflicts
for (sample_name in names(seurat_list)) {
  obj <- seurat_list[[sample_name]]
  bc <- colnames(obj)
  need_prefix <- !grepl(paste0("^", sample_name, "_"), bc)
  if (any(need_prefix)) {
    bc[need_prefix] <- paste0(sample_name, "_", bc[need_prefix])
    colnames(obj) <- bc
    rownames(obj@meta.data) <- bc
  }
  seurat_list[[sample_name]] <- obj
}

# Harmonize gene sets across all samples
all_genes <- unique(unlist(lapply(seurat_list, rownames)))
for (nm in names(seurat_list)) {
  obj <- seurat_list[[nm]]
  missing <- setdiff(all_genes, rownames(obj))
  if (length(missing) > 0) {
    zero <- matrix(
      0, nrow = length(missing), ncol = ncol(obj),
      dimnames = list(missing, colnames(obj))
    )
    counts_mat <- GetAssayData(obj, slot = "counts")
    obj <- CreateSeuratObject(counts = rbind(counts_mat, zero), meta.data = obj@meta.data)
  }
  seurat_list[[nm]] <- obj
}

# Merge all Seurat objects into a single combined object
combined_seurat_object <- Reduce(merge, seurat_list)
cat("All samples merged. Object dimensions:", dim(combined_seurat_object), "\n")

# Remove cells with zero total RNA counts
combined_seurat_object <- subset(combined_seurat_object, subset = nCount_RNA > 0 & is.finite(nCount_RNA))

# Normalize and scale the combined Seurat object using SCTransform
combined_seurat_object <- SCTransform(
  combined_seurat_object,
  vars.to.regress = NULL,
  verbose = FALSE,
  return.only.var.genes = FALSE,
  variable.features.n = 2000
)
cat("SCTransform completed. Final dimensions:", dim(combined_seurat_object), "\n")

# Save the final Seurat object to disk
saveRDS(combined_seurat_object, file = "combined_seurat_object.rds")

# Reload the Seurat object if needed
combined_seurat_object <- readRDS(file = "combined_seurat_object.rds")

set.seed(1024)

# Filter cells to retain only those belonging to the target species
target_species <- "Phocaeicola vulgatus"

if (!"species_info" %in% colnames(combined_seurat_object@meta.data)) {
  stop("species_info column not found in meta.data")
}

species_cells <- rownames(combined_seurat_object@meta.data)[
  combined_seurat_object@meta.data$species_info == target_species
]

if (length(species_cells) == 0) {
  stop(paste("No cells found for species:", target_species))
}

combined_seurat_object <- subset(combined_seurat_object, cells = species_cells)

cat("Cells retained after species filtering:",
    ncol(combined_seurat_object), "\n")

# Perform PCA on the filtered cells using highly variable genes
combined_seurat_object <- RunPCA(
  combined_seurat_object, 
  features = VariableFeatures(combined_seurat_object)
)

# Define a function to automatically detect the elbow point in PCA
find_pca_elbow <- function(
    seurat_obj, 
    reduction = "pca", 
    ndims = 50,
    axis_title_size = 25, 
    axis_text_size = 25, 
    axis_font = "sans",
    title_size = 30, 
    title_font = "sans"
) {
  
  # Extract the standard deviations of principal components
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop("Specified reduction not found")
  }
  pca_sd <- seurat_obj[[reduction]]@stdev[1:ndims]
  
  # Normalize PC indices and standard deviations for elbow calculation
  x <- 1:ndims
  y <- pca_sd
  x_norm <- (x - min(x)) / (max(x) - min(x))
  y_norm <- (y - min(y)) / (max(y) - min(y))
  
  # Compute perpendicular distance from diagonal to identify the elbow
  distance <- abs(y_norm - (1 - x_norm)) / sqrt(2)
  elbow_point <- which.max(distance)
  
  # Visualize PCA elbow plot with detected point highlighted
  plot_df <- data.frame(PC = x, SD = y)
  p <- ggplot(plot_df, aes(x = PC, y = SD)) +
    geom_point(size = 3, color = "#2E86AB") + 
    geom_line(color = "#2E86AB") +
    geom_vline(xintercept = elbow_point, linetype = "dashed", color = "red") +
    annotate("text", x = elbow_point + 1, y = max(y), 
             label = paste0("Elbow = PC", elbow_point), 
             color = "red", size = 8, hjust = 0) +
    labs(
      title = "PCA Elbow Plot",
      x = "Principal Component",
      y = "Standard Deviation"
    ) +
    theme_bw(base_size = 20) +
    theme(
      plot.title = element_text(
        family = title_font, size = title_size, face = "bold", hjust = 0.5
      ),
      axis.title = element_text(
        family = axis_font, size = axis_title_size, face = "bold"
      ),
      axis.text = element_text(
        family = axis_font, size = axis_text_size
      ),
      panel.grid.major = element_line(color = "grey80", linetype = "dotted"),
      panel.grid.minor = element_blank()
    )
  
  print(p)
  
  cat("Suggested number of PCs:", elbow_point, "\n")
  
  return(elbow_point)
}

# Detect the optimal number of PCs for downstream analysis
elbow_pc <- find_pca_elbow(combined_seurat_object)

# Perform UMAP dimensionality reduction and cluster cells
combined_seurat_object <- combined_seurat_object %>%
  RunUMAP(reduction = "pca", dims = 1:elbow_pc) %>%
  FindNeighbors(reduction = "pca", dims = 1:elbow_pc) %>%
  FindClusters(resolution = 0.2)

cat("UMAP and clustering completed. Object dimensions:",
    dim(combined_seurat_object), "\n")

# Export the clustering results and metadata for Phocaeicola vulgatus
Idents(combined_seurat_object) <- "seurat_clusters"

md_species <- combined_seurat_object@meta.data

if (!"seurat_clusters" %in% colnames(md_species)) {
  stop("seurat_clusters not found in metadata. Please run FindClusters first.")
}

# Create a tidy metadata table with key information for each cell
output_df_species <- data.frame(
  barcode   = rownames(md_species),
  cluster   = as.character(md_species$seurat_clusters),
  sample    = md_species$orig.ident,
  astronaut = md_species$astronaut,
  timepoint = md_species$timepoint,
  stage     = md_species$stage,
  stringsAsFactors = FALSE
)

# Save metadata table to a TSV file
write.table(
  output_df_species,
  file = "Phocaeicola_vulgatus_metadata.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Metadata for Phocaeicola vulgatus exported successfully\n")

# Define colors for each cluster and extract UMAP coordinates
clusters <- levels(Idents(combined_seurat_object))
n_clusters <- length(clusters)

cluster_colors <- setNames(
  colorRampPalette(brewer.pal(12, "Paired"))(n_clusters),
  clusters
)

umap_coords <- Embeddings(combined_seurat_object, "umap") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell")
umap_coords$cluster <- Idents(combined_seurat_object)[umap_coords$cell]

# Calculate ranges and define arrow parameters for mini-axis
x_min <- min(umap_coords$umap_1); x_max <- max(umap_coords$umap_1)
y_min <- min(umap_coords$umap_2); y_max <- max(umap_coords$umap_2)
x_range <- x_max - x_min; y_range <- y_max - y_min
axis_len_x <- 0.2 * x_range; axis_len_y <- 0.2 * y_range
axis_x0 <- x_min - 0.10 * x_range; axis_y0 <- y_min - 0.10 * y_range
axis_arrow <- arrow(type = "closed", length = unit(6, "mm"))

# Plot UMAP with cells colored by cluster
p1 <- DimPlot(
  combined_seurat_object, reduction = "umap", label = FALSE,
  cols = cluster_colors, pt.size = 0.1, raster = FALSE
) +
  labs(title = "UMAP by Cluster") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"),
    axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
    legend.text = element_text(size = 20, family = "Arial"), legend.key.height = unit(2, "lines"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 90, 90, unit = "pt")
  ) +
  coord_cartesian(clip = "off") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 10))) +
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0,
           linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y,
           linewidth = 2, arrow = axis_arrow) +
  annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range,
           label = "UMAP 1", size = 8, family = "Arial", vjust = 1) +
  annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2,
           label = "UMAP 2", size = 8, family = "Arial", angle = 90, hjust = 0.5)

# Save UMAP plot colored by cluster
ggsave("Figure5a_cluster.pdf", plot = p1, device = cairo_pdf, width = 10, height = 10, units = "in")

# Define colors for each astronaut
astronaut_colors <- c("Astronaut 1" = "#FF9A8B", "Astronaut 2" = "#FFE0AC", "Astronaut 3" = "#BBDED6")

# Add jitter to UMAP coordinates to improve visualization by astronaut
umap_coords <- Embeddings(combined_seurat_object, "umap") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  mutate(umap_1 = umap_1 + runif(n(), -0.02, 0.02),
         umap_2 = umap_2 + runif(n(), -0.02, 0.02))

# Recalculate axis ranges for astronaut plot
x_min <- min(umap_coords$umap_1); x_max <- max(umap_coords$umap_1)
y_min <- min(umap_coords$umap_2); y_max <- max(umap_coords$umap_2)
x_range <- x_max - x_min; y_range <- y_max - y_min
axis_len_x <- 0.2 * x_range; axis_len_y <- 0.2 * y_range
axis_x0 <- x_min - 0.10 * x_range; axis_y0 <- y_min - 0.10 * y_range

# Plot UMAP with cells colored by astronaut
p2 <- DimPlot(
  combined_seurat_object, group.by = "astronaut", reduction = "umap",
  label = FALSE, cols = astronaut_colors, pt.size = 0.1, raster = FALSE, alpha = 1
) +
  scale_color_manual(values = astronaut_colors, labels = names(astronaut_colors)) +
  labs(title = "UMAP by Astronaut") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"),
        axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size = 20, family = "Arial"), legend.key.height = unit(2, "lines"),
        legend.key.width = unit(1.5, "lines"), panel.grid = element_blank(),
        plot.margin = margin(10, 10, 90, 90, unit = "pt"),
        legend.position = "right", legend.direction = "vertical") +
  coord_cartesian(clip = "off") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 10))) +
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0,
           linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y,
           linewidth = 2, arrow = axis_arrow) +
  annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range,
           label = "UMAP 1", size = 8, family = "Arial") +
  annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2,
           label = "UMAP 2", size = 8, family = "Arial", angle = 90)

# Save UMAP plot colored by astronaut
ggsave("Figure5a_astronaut.pdf", plot = p2, device = cairo_pdf, width = 11.5, height = 10, units = "in")

# Define colors for each timepoint and set factor levels
time_colors <- c(
  "L-60"  = "#F4F1DE",
  "L-30"  = "#E6D8AD",
  "FD30"  = "#5F6F75",
  "FD90"  = "#81B295",
  "FD150" = "#A3BFA8",
  "R+1"   = "#9A031E",
  "R+7"   = "#BA3A26",
  "R+14"  = "#E07A5F"
)
combined_seurat_object$timepoint <- factor(
  combined_seurat_object$timepoint,
  levels = names(time_colors)
)

# Plot UMAP colored by timepoint
p3 <- DimPlot(
  combined_seurat_object, group.by = "timepoint", reduction = "umap",
  label = FALSE, cols = time_colors, pt.size = 0.1, raster = FALSE, alpha = 1
) +
  scale_color_manual(values = time_colors, labels = names(time_colors)) +
  labs(title = "UMAP by Timepoint") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"),
    axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
    legend.text = element_text(size = 20, family = "Arial"),
    legend.key.height = unit(2, "lines"), legend.key.width = unit(1.5, "lines"),
    panel.grid = element_blank(), panel.border = element_blank(),
    legend.position = "right", legend.direction = "vertical",
    plot.margin = margin(10, 10, 90, 90, unit = "pt")
  ) +
  coord_cartesian(clip = "off") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 10))) +
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0,
           linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y,
           linewidth = 2, arrow = axis_arrow) +
  annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range,
           label = "UMAP 1", size = 8, family = "Arial") +
  annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2,
           label = "UMAP 2", size = 8, family = "Arial", angle = 90)

# Save UMAP plot colored by timepoint
ggsave("Figure5a_timepoint.pdf", plot = p3, device = cairo_pdf, width = 10.5, height = 10, units = "in")

# Define colors for each stage
stage_colors <- c("L" = "#FF6F61", "FD" = "#6B5B95", "R" = "#88B04B")

# Plot UMAP colored by stage
p4 <- DimPlot(
  combined_seurat_object, group.by = "stage", reduction = "umap",
  label = FALSE, cols = stage_colors, pt.size = 0.1, raster = FALSE, alpha = 1
) +
  scale_color_manual(values = stage_colors, labels = c("Stage L", "Stage FD", "Stage R")) +
  labs(title = "UMAP by Stage") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"),
    axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
    legend.text = element_text(size = 20, family = "Arial"),
    legend.key.height = unit(2, "lines"), legend.key.width = unit(1.5, "lines"),
    panel.grid = element_blank(), panel.border = element_blank(),
    legend.position = "right", legend.direction = "vertical",
    plot.margin = margin(10, 10, 90, 90, unit = "pt")
  ) +
  coord_cartesian(clip = "off") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 10))) +
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0,
           linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y,
           linewidth = 2, arrow = axis_arrow) +
  annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range,
           label = "UMAP 1", size = 8, family = "Arial") +
  annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2,
           label = "UMAP 2", size = 8, family = "Arial", angle = 90)

# Save UMAP plot colored by stage
ggsave("Figure5a_stage.pdf", plot = p4, device = cairo_pdf, width = 11, height = 10, units = "in")

# Convert the Seurat object into a Monocle3 cell_data_set for trajectory analysis
cds <- as.cell_data_set(combined_seurat_object)
cds@clusters$UMAP$clusters <- Idents(combined_seurat_object)
cds@colData@listData$seurat_clusters <- Idents(combined_seurat_object)

# Recalculate size factors for accurate normalization in Monocle3
cds <- estimate_size_factors(cds)

# Define a function to automatically determine the optimal number of principal components
find_cds_num_dim <- function(cds, max_dim = 50, plot = TRUE) {
  expr_mat <- counts(cds)
  expr_mat_log <- log1p(expr_mat)
  gene_var <- Matrix::rowMeans(expr_mat_log^2) - (Matrix::rowMeans(expr_mat_log))^2
  expr_mat_log_filt <- expr_mat_log[gene_var > 0, ]
  cat("Retained genes with non-zero variance:", nrow(expr_mat_log_filt), "\n")
  
  n_pcs <- min(max_dim, ncol(expr_mat_log_filt))
  pca_res <- prcomp_irlba(t(expr_mat_log_filt), n = n_pcs)
  
  x <- 1:length(pca_res$sdev)
  y <- pca_res$sdev
  x_norm <- (x - min(x)) / (max(x) - min(x))
  y_norm <- (y - min(y)) / (max(y) - min(y))
  distance <- abs(y_norm - (1 - x_norm)) / sqrt(2)
  elbow_point <- which.max(distance)
  
  if (plot) {
    plot_df <- data.frame(PC = x, SD = y)
    p <- ggplot(plot_df, aes(x = PC, y = SD)) +
      geom_point(size = 3, color = "#2E86AB") +
      geom_line(color = "#2E86AB") +
      geom_vline(xintercept = elbow_point, linetype = "dashed", color = "red") +
      annotate("text", x = elbow_point + 1, y = max(y),
               label = paste0("Elbow = PC", elbow_point), color = "red", size = 5) +
      labs(title = "Monocle3 num_dim Elbow Plot", x = "Principal Component", y = "Standard Deviation") +
      theme_bw()
    print(p)
  }
  cat("Suggested num_dim:", elbow_point, "\n")
  return(list(num_dim = elbow_point, elbow_plot = p))
}

# Identify the optimal number of principal components for downstream analysis
num_dim_res <- find_cds_num_dim(cds, max_dim = 50, plot = TRUE)
num_dim_opt <- num_dim_res$num_dim

# Preprocess the CDS and perform dimensionality reduction for trajectory inference
cds <- preprocess_cds(cds, num_dim = num_dim_opt)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds, resolution = 1e-3)
cds <- learn_graph(cds, use_partition = FALSE)
cds <- order_cells(cds)

# Extract UMAP coordinates and calculate ranges for axis annotations
umap_mat <- reducedDims(cds)$UMAP
umap_df <- as.data.frame(umap_mat)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
x_min <- min(umap_df$UMAP_1)
x_max <- max(umap_df$UMAP_1)
y_min <- min(umap_df$UMAP_2)
y_max <- max(umap_df$UMAP_2)
x_range <- x_max - x_min
y_range <- y_max - y_min

# Define parameters for miniature axis arrows in the plots
axis_len_x <- 0.2 * x_range
axis_len_y <- 0.2 * y_range
axis_x0 <- x_min - 0.05 * x_range
axis_y0 <- y_min - 0.05 * y_range
axis_arrow <- arrow(type = "closed", length = unit(6, "mm"))

# Plot trajectory colored by cluster identity
p5 <- plot_cells(cds, color_cells_by = "seurat_clusters",
                 label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE) +
  geom_point(aes(color = seurat_clusters), size = 0.2) +
  scale_color_manual(values = cluster_colors) +
  labs(title = "Trajectory Plot by Cluster") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"),
        axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size = 20, family = "Arial"), legend.title = element_blank(),
        legend.key.height = unit(2, "lines"), panel.grid = element_blank(),
        plot.margin = margin(10, 10, 30, 30, unit = "pt")) +
  coord_cartesian(clip = "off") +
  guides(color = guide_legend(override.aes = list(size = 10))) +
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0,
           linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y,
           linewidth = 2, arrow = axis_arrow) +
  annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range,
           label = "UMAP 1", size = 8, family = "Arial", vjust = 1) +
  annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2,
           label = "UMAP 2", size = 8, family = "Arial", angle = 90, hjust = 0.5)
ggsave("Figure5c_cluster.pdf", plot = p5, device = cairo_pdf, width = 10, height = 10, units = "in")

# Plot trajectory colored by pseudotime progression
p6 <- plot_cells(cds, color_cells_by = "pseudotime",
                 label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE) +
  geom_point(aes(color = pseudotime(cds)), size = 0.2) +
  labs(title = "Trajectory Plot by Pseudotime") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"),
        axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size = 20, family = "Arial"), legend.title = element_blank(),
        legend.key.height = unit(2, "lines"), panel.grid = element_blank(),
        plot.margin = margin(10, 10, 30, 30, unit = "pt")) +
  coord_cartesian(clip = "off") +
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0,
           linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y,
           linewidth = 2, arrow = axis_arrow) +
  annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range,
           label = "UMAP 1", size = 8, family = "Arial", vjust = 1) +
  annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2,
           label = "UMAP 2", size = 8, family = "Arial", angle = 90, hjust = 0.5)
ggsave("Figure5c_pseudotime.pdf", plot = p6, device = cairo_pdf, width = 10.2, height = 10, units = "in")

# Plot trajectory colored by experimental timepoints
cds@colData$timepoint <- factor(cds@colData$timepoint,
                                levels = c("L-60", "L-30", "FD30", "FD90", "FD150", "R+1", "R+7", "R+14"))
p7 <- plot_cells(cds, color_cells_by = "timepoint",
                 label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE,
                 cell_size = 1, cell_stroke = 0) +
  geom_point(aes(color = timepoint), size = 0.2) +
  scale_color_manual(values = time_colors) +
  labs(title = "Trajectory Plot by Timepoint") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"),
        axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size = 20, family = "Arial"), legend.title = element_blank(),
        legend.key.height = unit(2, "lines"), panel.grid = element_blank(),
        plot.margin = margin(10, 10, 30, 30, unit = "pt")) +
  coord_cartesian(clip = "off") +
  guides(color = guide_legend(override.aes = list(size = 10))) +
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0,
           linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y,
           linewidth = 2, arrow = axis_arrow) +
  annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range,
           label = "UMAP 1", size = 8, family = "Arial", vjust = 1) +
  annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2,
           label = "UMAP 2", size = 8, family = "Arial", angle = 90, hjust = 0.5)
ggsave("Figure5c_timepoint.pdf", plot = p7, device = cairo_pdf, width = 10.5, height = 10, units = "in")

# Plot trajectory colored by experimental stages
cds@colData$stage <- factor(cds@colData$stage, levels = c("L", "FD", "R"))
p8 <- plot_cells(cds, color_cells_by = "stage",
                 label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE,
                 cell_size = 1, cell_stroke = 0) +
  geom_point(aes(color = stage), size = 0.2) +
  scale_color_manual(values = stage_colors,
                     breaks = c("L", "FD", "R"),
                     labels = c("Stage L", "Stage FD", "Stage R")) +
  labs(title = "Trajectory Plot by Stage") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"),
        axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
        legend.text = element_text(size = 20, family = "Arial"), legend.title = element_blank(),
        legend.key.height = unit(2, "lines"), panel.grid = element_blank(),
        plot.margin = margin(10, 10, 30, 30, unit = "pt")) +
  coord_cartesian(clip = "off") +
  guides(color = guide_legend(override.aes = list(size = 10))) +
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0,
           linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y,
           linewidth = 2, arrow = axis_arrow) +
  annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range,
           label = "UMAP 1", size = 8, family = "Arial", vjust = 1) +
  annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2,
           label = "UMAP 2", size = 8, family = "Arial", angle = 90, hjust = 0.5)
ggsave("Figure5c_stage.pdf", plot = p8, device = cairo_pdf, width = 11, height = 10, units = "in")

# Save pseudotime values for all cells
pseudotime_df <- data.frame(cell = rownames(colData(cds)), pseudotime = pseudotime(cds))
write.csv(pseudotime_df, "species_pseudotime.csv", row.names = FALSE)

# Identify genes whose expression changes significantly along pseudotime
deg_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 10)
deg_genes_sig <- deg_genes %>% filter(q_value < 0.05)
write.csv(deg_genes_sig, "species_pseudotime_deg.csv", row.names = TRUE)

# Export expression matrix of top DEGs with corresponding pseudotime
expr_mat <- counts(cds)
top_deg <- rownames(deg_genes_sig)
expr_mat_sub <- expr_mat[top_deg, ]
expr_df_t <- as.data.frame(t(expr_mat_sub))
expr_df_t$pseudotime <- pseudotime(cds)[colnames(expr_mat_sub)]
write.csv(expr_df_t, "species_pseudotime_deg_expression.csv", quote = FALSE)