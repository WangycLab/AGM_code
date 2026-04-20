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
library(forcats)

# Set the working directory to the current folder
setwd(getwd())

# Define data path based on current directory
data_path <- gsub("Space_Scripts[/\\\\]Figure2", "Store_Matrix", getwd())

# List top-level directories for each sample
name_list <- list.dirs(data_path, full.names = FALSE, recursive = FALSE)
seurat_list <- list()

# Initialize data frames to track sample statistics at different steps
initial_stats <- data.frame(Sample = character(), Genes = integer(), Cells = integer(), stringsAsFactors = FALSE)
after_filter_stats <- data.frame(Sample = character(), Genes = integer(), Cells = integer(), stringsAsFactors = FALSE)
after_seurat_stats <- data.frame(Sample = character(), Genes = integer(), Cells = integer(), stringsAsFactors = FALSE)

# Prepare lists to store cell-level non-zero gene counts
expr_count_raw <- list() 
expr_count_filtered <- list() 
gene_number_list <- list()  

# Define a t-test based function to filter cells
pass_ttest <- function(fraction_total_reads, fraction_total_reads2, fraction_total_reads3) {
  values <- c(fraction_total_reads2, fraction_total_reads3)
  values <- as.numeric(values)
  fraction_total_reads <- as.numeric(fraction_total_reads)
  
  values <- values[!is.na(values)]
  if(length(values) < 2 || var(values) == 0 || is.na(fraction_total_reads)) {
    return(FALSE)
  }
  
  test <- t.test(values, mu = fraction_total_reads)
  return(test$p.value < 0.05 & fraction_total_reads > max(values))
}

# Process each sample directory
for (sample_name in name_list) {
  cat("Processing sample:", sample_name, "\n")
  
  # Read species- and genus-level taxonomy reports
  report_S_file <- file.path(data_path, paste0(sample_name, "_sc_taxonomy_greedy_S.report"))
  report_G_file <- file.path(data_path, paste0(sample_name, "_sc_taxonomy_greedy_G.report"))
  
  if(!file.exists(report_S_file) || !file.exists(report_G_file)) {
    warning("Missing report files, skipping sample:", sample_name)
    next
  }
  
  bacteria_info <- read.delim(report_S_file, header = TRUE, sep = "\t")
  colnames(bacteria_info) <- c("BC", "name", "taxonomy_id", "taxonomy_lvl", "reads", "all_reads")
  
  bacteria_info_G <- read.delim(report_G_file, header = TRUE, sep = "\t")
  colnames(bacteria_info_G) <- c("BC", "name", "taxonomy_id", "taxonomy_lvl", "reads", "all_reads")
  
  # Read combined taxonomy report
  taxonomy_file <- file.path(data_path, paste0(sample_name, "_sc_taxonomy.report"))
  if(!file.exists(taxonomy_file)) {
    warning("Missing taxonomy.report file, skipping sample:", sample_name)
    next
  }
  tax_df <- read.delim(taxonomy_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Load the expression matrix
  sample_dir <- file.path(data_path, sample_name)
  if(!dir.exists(sample_dir)) {
    warning("10X data directory not found, skipping sample:", sample_name)
    next
  }
  
  expr_mat <- tryCatch(Read10X(data.dir = sample_dir), error = function(e) {
    warning("Failed to read 10X data, skipping sample:", sample_name)
    return(NULL)
  })
  
  if(is.null(expr_mat) || ncol(expr_mat) == 0) {
    warning("Expression matrix is empty, skipping sample:", sample_name)
    next
  }
  
  # Record initial sample statistics
  initial_stats <- rbind(initial_stats, data.frame(
    Sample = sample_name,
    Genes = nrow(expr_mat),
    Cells = ncol(expr_mat)
  ))
  
  # Count non-zero genes per cell before filtering
  nonzero_counts_raw <- Matrix::colSums(expr_mat != 0)
  expr_count_raw[[sample_name]] <- data.frame(
    cell_id = seq_along(nonzero_counts_raw),
    count   = nonzero_counts_raw,
    stringsAsFactors = FALSE
  )
  colnames(expr_count_raw[[sample_name]])[2] <- sample_name
  
  # Filter cells based on read fractions and t-test
  if(!"fraction_total_reads" %in% colnames(tax_df) || !"barcode" %in% colnames(tax_df)) {
    warning("Missing fraction_total_reads or barcode column, skipping sample:", sample_name)
    next
  }
  
  filtered_tax_df <- tax_df[tax_df$fraction_total_reads > 0.5, ]
  filtered_tax_df <- filtered_tax_df[
    apply(filtered_tax_df[, c("fraction_total_reads", "fraction_total_reads2", "fraction_total_reads3")], 1, function(x) {
      pass_ttest(x[1], x[2], x[3])
    }), ]
  
  keep_cells <- intersect(colnames(expr_mat), filtered_tax_df$barcode)
  expr_mat_filtered <- expr_mat[, keep_cells, drop = FALSE]
  
  # Select top cells (all kept here)
  nonzero_counts <- Matrix::colSums(expr_mat_filtered != 0)
  top_cells <- colnames(expr_mat_filtered)
  expr_mat_top <- expr_mat_filtered[, top_cells, drop = FALSE]
  
  # Record statistics after filtering
  after_filter_stats <- rbind(after_filter_stats, data.frame(
    Sample = sample_name,
    Genes = nrow(expr_mat_top),
    Cells = ncol(expr_mat_top)
  ))
  
  if(ncol(expr_mat_top) == 0) {
    warning("No cells left after filtering, skipping sample:", sample_name)
    next
  }
  
  # Create Seurat object and annotate metadata
  sample_obj <- CreateSeuratObject(counts = expr_mat_top, min.cells = 10, min.features = 10)
  BC_info <- bacteria_info[match(colnames(sample_obj), bacteria_info$BC), ]
  BC_info_G <- bacteria_info_G[match(colnames(sample_obj), bacteria_info_G$BC), ]
  
  sample_obj@meta.data$species_info <- BC_info$name
  sample_obj@meta.data$genus_info <- BC_info_G$name
  sample_obj@meta.data$orig.ident <- sample_name
  
  # Add gene number per cell
  gene_number <- Matrix::colSums(GetAssayData(sample_obj, slot = "counts") != 0)
  sample_obj@meta.data$gene_number <- gene_number[colnames(sample_obj)]
  
  # Store cell gene counts
  gene_number_list[[sample_name]] <- data.frame(
    cell_id = colnames(sample_obj),
    gene_number = sample_obj@meta.data$gene_number,
    sample = sample_name,
    stringsAsFactors = FALSE
  )
  
  # Record statistics after Seurat object creation
  after_seurat_stats <- rbind(after_seurat_stats, data.frame(
    Sample = sample_name,
    Genes = nrow(GetAssayData(sample_obj, slot = "counts")),
    Cells = ncol(sample_obj)
  ))
  
  # Store filtered non-zero gene counts
  nonzero_counts_filtered <- gene_number
  expr_count_filtered[[sample_name]] <- data.frame(
    cell_id = seq_along(nonzero_counts_filtered),
    count   = nonzero_counts_filtered,
    stringsAsFactors = FALSE
  )
  colnames(expr_count_filtered[[sample_name]])[2] <- sample_name
  
  # Add Seurat object to the list
  seurat_list[[sample_name]] <- sample_obj
  cat("Sample", sample_name, "processed, number of cells:", ncol(sample_obj), "\n")
}

# Save CSVs for sample-level statistics at each processing step
write.csv(initial_stats, "initial_stats.csv", row.names = FALSE)
write.csv(after_filter_stats, "after_filter_stats.csv", row.names = FALSE)
write.csv(after_seurat_stats, "after_seurat_stats.csv", row.names = FALSE)

# Make cell barcodes unique by adding sample prefix
for (sample_name in names(seurat_list)) {
  sample_obj <- seurat_list[[sample_name]]
  old_barcodes <- colnames(sample_obj)
  new_barcodes <- paste0(sample_name, "_", old_barcodes)
  colnames(sample_obj) <- new_barcodes
  rownames(sample_obj@meta.data) <- new_barcodes
  seurat_list[[sample_name]] <- sample_obj
}

# Merge all samples using the union of genes
all_genes <- unique(unlist(lapply(seurat_list, function(x) rownames(x))))

# Fill missing genes with zeros for each sample
for (sample_name in names(seurat_list)) {
  sample_obj <- seurat_list[[sample_name]]
  missing_genes <- setdiff(all_genes, rownames(sample_obj))
  
  if (length(missing_genes) > 0) {
    zero_mat <- matrix(0, nrow = length(missing_genes), ncol = ncol(sample_obj))
    rownames(zero_mat) <- missing_genes
    colnames(zero_mat) <- colnames(sample_obj)
    expr_mat <- rbind(GetAssayData(sample_obj, slot = "counts"), zero_mat)
    sample_obj <- CreateSeuratObject(counts = expr_mat, meta.data = sample_obj@meta.data)
  }
  
  seurat_list[[sample_name]] <- sample_obj
}

# Merge all Seurat objects into one combined object
if (length(seurat_list) > 1) {
  combined_seurat_object <- Reduce(function(x, y) merge(x, y), seurat_list)
} else if (length(seurat_list) == 1) {
  combined_seurat_object <- seurat_list[[1]]
} else {
  stop("No valid samples loaded")
}

cat("All samples merged, combined object dimensions:", dim(combined_seurat_object), "\n")

# Remove cells with zero or non-finite RNA counts
combined_seurat_object <- subset(
  combined_seurat_object,
  subset = nCount_RNA > 0 & is.finite(nCount_RNA)
)

# Normalize the Seurat object using SCTransform and prepare for PCA
combined_seurat_object <- SCTransform(
  combined_seurat_object, 
  vars.to.regress = NULL, 
  verbose = FALSE,
  return.only.var.genes = FALSE, 
  variable.features.n = 2000
)

# Run PCA on variable features to reduce dimensionality
combined_seurat_object <- RunPCA(combined_seurat_object, 
                                 features = VariableFeatures(combined_seurat_object))

# Inspect the variance explained by each principal component using an elbow plot
ElbowPlot(combined_seurat_object, ndims = 50, reduction = "pca")

# Function to automatically detect the optimal PCA elbow point
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
  
  # Ensure the requested reduction exists and extract standard deviations
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop("Specified reduction does not exist!")
  }
  pca_sd <- seurat_obj[[reduction]]@stdev[1:ndims]
  
  # Normalize PC indices and standard deviations for distance calculation
  x <- 1:ndims
  y <- pca_sd
  x_norm <- (x - min(x)) / (max(x) - min(x))
  y_norm <- (y - min(y)) / (max(y) - min(y))
  
  # Compute distance from each PC point to the diagonal line
  distance <- abs(y_norm - (1 - x_norm)) / sqrt(2)
  
  # Identify elbow point as the maximum distance from the diagonal
  elbow_point <- which.max(distance)
  
  # Prepare a dataframe for plotting
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
      plot.title = element_text(family = title_font, size = title_size, face = "bold", hjust = 0.5),
      axis.title = element_text(family = axis_font, size = axis_title_size, face = "bold"),
      axis.text = element_text(family = axis_font, size = axis_text_size),
      panel.grid.major = element_line(color = "grey80", linetype = "dotted"),
      panel.grid.minor = element_blank()
    )
  
  # Display the elbow plot for visual inspection
  print(p)  
  
  cat("Recommended number of PCs (Elbow Point):", elbow_point, "\n")
  
  return(elbow_point)
}

# Determine the number of principal components to use for downstream analysis
elbow_pc <- find_pca_elbow(combined_seurat_object)

# Apply Harmony for batch correction across different sample origins
set.seed(1024)
combined_seurat_object <- RunHarmony(
  object = combined_seurat_object,
  reduction = "pca",
  group.by.vars = "orig.ident", 
  reduction.save = "harmony",
  assay.use = "SCT" 
)

# Use Harmony embeddings to perform UMAP, find neighbors, and cluster cells
combined_seurat_object <- combined_seurat_object %>%
  RunUMAP(reduction = "harmony", dims = 1:elbow_pc) %>%
  FindNeighbors(reduction = "harmony", dims = 1:elbow_pc) %>%
  FindClusters(resolution = 0.3)

# Save the integrated Seurat object for future use
saveRDS(
  combined_seurat_object,
  file = "combined_seurat_object.rds"
)

# Load the saved Seurat object
combined_seurat_object <- readRDS("combined_seurat_object.rds")

# Extract UMAP coordinates for plotting additional annotations
umap_df <- Embeddings(combined_seurat_object, "umap") %>%
  as.data.frame()

# Compute ranges for the mini axis to position it outside the plot area
x_min <- min(umap_df$umap_1)
x_max <- max(umap_df$umap_1)
y_min <- min(umap_df$umap_2)
y_max <- max(umap_df$umap_2)
x_range <- x_max - x_min
y_range <- y_max - y_min

# Set positions and length for custom mini axes outside the plot
axis_x0 <- x_min - 0.15 * x_range
axis_y0 <- y_min - 0.15 * y_range
axis_len_x <- 0.25 * x_range
axis_len_y <- 0.25 * y_range

# Define arrow style for the mini axes
axis_arrow <- arrow(type = "closed", length = unit(0.25, "cm"))

# Assign distinct colors for each cluster
clusters <- levels(Idents(combined_seurat_object))
n_clusters <- length(clusters)
cluster_colors <- setNames(colorRampPalette(brewer.pal(12, "Paired"))(n_clusters), clusters)

# Plot UMAP colored by clusters and add mini axes
p1 <- DimPlot(
  combined_seurat_object,
  reduction = "umap",
  label = FALSE,
  cols = cluster_colors,
  pt.size = 0.1,
  raster = FALSE
) +
  labs(title = "UMAP by Cluster") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.text = element_text(size = 20, family = "Arial"),
    legend.title = element_blank(),
    legend.key.height = unit(1.5, "lines"),
    legend.key.width = unit(1.5, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    legend.position = "right",
    legend.direction = "vertical",
    plot.margin = margin(10, 10, 120, 120, unit = "pt")
  ) +
  coord_cartesian(clip = "off") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 8))) +
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0, linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y, linewidth = 2, arrow = axis_arrow) +
  annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.04 * y_range, label = "UMAP 1", size = 8, family = "Arial", vjust = 1) +
  annotate("text", x = axis_x0 - 0.04 * x_range, y = axis_y0 + axis_len_y / 2, label = "UMAP 2", size = 8, family = "Arial", angle = 90, hjust = 0.5)

# Save the cluster-colored UMAP plot as PDF
ggsave(filename = "Figure2i_cluster.pdf", plot = p1, width = 9, height = 8, units = "in", device = cairo_pdf)

# Set the order of samples for consistent coloring in plots
sample_order <- c("Con", "2mo", "4mo", "6mo")
combined_seurat_object$orig.ident <- factor(combined_seurat_object$orig.ident, levels = sample_order)

# Define colors for each sample group
n_samples <- length(sample_order)
sample_colors <- setNames(colorRampPalette(brewer.pal(9, "Set3"))(n_samples), sample_order)

# Shuffle cells to randomize plotting order for better visual clarity
combined_seurat_object <- combined_seurat_object[, sample(colnames(combined_seurat_object))]

# Plot UMAP colored by sample with mini axes
p2 <- DimPlot(
  combined_seurat_object,
  reduction = "umap",
  group.by = "orig.ident",
  pt.size = 0.1,
  cols = sample_colors,
  label = FALSE,
  raster = FALSE,
  alpha = 0.5
) +
  labs(title = "UMAP by Sample") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.text = element_text(size = 20, family = "Arial"),
    legend.title = element_blank(),
    legend.key.height = unit(2, "lines"),
    legend.key.width = unit(1.5, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    legend.position = "right",
    legend.direction = "vertical",
    plot.margin = margin(10, 10, 120, 120, unit = "pt")
  ) +
  coord_cartesian(clip = "off") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 8))) +
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0, linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y, linewidth = 2, arrow = axis_arrow) +
  annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.04 * y_range, label = "UMAP 1", size = 8, family = "Arial", vjust = 1) +
  annotate("text", x = axis_x0 - 0.04 * x_range, y = axis_y0 + axis_len_y / 2, label = "UMAP 2", size = 8, family = "Arial", angle = 90, hjust = 0.5)

# Save the sample-colored UMAP plot as PDF
ggsave(filename = "Figure2i_sample.pdf", plot = p2, width = 9.3, height = 8, units = "in", device = cairo_pdf)

# Set cluster identities to export metadata
Idents(combined_seurat_object) <- "seurat_clusters"

# Extract cell metadata including cluster, sample, and gene information
md <- combined_seurat_object@meta.data

# Ensure clustering information exists
if (!"seurat_clusters" %in% colnames(md)) {
  stop("No seurat_clusters found, please run FindClusters() first.")
}

# Add a clear sample label column
md$sample_label <- md$orig.ident

# Prepare a dataframe to export cell-level metadata
output_df <- data.frame(
  barcode     = rownames(md),
  cluster     = as.character(md$seurat_clusters),
  sample      = md$sample_label, 
  species     = md$species_info,
  genus       = md$genus_info,
  gene_number = md$gene_number,
  stringsAsFactors = FALSE
)

# Export the metadata to a tab-separated file
write.table(output_df,
            file = "cell_metadata.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)