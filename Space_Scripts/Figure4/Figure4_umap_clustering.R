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

# Set the working directory to the current folder
current_dir <- trimws(getwd())
setwd(current_dir)
cat("Working directory set to:", current_dir, "\n")

# Determine the project root directory based on folder structure
parts <- strsplit(current_dir, .Platform$file.sep)[[1]]
if (length(parts) >= 2 && all(tail(parts, 2) == c("Space_Scripts", "Figure4"))) {
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

combined_seurat_object <- RunPCA(combined_seurat_object, features = VariableFeatures(combined_seurat_object))
ElbowPlot(combined_seurat_object, ndims = 50, reduction = "pca") 

# Function to automatically find PCA elbow point
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
  
  # Get PCA standard deviation
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop("Specified reduction does not exist")
  }
  pca_sd <- seurat_obj[[reduction]]@stdev[1:ndims]
  
  # Normalize coordinates
  x <- 1:ndims
  y <- pca_sd
  x_norm <- (x - min(x)) / (max(x) - min(x))
  y_norm <- (y - min(y)) / (max(y) - min(y))
  
  # Calculate distance to diagonal
  distance <- abs(y_norm - (1 - x_norm)) / sqrt(2)
  
  # Elbow point is the max distance
  elbow_point <- which.max(distance)
  
  # Plot elbow
  plot_df <- data.frame(PC = x, SD = y)
  p <- ggplot(plot_df, aes(x = PC, y = SD)) +
    geom_point(size = 3, color = "#2E86AB") + 
    geom_line(color = "#2E86AB") +
    geom_vline(xintercept = elbow_point, linetype = "dashed", color = "red") +
    annotate("text", x = elbow_point + 1, y = max(y), 
             label = paste0("Elbow = PC", elbow_point), 
             color = "red", size = 8, hjust = 0) +
    labs(title = "PCA Elbow Plot", x = "Principal Component", y = "Standard Deviation") +
    theme_bw(base_size = 20) +
    theme(
      plot.title = element_text(family = title_font, size = title_size, face = "bold", hjust = 0.5),
      axis.title = element_text(family = axis_font, size = axis_title_size, face = "bold"),
      axis.text = element_text(family = axis_font, size = axis_text_size),
      panel.grid.major = element_line(color = "grey80", linetype = "dotted"),
      panel.grid.minor = element_blank()
    )
  
  # Show plot in RStudio
  print(p)  
  
  # Print suggested PC
  cat("Suggested number of PCs:", elbow_point, "\n")
  
  return(elbow_point)
}

# Run PCA elbow detection
elbow_pc <- find_pca_elbow(combined_seurat_object)

# Run Harmony batch correction
set.seed(1024)
combined_seurat_object <- RunHarmony(
  object = combined_seurat_object,
  reduction = "pca",
  group.by.vars = "orig.ident",
  reduction.save = "harmony",
  assay.use = "SCT"
)

# Run downstream analysis with Harmony
combined_seurat_object <- combined_seurat_object %>%
  RunUMAP(reduction = "harmony", dims = 1:elbow_pc) %>%
  FindNeighbors(reduction = "harmony", dims = 1:elbow_pc) %>%
  FindClusters(resolution = 0.2)

# Set cluster identities
Idents(combined_seurat_object) <- "seurat_clusters"

# Export cell metadata including stage column
md <- combined_seurat_object@meta.data
if (!"seurat_clusters" %in% colnames(md)) stop("No seurat_clusters column found")

output_df <- data.frame(
  barcode   = rownames(md),
  cluster   = as.character(md$seurat_clusters),
  sample    = md$orig.ident,
  astronaut = md$astronaut,
  timepoint = md$timepoint,
  stage     = md$stage,          
  species   = md$species_info,
  genus     = md$genus_info,
  gene_number = md$gene_number,
  stringsAsFactors = FALSE
)

write.table(output_df, file = "cell_metadata.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Prepare for differential expression analysis
combined_seurat_object <- PrepSCTFindMarkers(combined_seurat_object)

# Find markers for all clusters
markers <- FindAllMarkers(
  combined_seurat_object,
  assay = "SCT",
  only.pos = TRUE,
  min.pct = 0.05,
  logfc.threshold = 0.25
)
markers <- markers %>% arrange(cluster, desc(avg_log2FC))
write.table(markers, file = "all_clusters_DEGs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Extract normalized expression matrix and cluster info
counts <- GetAssayData(combined_seurat_object, assay = "SCT", slot = "data")
clusters <- combined_seurat_object$seurat_clusters

# Prepare gene list
markers$gene <- gsub("_", "-", markers$gene)
genes <- unique(markers$gene)
genes <- genes[genes %in% rownames(counts)]

# Prepare cluster names
cluster_ids <- sort(unique(clusters))
cluster_names <- paste0("Cluster_", cluster_ids)

# Initialize result matrices
fraction_mat <- matrix(0, nrow = length(genes), ncol = length(cluster_ids),
                       dimnames = list(genes, cluster_names))
mean_mat <- matrix(0, nrow = length(genes), ncol = length(cluster_ids),
                   dimnames = list(genes, cluster_names))

# Compute fraction and mean expression per cluster
for (cid in cluster_ids) {
  cluster_cells <- names(clusters[clusters == cid])
  sub_counts <- counts[genes, cluster_cells, drop = FALSE]
  fraction <- rowSums(sub_counts > 0) / length(cluster_cells)
  mean_expr <- rowMeans(sub_counts)
  fraction_mat[, paste0("Cluster_", cid)] <- fraction
  mean_mat[, paste0("Cluster_", cid)] <- mean_expr
}

# Save fraction and mean matrices
write.csv(fraction_mat, "gene_fraction_by_cluster.csv", quote = FALSE)
write.csv(mean_mat, "gene_mean_expr_by_cluster.csv", quote = FALSE)

# Save Seurat object
saveRDS(combined_seurat_object, file = "combined_seurat_object.rds")
combined_seurat_object <- readRDS(file = "combined_seurat_object.rds")

# Prepare cluster color palette for UMAP visualization
clusters <- levels(Idents(combined_seurat_object))
n_clusters <- length(clusters)

cluster_colors <- setNames(
  colorRampPalette(brewer.pal(12, "Paired"))(n_clusters),
  clusters
)

# Extract UMAP coordinates and cluster identities
umap_coords <- Embeddings(combined_seurat_object, "umap") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell")

umap_coords$cluster <- Idents(combined_seurat_object)[umap_coords$cell]

# Calculate UMAP ranges for positioning mini-axis
x_min <- min(umap_coords$umap_1)
x_max <- max(umap_coords$umap_1)
y_min <- min(umap_coords$umap_2)
y_max <- max(umap_coords$umap_2)
x_range <- x_max - x_min
y_range <- y_max - y_min

axis_len_x <- 0.2 * x_range
axis_len_y <- 0.2 * y_range
axis_x0 <- x_min - 0.05 * x_range
axis_y0 <- y_min - 0.05 * y_range

axis_arrow <- arrow(type = "closed", length = unit(5, "mm"))

# Generate UMAP plot colored by cluster
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
    plot.title.position = "panel",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.text = element_text(size = 20, family = "Arial"),
    legend.key.height = unit(2, "lines"),
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 120, 120, unit = "pt")
  ) +
  guides(color = guide_legend(
    ncol = 2,
    override.aes = list(size = 8)
  )) +
  
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x,
           y = axis_y0, yend = axis_y0,
           linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0,
           y = axis_y0, yend = axis_y0 + axis_len_y,
           linewidth = 2, arrow = axis_arrow) +
  
  annotate("text", x = axis_x0 + axis_len_x / 2,
           y = axis_y0 - 0.05 * y_range,
           label = "UMAP 1", size = 8, family = "Arial", vjust = 1) +
  annotate("text", x = axis_x0 - 0.05 * x_range,
           y = axis_y0 + axis_len_y / 2,
           label = "UMAP 2", size = 8, family = "Arial", angle = 90, hjust = 0.5)

p1

# Export the cluster-based UMAP figure as a PDF
ggsave(
  filename = "Figure4c_cluster.pdf",
  plot = p1,
  device = cairo_pdf,
  width = 11,
  height = 10,
  units = "in",
  limitsize = FALSE
)

# Define astronaut-specific colors for UMAP visualization
astronaut_colors <- c(
  "Astronaut 1" = "#FF9A8B",
  "Astronaut 2" = "#FFE0AC",
  "Astronaut 3" = "#BBDED6"
)

# Extract UMAP coordinates for astronaut visualization
umap_coords <- Embeddings(combined_seurat_object, "umap") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell")

# Calculate UMAP ranges for mini-axis placement
x_min <- min(umap_coords$umap_1)
x_max <- max(umap_coords$umap_1)
y_min <- min(umap_coords$umap_2)
y_max <- max(umap_coords$umap_2)
x_range <- x_max - x_min
y_range <- y_max - y_min

axis_len_x <- 0.2 * x_range
axis_len_y <- 0.2 * y_range
axis_x0 <- x_min - 0.05 * x_range
axis_y0 <- y_min - 0.05 * y_range
axis_arrow <- arrow(type = "closed", length = unit(5, "mm"))

umap_coords <- umap_coords[sample(nrow(umap_coords)), ]

# Generate UMAP plot colored by astronaut identity
p2 <- DimPlot(
  combined_seurat_object,
  group.by = "astronaut",
  reduction = "umap",
  label = FALSE,
  cols = astronaut_colors,
  pt.size = 0.1,
  raster = FALSE,
  alpha = 1
) +
  labs(title = "UMAP by Astronaut") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"),
    plot.title.position = "panel",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.text = element_text(size = 20, family = "Arial"),
    legend.key.height = unit(2, "lines"),
    legend.key.width = unit(1.5, "lines"),
    legend.position = "right",
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 120, 120, unit = "pt")
  ) +
  guides(color = guide_legend(
    ncol = 1,
    override.aes = list(size = 8)
  )) +
  
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x,
           y = axis_y0, yend = axis_y0,
           linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0,
           y = axis_y0, yend = axis_y0 + axis_len_y,
           linewidth = 2, arrow = axis_arrow) +
  
  annotate("text", x = axis_x0 + axis_len_x / 2,
           y = axis_y0 - 0.05 * y_range,
           label = "UMAP 1", size = 8, family = "Arial", vjust = 1) +
  annotate("text", x = axis_x0 - 0.05 * x_range,
           y = axis_y0 + axis_len_y / 2,
           label = "UMAP 2", size = 8, family = "Arial", angle = 90, hjust = 0.5)

p2

# Export the astronaut-based UMAP figure as a PDF
ggsave(
  filename = "FigureS6_astronaut.pdf",
  plot = p2,
  device = cairo_pdf,
  width = 11.5,
  height = 10,
  units = "in",
  limitsize = FALSE
)

# Define timepoint colors for UMAP visualization
time_colors <- c(
  "L-60" = "#F4F1DE",
  "L-30" = "#E6D8AD",
  "FD30" = "#5F6F75",
  "FD90" = "#81B295",
  "FD150" = "#A3BFA8",
  "R+1" = "#9A031E",
  "R+7" = "#BA3A26",
  "R+14" = "#E07A5F"
)

combined_seurat_object$timepoint <- factor(
  combined_seurat_object$timepoint,
  levels = names(time_colors)
)

# Extract UMAP coordinates and randomize plotting order
umap_coords <- Embeddings(combined_seurat_object, "umap") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell")

umap_coords$timepoint <- combined_seurat_object$timepoint[umap_coords$cell]

umap_coords <- umap_coords[sample(nrow(umap_coords)), ]

# Calculate UMAP ranges for mini-axis positioning
x_min <- min(umap_coords$umap_1)
x_max <- max(umap_coords$umap_1)
y_min <- min(umap_coords$umap_2)
y_max <- max(umap_coords$umap_2)
x_range <- x_max - x_min
y_range <- y_max - y_min

axis_len_x <- 0.2 * x_range
axis_len_y <- 0.2 * y_range
axis_x0 <- x_min - 0.05 * x_range
axis_y0 <- y_min - 0.05 * y_range
axis_arrow <- arrow(type = "closed", length = unit(5, "mm"))

# Generate UMAP plot colored by timepoint
p3 <- ggplot(umap_coords, aes(x = umap_1, y = umap_2, color = timepoint)) +
  geom_point(size = 0.1, alpha = 1) +
  scale_color_manual(values = time_colors) +
  labs(title = "UMAP by Timepoint", color = NULL) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"),
    plot.title.position = "panel",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.text = element_text(size = 20, family = "Arial"),
    legend.key.height = unit(2, "lines"),
    legend.key.width = unit(1.5, "lines"),
    legend.position = "right",
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 120, 120, unit = "pt")
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 8))) +
  
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x,
           y = axis_y0, yend = axis_y0,
           linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0,
           y = axis_y0, yend = axis_y0 + axis_len_y,
           linewidth = 2, arrow = axis_arrow) +
  annotate("text", x = axis_x0 + axis_len_x / 2,
           y = axis_y0 - 0.05 * y_range,
           label = "UMAP 1", size = 8, family = "Arial", vjust = 1) +
  annotate("text", x = axis_x0 - 0.05 * x_range,
           y = axis_y0 + axis_len_y / 2,
           label = "UMAP 2", size = 8, family = "Arial", angle = 90, hjust = 0.5)

p3

# Export the timepoint-based UMAP figure as a PDF
ggsave(
  filename = "FigureS6_timepoint.pdf",
  plot = p3,
  device = cairo_pdf,
  width = 10.8,
  height = 10,
  units = "in",
  limitsize = FALSE
)

# Calculate species abundance and identify the top 20 species
species_abundance <- combined_seurat_object@meta.data %>%
  group_by(species_info) %>%
  summarise(cell_count = n()) %>%
  arrange(desc(cell_count))

top_species <- head(species_abundance$species_info, 20)

combined_seurat_object$top_species <- ifelse(
  combined_seurat_object$species_info %in% top_species,
  combined_seurat_object$species_info,
  "Others"
)

species_levels <- c(top_species, "Others")
combined_seurat_object$top_species <- factor(
  combined_seurat_object$top_species,
  levels = species_levels
)

# Generate color palette for top species and Others
top_colors <- colorRampPalette(brewer.pal(12, "Set3"))(20)
species_colors <- setNames(
  c(top_colors, "grey90"),
  species_levels
)

# Extract UMAP coordinates and shuffle cell order
umap_coords <- Embeddings(combined_seurat_object, "umap") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell")

umap_coords$species <- combined_seurat_object$top_species[umap_coords$cell]
umap_coords <- umap_coords[sample(nrow(umap_coords)), ]

# Compute coordinate ranges for mini-axis drawing
x_min <- min(umap_coords$umap_1)
x_max <- max(umap_coords$umap_1)
y_min <- min(umap_coords$umap_2)
y_max <- max(umap_coords$umap_2)
x_range <- x_max - x_min
y_range <- y_max - y_min

axis_len_x <- 0.2 * x_range
axis_len_y <- 0.2 * y_range
axis_x0 <- x_min - 0.05 * x_range
axis_y0 <- y_min - 0.05 * y_range
axis_arrow <- arrow(type = "closed", length = unit(5, "mm"))

# Draw UMAP visualization colored by top species
p4 <- ggplot(umap_coords, aes(x = umap_1, y = umap_2, color = species)) +
  geom_point(size = 0.1, alpha = 1) +
  scale_color_manual(values = species_colors) +
  labs(title = "UMAP by Top 20 Species", color = NULL) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"),
    plot.title.position = "panel",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.text = element_text(size = 15, face = "italic"),
    legend.key.height = unit(1.3, "lines"),
    legend.key.width = unit(1.3, "lines"),
    legend.position = "right",
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 120, 120, unit = "pt")
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 6))) +
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x,
           y = axis_y0, yend = axis_y0,
           linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0,
           y = axis_y0, yend = axis_y0 + axis_len_y,
           linewidth = 2, arrow = axis_arrow) +
  annotate("text", x = axis_x0 + axis_len_x / 2,
           y = axis_y0 - 0.05 * y_range,
           label = "UMAP 1", size = 8, family = "Arial", vjust = 1) +
  annotate("text", x = axis_x0 - 0.05 * x_range,
           y = axis_y0 + axis_len_y / 2,
           label = "UMAP 2", size = 8, family = "Arial", angle = 90, hjust = 0.5)

p4

# Save UMAP figure as editable vector PDF
ggsave(
  filename = "FigureS6_species.pdf",
  plot = p4,
  device = cairo_pdf,
  width = 13,
  height = 10,
  units = "in",
  limitsize = FALSE
)

# Recalculate abundance to obtain top 20 species list
species_abundance <- combined_seurat_object@meta.data %>%
  group_by(species_info) %>%
  summarise(cell_count = n()) %>%
  arrange(desc(cell_count))

top20_species <- head(species_abundance$species_info, 20)

# Define color mapping for species highlighting
top_colors <- colorRampPalette(brewer.pal(12, "Set3"))(20)
species_grouped_levels <- c(top20_species, "Others")
species_grouped_colors <- setNames(c(top_colors, "grey90"), species_grouped_levels)

# Create directory to store species-specific UMAP plots
output_dir <- "UMAP_top20_species"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Generate individual UMAP plots highlighting each species
for (sp in top20_species) {
  
  umap_df <- Embeddings(combined_seurat_object, "umap") %>% as.data.frame()
  colnames(umap_df)[1:2] <- c("UMAP1", "UMAP2")
  umap_df$species_info <- combined_seurat_object$species_info
  
  umap_df$highlight <- ifelse(umap_df$species_info == sp, sp, "Others")
  umap_df <- umap_df[sample(nrow(umap_df)), ]
  
  plot_colors <- c()
  plot_colors[sp] <- species_grouped_colors[sp]
  plot_colors["Others"] <- "grey90"
  
  p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
    geom_point(data = subset(umap_df, highlight == "Others"),
               aes(color = highlight), size = 0.1, alpha = 0.3) +
    geom_point(data = subset(umap_df, highlight == sp),
               aes(color = highlight), size = 0.1, alpha = 1) +
    scale_color_manual(values = plot_colors) +
    ggtitle(sp) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "italic", size = 20, hjust = 0.5),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_blank()
    )
  
  ggsave(file.path(output_dir, paste0("UMAP_", sp, ".jpg")),
         plot = p, width = 5, height = 5, dpi = 300)
}

# Calculate genus abundance and select the top 15 genera
genus_abundance <- combined_seurat_object@meta.data %>%
  group_by(genus_info) %>%
  summarise(cell_count = n()) %>%
  arrange(desc(cell_count))

top15_genus <- head(genus_abundance$genus_info, 15)

# Create a grouped genus column where non-top genera are labeled as Others
combined_seurat_object$top15_genus <- ifelse(
  combined_seurat_object$genus_info %in% top15_genus,
  combined_seurat_object$genus_info, "Others"
)

combined_seurat_object$top15_genus <- factor(
  combined_seurat_object$top15_genus,
  levels = c(top15_genus, "Others")
)

# Generate color palette for the top genera and Others
top_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(15)
genus_grouped_colors <- setNames(c(top_colors, "grey90"), levels(combined_seurat_object$top15_genus))

# Extract UMAP coordinates and compute ranges for mini-axis
umap_coords <- Embeddings(combined_seurat_object, "umap") %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("cell")

colnames(umap_coords)[2:3] <- c("UMAP_1", "UMAP_2")
umap_coords$genus <- combined_seurat_object$top15_genus[umap_coords$cell]
umap_coords <- umap_coords[sample(nrow(umap_coords)), ]

x_min <- min(umap_coords$UMAP_1)
x_max <- max(umap_coords$UMAP_1)
y_min <- min(umap_coords$UMAP_2)
y_max <- max(umap_coords$UMAP_2)

x_range <- x_max - x_min
y_range <- y_max - y_min

axis_len_x <- 0.2 * x_range
axis_len_y <- 0.2 * y_range
axis_x0 <- x_min - 0.05 * x_range
axis_y0 <- y_min - 0.05 * y_range
axis_arrow <- arrow(type = "closed", length = unit(5, "mm"))

# Generate UMAP plot colored by the top 15 genera
p5 <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = genus), size = 0.1, alpha = 1) +
  scale_color_manual(values = genus_grouped_colors) +
  labs(title = "UMAP by Top 15 Genus") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"),
    plot.title.position = "panel",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.text = element_text(size = 15, face = "italic"),
    legend.key.height = unit(1.5, "lines"),
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 120, 120, unit = "pt")
  ) +
  guides(color = guide_legend(
    ncol = 1,
    override.aes = list(size = 6),
    title = NULL
  )) +
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x,
           y = axis_y0, yend = axis_y0,
           linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0,
           y = axis_y0, yend = axis_y0 + axis_len_y,
           linewidth = 2, arrow = axis_arrow) +
  annotate("text", x = axis_x0 + axis_len_x / 2,
           y = axis_y0 - 0.05 * y_range,
           label = "UMAP 1", size = 8, family = "Arial", vjust = 1) +
  annotate("text", x = axis_x0 - 0.05 * x_range,
           y = axis_y0 + axis_len_y / 2,
           label = "UMAP 2", size = 8, family = "Arial", angle = 90, hjust = 0.5)

p5

# Save UMAP figure as editable vector PDF
ggsave("FigureS6_genus.pdf",
       plot = p5,
       device = cairo_pdf,
       width = 11.5,
       height = 10,
       units = "in",
       limitsize = FALSE)

# Recalculate abundance to retrieve the top 15 genera
genus_abundance <- combined_seurat_object@meta.data %>%
  group_by(genus_info) %>%
  summarise(cell_count = n()) %>%
  arrange(desc(cell_count))

top15_genus <- head(genus_abundance$genus_info, 15)

# Define color mapping for genus highlighting
top_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(15)
genus_grouped_levels <- c(top15_genus, "Others")
genus_grouped_colors <- setNames(c(top_colors, "grey90"), genus_grouped_levels)

# Create directory for genus-specific UMAP plots
output_dir <- "UMAP_top15_genus"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Generate individual UMAP plots highlighting each genus
for (gn in top15_genus) {
  
  umap_df <- Embeddings(combined_seurat_object, "umap") %>% as.data.frame()
  colnames(umap_df)[1:2] <- c("UMAP1", "UMAP2")
  umap_df$genus_info <- combined_seurat_object$genus_info
  umap_df$highlight <- ifelse(umap_df$genus_info == gn, gn, "Others")
  
  plot_colors <- c()
  plot_colors[gn] <- genus_grouped_colors[gn]
  plot_colors["Others"] <- "grey90"
  
  p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
    geom_point(data = subset(umap_df, highlight == "Others"),
               aes(color = highlight), size = 0.1, alpha = 0.3) +
    geom_point(data = subset(umap_df, highlight == gn),
               aes(color = highlight), size = 0.1, alpha = 1) +
    scale_color_manual(values = plot_colors) +
    ggtitle(gn) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "italic", size = 20, hjust = 0.5),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_blank()
    )
  
  ggsave(file.path(output_dir, paste0("UMAP_", gn, ".jpg")),
         plot = p, width = 5, height = 5, dpi = 300)
}