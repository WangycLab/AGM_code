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
library(BPCells)

# Set project paths
current_dir <- trimws(getwd())
setwd(current_dir)
cat("Working directory set to:", current_dir, "\n")

parts <- strsplit(current_dir, .Platform$file.sep)[[1]]
if (length(parts) >= 2 && all(tail(parts, 2) == c("Space_Scripts", "Figure5"))) project_root <- paste(parts[1:(length(parts) - 2)], collapse = .Platform$file.sep) else project_root <- current_dir
cat("Project root:", project_root, "\n")

data_path <- file.path(project_root, "Space_Matrix")
cat("Data path:", data_path, "\n")

GLOBAL_PURITY_CUTOFF <- 0.6
name_list <- list.dirs(data_path, full.names = FALSE, recursive = FALSE)
seurat_list <- list()

initial_stats <- data.frame(Sample = character(), Genes = integer(), Cells = integer())
after_filter_stats <- data.frame(Sample = character(), Genes = integer(), Cells = integer())
after_seurat_stats <- data.frame(Sample = character(), Genes = integer(), Cells = integer())
expr_count_raw <- list()
expr_count_filtered <- list()
gene_number_list <- list()

# Cell filtering
pass_ttest <- function(fraction_total_reads, fraction_total_reads2, fraction_total_reads3) {
  values <- as.numeric(c(fraction_total_reads2, fraction_total_reads3))
  fraction_total_reads <- as.numeric(fraction_total_reads)
  values <- values[!is.na(values)]
  if (length(values) < 2 || var(values) == 0 || is.na(fraction_total_reads)) return(FALSE)
  test <- t.test(values, mu = fraction_total_reads)
  test$p.value <= 0.05 && fraction_total_reads > max(values)
}

# Process samples
for (sample_name in name_list) {
  cat("Processing sample:", sample_name, "\n")
  
  report_S_file <- file.path(data_path, paste0(sample_name, "_sc_taxonomy.report"))
  report_G_file <- file.path(data_path, paste0(sample_name, "_sc_taxonomy_G.report"))
  
  if (!file.exists(report_S_file) || !file.exists(report_G_file)) {
    cat("Missing files for:", sample_name, "\n")
    next
  }
  
  bacteria_info <- read.delim(report_S_file, header = TRUE)
  colnames(bacteria_info) <- c("BC", "name", "taxonomy_id", "taxonomy_lvl", "reads", "all_reads")
  
  bacteria_info_G <- read.delim(report_G_file, header = TRUE)
  colnames(bacteria_info_G) <- c("BC", "name", "taxonomy_id", "taxonomy_lvl", "reads", "all_reads")
  
  tax_df <- read.delim(report_S_file, header = TRUE, stringsAsFactors = FALSE)
  sample_dir <- file.path(data_path, sample_name)
  
  expr_mat <- tryCatch(Read10X(sample_dir), error = function(e) NULL)
  if (is.null(expr_mat) || ncol(expr_mat) == 0) {
    cat("Empty matrix:", sample_name, "\n")
    next
  }
  
  # Initial statistics
  initial_stats <- rbind(initial_stats, data.frame(Sample = sample_name, Genes = nrow(expr_mat), Cells = ncol(expr_mat)))
  expr_count_raw[[sample_name]] <- data.frame(cell_id = colnames(expr_mat), count = Matrix::colSums(expr_mat != 0))
  
  # Filter cells
  filtered_tax_df <- tax_df[tax_df$fraction_total_reads >= GLOBAL_PURITY_CUTOFF, ]
  filtered_tax_df <- filtered_tax_df[apply(filtered_tax_df[, c("fraction_total_reads", "fraction_total_reads2", "fraction_total_reads3")], 1, function(x) pass_ttest(x[1], x[2], x[3])), ]
  
  keep_cells <- intersect(colnames(expr_mat), filtered_tax_df$barcode)
  expr_mat <- expr_mat[, keep_cells, drop = FALSE]
  
  after_filter_stats <- rbind(after_filter_stats, data.frame(Sample = sample_name, Genes = nrow(expr_mat), Cells = ncol(expr_mat)))
  
  if (ncol(expr_mat) == 0) {
    cat("No cells left after filtering:", sample_name, "\n")
    next
  }
  
  # Create Seurat object
  sample_obj <- CreateSeuratObject(expr_mat, min.cells = 30, min.features = 30)
  
  BC_info <- bacteria_info[match(colnames(sample_obj), bacteria_info$BC), ]
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
  
  # Per-cell statistics
  sample_obj$gene_number <- Matrix::colSums(GetAssayData(sample_obj, slot = "counts") != 0)
  sample_obj$total_counts <- Matrix::colSums(GetAssayData(sample_obj, slot = "counts"))
  
  gene_number_list[[sample_name]] <- data.frame(cell_id = colnames(sample_obj), gene_number = sample_obj$gene_number, total_counts = sample_obj$total_counts, sample = sample_name)
  
  after_seurat_stats <- rbind(after_seurat_stats, data.frame(Sample = sample_name, Genes = nrow(sample_obj), Cells = ncol(sample_obj)))
  expr_count_filtered[[sample_name]] <- data.frame(cell_id = colnames(sample_obj), count = sample_obj$gene_number)
  
  seurat_list[[sample_name]] <- sample_obj
  
  cat("Finished:", sample_name, "| global purity cutoff =", GLOBAL_PURITY_CUTOFF, "| cells =", ncol(sample_obj), "\n")
}

# Save summary statistics
write.csv(initial_stats, "initial_stats.csv", row.names = FALSE)
write.csv(after_filter_stats, "after_filter_stats.csv", row.names = FALSE)
write.csv(after_seurat_stats, "after_seurat_stats.csv", row.names = FALSE)

# Add sample prefixes to barcodes
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

# Harmonize gene sets
all_genes <- unique(unlist(lapply(seurat_list, rownames)))
for (nm in names(seurat_list)) {
  obj <- seurat_list[[nm]]
  missing <- setdiff(all_genes, rownames(obj))
  if (length(missing) > 0) {
    zero <- matrix(0, nrow = length(missing), ncol = ncol(obj), dimnames = list(missing, colnames(obj)))
    counts_mat <- GetAssayData(obj, slot = "counts")
    obj <- CreateSeuratObject(counts = rbind(counts_mat, zero), meta.data = obj@meta.data)
  }
  seurat_list[[nm]] <- obj
}

# Merge samples
combined_seurat_object <- Reduce(merge, seurat_list)
cat("All samples merged. Object dimensions:", dim(combined_seurat_object), "\n")

combined_seurat_object <- subset(combined_seurat_object, subset = nCount_RNA > 0 & is.finite(nCount_RNA))

# Export cell metadata
md <- combined_seurat_object@meta.data
output_df <- data.frame(barcode = rownames(md), sample = md$orig.ident, astronaut = md$astronaut, timepoint = md$timepoint, stage = md$stage, species = md$species_info, genus = md$genus_info, gene_number = md$gene_number, UMI = md$total_counts, stringsAsFactors = FALSE)
write.table("cell_metadata0.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Filter species with at least 500 cells
species_table <- table(combined_seurat_object$species_info)
keep_species <- names(species_table[species_table >= 500])
cat("Kept species (>=500 cells):\n")
print(keep_species)

combined_seurat_object <- subset(combined_seurat_object, subset = species_info %in% keep_species)

cat("After filtering:\n")
cat("Number of cells:", ncol(combined_seurat_object), "\n")
cat("Number of species:", length(unique(combined_seurat_object$species_info)), "\n")

# SCTransform and PCA
combined_seurat_object <- SCTransform(combined_seurat_object, vars.to.regress = NULL, verbose = FALSE, return.only.var.genes = FALSE, variable.features.n = 2000)
cat("SCTransform completed. Final dimensions:", dim(combined_seurat_object), "\n")

# Save Seurat object
saveRDS(combined_seurat_object, file = file.path(project_root, "combined_seurat_object.rds"))

combined_seurat_object <- readRDS(file = file.path(project_root, "combined_seurat_object.rds"))

set.seed(1024)

# Filter target species
target_species <- "Phocaeicola vulgatus"
if (!"species_info" %in% colnames(combined_seurat_object@meta.data)) stop("species_info column not found in meta.data")
species_cells <- rownames(combined_seurat_object@meta.data)[combined_seurat_object@meta.data$species_info == target_species]
if (length(species_cells) == 0) stop(paste("No cells found for species:", target_species))
combined_seurat_object <- subset(combined_seurat_object, cells = species_cells)

# PCA
combined_seurat_object <- RunPCA(combined_seurat_object, features = VariableFeatures(combined_seurat_object))

# Detect PCA elbow
find_pca_elbow <- function(seurat_obj, reduction = "pca", ndims = 50, axis_title_size = 25, axis_text_size = 25, axis_font = "sans", title_size = 30, title_font = "sans") {
  if (!reduction %in% names(seurat_obj@reductions)) stop("Specified reduction not found")
  pca_sd <- seurat_obj[[reduction]]@stdev[1:ndims]
  x <- 1:ndims; y <- pca_sd
  x_norm <- (x - min(x)) / (max(x) - min(x)); y_norm <- (y - min(y)) / (max(y) - min(y))
  distance <- abs(y_norm - (1 - x_norm)) / sqrt(2)
  elbow_point <- which.max(distance)
  plot_df <- data.frame(PC = x, SD = y)
  p <- ggplot(plot_df, aes(PC, SD)) +
    geom_point(size = 3, color = "#2E86AB") +
    geom_line(color = "#2E86AB") +
    geom_vline(xintercept = elbow_point, linetype = "dashed", color = "red") +
    annotate("text", x = elbow_point + 1, y = max(y), label = paste0("Elbow = PC", elbow_point), color = "red", size = 8, hjust = 0) +
    labs(title = "PCA Elbow Plot", x = "Principal Component", y = "Standard Deviation") +
    theme_bw(base_size = 20) +
    theme(
      plot.title = element_text(family = title_font, size = title_size, face = "bold", hjust = 0.5),
      axis.title = element_text(family = axis_font, size = axis_title_size, face = "bold"),
      axis.text = element_text(family = axis_font, size = axis_text_size),
      panel.grid.major = element_line(color = "grey80", linetype = "dotted"),
      panel.grid.minor = element_blank()
    )
  print(p)
  cat("Suggested number of PCs:", elbow_point, "\n")
  elbow_point
}

elbow_pc <- find_pca_elbow(combined_seurat_object)

# UMAP and clustering
combined_seurat_object <- combined_seurat_object %>%
  RunUMAP(reduction = "pca", dims = 1:elbow_pc) %>%
  FindNeighbors(reduction = "pca", dims = 1:elbow_pc) %>%
  FindClusters(resolution = 0.1)

# Astronaut colors
astronaut_colors <- c("Astronaut 1" = "#FF9A8B", "Astronaut 2" = "#FFE0AC", "Astronaut 3" = "#BBDED6")

# Timepoint colors
time_colors <- c("L-60" = "#F4F1DE", "L-30" = "#E6D8AD", "FD30" = "#5F6F75", "FD90" = "#81B295", "FD150" = "#A3BFA8", "R+1" = "#9A031E", "R+7" = "#BA3A26", "R+14" = "#E07A5F")

# Stage colors
stage_colors <- c("L" = "#FF6F61", "FD" = "#6B5B95", "R" = "#88B04B")

# Convert Seurat object to Monocle3 CDS
cds <- as.cell_data_set(combined_seurat_object)
cds@clusters$UMAP$clusters <- Idents(combined_seurat_object)
cds@colData@listData$seurat_clusters <- Idents(combined_seurat_object)

# Convert counts matrix to BPCells and preprocess
cds <- convert_counts_matrix(cds, matrix_control = list(matrix_class = "BPCells", matrix_path = "Monocle3_BPCells"))
cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)

# Extract PCA variance
get_pc_variance <- function(cds) {
  if ("prop_var_expl" %in% names(cds@reduce_dim_aux)) return(cds@reduce_dim_aux$prop_var_expl)
  if ("PCA" %in% names(cds@reduce_dim_aux)) {
    pca_aux <- cds@reduce_dim_aux$PCA
    if ("prop_var_expl" %in% names(pca_aux)) return(pca_aux$prop_var_expl)
    if ("model" %in% names(pca_aux) && !is.null(pca_aux$model$sdev)) {
      var <- pca_aux$model$sdev^2
      return(var / sum(var))
    }
  }
  rd <- SingleCellExperiment::reducedDims(cds)
  if ("PCA" %in% names(rd)) {
    var <- apply(rd$PCA, 2, var)
    return(var / sum(var))
  }
  stop("Cannot locate PCA variance information.")
}

# Detect PCA elbow
find_elbow <- function(var_explained, plot = TRUE) {
  x <- seq_along(var_explained)
  x_norm <- (x - min(x)) / (max(x) - min(x))
  y_norm <- (var_explained - min(var_explained)) / (max(var_explained) - min(var_explained))
  elbow <- which.max(abs(y_norm - (1 - x_norm)) / sqrt(2))
  if (plot) {
    plot(x, var_explained, type = "b", pch = 16, xlab = "Principal Components", ylab = "Variance Explained", main = "Elbow Plot")
    abline(v = elbow, col = "red", lty = 2, lwd = 2)
    text(elbow, max(var_explained), labels = paste0("PC ", elbow), pos = 4, col = "red")
  }
  elbow
}

var_explained <- get_pc_variance(cds)
num_dim_opt <- find_elbow(var_explained)
cds <- preprocess_cds(cds, num_dim = num_dim_opt)

# UMAP, clustering, and trajectory
cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds)
cds <- learn_graph(cds, use_partition = FALSE)
cds <- order_cells(cds)

# Randomize cell order
set.seed(1024)
cds <- cds[, sample(colnames(cds), ncol(cds), replace = FALSE)]

# Extract UMAP coordinates and define mini-axis
umap_df <- as.data.frame(reducedDims(cds)$UMAP)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
x_range <- diff(range(umap_df$UMAP_1)); y_range <- diff(range(umap_df$UMAP_2))
axis_len_x <- 0.2 * x_range; axis_len_y <- 0.2 * y_range
axis_x0 <- min(umap_df$UMAP_1) - 0.05 * x_range; axis_y0 <- min(umap_df$UMAP_2) - 0.05 * y_range
axis_arrow <- arrow(type = "closed", length = unit(6, "mm"))

# Add mini-axis
add_mini_axis <- function(p) {
  p + coord_cartesian(clip = "off") +
    annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0, linewidth = 2, arrow = axis_arrow) +
    annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y, linewidth = 2, arrow = axis_arrow) +
    annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range, label = "UMAP 1", size = 8, family = "Arial", vjust = 1) +
    annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2, label = "UMAP 2", size = 8, family = "Arial", angle = 90, hjust = 0.5)
}

# Common trajectory plot settings
trajectory_theme <- theme_minimal() + theme(
  plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"),
  axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
  legend.text = element_text(size = 20, family = "Arial"), legend.title = element_blank(),
  legend.key.height = unit(2, "lines"), panel.grid = element_blank(),
  plot.margin = margin(10, 10, 30, 30, unit = "pt")
)

# Plot and save trajectory
plot_trajectory <- function(color_by, title, file, width, scale = NULL, labels = NULL) {
  p <- plot_cells(cds, color_cells_by = color_by, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, cell_size = 0.8, alpha = 0.7) +
    labs(title = title) + trajectory_theme +
    guides(color = guide_legend(override.aes = list(size = 10)))
  p <- add_mini_axis(p)
  if (!is.null(scale)) p <- p + scale_color_manual(values = scale, labels = if (is.null(labels)) names(scale) else labels)
  ggsave(file, p, device = cairo_pdf, width = width, height = 10, units = "in")
  p
}

# Pseudotime trajectory
p1 <- plot_trajectory("pseudotime", "Trajectory Plot by Pseudotime", "FigureS5a_pseudotime.pdf", 10.2)

# Timepoint trajectory
time_levels <- c("L-60", "L-30", "FD30", "FD90", "FD150", "R+1", "R+7", "R+14")
cds@colData$timepoint <- factor(cds@colData$timepoint, levels = time_levels)
p2 <- plot_trajectory("timepoint", "Trajectory Plot by Timepoint", "Figure5a_timepoint.pdf", 10.5, time_colors)

# Stage trajectory
cds@colData$stage <- factor(cds@colData$stage, levels = c("L", "FD", "R"))
p3 <- plot_trajectory("stage", "Trajectory Plot by Stage", "FigureS5b_stage.pdf", 11, stage_colors, c("Stage L", "Stage FD", "Stage R"))

# Save pseudotime values
pseudotime_df <- data.frame(cell = rownames(colData(cds)), pseudotime = pseudotime(cds))
write.csv(pseudotime_df, "species_pseudotime.csv", row.names = FALSE)

# Save and reload objects
saveRDS(combined_seurat_object, file = file.path(project_root, "combined_seurat_object_PV.rds"))
saveRDS(cds, file = file.path(project_root, "monocle3_cds_PV.rds"))
combined_seurat_object <- readRDS(file = file.path(project_root, "combined_seurat_object_PV.rds"))
cds <- readRDS(file = file.path(project_root, "monocle3_cds_PV.rds"))

# Identify pseudotime-associated genes
deg_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 10)
deg_genes_sig <- deg_genes %>% filter(q_value <= 0.05)
write.csv(deg_genes_sig, "species_pseudotime_deg.csv", row.names = TRUE)

# Export top DEG expression matrix with pseudotime
expr_mat_sub <- as(counts(cds)[rownames(deg_genes_sig), ], "dgCMatrix")
expr_df_t <- as.data.frame(as.matrix(t(expr_mat_sub)))
expr_df_t$pseudotime <- pseudotime(cds)
write.csv(expr_df_t, "species_pseudotime_deg_expression.csv", quote = FALSE)