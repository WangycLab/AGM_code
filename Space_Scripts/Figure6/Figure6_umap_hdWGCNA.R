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

# Set project paths
current_dir <- trimws(getwd())
setwd(current_dir)
cat("Working directory set to:", current_dir, "\n")

parts <- strsplit(current_dir, .Platform$file.sep)[[1]]
if (length(parts) >= 2 && all(tail(parts, 2) == c("Space_Scripts", "Figure6"))) project_root <- paste(parts[1:(length(parts) - 2)], collapse = .Platform$file.sep) else project_root <- current_dir
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

# Export cell metadata
Idents(combined_seurat_object) <- "seurat_clusters"
md <- combined_seurat_object@meta.data
if (!"seurat_clusters" %in% colnames(md)) stop("No seurat_clusters column found")

output_df <- data.frame(
  barcode = rownames(md),
  cluster = as.character(md$seurat_clusters),
  sample = md$orig.ident,
  astronaut = md$astronaut,
  timepoint = md$timepoint,
  stage = md$stage,
  gene_number = md$gene_number,
  UMI = md$total_counts,
  stringsAsFactors = FALSE
)

write.table(output_df, file = "cell_metadata.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Find cluster markers
combined_seurat_object <- PrepSCTFindMarkers(combined_seurat_object)
markers <- FindAllMarkers(combined_seurat_object, assay = "SCT", only.pos = TRUE, logfc.threshold = 0.25)

# Calculate the number of expressing cells for each marker
expression_matrix <- GetAssayData(combined_seurat_object, assay = "SCT", slot = "data")
cluster_ids <- Idents(combined_seurat_object)
cells_expressing <- c()

for (i in seq_len(nrow(markers))) {
  gene <- markers[i, "gene"]; cluster <- markers[i, "cluster"]
  cells_in_cluster <- names(cluster_ids)[cluster_ids == cluster]
  cells_expressing <- c(cells_expressing, sum(expression_matrix[gene, cells_in_cluster] > 0))
}

markers$cells_expressing <- cells_expressing
markers <- markers %>% arrange(cluster, desc(avg_log2FC))
write.table(markers, file = "all_clusters_DEGs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Calculate gene expression fraction and mean expression
counts <- GetAssayData(combined_seurat_object, assay = "SCT", slot = "data")
clusters <- combined_seurat_object$seurat_clusters
markers$gene <- gsub("_", "-", markers$gene)
genes <- unique(markers$gene); genes <- genes[genes %in% rownames(counts)]
cluster_ids <- sort(unique(clusters)); cluster_names <- paste0("Cluster_", cluster_ids)

fraction_mat <- matrix(0, nrow = length(genes), ncol = length(cluster_ids), dimnames = list(genes, cluster_names))
mean_mat <- matrix(0, nrow = length(genes), ncol = length(cluster_ids), dimnames = list(genes, cluster_names))

for (cid in cluster_ids) {
  cluster_cells <- names(clusters[clusters == cid]); sub_counts <- counts[genes, cluster_cells, drop = FALSE]
  fraction_mat[, paste0("Cluster_", cid)] <- rowSums(sub_counts > 0) / length(cluster_cells)
  mean_mat[, paste0("Cluster_", cid)] <- rowMeans(sub_counts)
}

write.csv(fraction_mat, "gene_fraction_by_cluster.csv", quote = FALSE)
write.csv(mean_mat, "gene_mean_expr_by_cluster.csv", quote = FALSE)

# Randomize cell order
random_cells <- sample(colnames(combined_seurat_object), ncol(combined_seurat_object), replace = FALSE)
combined_seurat_object <- combined_seurat_object[, random_cells]

Idents(combined_seurat_object) <- "seurat_clusters"

cluster_annotation <- c(
  "0" = "Cofactor metabolism|0",
  "1" = "Redox homeostasis|1",
  "2" = "Uncharacterized function|2",
  "3" = "Nutrient transport|3",
  "4" = "Mobile genetic elements|4",
  "5" = "Protein stress response|5"
)

# Annotate and export metadata
md <- combined_seurat_object@meta.data
if (!"seurat_clusters" %in% colnames(md)) stop("No seurat_clusters column found")
md$cluster_annotation <- cluster_annotation[as.character(md$seurat_clusters)]
table(md$cluster_annotation)

output_df <- data.frame(
  barcode = rownames(md), cluster = md$cluster_annotation, sample = md$orig.ident,
  astronaut = md$astronaut, timepoint = md$timepoint, stage = md$stage,
  gene_number = md$gene_number, UMI = md$total_counts, stringsAsFactors = FALSE
)
write.table(output_df, "cell_metadata_annotated.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Functional cluster UMAP
clusters <- levels(Idents(combined_seurat_object))
cluster_colors <- setNames(colorRampPalette(brewer.pal(12, "Paired"))(length(clusters)), clusters)

umap_coords <- Embeddings(combined_seurat_object, "umap") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell")

umap_coords$cluster <- as.character(Idents(combined_seurat_object)[umap_coords$cell])
umap_coords$cluster_label <- factor(cluster_annotation[umap_coords$cluster], levels = cluster_annotation)
umap_coords_random <- umap_coords[sample(nrow(umap_coords)), ]

x_min <- min(umap_coords$umap_1); x_max <- max(umap_coords$umap_1)
y_min <- min(umap_coords$umap_2); y_max <- max(umap_coords$umap_2)
x_range <- x_max - x_min; y_range <- y_max - y_min
axis_len_x <- 0.2 * x_range; axis_len_y <- 0.2 * y_range
axis_x0 <- x_min - 0.05 * x_range; axis_y0 <- y_min - 0.05 * y_range
axis_arrow <- arrow(type = "closed", length = unit(5, "mm"))

p1 <- ggplot(umap_coords_random, aes(umap_1, umap_2, color = cluster_label)) +
  geom_point(size = 0.1, stroke = 0) +
  scale_color_manual(values = setNames(cluster_colors, cluster_annotation)) +
  labs(title = "UMAP by Functional Cluster", color = NULL) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"),
    plot.title.position = "panel", axis.text = element_blank(), axis.ticks = element_blank(),
    axis.title = element_blank(), legend.text = element_text(size = 20, family = "Arial"),
    legend.key.height = unit(2, "lines"), panel.grid = element_blank(),
    plot.margin = margin(20, 20, 120, 120, unit = "pt")
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 8), title = NULL)) +
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0, linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y, linewidth = 2, arrow = axis_arrow) +
  annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range, label = "UMAP 1", size = 8, family = "Arial", vjust = 1) +
  annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2, label = "UMAP 2", size = 8, family = "Arial", angle = 90, hjust = 0.5)

ggsave("Figure6a.pdf", p1, device = cairo_pdf, width = 14, height = 10, units = "in", limitsize = FALSE)

# Gene expression UMAPs
selected_timepoint <- "FD30"
gene_list <- c("BVU-RS02475", "MGYG000001072-00455", "MGYG000001364-02154", "MGYG000003681-01367")
expr_threshold <- 0

umap_df <- Embeddings(combined_seurat_object, "umap") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell")
colnames(umap_df)[2:3] <- c("UMAP_1", "UMAP_2")

expr_mat <- GetAssayData(combined_seurat_object, assay = "SCT", slot = "data")
fd_cells <- rownames(combined_seurat_object@meta.data)[combined_seurat_object$timepoint == selected_timepoint]

for (selected_gene in gene_list) {
  if (!(selected_gene %in% rownames(expr_mat))) {
    warning(paste("Gene", selected_gene, "not found in expression matrix! Skipping."))
    next
  }
  
  # Define cell groups
  plot_df <- umap_df
  plot_df$ColorGroup <- "Other Timepoints"
  plot_df$ColorGroup[plot_df$cell %in% fd_cells] <- "Timepoint FD30"
  high_cells <- colnames(expr_mat)[expr_mat[selected_gene, ] > expr_threshold]
  plot_df$ColorGroup[plot_df$cell %in% intersect(fd_cells, high_cells)] <- "Selected Gene"
  plot_df$ColorGroup <- factor(plot_df$ColorGroup, levels = c("Other Timepoints", "Timepoint FD30", "Selected Gene"))
  
  # Plot colors
  color_map <- c(
    "Other Timepoints" = "#E0E0E0",
    "Timepoint FD30" = "#9E9E9E",
    "Selected Gene" = "#B22222"
  )
  
  x_min <- min(plot_df$UMAP_1); x_max <- max(plot_df$UMAP_1)
  y_min <- min(plot_df$UMAP_2); y_max <- max(plot_df$UMAP_2)
  x_range <- x_max - x_min; y_range <- y_max - y_min
  axis_len_x <- 0.2 * x_range; axis_len_y <- 0.2 * y_range
  axis_x0 <- x_min - 0.10 * x_range; axis_y0 <- y_min - 0.10 * y_range
  axis_arrow <- arrow(type = "closed", length = unit(6, "mm"))
  
  # Plot UMAP with explicit layer order
  p <- ggplot() +
    geom_point(data = subset(plot_df, ColorGroup == "Other Timepoints"),
               aes(UMAP_1, UMAP_2), color = color_map["Other Timepoints"], size = 1) +
    geom_point(data = subset(plot_df, ColorGroup == "Timepoint FD30"),
               aes(UMAP_1, UMAP_2), color = color_map["Timepoint FD30"], size = 1) +
    geom_point(data = subset(plot_df, ColorGroup == "Selected Gene"),
               aes(UMAP_1, UMAP_2), color = color_map["Selected Gene"], size = 1.2) +
    annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0, linewidth = 1.5, arrow = axis_arrow) +
    annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y, linewidth = 1.5, arrow = axis_arrow) +
    annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range, label = "UMAP 1", size = 6, vjust = 1) +
    annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2, label = "UMAP 2", size = 6, angle = 90, hjust = 0.5) +
    ggtitle(bquote(italic(.(selected_gene)))) +
    guides(color = "none") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 30),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      legend.position = "none"
    )
  
  ggsave(paste0("Figure6g_", selected_gene, ".pdf"), p, width = 8, height = 8, dpi = 300)
}

# Setup
theme_set(theme_cowplot())
enableWGCNAThreads(nThreads = 16)

combined_seurat_object <- SetupForWGCNA(
  combined_seurat_object, gene_select = "fraction", fraction = 0.001, wgcna_name = "combined"
)
selected_genes <- GetWGCNAGenes(combined_seurat_object, wgcna_name = "combined")
cat("Number of selected genes:", length(selected_genes), "\n")

# Metacells and expression matrix
combined_seurat_object <- MetacellsByGroups(
  seurat_obj = combined_seurat_object, reduction = "pca", k = 25, max_shared = 10, ident.group = "seurat_clusters"
)
combined_seurat_object <- NormalizeMetacells(combined_seurat_object)
combined_seurat_object <- subset(combined_seurat_object, features = selected_genes)
combined_seurat_object <- SetDatExpr(combined_seurat_object, assay = "SCT", layer = "data")

# Network construction
combined_seurat_object <- TestSoftPowers(combined_seurat_object, networkType = "signed")
plot_list <- PlotSoftPowers(combined_seurat_object)

withr::with_pdf("FigureS7a.pdf", width = 10, height = 9, code = {
  p_all <- wrap_plots(plot_list, ncol = 2)
  print(p_all)
})

power_table <- GetPowerTable(combined_seurat_object)
head(power_table)

combined_seurat_object <- ConstructNetwork(
  combined_seurat_object, tom_name = "combined_network", overwrite_tom = TRUE,
  minModuleSize = 30, mergeCutHeight = 0.25
)

# Export network
load("TOM/combined_network_TOM.rda")

gene_names <- combined_seurat_object@misc$combined$wgcna_genes
TOM_mat <- as.matrix(consTomDS)
rownames(TOM_mat) <- gene_names; colnames(TOM_mat) <- gene_names
TOM_mat[lower.tri(TOM_mat, diag = TRUE)] <- NA

edge_df <- as.data.frame(as.table(TOM_mat))
colnames(edge_df) <- c("Gene1", "Gene2", "Weight")
edge_df <- subset(edge_df, !is.na(Weight) & Weight >= 0)

write.csv(edge_df, "combined_WGCNA_edges.csv", row.names = FALSE, quote = FALSE)
write.csv(GetModules(combined_seurat_object), "combined_WGCNA_nodes.csv", row.names = FALSE, quote = FALSE)

withr::with_pdf("FigureS7b.pdf", width = 8, height = 5, code = {
  PlotDendrogram(combined_seurat_object, main = "Gene dendrogram")
})

# Module analysis
combined_seurat_object <- ModuleEigengenes(combined_seurat_object)
combined_seurat_object <- ModuleConnectivity(combined_seurat_object)
PlotKMEs(combined_seurat_object, ncol = 5)

modules <- GetModules(combined_seurat_object) %>% subset(module != "grey")

hub_df <- GetHubGenes(combined_seurat_object, n_hubs = 10)
write.csv(hub_df, "hub_genes_top10.csv", row.names = FALSE, quote = FALSE)

hub_df <- GetHubGenes(combined_seurat_object, n_hubs = 50)
write.csv(hub_df, "hub_genes_top50.csv", row.names = FALSE, quote = FALSE)

hub_df <- GetHubGenes(combined_seurat_object, n_hubs = Inf)
write.csv(hub_df, "hub_genes_all.csv", row.names = FALSE, quote = FALSE)

# Hub gene expression
hub_genes <- read.csv("hub_genes_top50.csv")
expr_matrix <- GetAssayData(combined_seurat_object, assay = "SCT", slot = "data")
hub_gene_expr_matrix <- expr_matrix[hub_genes$gene_name, , drop = FALSE]
write.csv(t(as.data.frame(hub_gene_expr_matrix)), "hub_genes_expression_matrix.csv", quote = FALSE)

# Summarize hub gene expression by timepoint
timepoints <- combined_seurat_object$timepoint
timepoint_ids <- sort(unique(timepoints))

fraction_mat <- matrix(
  0, nrow = length(hub_genes$gene_name), ncol = length(timepoint_ids),
  dimnames = list(hub_genes$gene_name, timepoint_ids)
)
mean_mat <- fraction_mat

for (tid in timepoint_ids) {
  timepoint_cells <- names(timepoints[timepoints == tid])
  sub_counts <- hub_gene_expr_matrix[, timepoint_cells, drop = FALSE]
  fraction_mat[, tid] <- rowSums(sub_counts > 0) / length(timepoint_cells)
  mean_mat[, tid] <- rowMeans(sub_counts)
}

write.csv(fraction_mat, "hub_gene_fraction_by_timepoint.csv", quote = FALSE)
write.csv(mean_mat, "hub_gene_mean_expr_by_timepoint.csv", quote = FALSE)

# Module scores and visualization
combined_seurat_object <- ModuleExprScore(combined_seurat_object, n_genes = 50, method = "UCell")

plot_list <- ModuleFeaturePlot(
  combined_seurat_object, features = "scores", order = "shuffle", ucell = TRUE
)
plot_list <- lapply(plot_list, function(p) p + geom_point(size = 0.01) + theme(plot.title = element_text(size = 30, face = "bold")))

# Module radar plots
combined_seurat_object$cluster <- as.character(combined_seurat_object$seurat_clusters)

withr::with_pdf("Figure6d.pdf", width = 8, height = 8, code = {
  p <- ModuleRadarPlot(combined_seurat_object, group.by = "timepoint", axis.label.size = 4, grid.label.size = 0)
  print(p)
})

withr::with_pdf("Figure6c.pdf", width = 8, height = 8, code = {
  p <- ModuleRadarPlot(combined_seurat_object, group.by = "seurat_clusters", axis.label.size = 4, grid.label.size = 0)
  print(p)
})

# Save and read Seurat object
saveRDS(combined_seurat_object, file.path(project_root, "combined_seurat_object_hdWGCNA.rds"))
combined_seurat_object <- readRDS(file.path(project_root, "combined_seurat_object_hdWGCNA.rds"))