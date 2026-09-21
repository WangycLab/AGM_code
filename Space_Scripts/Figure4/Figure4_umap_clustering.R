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

# Set project paths
current_dir <- trimws(getwd())
setwd(current_dir)
cat("Working directory set to:", current_dir, "\n")

parts <- strsplit(current_dir, .Platform$file.sep)[[1]]
if (length(parts) >= 2 && all(tail(parts, 2) == c("Space_Scripts", "Figure4"))) project_root <- paste(parts[1:(length(parts) - 2)], collapse = .Platform$file.sep) else project_root <- current_dir
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

combined_seurat_object <- RunPCA(combined_seurat_object, features = VariableFeatures(combined_seurat_object))
ElbowPlot(combined_seurat_object, ndims = 50, reduction = "pca")

# Automatically detect the PCA elbow point
find_pca_elbow <- function(seurat_obj, reduction = "pca", ndims = 50, axis_title_size = 25, axis_text_size = 25, axis_font = "sans", title_size = 30, title_font = "sans") {
  if (!reduction %in% names(seurat_obj@reductions)) stop("Specified reduction does not exist")
  
  pca_sd <- seurat_obj[[reduction]]@stdev[1:ndims]
  x <- 1:ndims
  y <- pca_sd
  x_norm <- (x - min(x)) / (max(x) - min(x))
  y_norm <- (y - min(y)) / (max(y) - min(y))
  distance <- abs(y_norm - (1 - x_norm)) / sqrt(2)
  elbow_point <- which.max(distance)
  
  plot_df <- data.frame(PC = x, SD = y)
  p <- ggplot(plot_df, aes(x = PC, y = SD)) +
    geom_point(size = 3, color = "#2E86AB") +
    geom_line(color = "#2E86AB") +
    geom_vline(xintercept = elbow_point, linetype = "dashed", color = "red") +
    annotate("text", x = elbow_point + 1, y = max(y), label = paste0("Elbow = PC", elbow_point), color = "red", size = 8, hjust = 0) +
    labs(title = "PCA Elbow Plot", x = "Principal Component", y = "Standard Deviation") +
    theme_bw(base_size = 20) +
    theme(plot.title = element_text(family = title_font, size = title_size, face = "bold", hjust = 0.5), axis.title = element_text(family = axis_font, size = axis_title_size, face = "bold"), axis.text = element_text(family = axis_font, size = axis_text_size), panel.grid.major = element_line(color = "grey80", linetype = "dotted"), panel.grid.minor = element_blank())
  
  print(p)
  cat("Suggested number of PCs:", elbow_point, "\n")
  return(elbow_point)
}

elbow_pc <- find_pca_elbow(combined_seurat_object)

# Harmony integration
set.seed(1024)
combined_seurat_object <- RunHarmony(object = combined_seurat_object, reduction = "pca", group.by.vars = "orig.ident", reduction.save = "harmony", assay.use = "SCT")

combined_seurat_object <- combined_seurat_object %>%
  RunUMAP(reduction = "harmony", dims = 1:elbow_pc) %>%
  FindNeighbors(reduction = "harmony", dims = 1:elbow_pc) %>%
  FindClusters(resolution = 0.1)

# Save Seurat object
saveRDS(combined_seurat_object, file = file.path(project_root, "combined_seurat_object.rds"))

combined_seurat_object <- readRDS(file = file.path(project_root, "combined_seurat_object.rds"))

# Set cluster identities
Idents(combined_seurat_object) <- "seurat_clusters"

# Export cell metadata
md <- combined_seurat_object@meta.data
if (!"seurat_clusters" %in% colnames(md)) stop("No seurat_clusters column found")

output_df <- data.frame(barcode = rownames(md), cluster = as.character(md$seurat_clusters), sample = md$orig.ident, astronaut = md$astronaut, timepoint = md$timepoint, stage = md$stage, species = md$species_info, genus = md$genus_info, gene_number = md$gene_number, UMI = md$total_counts, stringsAsFactors = FALSE)
write.table(output_df, file = file.path("cell_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Find cluster markers
combined_seurat_object <- PrepSCTFindMarkers(combined_seurat_object)
markers <- FindAllMarkers(combined_seurat_object, assay = "SCT", only.pos = TRUE, logfc.threshold = 0.25)

# Calculate the number of expressing cells
expression_matrix <- GetAssayData(combined_seurat_object, assay = "SCT", slot = "data")
cluster_ids <- Idents(combined_seurat_object)
cells_expressing <- c()

for (i in 1:nrow(markers)) {
  gene <- markers[i, "gene"]
  cluster <- markers[i, "cluster"]
  cells_in_cluster <- names(cluster_ids)[cluster_ids == cluster]
  expr_values <- expression_matrix[gene, cells_in_cluster]
  cells_expressing <- c(cells_expressing, sum(expr_values > 0))
}

markers$cells_expressing <- cells_expressing
markers <- markers %>% arrange(cluster, desc(avg_log2FC))
write.table(markers, file = file.path("all_clusters_DEGs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Calculate gene expression statistics by cluster
counts <- GetAssayData(combined_seurat_object, assay = "SCT", slot = "data")
clusters <- combined_seurat_object$seurat_clusters

markers$gene <- gsub("_", "-", markers$gene)
genes <- unique(markers$gene)
genes <- genes[genes %in% rownames(counts)]

cluster_ids <- sort(unique(clusters))
cluster_names <- paste0("Cluster_", cluster_ids)

fraction_mat <- matrix(0, nrow = length(genes), ncol = length(cluster_ids), dimnames = list(genes, cluster_names))
mean_mat <- matrix(0, nrow = length(genes), ncol = length(cluster_ids), dimnames = list(genes, cluster_names))

for (cid in cluster_ids) {
  cluster_cells <- names(clusters[clusters == cid])
  sub_counts <- counts[genes, cluster_cells, drop = FALSE]
  fraction_mat[, paste0("Cluster_", cid)] <- rowSums(sub_counts > 0) / length(cluster_cells)
  mean_mat[, paste0("Cluster_", cid)] <- rowMeans(sub_counts)
}

# Save expression matrices
write.csv(fraction_mat, file.path("gene_fraction_by_cluster.csv"), quote = FALSE)
write.csv(mean_mat, file.path("gene_mean_expr_by_cluster.csv"), quote = FALSE)

# Set cluster identities
Idents(combined_seurat_object) <- "seurat_clusters"

# Define functional annotations
cluster_annotation <- c(
  "0" = "Oxidative stress response|0", "1" = "Propionate production|1", "2" = "Cofactor biosynthesis|2",
  "3" = "Central carbon metabolism|3", "4" = "SCFA metabolism|4", "5" = "Oxidative stress adaptation|5",
  "6" = "Cell envelope remodeling|6", "7" = "Host interface adaptation|7", "8" = "Environmental sensing|8",
  "9" = "Anaerobic metabolism|9", "10" = "Envelope and motility|10", "11" = "Metabolic flexibility|11",
  "12" = "Propionate utilization|12", "13" = "Anaerobic energy metabolism|13", "14" = "Organic acid metabolism|14"
)

# Annotate and export metadata
md <- combined_seurat_object@meta.data
if (!"seurat_clusters" %in% colnames(md)) stop("No seurat_clusters column found")
md$cluster_annotation <- cluster_annotation[as.character(md$seurat_clusters)]
table(md$cluster_annotation)

output_df <- data.frame(barcode = rownames(md), cluster = md$cluster_annotation, sample = md$orig.ident, astronaut = md$astronaut, timepoint = md$timepoint, stage = md$stage, species = md$species_info, genus = md$genus_info, gene_number = md$gene_number, UMI = md$total_counts, stringsAsFactors = FALSE)
write.table(output_df, file = file.path("cell_metadata_annotated.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("Annotated metadata exported successfully!\n")

# Prepare colors
clusters <- levels(Idents(combined_seurat_object))
n_clusters <- length(clusters)
cluster_colors <- setNames(colorRampPalette(brewer.pal(12, "Paired"))(n_clusters), clusters)

# Extract UMAP coordinates
umap_coords <- Embeddings(combined_seurat_object, "umap") %>% as.data.frame() %>% tibble::rownames_to_column("cell")
umap_coords$cluster <- as.character(Idents(combined_seurat_object)[umap_coords$cell])
umap_coords$cluster_label <- cluster_annotation[umap_coords$cluster]
umap_coords$cluster_label <- factor(umap_coords$cluster_label, levels = cluster_annotation)
umap_coords_random <- umap_coords[sample(nrow(umap_coords)), ]

# Calculate mini-axis position
x_min <- min(umap_coords$umap_1); x_max <- max(umap_coords$umap_1)
y_min <- min(umap_coords$umap_2); y_max <- max(umap_coords$umap_2)
x_range <- x_max - x_min; y_range <- y_max - y_min
axis_len_x <- 0.2 * x_range; axis_len_y <- 0.2 * y_range
axis_x0 <- x_min - 0.05 * x_range; axis_y0 <- y_min - 0.05 * y_range
axis_arrow <- arrow(type = "closed", length = unit(5, "mm"))

# UMAP plot
p1 <- ggplot() +
  geom_point(data = umap_coords_random, aes(x = umap_1, y = umap_2, color = cluster_label), size = 0.1, stroke = 0, raster = FALSE) +
  scale_color_manual(values = setNames(cluster_colors, cluster_annotation)) +
  labs(title = "UMAP by Functional Cluster", color = NULL) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"),
    plot.title.position = "panel",
    axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
    legend.text = element_text(size = 20, family = "Arial"),
    legend.key.height = unit(2, "lines"),
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 120, 120, unit = "pt")
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 8), title = NULL)) +
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0, linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y, linewidth = 2, arrow = axis_arrow) +
  annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range, label = "UMAP 1", size = 8, family = "Arial", vjust = 1) +
  annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2, label = "UMAP 2", size = 8, family = "Arial", angle = 90, hjust = 0.5)

# Export PDF
ggsave(filename = file.path("Figure4b_cluster.pdf"), plot = p1, device = cairo_pdf, width = 15, height = 10, units = "in", limitsize = FALSE)

# Astronaut colors
astronaut_colors <- c("Astronaut 1" = "#FF9A8B", "Astronaut 2" = "#FFE0AC", "Astronaut 3" = "#BBDED6")

# UMAP by astronaut
umap_coords <- Embeddings(combined_seurat_object, "umap") %>% as.data.frame() %>% tibble::rownames_to_column("cell")
umap_coords$astronaut <- combined_seurat_object$astronaut[umap_coords$cell]
umap_random <- umap_coords[sample(nrow(umap_coords)), ]

x_min <- min(umap_coords$umap_1); x_max <- max(umap_coords$umap_1)
y_min <- min(umap_coords$umap_2); y_max <- max(umap_coords$umap_2)
x_range <- x_max - x_min; y_range <- y_max - y_min
axis_len_x <- 0.2 * x_range; axis_len_y <- 0.2 * y_range
axis_x0 <- x_min - 0.05 * x_range; axis_y0 <- y_min - 0.05 * y_range
axis_arrow <- arrow(type = "closed", length = unit(5, "mm"))

p2 <- ggplot() +
  geom_point(data = umap_random, aes(x = umap_1, y = umap_2, color = astronaut), size = 0.1, alpha = 1, stroke = 0, raster = FALSE) +
  scale_color_manual(values = astronaut_colors) +
  labs(title = "UMAP by Astronaut", color = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"), plot.title.position = "panel", axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), legend.text = element_text(size = 20, family = "Arial"), legend.key.height = unit(2, "lines"), legend.key.width = unit(1.5, "lines"), legend.position = "right", panel.grid = element_blank(), plot.margin = margin(20, 20, 120, 120, unit = "pt")) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 8), title = NULL)) +
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0, linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y, linewidth = 2, arrow = axis_arrow) +
  annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range, label = "UMAP 1", size = 8, family = "Arial", vjust = 1) +
  annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2, label = "UMAP 2", size = 8, family = "Arial", angle = 90, hjust = 0.5)

ggsave(filename = file.path("FigureS3f_astronaut.pdf"), plot = p2, device = cairo_pdf, width = 11.5, height = 10, units = "in", limitsize = FALSE)

# Timepoint colors
time_colors <- c("L-60" = "#F4F1DE", "L-30" = "#E6D8AD", "FD30" = "#5F6F75", "FD90" = "#81B295", "FD150" = "#A3BFA8", "R+1" = "#9A031E", "R+7" = "#BA3A26", "R+14" = "#E07A5F")
combined_seurat_object$timepoint <- factor(combined_seurat_object$timepoint, levels = names(time_colors))

# UMAP by timepoint
umap_coords <- Embeddings(combined_seurat_object, "umap") %>% as.data.frame() %>% tibble::rownames_to_column("cell")
umap_coords$timepoint <- combined_seurat_object$timepoint[umap_coords$cell]
umap_coords <- umap_coords[sample(nrow(umap_coords)), ]

x_min <- min(umap_coords$umap_1); x_max <- max(umap_coords$umap_1)
y_min <- min(umap_coords$umap_2); y_max <- max(umap_coords$umap_2)
x_range <- x_max - x_min; y_range <- y_max - y_min
axis_len_x <- 0.2 * x_range; axis_len_y <- 0.2 * y_range
axis_x0 <- x_min - 0.05 * x_range; axis_y0 <- y_min - 0.05 * y_range
axis_arrow <- arrow(type = "closed", length = unit(5, "mm"))

p3 <- ggplot(umap_coords, aes(x = umap_1, y = umap_2, color = timepoint)) +
  geom_point(size = 0.1, alpha = 1, stroke = 0) +
  scale_color_manual(values = time_colors) +
  labs(title = "UMAP by Timepoint", color = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"), plot.title.position = "panel", axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), legend.text = element_text(size = 20, family = "Arial"), legend.key.height = unit(2, "lines"), legend.key.width = unit(1.5, "lines"), legend.position = "right", panel.grid = element_blank(), plot.margin = margin(20, 20, 120, 120, unit = "pt")) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 8), title = NULL)) +
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0, linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y, linewidth = 2, arrow = axis_arrow) +
  annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range, label = "UMAP 1", size = 8, family = "Arial", vjust = 1) +
  annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2, label = "UMAP 2", size = 8, family = "Arial", angle = 90, hjust = 0.5)

ggsave(filename = file.path("FigureS3g_timepoint.pdf"), plot = p3, device = cairo_pdf, width = 10.8, height = 10, units = "in", limitsize = FALSE)

# Species colors
species_colors <- c(
  "Agathobacter faecis" = "#8A97A8", "Agathobacter rectalis" = "#B4C2D3",
  "Bacteroides eggerthii" = "#B4D1AC", "Bacteroides fragilis" = "#D8A27E", "Bacteroides intestinalis" = "#E8C8B0",
  "Bacteroides ovatus" = "#B8A878", "Bacteroides stercoris" = "#CBB987", "Bacteroides thetaiotaomicron" = "#D8C89A",
  "Bacteroides uniformis" = "#A89B68", "Bacteroides xylanisolvens" = "#E2D2B0",
  "CAG-81 sp900066785" = "#C8C2BC", "Clostridium_Q sp003024715" = "#C898A8",
  "Enterocloster bolteae" = "#BCA8C8", "Enterocloster sp000431375" = "#CDB6D6", "Enterocloster sp001517625" = "#A888B8",
  "Escherichia coli_D" = "#9674A8", "Faecalibacterium prausnitzii_C" = "#7FA6D8",
  "Fusicatenibacter saccharivorans" = "#D0B8D8", "Fusobacterium_A mortiferum" = "#8F8F78", "Fusobacterium_A varium" = "#B1B194",
  "Megamonas funiformis" = "#80A8A8", "Parabacteroides distasonis" = "#A8C8C8", "Parabacteroides merdae" = "#D88888",
  "Parasutterella excrementihominis" = "#F0B8B8", "Phascolarctobacterium faecium" = "#C9A6D8",
  "Phocaeicola dorei" = "#E8D090", "Phocaeicola massiliensis" = "#D8C4A0", "Phocaeicola vulgatus" = "#A88C78",
  "Prevotella stercorea" = "#E0C8B8", "Roseburia intestinalis" = "#8298B0", "Roseburia sp900552665" = "#B0C8D8",
  "Sutterella wadsworthensis" = "#98BC98", "UBA7182 sp003480725" = "#C48898"
)
species_colors <- c(species_colors, "Others" = "#BDBDBD")

# Select top 20 species
species_abundance <- combined_seurat_object@meta.data %>% group_by(species_info) %>% summarise(cell_count = n(), .groups = "drop") %>% arrange(desc(cell_count))
all_valid_species <- species_abundance$species_info[species_abundance$species_info %in% names(species_colors)]
all_valid_species <- setdiff(all_valid_species, "Others")
top20_species <- head(all_valid_species, 20)

combined_seurat_object$top_species <- ifelse(combined_seurat_object$species_info %in% top20_species, combined_seurat_object$species_info, "Others")
species_levels <- c(top20_species, "Others")
combined_seurat_object$top_species <- factor(combined_seurat_object$top_species, levels = species_levels)
species_colors_use <- species_colors[species_levels]

# UMAP by top 20 species
umap_coords <- Embeddings(combined_seurat_object, "umap") %>% as.data.frame() %>% tibble::rownames_to_column("cell")
umap_coords$species <- combined_seurat_object$top_species[umap_coords$cell]
umap_coords <- umap_coords[sample(nrow(umap_coords)), ]

x_min <- min(umap_coords$umap_1); x_max <- max(umap_coords$umap_1)
y_min <- min(umap_coords$umap_2); y_max <- max(umap_coords$umap_2)
x_range <- x_max - x_min; y_range <- y_max - y_min
axis_len_x <- 0.2 * x_range; axis_len_y <- 0.2 * y_range
axis_x0 <- x_min - 0.05 * x_range; axis_y0 <- y_min - 0.05 * y_range
axis_arrow <- arrow(type = "closed", length = unit(5, "mm"))

p4 <- ggplot(umap_coords, aes(x = umap_1, y = umap_2, color = species)) +
  geom_point(size = 0.1, alpha = 1, stroke = 0) +
  scale_color_manual(values = species_colors_use) +
  labs(title = "UMAP by Species", color = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"), plot.title.position = "panel", axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), legend.text = element_text(size = 15, face = "italic"), legend.key.height = unit(1.3, "lines"), legend.key.width = unit(1.3, "lines"), legend.position = "right", panel.grid = element_blank(), plot.margin = margin(20, 20, 120, 120, unit = "pt")) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 6), title = NULL)) +
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0, linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y, linewidth = 2, arrow = axis_arrow) +
  annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range, label = "UMAP 1", size = 8, family = "Arial", vjust = 1) +
  annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2, label = "UMAP 2", size = 8, family = "Arial", angle = 90, hjust = 0.5)

ggsave(filename = file.path("FigureS6e_species.pdf"), plot = p4, device = cairo_pdf, width = 13, height = 10, units = "in", limitsize = FALSE)

# Select top 15 genera
genus_abundance <- combined_seurat_object@meta.data %>% group_by(genus_info) %>% summarise(cell_count = n(), .groups = "drop") %>% arrange(desc(cell_count))
top15_genus <- head(genus_abundance$genus_info, 15)

combined_seurat_object$top15_genus <- ifelse(combined_seurat_object$genus_info %in% top15_genus, combined_seurat_object$genus_info, "Others")
combined_seurat_object$top15_genus <- factor(combined_seurat_object$top15_genus, levels = c(top15_genus, "Others"))

top_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(15)
genus_grouped_colors <- setNames(c(top_colors, "#BDBDBD"), levels(combined_seurat_object$top15_genus))

# UMAP by top 15 genera
umap_coords <- Embeddings(combined_seurat_object, "umap") %>% as.data.frame() %>% tibble::rownames_to_column("cell")
colnames(umap_coords)[2:3] <- c("UMAP_1", "UMAP_2")
umap_coords$genus <- combined_seurat_object$top15_genus[umap_coords$cell]
umap_coords <- umap_coords[sample(nrow(umap_coords)), ]

x_min <- min(umap_coords$UMAP_1); x_max <- max(umap_coords$UMAP_1)
y_min <- min(umap_coords$UMAP_2); y_max <- max(umap_coords$UMAP_2)
x_range <- x_max - x_min; y_range <- y_max - y_min
axis_len_x <- 0.2 * x_range; axis_len_y <- 0.2 * y_range
axis_x0 <- x_min - 0.05 * x_range; axis_y0 <- y_min - 0.05 * y_range
axis_arrow <- arrow(type = "closed", length = unit(5, "mm"))

p5 <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = genus), size = 0.1, alpha = 1, stroke = 0) +
  scale_color_manual(values = genus_grouped_colors) +
  labs(title = "UMAP by Top 15 Genus", color = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"), plot.title.position = "panel", axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), legend.text = element_text(size = 15, face = "italic"), legend.key.height = unit(1.5, "lines"), legend.position = "right", panel.grid = element_blank(), plot.margin = margin(20, 20, 120, 120, unit = "pt")) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 6), title = NULL)) +
  annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0, linewidth = 2, arrow = axis_arrow) +
  annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y, linewidth = 2, arrow = axis_arrow) +
  annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range, label = "UMAP 1", size = 8, family = "Arial", vjust = 1) +
  annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2, label = "UMAP 2", size = 8, family = "Arial", angle = 90, hjust = 0.5)

ggsave(filename = file.path("FigureS3d_genus.pdf"), plot = p5, device = cairo_pdf, width = 11.5, height = 10, units = "in", limitsize = FALSE)

# Define colors
astronaut_colors <- c("Astronaut 1" = "#FF9A8B", "Astronaut 2" = "#FFE0AC", "Astronaut 3" = "#BBDED6")

# UMAP split by astronaut
p2_split <- DimPlot(combined_seurat_object, group.by = "astronaut", split.by = "astronaut", reduction = "umap", label = FALSE, cols = astronaut_colors, pt.size = 0.1, alpha = 0.5, raster = FALSE, shuffle = TRUE, seed = 123) +
  scale_color_manual(values = astronaut_colors, breaks = names(astronaut_colors), labels = names(astronaut_colors)) +
  labs(title = "UMAP distribution across astronauts") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), legend.text = element_text(size = 20, family = "Arial"), legend.key.height = unit(2, "lines"), legend.key.width = unit(1.5, "lines"), panel.grid = element_blank(), panel.border = element_blank(), legend.position = "right", strip.text = element_text(size = 30, family = "Arial"), plot.margin = margin(10, 10, 90, 90, unit = "pt")) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 10)))

ggsave("FigureS3f_astronaut_split.pdf", plot = p2_split, device = cairo_pdf, width = 15, height = 6, units = "in")

# Define timepoint colors and order
time_colors <- c("L-60" = "#F4F1DE", "L-30" = "#E6D8AD", "FD30" = "#5F6F75", "FD90" = "#81B295", "FD150" = "#A3BFA8", "R+1" = "#9A031E", "R+7" = "#BA3A26", "R+14" = "#E07A5F")
combined_seurat_object$timepoint <- factor(combined_seurat_object$timepoint, levels = names(time_colors))

# UMAP split by timepoint
p3_split <- DimPlot(combined_seurat_object, group.by = "timepoint", split.by = "timepoint", reduction = "umap", label = FALSE, cols = time_colors, pt.size = 0.1, alpha = 0.5, raster = FALSE, shuffle = TRUE, seed = 123, ncol = 4) +
  scale_color_manual(values = time_colors, breaks = names(time_colors), labels = names(time_colors)) +
  labs(title = "UMAP distribution across timepoints") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 40, family = "Arial"), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), legend.text = element_text(size = 20, family = "Arial"), legend.key.height = unit(2, "lines"), legend.key.width = unit(1.5, "lines"), panel.grid = element_blank(), panel.border = element_blank(), legend.position = "right", strip.text = element_text(size = 25, family = "Arial"), plot.margin = margin(10, 10, 90, 90, unit = "pt")) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 10)))

ggsave("FigureS3g_timepoint_split.pdf", plot = p3_split, device = cairo_pdf, width = 24, height = 12, units = "in")