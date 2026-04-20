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
if (length(parts) >= 2 && all(tail(parts, 2) == c("Space_Scripts", "Figure6"))) {
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

saveRDS(
  combined_seurat_object,
  file = "combined_seurat_object.rds"
)

combined_seurat_object <- readRDS(
  file = "combined_seurat_object.rds"
)

set.seed(1024)

# Filter cells to include only the target species for PCA analysis
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

cat("Cells retained after species filtering:", ncol(combined_seurat_object), "\n")

# Perform PCA on the filtered Seurat object
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
  # Extract PCA standard deviations
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop("Specified reduction not found")
  }
  pca_sd <- seurat_obj[[reduction]]@stdev[1:ndims]
  
  # Normalize coordinates for elbow detection
  x <- 1:ndims
  y <- pca_sd
  x_norm <- (x - min(x)) / (max(x) - min(x))
  y_norm <- (y - min(y)) / (max(y) - min(y))
  
  # Compute perpendicular distance from the diagonal
  distance <- abs(y_norm - (1 - x_norm)) / sqrt(2)
  
  # Determine elbow point
  elbow_point <- which.max(distance)
  
  # Visualize PCA elbow plot
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
  
  print(p)
  cat("Suggested number of PCs:", elbow_point, "\n")
  return(elbow_point)
}

# Identify the optimal number of principal components
elbow_pc <- find_pca_elbow(combined_seurat_object)

# Run UMAP dimensional reduction and clustering on selected PCs
combined_seurat_object <- combined_seurat_object %>%
  RunUMAP(reduction = "pca", dims = 1:elbow_pc) %>%
  FindNeighbors(reduction = "pca", dims = 1:elbow_pc) %>%
  FindClusters(resolution = 0.2)

cat("UMAP and clustering completed. Object dimensions:", dim(combined_seurat_object), "\n")

# Prepare Seurat object for hdWGCNA analysis
theme_set(theme_cowplot())
enableWGCNAThreads(nThreads = 16)

# Select highly variable genes and construct metacells
combined_seurat_object <- SetupForWGCNA(
  combined_seurat_object,
  gene_select = "fraction",
  fraction = 0.001,
  wgcna_name = "combined"
)

selected_genes <- GetWGCNAGenes(combined_seurat_object, wgcna_name = "combined")
cat("Number of selected genes:", length(selected_genes), "\n")

# Build metacells based on cluster identity
combined_seurat_object <- MetacellsByGroups(
  seurat_obj = combined_seurat_object,
  reduction = "pca",
  k = 25,
  max_shared = 10,
  ident.group = "seurat_clusters"
)

# Normalize metacell expression
combined_seurat_object <- NormalizeMetacells(combined_seurat_object)

# Subset Seurat object to selected genes for WGCNA
combined_seurat_object <- subset(
  combined_seurat_object,
  features = selected_genes
)

# Set the expression matrix for WGCNA analysis
combined_seurat_object <- SetDatExpr(
  combined_seurat_object,
  assay = "SCT",
  layer = "data"
)

# Test candidate soft-threshold powers
combined_seurat_object <- TestSoftPowers(
  combined_seurat_object,
  networkType = "signed"
)

# Visualize diagnostics for soft-thresholding powers
plot_list <- PlotSoftPowers(combined_seurat_object)
withr::with_pdf(
  file = "FigureS8a.pdf",
  width = 10, 
  height = 9,
  code = {
    p_all <- wrap_plots(plot_list, ncol = 2)
    print(p_all)
  }
)

# Inspect recommended power table
power_table <- GetPowerTable(combined_seurat_object)
head(power_table)

# Construct co-expression network and detect modules
combined_seurat_object <- ConstructNetwork(
  combined_seurat_object,
  tom_name = "combined_network",
  overwrite_tom = TRUE,
  minModuleSize = 30,
  mergeCutHeight = 0.25
)

# Load TOM matrix for downstream analysis
tom_file <- "TOM/combined_network_TOM.rda"
load(tom_file)

# Extract gene names and build edge list for network
gene_names <- combined_seurat_object@misc$combined$wgcna_genes
TOM_mat <- as.matrix(consTomDS)
rownames(TOM_mat) <- gene_names
colnames(TOM_mat) <- gene_names
TOM_mat[lower.tri(TOM_mat, diag = TRUE)] <- NA

edge_df <- as.data.frame(as.table(TOM_mat))
colnames(edge_df) <- c("Gene1", "Gene2", "Weight")
edge_df <- edge_df[!is.na(edge_df$Weight), ]
edge_df <- subset(edge_df, Weight >= 0)

# Save network edges and node/module annotations
write.csv(edge_df, "combined_WGCNA_edges.csv", row.names = FALSE, quote = FALSE)
node_df <- GetModules(combined_seurat_object)
write.csv(node_df, "combined_WGCNA_nodes.csv", row.names = FALSE, quote = FALSE)

# Save dendrogram of genes
withr::with_pdf(
  file = "FigureS8b.pdf",
  width = 8,
  height = 5,
  code = {
    PlotDendrogram(combined_seurat_object, main = "Gene dendrogram")
  }
)

# Compute module eigengenes and gene connectivity
combined_seurat_object <- ModuleEigengenes(combined_seurat_object)
combined_seurat_object <- ModuleConnectivity(combined_seurat_object)

# Visualize module eigengene correlations
PlotKMEs(combined_seurat_object, ncol = 5)

# Identify hub genes in each module and save top hubs
modules <- GetModules(combined_seurat_object) %>% subset(module != "grey")
hub_df <- GetHubGenes(combined_seurat_object, n_hubs = 10)
write.csv(hub_df, "hub_genes_top10.csv", row.names = FALSE, quote = FALSE)
hub_df <- GetHubGenes(combined_seurat_object, n_hubs = 50)
write.csv(hub_df, "hub_genes_top50.csv", row.names = FALSE, quote = FALSE)

# Export hub gene expression matrix transposed for downstream analysis
hub_genes <- read.csv("hub_genes_top50.csv")
expr_matrix <- GetAssayData(combined_seurat_object, assay = "SCT", slot = "data")
hub_gene_expr_matrix <- expr_matrix[hub_genes$gene_name, , drop = FALSE]
hub_gene_expr_df_transposed <- t(as.data.frame(hub_gene_expr_matrix))
write.csv(hub_gene_expr_df_transposed, "hub_genes_expression_matrix.csv", row.names = TRUE, quote = FALSE)

# Summarize hub gene expression by timepoint
timepoints <- combined_seurat_object$timepoint
timepoint_ids <- sort(unique(timepoints))
fraction_mat <- matrix(0, nrow = length(hub_genes$gene_name), ncol = length(timepoint_ids),
                       dimnames = list(hub_genes$gene_name, paste0(timepoint_ids)))
mean_mat <- fraction_mat

for (tid in timepoint_ids) {
  timepoint_cells <- names(timepoints[timepoints == tid])
  sub_counts <- hub_gene_expr_matrix[, timepoint_cells, drop = FALSE]
  fraction_mat[, paste0(tid)] <- rowSums(sub_counts > 0) / length(timepoint_cells)
  mean_mat[, paste0(tid)] <- rowMeans(sub_counts)
}

write.csv(fraction_mat, "hub_gene_fraction_by_timepoint.csv", quote = FALSE)
write.csv(mean_mat, "hub_gene_mean_expr_by_timepoint.csv", quote = FALSE)

# Compute module expression scores and visualize patterns
combined_seurat_object <- ModuleExprScore(combined_seurat_object, n_genes = 50, method = "UCell")
plot_list <- ModuleFeaturePlot(combined_seurat_object, features = "scores", order = "shuffle", ucell = TRUE)
plot_list <- lapply(plot_list, function(p) p + geom_point(size = 0.01) + theme(plot.title = element_text(size = 30, face = "bold")))

withr::with_pdf(file = "Figure6a.pdf", width = 12, height = 10, code = {
  p_all <- wrap_plots(plot_list, ncol = 3)
  print(p_all)
})

# Visualize module radar plots grouped by timepoint
combined_seurat_object$cluster <- as.character(combined_seurat_object$seurat_clusters)
withr::with_pdf(file = "Figure6b.pdf", width = 8, height = 8, code = {
  p <- ModuleRadarPlot(combined_seurat_object, group.by = "timepoint", axis.label.size = 4, grid.label.size = 0)
  print(p)
})

# Save the final hdWGCNA Seurat object
saveRDS(combined_seurat_object, file = "combined_seurat_object_hdWGCNA.rds")
cat("hdWGCNA analysis completed and saved\n")

# set selected timepoints and gene list
selected_timepoint <- "FD30"
gene_list <- c("BVU-RS05735", "BVU-RS06335", "MGYG000001306-00848", "MGYG000002478-01902")
expr_threshold <- 0 

# loop over all genes in the list
for (selected_gene in gene_list) {
  
  # extract UMAP coordinates
  umap_df <- Embeddings(combined_seurat_object, "umap") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell")
  colnames(umap_df)[2:3] <- c("UMAP_1", "UMAP_2")
  
  # get gene expression matrix
  expr_mat <- GetAssayData(combined_seurat_object, assay = "SCT", slot = "data")
  if (!(selected_gene %in% rownames(expr_mat))) {
    warning(paste("Gene", selected_gene, "not found in expression matrix! Skipping."))
    next
  }
  
  # build color groups
  umap_df$ColorGroup <- "Other Timepoints"
  fd_cells <- rownames(combined_seurat_object@meta.data)[combined_seurat_object$timepoint == selected_timepoint]
  umap_df$ColorGroup[umap_df$cell %in% fd_cells] <- "Timepoint FD30"
  high_cells <- colnames(expr_mat)[expr_mat[selected_gene, ] > expr_threshold]
  umap_df$ColorGroup[umap_df$cell %in% intersect(fd_cells, high_cells)] <- "Selected Gene"
  
  # define color mapping
  color_map <- c(
    "Other Timepoints" = "#E0E0E0",
    "Timepoint FD30"   = "#A9A9A9",
    "Selected Gene"    = "#B22222"
  )
  umap_df$ColorGroup <- factor(umap_df$ColorGroup,
                               levels = c("Other Timepoints", "Timepoint FD30", "Selected Gene"))
  
  # compute UMAP range for mini axes
  x_min <- min(umap_df$UMAP_1); x_max <- max(umap_df$UMAP_1)
  y_min <- min(umap_df$UMAP_2); y_max <- max(umap_df$UMAP_2)
  x_range <- x_max - x_min; y_range <- y_max - y_min
  axis_len_x <- 0.2 * x_range; axis_len_y <- 0.2 * y_range
  axis_x0 <- x_min - 0.10 * x_range; axis_y0 <- y_min - 0.10 * y_range
  axis_arrow <- arrow(type = "closed", length = unit(6, "mm"))
  
  # split cells by color group
  gray_cells <- umap_df %>% filter(ColorGroup == "Other Timepoints")
  darkgray_cells <- umap_df %>% filter(ColorGroup == "Timepoint FD30")
  red_cells <- umap_df %>% filter(ColorGroup == "Selected Gene")
  
  # plot UMAP with layered points
  p <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = ColorGroup)) +
    geom_point(data = gray_cells, size = 1) +
    geom_point(data = darkgray_cells, size = 1) +
    geom_point(data = red_cells, size = 1.2, shape = 16) +
    annotate("segment", x = axis_x0, xend = axis_x0 + axis_len_x, y = axis_y0, yend = axis_y0, linewidth = 1.5, arrow = axis_arrow) +
    annotate("segment", x = axis_x0, xend = axis_x0, y = axis_y0, yend = axis_y0 + axis_len_y, linewidth = 1.5, arrow = axis_arrow) +
    annotate("text", x = axis_x0 + axis_len_x / 2, y = axis_y0 - 0.05 * y_range, label = "UMAP 1", size = 6, vjust = 1) +
    annotate("text", x = axis_x0 - 0.05 * x_range, y = axis_y0 + axis_len_y / 2, label = "UMAP 2", size = 6, angle = 90, hjust = 0.5) +
    ggtitle(bquote(italic(.(selected_gene)))) +
    scale_color_manual(values = color_map, drop = FALSE) +
    guides(color = guide_legend(title = NULL, override.aes = list(size = 6), ncol = 1)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 30),
          axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
          panel.grid = element_blank(), panel.border = element_blank(),
          legend.position = "right", legend.text = element_text(size = 20))
  
  # save plot as PDF
  pdf_file <- paste0("Figure6d_", selected_gene, ".pdf")
  ggsave(pdf_file, plot = p, width = 8.5, height = 6, dpi = 300)
  
  print(p)
}
