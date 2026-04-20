# ===============================================
# Author : ZHENG XINGHAI
# Date   : 2025-11-30
# Project: Astronaut mscRNA-seq Analysis
# ===============================================

library(limma)
library(ggplot2)
library(ggrepel)

# Set the working directory to the location of the current script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the scaled species abundance matrix
df <- read.csv("species_abundance_profile_scaled.csv", row.names = 1)

# Define column names corresponding to L and FD stages
L_columns <- c("L_1", "L_2", "L_3")
FD_columns <- c("FD_1", "FD_2", "FD_3")

# Create sample metadata describing the experimental condition of each column
sample_info <- data.frame(
  condition = factor(
    c(rep("L", length(L_columns)), rep("FD", length(FD_columns)))
  )
)

# Construct the design matrix for the linear model
design <- model.matrix(~ 0 + condition, data = sample_info)
colnames(design) <- c("L", "FD")

# Fit the linear model using the limma framework
fit <- lmFit(df[, c(L_columns, FD_columns)], design)
contrast.matrix <- makeContrasts(FDvsL = FD - L, levels = design)

# Apply contrast fitting and empirical Bayes moderation
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract differential abundance results for all species
res_FD_L <- topTable(
  fit2,
  coef = "FDvsL",
  adjust = "fdr",
  number = nrow(df)
)

# Save the differential analysis results to a CSV file
write.csv(
  as.data.frame(res_FD_L),
  file = "FD_vs_L_results.csv"
)

# Prepare a data frame for volcano plot visualization
res_df <- res_FD_L
res_df$Species <- rownames(res_df)

# Assign significance categories based on P-value and log fold change
res_df$group <- "Not_sig"
res_df$group[res_df$P.Value < 0.05 & res_df$logFC > 0] <- "FD_down"
res_df$group[res_df$P.Value < 0.05 & res_df$logFC < 0] <- "FD_up"

# Generate shortened species labels for cleaner annotation
res_df$label_short <- sapply(res_df$Species, function(x) {
  parts <- strsplit(x, " ")[[1]]
  if (grepl("UBA", x)) {
    return(parts[1])
  } else {
    first_part <- paste0(toupper(substr(parts[1], 1, 1)), ". ")
    if (length(parts) > 1) {
      rest <- paste(parts[2:length(parts)], collapse = " ")
      return(paste0(first_part, rest))
    } else {
      return(first_part)
    }
  }
})

# Create a volcano plot to visualize differential species abundance
volcano_plot <- ggplot(
  res_df,
  aes(x = -logFC, y = -log10(P.Value))
) +
  geom_point(aes(color = group), size = 3, alpha = 0.8) +
  geom_text_repel(
    data = subset(res_df, group != "Not_sig"),
    aes(label = label_short),
    size = 6,
    fontface = "italic",
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf,
    segment.color = "grey50"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "darkred",
    size = 1.2
  ) +
  annotate(
    "text",
    x = 0,
    y = -log10(0.05) - 0.3,
    label = "P-value = 0.05",
    color = "darkred",
    fontface = "bold",
    size = 6
  ) +
  scale_color_manual(
    values = c(
      "FD_up" = "red",
      "FD_down" = "blue",
      "Not_sig" = "grey"
    )
  ) +
  theme_minimal() +
  labs(
    title = "FD vs L",
    x = "log2(Fold Change)",
    y = "-log10(P-value)"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 25),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18)
  )

# Export the volcano plot as a PDF figure
ggsave(
  filename = "Figure3f.pdf",
  plot = volcano_plot,
  device = "pdf",
  width = 8,
  height = 5
)