# ===============================================
# Author : ZHENG XINGHAI
# Date   : 2025-11-30
# Project: Astronaut mscRNA-seq Analysis
# ===============================================

library(gridExtra)
library(ggpubr)
library(patchwork)

# Set working directory to the script location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Function to retain valid species names
valid_index <- function(name) {
  if (length(unlist(str_extract_all(name, "[A-Z]"))) != 1) return(FALSE)
  if (str_detect(name, "[-0-9]")) return(FALSE)
  return(TRUE)
}

# Load species abundance tables
df  <- read_csv("Meta_species.csv")
df0 <- read_csv("MscRNA_species.csv")

# Use first column as row names
df  <- df  %>% column_to_rownames(var = colnames(df)[1])
df0 <- df0 %>% column_to_rownames(var = colnames(df0)[1])

# Standardize species naming
rownames(df0) <- gsub(" ", "_", rownames(df0))

# Keep only valid species
df  <- df[ sapply(rownames(df), valid_index), ]
df0 <- df0[ sapply(rownames(df0), valid_index), ]

# Align samples shared between datasets
keep_cols <- intersect(colnames(df), colnames(df0))
df  <- df[, keep_cols]
df0 <- df0[, keep_cols]

# Compute sample-wise Pearson correlations
pearson_result <- data.frame()
for (sample in colnames(df0)) {
  if (!(sample %in% colnames(df))) next
  sp0 <- rownames(df0)[df0[[sample]] > 0 & df0[[sample]] < 1]
  sp1 <- rownames(df)[ df[[sample]]  > 0 & df[[sample]]  < 1]
  inter <- intersect(sp0, sp1)
  if (length(inter) < 3) next
  v0 <- as.numeric(df0[inter, sample])
  v1 <- as.numeric(df[inter, sample])
  cor_test <- cor.test(v0, v1, method = "pearson")
  pearson_result <- rbind(
    pearson_result,
    data.frame(Sample = sample, R = as.numeric(cor_test$estimate), P = cor_test$p.value)
  )
}

# Prepare plotting data
volcano_df <- pearson_result %>%
  mutate(
    logP = -log10(P),
    Time = factor(sub("_.*", "", Sample), levels = c("L-30", "FD30", "FD90", "FD150", "R+1", "R+7")),
    Subject = factor(sub(".*_", "", Sample), levels = c("1", "2", "3"))
  )

# Define color and shape mappings
time_cols <- c("L-30"="#E6D8AD","FD30"="#5F6F75","FD90"="#81B295","FD150"="#A3BFA8","R+1"="#9A031E","R+7"="#BA3A26")
subject_shapes <- c("1"=21,"2"=22,"3"=24)

# Generate volcano-style correlation plot
p_volcano <- ggplot(volcano_df, aes(x=R, y=logP, fill=Time, shape=Subject)) +
  geom_point(size=4.5, stroke=0.6, colour="black") +
  scale_fill_manual(values=time_cols) +
  scale_shape_manual(values=subject_shapes) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", linewidth=0.5, color="grey40") +
  labs(x="Pearson correlation (R)", y=expression(-log[10](italic(P))), fill="Timepoint", shape="Astronaut") +
  theme_classic(base_size=15) +
  theme(
    axis.line = element_line(linewidth=0.8),
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size=18),
    axis.text.x  = element_text(size=15),
    axis.text.y  = element_text(size=15),
    legend.position="right"
  ) +
  guides(
    fill = guide_legend(override.aes=list(shape=21, colour="black", size=5)),
    shape = guide_legend(override.aes=list(fill="grey80", colour="black", size=5))
  )

# Export figure
ggsave("Figure3b.pdf", p_volcano, width=5.6, height=4)