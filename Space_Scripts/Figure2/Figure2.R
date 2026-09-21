# ===============================================
# Author : DING LIGUO
# Date   : 2025-11-30
# Project: Astronaut mscRNA-seq Analysis
# ===============================================

library(Seurat)
library(dplyr)
library(tibble)
library(edgeR)
library(data.table)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(tidyr)
library(patchwork)

setwd(getwd())
sample_list <- c("Con","2mo","4mo","6mo")
matrix_dir <- file.path(dirname(dirname(getwd())),"Store_Matrix")
sample_cols <- c("Con"="#59A14F","2mo"="#4E79A7","4mo"="#F28E2B","6mo"="#E15759")

## Figure 2e
plot_df <- data.frame()
for(s in sample_list){
  cat("Processing:",s,"\n")
  obj <- CreateSeuratObject(Read10X(file.path(matrix_dir,s)),min.cells=10,min.features=10)
  plot_df <- rbind(plot_df,data.frame(sample=s,nFeature_RNA=obj$nFeature_RNA,nCount_RNA=obj$nCount_RNA))
}
plot_df <- plot_df %>% filter(nFeature_RNA<=1000)
plot_df$sample <- factor(plot_df$sample,levels=sample_list)

p2e <- ggplot(plot_df,aes(sample,nFeature_RNA,fill=sample)) +
  geom_violin(scale="width",trim=TRUE,color=NA,alpha=0.95) +
  geom_boxplot(width=0.12,outlier.shape=NA,fill="white",color="black",linewidth=0.3) +
  stat_summary(fun=median,geom="point",size=1.5,color="black") +
  scale_fill_manual(values=sample_cols) +
  scale_y_continuous(trans="log10") +
  labs(x=NULL,y="Detected genes per cell") +
  theme_classic(base_size=14) +
  theme(legend.position="none",axis.text.x=element_text(size=12,color="black"),axis.text.y=element_text(color="black"),axis.title.y=element_text(face="bold"))

ggsave("Figure2e.pdf",p2e,width=4,height=4.5,units="in",useDingbats=FALSE)

## Figure 2f
groups <- sample_list[-1]
con_df <- fread("Con_species.tsv") %>% select(name,fraction_total_reads) %>% rename(Con=fraction_total_reads)

species_all <- lapply(groups,function(g){
  tmp <- fread(paste0(g,"_species.tsv")) %>% select(name,fraction_total_reads) %>% rename(!!g:=fraction_total_reads)
  merged <- full_join(con_df,tmp,by="name")
  merged[is.na(merged)] <- 0
  merged$group <- g
  merged$group_prop <- merged[[g]]
  merged
})
names(species_all) <- groups
species_long <- bind_rows(species_all)

cor_label <- species_long %>%
  group_by(group) %>%
  summarise(R=cor(Con,group_prop,method="pearson"),P=cor.test(Con,group_prop,method="pearson")$p.value,.groups="drop")

p2f <- ggplot(species_long,aes(Con,group_prop)) +
  geom_point(aes(fill=group),shape=21,color="black",size=3.5,alpha=0.9,stroke=0.4) +
  geom_smooth(aes(color=group),method="lm",se=FALSE,linewidth=1) +
  scale_fill_manual(values=c("2mo"="#4E79A7","4mo"="#F28E2B","6mo"="#E15759")) +
  scale_color_manual(values=c("2mo"="#4E79A7","4mo"="#F28E2B","6mo"="#E15759")) +
  annotate("text",x=Inf,y=Inf,hjust=1.05,vjust=1.2,label=paste0("2mo: R = ",round(cor_label$R[1],3),", P = ",signif(cor_label$P[1],3)),color="#4E79A7",size=4.5) +
  annotate("text",x=Inf,y=Inf,hjust=1.05,vjust=2.8,label=paste0("4mo: R = ",round(cor_label$R[2],3),", P = ",signif(cor_label$P[2],3)),color="#F28E2B",size=4.5) +
  annotate("text",x=Inf,y=Inf,hjust=1.05,vjust=4.4,label=paste0("6mo: R = ",round(cor_label$R[3],3),", P = ",signif(cor_label$P[3],3)),color="#E15759",size=4.5) +
  labs(x="Con species abundance",y="Preserved sample abundance",title="Species abundance correlation") +
  theme_classic(base_size=15) +
  theme(plot.title=element_text(face="bold",hjust=0.5),axis.title=element_text(face="bold"),axis.text=element_text(color="black"),legend.title=element_blank(),legend.position=c(0.8,0.2))

ggsave("Figure2f.pdf",p2f,width=5,height=5,units="in",useDingbats=FALSE)

## Figure 2g
sample_paths <- setNames(file.path(matrix_dir,sample_list),sample_list)
groups <- sample_list[-1]

get_logCPM <- function(path){
  counts_mat <- CreateSeuratObject(Read10X(path,gene.column=1),min.cells=1,min.features=1)@assays$RNA$counts
  tibble(gene=rownames(counts_mat),expr=log(rowSums(counts_mat)/sum(counts_mat)*1e6+1))
}

expr_list <- lapply(sample_paths,get_logCPM)
names(expr_list) <- sample_list
con_df <- expr_list$Con
top95_genes <- function(x) names(x)[x>=quantile(x,0.05)]

expr_all <- lapply(groups,function(g){
  merged <- inner_join(con_df,expr_list[[g]],by="gene")
  colnames(merged) <- c("gene","Con","group_expr")
  shared_genes <- intersect(top95_genes(setNames(merged$Con,merged$gene)),top95_genes(setNames(merged$group_expr,merged$gene)))
  merged <- merged %>% filter(gene %in% shared_genes)
  merged$group <- g
  merged
})
names(expr_all) <- groups
expr_all <- bind_rows(expr_all)

gene_cor <- bind_rows(lapply(groups,function(g){
  df <- expr_all %>% filter(group==g)
  cor_res <- cor.test(df$Con,df$group_expr,method="pearson")
  data.frame(group=g,R=cor_res$estimate,P=cor_res$p.value,n=nrow(df))
}))
gene_cor

p2g <- ggplot(expr_all,aes(Con,group_expr)) +
  geom_point(aes(color=group),alpha=0.35,size=1) +
  geom_smooth(aes(color=group),method="lm",se=FALSE,linewidth=0.9) +
  scale_color_manual(values=c("2mo"="#4E79A7","4mo"="#F28E2B","6mo"="#E15759")) +
  annotate("text",x=Inf,y=Inf,hjust=1.05,vjust=1.2,label=paste0("2mo: R = ",round(gene_cor$R[1],3),", P = ",signif(gene_cor$P[1],3)),color="#4E79A7",size=4.3) +
  annotate("text",x=Inf,y=Inf,hjust=1.05,vjust=2.7,label=paste0("4mo: R = ",round(gene_cor$R[2],3),", P = ",signif(gene_cor$P[2],3)),color="#F28E2B",size=4.3) +
  annotate("text",x=Inf,y=Inf,hjust=1.45,vjust=4.2,label=paste0("6mo: R = ",round(gene_cor$R[3],3),", P = ",signif(gene_cor$P[3],3)),color="#E15759",size=4.3) +
  labs(x="Con gene expression log(CPM+1)",y="Preserved sample expression log(CPM+1)",title="Gene expression correlation") +
  theme_classic(base_size=15) +
  theme(plot.title=element_text(face="bold",hjust=0.5,size=16),axis.title=element_text(face="bold",size=14),axis.text=element_text(color="black",size=12),axis.line=element_line(linewidth=0.8),legend.position="none") +
  coord_fixed(ratio=1)

ggsave("Figure2g.pdf",p2g,width=5,height=5,units="in",useDingbats=FALSE)

## Figure 2h
metaG <- fread("Fig2_metaG.txt",sep="\t",header=TRUE)
colnames(metaG) <- c("species","Meta_G","Con","2mo","4mo","6mo")
metaG_filtered <- metaG[Meta_G>0.01]
species_order <- metaG_filtered[order(-Meta_G)]$species
metaG_filtered$species <- factor(metaG_filtered$species,levels=species_order)
genus_colors <- setNames(colorRampPalette(brewer.pal(12,"Set3"))(length(species_order)),species_order)

plot_cor <- function(df,ycol,title){
  cor_res <- cor.test(df$Meta_G,df[[ycol]],method="pearson")
  ggplot(df,aes(Meta_G,.data[[ycol]],fill=species)) +
    geom_point(size=4,shape=21,color="black",stroke=0.3) +
    geom_smooth(method="lm",se=TRUE,aes(group=1),color=NA,fill="grey70",alpha=0.25) +
    scale_fill_manual(values=genus_colors) +
    labs(title=title,x="Meta-G relative abundance",y="mscRNA-seq relative abundance") +
    annotate("text",x=max(df$Meta_G)*0.05,y=max(df[[ycol]])*0.95,hjust=0,label=paste0("R = ",sprintf("%.3f",cor_res$estimate),"\nP = ",format(cor_res$p.value,scientific=TRUE,digits=2)),size=5) +
    theme_pubr() +
    theme(legend.position="none",plot.title=element_text(hjust=0.5,face="bold"),axis.title=element_text(face="bold"))
}

p2h_Con <- plot_cor(metaG_filtered,"Con","Con vs Meta-G")
p2h_2mo <- plot_cor(metaG_filtered,"2mo","2mo vs Meta-G")
p2h_4mo <- plot_cor(metaG_filtered,"4mo","4mo vs Meta-G")
p2h_6mo <- plot_cor(metaG_filtered,"6mo","6mo vs Meta-G")

p2h <- (p2h_Con+p2h_2mo)/(p2h_4mo+p2h_6mo)
ggsave("Figure2h.pdf",p2h,width=8,height=8,units="in",useDingbats=FALSE)