# AGM_code

For questions regarding the dataset or analysis, please contact:

- **Yongcheng Wang** - [yongcheng@zju.edu.cn](mailto:yongcheng@zju.edu.cn)
- **Xinghai Zheng** - [xhzheng@zju.edu.cn](mailto:xhzheng@zju.edu.cn)

## Data and Script Preparation

The **Space_Scripts** folder contains the analysis scripts for the **spacecraft astronaut gut microbiome single-cell RNA sequencing (mscRNA-seq) data**.

The corresponding **expression matrices**, **species annotation files**, and **gene annotation file** for the dataset are available at the following link:

[https://doi.org/10.6084/m9.figshare.33395212](https://doi.org/10.6084/m9.figshare.33395212)

Please download and extract the dataset to obtain the **Space_Matrix** and **Store_Matrix** folders, as well as the file **reference_genome_gene_annotation.tsv**.

For subsequent analyses, ensure that the **Space_Matrix** folder, **Store_Matrix** folder, **reference_genome_gene_annotation.tsv** file, and the **Space_Scripts** folder are all placed in the same directory path, as shown below:

```
project_directory/
├── Space_Scripts
│   ├── Figure2
│   ├── Figure3
│   ├── Figure4
│   ├── Figure5
│   └── Figure6
├── Space_Matrix
├── Store_Matrix
└── reference_genome_gene_annotation.tsv
```

Once the folders are organized as described, you can directly execute the analysis scripts.

---

## Software and Version Information

The major software and packages used in the analysis include:

| Software / Package | Version | Main application |
|---|---|---|
| **Seurat** | **v5.3.0** | Single-cell expression matrix processing, normalization, dimensionality reduction, clustering, and downstream single-cell analysis |
| **Harmony** | **v1.2.3** | Batch-effect correction and integration of single-cell datasets |
| **monocle3** | **v1.4.26** | Trajectory inference and pseudotime analysis |
| **hdWGCNA** | **v0.4.6** | Single-cell co-expression network construction and module analysis |
| **WGCNA** | **v1.74** | Weighted gene co-expression network analysis |
| **BPCells** | **v0.3.1** | Efficient processing and storage of large single-cell expression matrices |
| **SeuratWrappers** | **v0.4.0** | Integration of additional single-cell analysis tools with Seurat |
| **edgeR** | **v4.0.16** | Differential expression and count-based statistical analysis |
| **limma** | **v3.58.1** | Differential expression and statistical modeling |
| **glmGamPoi** | **v1.14.3** | Gamma-Poisson generalized linear modeling for single-cell count data |
| **Matrix** | **v1.6.5** | Sparse matrix representation and manipulation of single-cell expression data |
| **irlba** | **v2.3.5.1** | Fast truncated singular value decomposition for dimensionality reduction |
| **dplyr** | **v1.1.4** | Data manipulation and transformation |
| **tibble** | **v3.3.0** | Data-frame and metadata handling |
| **data.table** | **v1.17.8** | Efficient processing of large tabular datasets |
| **tidyr** | **v1.3.1** | Data reshaping and tidying |
| **tidyverse** | **v2.0.0** | Integrated data manipulation and analysis workflows |
| **ggplot2** | **v4.0.0** | Data visualization and figure generation |
| **ggpubr** | **v0.6.2** | Statistical visualization and publication-quality figure generation |
| **ggrepel** | **v0.9.6** | Non-overlapping text labels in plots |
| **ggalluvial** | **v0.12.6** | Alluvial and flow diagram visualization |
| **patchwork** | **v1.3.2** | Combination and arrangement of multiple plots |
| **cowplot** | **v1.2.0** | Plot arrangement and figure composition |
| **RColorBrewer** | **v1.1.3** | Color palettes for data visualization |
| **ggsci** | **v3.2.0** | Scientific journal-style color palettes |
| **viridis** | **v0.6.5** | Color palettes for quantitative data visualization |

Custom R scripts and analysis pipelines used for data preprocessing, quality control, clustering, differential expression analysis, trajectory inference, co-expression network analysis, and figure generation are provided in the corresponding code repository.

---

## Bioinformatics Analysis Workflow

The scripts used to generate each figure are organized into separate folders (`Figure 2`–`Figure 6`) within the **Space_Scripts** directory.  
The execution order of the scripts for reproducing each figure is described below.

---

## Figure 2

Navigate to the folder:

```
Space_Scripts/Figure2
```

To generate **Figure 2e**, **Figure 2f**, **Figure 2g**, and **Figure 2h**, run:

```
Figure2.R
```

---

## Figure 3

Navigate to the folder:

```
Space_Scripts/Figure3
```

To generate **Figure 3b** and **Figure 3c**, run the following scripts in the following order:

```
Figure3_qc_primary.py
Figure3_qc_secondary.py
Figure3bc_similarity.py
```

To generate **Figure 3d**, run:

```
Figure3d_abundance.R
```

To generate **Figure 3e**, run the following scripts in sequence:

```
Figure3ef_stage_scaled.py
Figure3e_volcano.R
```

To generate **Figure 3f**, run:

```
Figure3f_stage_heatmap.py
```

---

## Figure 4

Navigate to the folder:

```
Space_Scripts/Figure4
```

To generate **Figure 4a**, run:

```
Figure4a.R
```

To generate **Figure 4b**, run the following scripts in sequence:

```
Figure4_umap_clustering.R
Figure4b_composition_umap.py
```

To generate **Figure 4c**, run:

```
Figure4c_cluster_feature.py
```

To generate **Figure 4d**, run:

```
Figure4d_genus.py
```

To generate **Figure 4e**, run:

```
Figure4e_astronaut.py
```

To generate **Figure 4f**, run:

```
Figure4f_stage.py
```

---

## Figure 5

Navigate to the folder:

```
Space_Scripts/Figure5
```

To generate **Figure 5a**, run:

```
Figure5_pseudotime.R
```

To generate **Figure 5b**, run:

```
Figure5b.py
```

To generate **Figure 5c**, run:

```
Figure5c.py
```

To generate **Figure 5d** and **Figure 5e**, run:

```
Figure5de.py
```

To generate **Figure 5f**, run the following scripts in sequence:

```
Figure5f_annotation.py
Figure5f_cog_enrichment.py
Figure5f_cog_heatmap.py
```

To generate **Figure 5g**, run the following scripts in sequence:

```
Figure5g_annotation.py
Figure5g_cog_enrichment.py
Figure5g_cog_heatmap.py
```

Here, only the code for *Phocaeicola vulgatus* is provided. The code for other species can be generated following the same workflow as that for *Phocaeicola vulgatus*.

---

## Figure 6

Navigate to the folder:

```
Space_Scripts/Figure6
```

To generate **Figure 6a**, **Figure 6c**, **Figure 6d**, and **Figure 6g**, run:

```
Figure6_umap_hdWGCNA.R
```

To generate **Figure 6b**, run:

```
Figure6b_timepoint.py
```

To generate **Figure 6e**, run the following scripts in sequence:

```
Figure6e_hub_annotation.py
Figure6e_hub_cog_enrichment.py
Figure6e_hub_cog_heatmap.py
```

To generate **Figure 6f**, run:

```
Figure6f_hub_timepoint_bubble.py
```
