# AGM_code

For questions regarding the dataset or analysis, please contact:

- **Yongcheng Wang** - [yongcheng@zju.edu.cn](mailto:yongcheng@zju.edu.cn)
- **Xinghai Zheng** - [xhzheng@zju.edu.cn](mailto:xhzheng@zju.edu.cn)

## Data and Script Preparation

The **Space_Scripts** folder contains the analysis scripts for the **spacecraft astronaut gut microbiome single-cell RNA sequencing (mscRNA-seq) data**.

The corresponding **expression matrices**, **species annotation files**, and **gene annotation file** for the dataset are available at the following link:

[https://ngdc.cncb.ac.cn/omix/preview/hSk6yG5k](https://ngdc.cncb.ac.cn/omix/preview/hSk6yG5k)

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

## Bioinformatics Analysis Workflow

The scripts used to generate each figure are organized into separate folders (`Figure 2`–`Figure 6`) within the **Space_Scripts** directory.  
The execution order of the scripts for reproducing each figure is described below.

---

## Figure 2

Navigate to the folder:

```
Space_Scripts/Figure2
```

To generate **Figure 2e** and **Figure 2f**, run the scripts in the following order:

```
Figure2ef_cell_qc_primary.py
Figure3ef_species_similarity.py
```

To generate **Figure 2i**, run:

```
Figure2gij_umap_clustering.R
```

To generate **Figure 2g** and **Figure 2j**, run the following scripts sequentially:

```
Figure2g_gene_counts.py
Figure2j_sankey_mapping.py
```

To generate **Figure 2h**, run the following scripts in sequence:

```
Figure2h_expression_profile.R
Figure2h_expression_similarity.py
```

---

## Figure 3

Navigate to the folder:

```
Space_Scripts/Figure3
```

To generate **Figure 3b**, run:

```
Figure3b_msc_vs_meta.R
```

To generate **Figure 3c** and **Figure 3d**, run the scripts in the following order:

```
Figure3_cell_qc_primary.py
Figure3_cell_qc_secondary.py
Figure3cd_similarity_mds_jsd.py
```

To generate **Figure 3e**, run:

```
Figure3e_abundance_profile.R
```

To generate **Figure 3f**, run the following scripts in sequence:

```
Figure3fg_stage_scaled.py
Figure3f_stage_volcano.R
```

To generate **Figure 3g**, run:

```
Figure3g_stage_heatmap.py
```

---

## Figure 4

Navigate to the folder:

```
Space_Scripts/Figure4
```

To generate **Figure 4a** and **Figure 4b**, run the following scripts in sequence:

```
Figure4_umap_clustering.R
Figure4ab_cell_gene_counts.py
```

To generate **Figure 4c**, run:

```
Figure4c_composition_umap.py
```

To generate **Figure 4d**, run:

```
Figure4d_species_feature.py
```

To generate **Figure 4e**, run the scripts in the following order:

```
Figure4e_cluster_annotation.py
Figure4e_cluster_similarity.py
Figure4e_cluster_feature.py
```

---

## Figure 5

Navigate to the folder:

```
Space_Scripts/Figure5
```

To generate **Figure 5a** and **Figure 5c**, run the following script:

```
Figure5_pseudotime.R
```

To generate **Figure 5b**, run the scripts in sequence:

```
Figure5b_cluster_timepoint.py
Figure5b_cluster_timepoint.R
```

To generate **Figure 5d**, run:

```
Figure5d_timepoint_pseudotime.py
```

To generate **Figure 5e**, run:

```
Figure5e_pseudotime_clustering.py
```

To generate **Figure 5f**, run the scripts in the following order:

```
Figure5f_cluster_annotation.py
Figure5f_cluster_cog_enrichment.py
Figure5f_cluster_cog_heatmap.py
```

---

## Figure 6

Navigate to the folder:

```
Space_Scripts/Figure6
```

To generate **Figure 6a**, **Figure 6b**, and **Figure 6d**, run:

```
Figure6_hdWGCNA.R
```

To generate **Figure 6c**, run the following scripts in sequence:

```
Figure6cf_hub_annotation.py
Figure6c_hub_timepoint_bubble.py
```

To generate **Figure 6f**, run the following scripts in sequence:

```
Figure6f_hub_cog_enrichment.py
Figure6f_hub_cog_heatmap.py
```

To generate **Figure 6e**, run:

```
Figure6e_edge_annotation.py
```

This will generate the base files for importing into Cytoscape to create **Figure 6e**.
