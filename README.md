# Single-Cell-RNA-seq-Analysis-of-Macrophage-Responses-to-Influenza-Infection-in-Mice

## Introduction

Through the usage of a single cell RNA sequencing data collected from Kazer et al, 2025, this study aimed to investigate the differences between macrophage recruitment across peak viral load at day 2 (D02) and myeloid recruitment at day 5 (D05) in mice.  

Influenza A virus is a highly contagious respiratory pathogen impacting almost 10% of the world's population each year (Javanian et al., 2021). The virus typically enters the body through respiratory epithelial cells when inhaled (Bouvier & Palese, 2008). Once the virus has entered the host, an immune response begins involving both innate and adaptive responses. Changes in gene expression across infection timepoints demonstrate shifts in immune activity and the recruitment of various immune cells as infection progresses. In the early phase of infection innate, immune cells including macrophages and monocytes are rapidly recruited to the site of infection as pathogens are detected (Latino & Gonzalez, 2021). Macrophages also work to upregulate antiviral interferon stimulated genes to further combat infection (Li et al., 2023).  

To analyze the dataset used, a computational pipeline was implemented using established bioinformatics tools chosen for their suitability for this project. Seurat was used for data processing, normalization, clustering, and visualization, as it is a widely used toolkit for scRNA-seq analysis (Hao et al., 2024). To correct for batch effects, Harmony was selected for data integration (Korsunsky et al., 2019). Harmony was selected over canonical correlation analysis (CCA) because it preserves biological data while being extremely efficient allowing it to run significantly better on the local hardware used for this analysis, with studies noted it requires 30 to 50 times less memory than MultiCCA (Korsunsky et al., 2019). Cell type annotation was performed using SingleR (Aran et al., 2019) with the MouseRNAseqData reference from the celldex package, which provides mouse cell type labels. 

Differential expression analysis between macrophages at D02 and D05 was performed using a pseudobulk approach with DESeq2 (Love et al., 2014). Functional enrichment was conducted through the use of Gene Set Enrichment Analysis (GSEA) and over-representation analysis (ORA) implemented in clusterProfiler (Yu et al., 2012), enabling biological interpretation of the transcriptional differences between timepoints.
  
## Methods

### Data collection 

The single cell RNA sequencing dataset was obtained from Kazer et al. 2025 (https://doi.org/10.1016/j.immuni.2024.06.005), which profiled nasal immune responses to primary influenza A virus infection across three nasal tissues, Respiratory Mucosa (RM), Olfactory Mucosa (OM), and Lateral Nasal Gland (LNG) at five different timepoints, Naive, 2, 5, 8, and 14 days after infection (dpi). The dataset used was in the form of a pre-processed Seurat object and loaded directly into R using readRDS(). All analyses were performed in R (v. 4.5.3) using Seurat (v. 5.0) (Hao et al., 2024).

### Quality control

Quality control was performed to remove low quality cells prior to analysis. Mitochondrial gene content was calculated for each cell using PercentageFeatureSet() with the pattern "^mt-" , for the mouse genome (Mus musculus). Cells were kept in the data if they expressed more than 200 genes and had less than 15% mitochondrial reads. Gene expression counts were normalized using NormalizeData() with the default log normalization method. Highly variable features were identified using FindVariableFeatures(), data was then scaled using ScaleData(), and principal component analysis was performed using RunPCA(). The number of significant principal components was assessed using an elbow plot generated with ElbowPlot(), and the top 30 PCs were used for downstream analysis.

### Batch correction and integration

To correct for batch effects across samples, data integration was performed using Harmony (Korsunsky et al., 2019) with the IntegrateLayers() function in Seurat with method = HarmonyIntegration, using the PCA data as input. Layers were rejoined after integration using JoinLayers() prior to clustering.

### Clustering and UMAP visualization

Cell clustering was performed using FindNeighbors() and FindClusters() in Seurat, using the harmony corrected data with dims = 1:30 and a resolution of 0.5. Uniform Manifold Approximation and Projection (UMAP) was performed using RunUMAP() for visualization. Clusters were visualized using DimPlot(), colored by cluster, tissue of origin, and timepoint.

### Cell cluster annotation

Automated cell type annotation was performed using SingleR (v. 2.1.0) (Aran et al., 2019) with the MouseRNAseqData reference from the celldex package. Annotation was performed at the cluster level, which assigned a single label per cluster based on the reference dataset. Predicted labels were mapped back to individual cells and stored in the celltype metadata column. Annotations were validated by examining feature plots of canonical cell type marker genes including Cd3d, Cd19, Adgre1, and Omp.

### Differential expression analysis

Differential expression analysis was performed to identify genes significantly differentially expressed in macrophages between day 2 at peak viral load and day 5 at the myeloid recruitment phase following influenza infection. Macrophage cells were first subset from the annotated dataset and filtered to include only D02 and D05 timepoints. A pseudobulk approach was used by aggregating gene expression counts by mouse ID and timepoint. Differential expression analysis was then conducted using DESeq2 via the FindMarkers() function in Seurat. Expression patterns of selected genes of interest were visualized using violin plots.

### Functional enrichment analysis

Functional enrichment was performed using the clusterProfiler (v. 4.16.0) package (Yu et al., 2012). The org.Mm.eg.db mouse annotation database was used for analysis. Gene Set Enrichment Analysis was performed using gseGO() on the ranked gene list for DESeq2 compared with Gene Ontology Biological Process terms and a p-value cutoff of 0.05. Results were visualized as a ridge plot using ridgeplot(). Over-representation analysis was performed on significantly upregulated genes at each timepoint using enrichGO(). ORA results were visualized using a dot plot.

## Results

<div align="center">

![fig1](https://github.com/Kalixie/Single-Cell-RNA-seq-Analysis-of-Macrophage-Responses-to-Influenza-Infection-in-Mice/blob/main/figures/qc_after_filter.png)

</div>

The majority of cells in the Seurat object displayed between 1,000 and 3,500 detected genes (nFeature_RNA), total molecules between 5,000 and 8,000 (nCount_RNA), and most mitochondrial transcripts appeared below 10%. The distributions were consistent across all timepoints as the data was likely pre-processed by the original authors. After filtering thresholds were applied the distributions were mostly unchanged suggesting that the data is of high quality.

After normalization, integration, and clustering, a UMAP plot revealed 36 distinct cell populations. Cells formed well separated clusters with minimal overlap between groups. Coloring cells by tissue of origin showed that all three tissues including Respiratory Mucosa (RM), Olfactory Mucosa (OM), and Lateral Nasal Gland (LNG) contributed cells across most clusters, with some tissue specific clustering visible such as in the later annotated neuronal cluster which was predominantly OM derived. Coloring by timepoint showed that all five timepoints were distributed across all clusters without strong separation by timepoint, indicating that Harmony integration properly corrected for batch effects and that cells clustered by cell type rather than sample.

Automated cell type annotation using SingleR with the MouseRNAseqData reference assigned labels to 13 identified clusters. Cell populations included Macrophages, Monocytes, T cells, NK cells, B cells, Granulocytes, Neurons, Fibroblasts, Epithelial cells, Endothelial cells, Erythrocytes, and Adipocytes. Separation was clear between most of the groups, although Fibroblasts had a spread out structure. Feature plots of marker genes specific to cell clusters confirmed the annotation assignments with Cd3d expression in the T cell cluster, Cd19 in the B cell cluster, Adgre1 in the Macrophage cluster, and Omp showed strong specific expression in the Neuron cluster consistent with the presence of olfactory sensory neurons.


## Discussion

Clustering revealed well separated cell populations, suggesting that transcriptional differences between cell types were clearly distinct between cells in the data. Although 36 clusters were identified pre-annotation, only 12 clusters were annotated automatically by SingleR, potentially influenced by the resolution. Despite this, annotated clusters were well separated and closely aligned with the results reference dataset paper by Kazer et al, with distinct groupings of Macrophages, Monocytes, T cells, NK cells, B cells, Granulocytes, Neurons, Fibroblasts, Epithelial cells, and Endothelial cells with variation likely due to filtering and annotation dataset differences between pipelines (Kazer et al., 2024). The presence of all five timepoints across clusters in the timepoint UMAP confirmed that Harmony integration properly removed sample level variation, and that clustering reflected biological cell identities rather than batch effects. Batch effects can lead to improper conclusions that are not based on biological data by clustering incorrectly, therefore correctness is required for proper downstream analysis (Yu et al., 2024). 

Validation of cell type assignments using canonical marker genes further supported the accuracy of the annotation. The annotation of neuronal cells within the olfactory mucosa is consistent with the presence of olfactory sensory neurons, further confirmed through the feature plot of Olfactory Marker Protein (Omp) identified in the Neuron cluster, which is present in mature olfactory sensory neurons (Gong, 2012). The expression patterns of Cd3d, Cd19, and Adgre1 correspond to T cells, B cells, and macrophages, respectively, further supporting the results of the annotation. 

Pseudobulk differential expression analysis between macrophages at D02 and D05 identified 5 genes significantly upregulated at D02 and 18 genes significantly upregulated at D05. The greater number of genes upregulated at D05 compared to D02 suggests a higher transcriptional response at the peak myeloid recruitment phase, results consistent with the more active antiviral response expected at this timepoint. Myeloid cells, such as monocytes and macrophages, produce inflammatory cytokines such as IL-12 which initiate an immune response cascade and lead to the recruitment of adaptive T-cells (Stegalmeier et al., 2019).


## References
