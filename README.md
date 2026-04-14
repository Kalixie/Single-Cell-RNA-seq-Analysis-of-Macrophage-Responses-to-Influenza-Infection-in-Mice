# Single-Cell-RNA-seq-Analysis-of-Macrophage-Responses-to-Influenza-Infection-in-Mice

## Introduction

Through the usage of single-cell RNA sequencing data collected from Kazer et al, 2024, this study aimed to investigate the differences between macrophage recruitment across peak viral load at day 2 (D02) and myeloid recruitment at day 5 (D05) in mice.  

Influenza A virus is a highly contagious respiratory pathogen impacting almost 10% of the world's population each year (Javanian et al., 2021). The virus typically enters the body through respiratory epithelial cells when inhaled (Bouvier & Palese, 2008). Once the virus has entered the host, an immune response begins involving both innate and adaptive responses. Changes in gene expression across infection time reflect shifts in immune activity and the recruitment of various immune cells as infection progresses. In the early phase of infection, innate immune cells, including macrophages and monocytes, are rapidly recruited to the site of infection as pathogens are detected (Latino & Gonzalez, 2021). Macrophages also work to upregulate antiviral interferon-stimulated genes to further combat infection (Li et al., 2023).  

To analyze the dataset used, a computational pipeline was implemented using established bioinformatics tools chosen for their suitability for this project. Seurat was used for data processing, normalization, clustering, and visualization, as it is a widely used toolkit for scRNA-seq analysis (Hao et al., 2024). To correct for batch effects, Harmony was selected for data integration (Korsunsky et al., 2019). Harmony was selected over canonical correlation analysis (CCA) because it preserves biological data while being extremely efficient, allowing it to run significantly better on the local hardware used for this analysis, with studies noting it requires 30 to 50 times less memory than MultiCCA (Korsunsky et al., 2019). Cell type annotation was performed using SingleR, which enabled automated assignment of cell identities by comparing each cell’s expression profile to specified reference datasets (Aran et al., 2019). The MouseRNAseqData reference from the celldex package was selected as it provides an annotation of mouse cell types for both immune and non-immune cells that may be found in samples.

Differential expression analysis between macrophages at D02 and D05 was performed using a pseudobulk approach with DESeq2 (Love et al., 2014). DESeq2 was selected over alternative methods such as the default Seurat differential expression test, as it models gene expression using a negative binomial distribution and accounts for variability between replicates, reducing false positives (Squair et al., 2021). Functional enrichment was conducted through the use of Gene Set Enrichment Analysis (GSEA) and over-representation analysis (ORA) implemented in clusterProfiler (Yu et al., 2012), enabling biological interpretation of the transcriptional differences between timepoints.
  
## Methods

### Data collection 

The single cell RNA sequencing dataset was obtained from Kazer et al. 2024 (https://doi.org/10.1016/j.immuni.2024.06.005), which profiled nasal immune responses to primary influenza A virus infection across three nasal tissues, Respiratory Mucosa (RM), Olfactory Mucosa (OM), and Lateral Nasal Gland (LNG) at five different timepoints, Naive, 2, 5, 8, and 14 days after infection (dpi). The dataset used was in the form of a pre-processed Seurat object and loaded directly into R using readRDS(). All analyses were performed in R (v. 4.5.3) using Seurat (v. 5.0) (Hao et al., 2024).

### Quality control

Quality control was performed to remove low-quality cells prior to analysis. Mitochondrial gene content was calculated for each cell using PercentageFeatureSet() with the pattern "^mt-" , for the mouse genome (Mus musculus). Cells were kept in the data if they expressed more than 200 genes and had less than 15% mitochondrial reads. Gene expression counts were normalized using NormalizeData() with the default log normalization method. Highly variable features were identified using FindVariableFeatures(), data was then scaled using ScaleData(), and principal component analysis was performed using RunPCA(). The number of significant principal components was assessed using an elbow plot generated with ElbowPlot(), and the top 30 PCs were used for downstream analysis.

### Batch correction and integration

To correct for batch effects across samples, data integration was performed using Harmony (v. 1.2.4) (Korsunsky et al., 2019) with the IntegrateLayers() function in Seurat with method = HarmonyIntegration, using the PCA data as input. Layers were rejoined after integration using JoinLayers() prior to clustering.

### Clustering and UMAP visualization

Cell clustering was performed using FindNeighbors() and FindClusters() in Seurat, using the harmony corrected data with dims = 1:30 and a resolution of 0.5. Uniform Manifold Approximation and Projection (UMAP) was performed using RunUMAP() for visualization. Clusters were visualized using DimPlot(), colored by cluster, tissue of origin, and timepoint.

### Cell cluster annotation

Automated cell type annotation was performed using SingleR (v. 2.1.0) (Aran et al., 2019) with the MouseRNAseqData reference from the celldex package (v. 1.18.0) (Aran et al., 2019). Annotation was performed at the cluster level, which assigned a single label per cluster based on the reference dataset. Predicted labels were mapped back to individual cells and stored in the celltype metadata column. Annotations were validated by examining feature plots of canonical cell type marker genes, including Cd3d, Cd19, Adgre1, and Omp.

### Differential expression analysis

Differential expression analysis was performed to identify genes significantly differentially expressed in macrophages between day 2 at peak viral load and day 5 at the myeloid recruitment phase following influenza infection. Macrophage cells were first subset from the annotated dataset and filtered to include only D02 and D05 timepoints. A pseudobulk approach was used by aggregating gene expression counts by mouse ID and time point. Differential expression analysis was then conducted using DESeq2 via the FindMarkers() function in Seurat. Expression patterns of selected genes of interest were visualized using violin plots.

### Functional enrichment analysis

Functional enrichment was performed using the clusterProfiler (v. 4.16.0) package (Yu et al., 2012). The org.Mm.eg.db mouse annotation database was used for analysis. Gene Set Enrichment Analysis was performed using gseGO() on the ranked gene list for DESeq2 (Love et al., 2014) compared with Gene Ontology Biological Process terms and a p-value cutoff of 0.05. Results were visualized as a ridge plot using ridgeplot(). Over-representation analysis was performed on significantly upregulated genes at each time point using enrichGO(). ORA results were visualized using a dot plot.

## Results

<div align="center">

![fig1](https://github.com/Kalixie/Single-Cell-RNA-seq-Analysis-of-Macrophage-Responses-to-Influenza-Infection-in-Mice/blob/main/figures/qc_after_filter.png)

</div>

##### Figure 1: Quality control violin plot after filtering. The violin plot shows the distribution of the number of detected genes (nFeature_RNA), total UMI counts (nCount_RNA), and mitochondrial transcript percentage (percent.mt) across all five timepoints: naive, D02, D05, D08, D14, after applying nFeature_RNA > 200 and percent.mt < 15% (bottom).

The majority of cells in the Seurat object displayed between 1,000 and 3,500 detected genes (nFeature_RNA), total molecules between 5,000 and 8,000 (nCount_RNA), and most mitochondrial transcripts appeared below 10%. The distributions were consistent across all timepoints, as the data was likely pre-processed by the original authors. After filtering, thresholds were applied, and the distributions were mostly unchanged, suggesting that the data is of high quality (Figure 1).

<div align="center">

![fig2](https://github.com/Kalixie/Single-Cell-RNA-seq-Analysis-of-Macrophage-Responses-to-Influenza-Infection-in-Mice/blob/main/figures/umap_clusters.png)

</div>

##### Figure 2: UMAP plot of all cells colored by cluster number. Each point represents a single cell following Uniform Manifold Approximation and Projection and Harmony integration. Cluster numbers are displayed for each cluster.

<div align="center">

![fig3](https://github.com/Kalixie/Single-Cell-RNA-seq-Analysis-of-Macrophage-Responses-to-Influenza-Infection-in-Mice/blob/main/figures/umap_tissue.png)

</div>

##### Figure 3: Figure 3: UMAP plot colored by tissue of origin. Cells are colored by the tissue from which they were collected, including Respiratory Mucosa (RM) in blue, Olfactory Mucosa (OM) in green, and Lateral Nasal Gland (LNG) in red. Overlap can be seen by the multicoloured sections throughout sections of the plot. 

<div align="center">

![fig4](https://github.com/Kalixie/Single-Cell-RNA-seq-Analysis-of-Macrophage-Responses-to-Influenza-Infection-in-Mice/blob/main/figures/umap_time.png)

</div>

##### Figure 4: UMAP plot colored by timepoint post-infection. Cells are colored by the timepoint at which they were collected, including Naive in magenta, D02 in red, D05 in yellow, D08 in green, and D14 in blue. The mixing of timepoints within clusters indicates successful batch correction by Harmony integration.

After normalization, integration, and clustering, a UMAP plot revealed 36 distinct cell populations. Cells formed well-separated clusters with minimal overlap between groups (Figure 2). Coloring cells by tissue of origin showed that all three tissues, including Respiratory Mucosa (RM), Olfactory Mucosa (OM), and Lateral Nasal Gland (LNG), contributed cells across most clusters, with some tissue-specific clustering visible, such as in the later annotated neuronal cluster, which was predominantly OM-derived (Figure 3). Coloring by timepoint showed that all five timepoints were distributed across all clusters without strong separation by timepoint, indicating that Harmony integration properly corrected for batch effects and that cells clustered by cell type rather than sample (Figure 4).

<div align="center">

![fig5](https://github.com/Kalixie/Single-Cell-RNA-seq-Analysis-of-Macrophage-Responses-to-Influenza-Infection-in-Mice/blob/main/figures/umap_celltype.png)

</div>

##### Figure 5: UMAP plot colored by automated cell type annotation. Cell type labels were assigned using SingleR with the MouseRNAseqData reference and annotated at the cluster level. Twelve distinct cell populations were identified and can be distinguished by colour and label.

<div align="center">

![fig6](https://github.com/Kalixie/Single-Cell-RNA-seq-Analysis-of-Macrophage-Responses-to-Influenza-Infection-in-Mice/blob/main/figures/featureplot_markers.png)

</div>

##### Figure 6: Feature plots of cell type marker genes. Expression of Cd3d (T cells), Cd19 (B cells), Adgre1 (macrophages), and Omp (olfactory sensory neurons) projected onto the UMAP embedding. Color intensity in purple reflects the expression level. Grey indicates no detected expression within the cluster.

Automated cell type annotation using SingleR with the MouseRNAseqData reference assigned labels to 12 identified clusters. Cell populations included Macrophages, Monocytes, T cells, NK cells, B cells, Granulocytes, Neurons, Fibroblasts, Epithelial cells, Endothelial cells, Erythrocytes, and Adipocytes (Figure 5). Separation was clear between most of the groups, although Fibroblasts had a spread-out structure. Feature plots of marker genes specific to cell clusters confirmed the annotation assignments with Cd3d expression in the T cell cluster, Cd19 in the B cell cluster, Adgre1 in the Macrophage cluster, and Omp showed strong specific expression in the Neuron cluster, consistent with the presence of olfactory sensory neurons (Figure 6).

<div align="center">

![fig7](https://github.com/Kalixie/Single-Cell-RNA-seq-Analysis-of-Macrophage-Responses-to-Influenza-Infection-in-Mice/blob/main/figures/gsea_ridgeplot.png)

</div>

##### Figure 7: Gene Set Enrichment Analysis (GSEA) ridge plot for macrophages comparing D02 vs D05. Gene sets from Gene Ontology Biological Process are shown ranked by enrichment score. Distributions shifted toward positive values indicate enrichment at D02, while distributions shifted toward negative values indicate enrichment at D05. Color indicates adjusted p-value.

<div align="center">

![fig8](https://github.com/Kalixie/Single-Cell-RNA-seq-Analysis-of-Macrophage-Responses-to-Influenza-Infection-in-Mice/blob/main/figures/ora_dotplot_D02.png)

</div>

Pseudobulk differential expression analysis using DESeq2 identified genes significantly changed in Macrophages between D02 and D05. Gene Set Enrichment Analysis of all ranked differentially expressed genes in macrophages between D02 and D05 revealed enrichment of multiple biological processes at D05 (Figure 7). All pathways displayed negative normalized enrichment scores and a large upregulation of genes at D05 relative to D02. 

##### Figure 8: Over-Representation Analysis (ORA) dot plot of genes significantly upregulated in macrophages at D02 relative to D05. Gene Ontology Biological Process terms are shown on the y-axis. Dot size shows the gene ratio, and color depicts the adjusted p-value.

<div align="center">

![fig9](https://github.com/Kalixie/Single-Cell-RNA-seq-Analysis-of-Macrophage-Responses-to-Influenza-Infection-in-Mice/blob/main/figures/ora_dotplot_D05.png)

</div>

##### Figure 9. Over-Representation Analysis (ORA) dot plot of genes significantly upregulated in macrophages at D05 relative to D02. Gene Ontology Biological Process terms are shown on the y-axis. Dot size shows the gene ratio and color depicts the adjusted p-value.

Overrepresentation analysis of significantly upregulated genes at D02 revealed enrichment of pathways related to actin filament organization, regulation of bone resorption, bone remodeling, regulation of tissue remodeling, and cellular pigmentation (Figure 8). ORA of significantly upregulated genes at D05 showed enrichment of response to virus as the most significant pathway, along with macrophage activation involved in immune response, and several nucleotide metabolic processes (Figure 9).

## Discussion

Clustering revealed well-separated cell populations, suggesting that transcriptional differences between cell types were clearly distinct between cells in the data. Although 36 clusters were identified pre-annotation, only 12 clusters were annotated automatically by SingleR, potentially influenced by the resolution. Despite this, annotated clusters were well separated and closely aligned with the results reference dataset paper by Kazer et al, with distinct groupings of Macrophages, Monocytes, T cells, NK cells, B cells, Granulocytes, Neurons, Fibroblasts, Epithelial cells, and Endothelial cells with variation likely due to filtering and annotation dataset differences between pipelines (Kazer et al., 2024). The presence of all five timepoints across clusters in the timepoint UMAP confirmed that Harmony integration properly removed sample-level variation and that clustering reflected biological cell identities rather than batch effects. Batch effects can lead to improper conclusions that are not based on biological data by clustering incorrectly, and high accuracy is required for proper downstream analysis (Yu et al., 2024). 

Validation of cell type assignments using canonical marker genes further supported the accuracy of the annotation. The annotation of neuronal cells within the olfactory mucosa is consistent with the presence of olfactory sensory neurons, further confirmed through the feature plot of Olfactory Marker Protein (Omp) identified in the Neuron cluster, which is present in mature olfactory sensory neurons (Gong, 2012). The expression patterns of Cd3d, Cd19, and Adgre1 correspond to T cells, B cells, and macrophages, respectively, further supporting the results of the annotation. 

Pseudobulk differential expression analysis between macrophages at D02 and D05 identified 5 genes significantly upregulated at D02 and 18 genes significantly upregulated at D05. The greater number of genes upregulated at D05 compared to D02 suggests a higher transcriptional response at the peak myeloid recruitment phase, results consistent with the increased active antiviral response expected at this time point. Myeloid cells, such as monocytes and macrophages, produce inflammatory cytokines such as IL-12, which initiate an immune response cascade and lead to the recruitment of adaptive T-cells (Stegalmeier et al., 2019).

GSEA revealed high levels of enrichment between Macrophage pathways when contrasting D02 with D05. When comparing both times, the resulting GSEA plot displayed negative log2FC changes throughout the plot, which suggests downregulation in D02 in comparison to D05, therefore the pathways were upregulated in D05. Antiviral and interferon-related pathways were notable at D05, including response to virus, defense response to virus, regulation of viral process, response to type I interferon, and microtubule transport-related pathways. These findings support previous studies noting that macrophages produce and respond to type 1 interferons as part of the innate immune response, leading to increased interferon-stimulated genes and the further recruitment of immune cells (Subramani et al., 2023).

Genes upregulated at D05 for over-representation analysis were enriched for response to virus, macrophage activation involved in immune response, and various metabolic processes, confirming that macrophages at peak myeloid recruitment are upregulated transcriptionally for antiviral defense. The enrichment of nucleotide metabolic processes, including pyrimidine deoxyribonucleotide and deoxyribonucleotide metabolic processes at D05 may also reflect increased cellular proliferation associated with macrophage activation (O'Neill et al., 2016). Genes upregulated at D02 were enriched for actin filament organization, bone resorption, and tissue and bone-related pathways. Actin transport systems are used in the trafficking of viruses, which potentially supports what is seen at D02 as infection reaches peak viral load (Jiang et al., 2025). Previous studies have noted that viral infections can impact bone homeostasis, potentially supporting the upregulation of bone resorption pathways (Caetano et al., 2024).

There are several limitations in this analysis that can be improved in future work. Cell type annotation was performed using the MouseRNAseqData general reference rather than the tissue label transfer approach used by Kazer et al. (2024), which may have resulted in the slightly varying annotation results seen. The analysis was restricted to macrophages at two timepoints, and future analyses could extend to other cell types and timepoints. Overall, the transcriptional changes observed in macrophages between D02 and D05 are generally biologically consistent with the known progression of innate immune responses during influenza infection.

## References

Aran, D., Looney, A. P., Liu, L., Wu, E., Fong, V., Hsu, A., Chak, S., Naikawadi, R. P., Wolters, P. J., Abate, A. R., Butte, A. J., & Bhattacharya, M. (2019). Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. Nature Immunology, 20(2), 163–172. https://doi.org/10.1038/s41590-018-0276-y

Bouvier, N. M., & Palese, P. (2008). The biology of influenza viruses. Vaccine, 26(Suppl 4), D49–D53. https://doi.org/10.1016/j.vaccine.2008.07.039

Caetano, C. C. S., Azamor, T., Meyer, N. M., Onwubueke, C., Calabrese, C. M., Calabrese, L. H., Visperas, A., Piuzzi, N. S., Husni, M. E., Foo, S. S., & Chen, W. (2024). Mechanistic insights into bone remodelling dysregulation by human viral pathogens. Nature Microbiology, 9(2), 322–335. https://doi.org/10.1038/s41564-023-01586-6

Gong, Q. (2012). Culture of mouse olfactory sensory neurons. Current Protocols in Neuroscience, 58(1), Unit 3.24. https://doi.org/10.1002/0471142301.ns0324s58

Hao, Y., Stuart, T., Kowalski, M. H., Choudhary, S., Hoffman, P., Hartman, A., Srivastava, A., Molla, G., Madad, S., Fernandez-Granda, C., & Satija, R. (2024). Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nature Biotechnology, 42(2), 293–304. https://doi.org/10.1038/s41587-023-01767-y

Javanian, M., Barary, M., Ghebrehewet, S., Koppolu, V., Vasigala, V., & Ebrahimpour, S. (2021). A brief review of influenza virus infection. Journal of Medical Virology, 93(8), 4638–4646. https://doi.org/10.1002/jmv.26990

Jiang, M., Zou, J., Jin, Y., Jiang, C., Tu, S., Chen, T., Guo, J., Cheng, Y., Jin, M., Chen, H., & Zhou, H. (2025). Adducin-1 facilitates influenza virus endosomal trafficking and uncoating by regulating branched actin dynamics and myosin IIB activity. Advanced Science, 12(28), e2417318. https://doi.org/10.1002/advs.202417318

Kazer, S. W., Match, C. M., Langan, E. M., Messou, M. A., LaSalle, T. J., O'Leary, E., Marbourg, J., Naughton, K., von Andrian, U. H., & Ordovas-Montanes, J. (2024). Primary nasal influenza infection rewires tissue-scale memory response dynamics. Immunity, 57(8), 1955–1974.e8. https://doi.org/10.1016/j.immuni.2024.06.005

Korsunsky, I., Millard, N., Fan, J., Slowikowski, K., Zhang, F., Wei, K., Baglaenko, Y., Brenner, M., Loh, P. R., & Raychaudhuri, S. (2019). Fast, sensitive and accurate integration of single-cell data with Harmony. Nature Methods, 16(12), 1289–1296. https://doi.org/10.1038/s41592-019-0619-0

Latino, I., & Gonzalez, S. F. (2021). Spatio-temporal profile of innate inflammatory cells and mediators during influenza virus infection. Current Opinion in Physiology, 19, 175–186. https://doi.org/10.1016/j.cophys.2020.10.008

Li, H., Wang, A., Zhang, Y., & Wei, F. (2023). Diverse roles of lung macrophages in the immune response to influenza A virus. Frontiers in Microbiology, 14, 1260543. https://doi.org/10.3389/fmicb.2023.1260543

Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15(12), 550. https://doi.org/10.1186/s13059-014-0550-8

O'Neill, L. A., & Pearce, E. J. (2016). Immunometabolism governs dendritic cell and macrophage function. Journal of Experimental Medicine, 213(1), 15–23. https://doi.org/10.1084/jem.20151570

Squair JW, Gautier M, Kathe C, Anderson MA, James ND, Hutson TH, Hudelle R, Qaiser T, Matson KJE, Barraud Q, Levine AJ, La Manno G, Skinnider MA, Courtine G. Confronting false discoveries in single-cell differential expression. Nat Commun. 2021 Sep 28;12(1):5692. doi: 10.1038/s41467-021-25960-2. PMID: 34584091; PMCID: PMC8479118.

Stegelmeier, A. A., van Vloten, J. P., Mould, R. C., Klafuric, E. M., Minott, J. A., Wootton, S. K., Bridle, B. W., & Karimi, K. (2019). Myeloid cells during viral infections and inflammation. Viruses, 11(2), 168. https://doi.org/10.3390/v11020168

Subramani, A., Hite, M. E. L., Garcia, S., Maxwell, J., Kondee, H., Millican, G. E., McClelland, E. E., Seipelt-Thiemann, R. L., & Nelson, D. E. (2023). Regulation of macrophage IFNγ-stimulated gene expression by the transcriptional coregulator CITED1. Journal of Cell Science, 136(1), jcs260529. https://doi.org/10.1242/jcs.260529

Yu, G., Wang, L. G., Han, Y., & He, Q. Y. (2012). clusterProfiler: An R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology, 16(5), 284–287. https://doi.org/10.1089/omi.2011.0118

Yu, Y., Mai, Y., Zheng, Y., & Shi, L. (2024). Assessing and mitigating batch effects in large-scale omics studies. Genome Biology, 25(1), 254. https://doi.org/10.1186/s13059-024-03401-9

