# Load Packages ----

library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(dplyr)
library(SingleR)
library(celldex)
library(clusterProfiler)
library(org.Mm.eg.db)

set.seed(23)

# Load data ----

seurat_obj <- readRDS("seurat_ass4.rds")

# seurat_obj <- readRDS("checkpoint_annotated.rds")

head(seurat_obj@meta.data)

# Map metadata to standard names

seurat_obj$tissue <- seurat_obj$organ_custom
seurat_obj$timepoint <- seurat_obj$time

table(seurat_obj$tissue)
table(seurat_obj$timepoint)

# Quality Control ----

# Calculate mtRNA percentage

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

png("images/qc_before_filter.png", width = 900, height = 500)
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
dev.off()

# Filter based on thresholds

seurat_obj <- subset(seurat_obj,
  subset = nFeature_RNA > 200 &
    percent.mt < 15
)

png("images/qc_after_filter.png", width = 900, height = 500)
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
dev.off()

# Normalization

seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj$biosample_id)

seurat_obj <- NormalizeData(seurat_obj)

# Variable features

seurat_obj <- FindVariableFeatures(seurat_obj)

# Scaling

seurat_obj <- ScaleData(seurat_obj)

# PCA

seurat_obj <- RunPCA(seurat_obj)

png("images/elbow_plot.png", width = 800, height = 500)
ElbowPlot(seurat_obj, ndims = 50)
dev.off()

# Integration using Harmony to correct batch effects ----

seurat_obj <- IntegrateLayers(
  object         = seurat_obj,
  method         = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction  = "harmony",
  verbose        = FALSE
)

seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])

saveRDS(seurat_obj, "checkpoint_integrated.rds")

# Clustering & UMAPs ----

seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:30)

seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, reduction = "harmony")

# UMAP of clusters

png("images/umap_clusters.png", width = 800, height = 600)
DimPlot(seurat_obj, reduction = "umap", label = TRUE)
dev.off()

# UMAP colored by timepoint and tissue to check integration

png("images/umap_time.png", width = 800, height = 600)
DimPlot(seurat_obj, reduction = "umap", group.by = "timepoint")
dev.off()

png("images/umap_tissue.png", width = 800, height = 600)
DimPlot(seurat_obj, reduction = "umap", group.by = "tissue")
dev.off()

saveRDS(seurat_obj, "checkpoint_clustered.rds")

# Automatic Annotation (SingleR) ----

ref <- celldex::MouseRNAseqData()

expr <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")

pred <- SingleR(
  test     = expr,
  ref      = ref,
  labels   = ref$label.main,
  clusters = Idents(seurat_obj)
)

label_map <- setNames(pred$labels, rownames(pred))

celltype_vector <- label_map[as.character(seurat_obj$seurat_clusters)]
names(celltype_vector) <- colnames(seurat_obj)

seurat_obj@meta.data$celltype <- celltype_vector

table(seurat_obj@meta.data$celltype)

# Plot

png("images/umap_celltype.png", width = 800, height = 600)
DimPlot(seurat_obj, group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend()
dev.off()

saveRDS(seurat_obj, "checkpoint_annotated.rds")

# Feature plot, genes of interest ----

# Markers to validate annotation

png("images/featureplot_markers.png", width = 900, height = 700)
FeaturePlot(seurat_obj, features = c("Cd3d", "Cd19", "Adgre1", "Omp"))
dev.off()

# Differntial Expression of D02 and D05 - Pseudobulk with DESeq2

# Subset to macrophages only first

Idents(seurat_obj) <- "celltype"
macro <- subset(seurat_obj, idents = "Macrophages")

# Only keep D02 and D05 timepoints

macro <- subset(macro, subset = timepoint %in% c("D02", "D05"))

# Pseudobulk aggregation

macro <- subset(macro, subset = mouse_id != "")
table(macro$mouse_id)

pseudo_macro <- AggregateExpression(macro,
  assays = "RNA",
  return.seurat = TRUE,
  group.by = c("mouse_id", "timepoint")
)

# Create identity column

Idents(pseudo_macro) <- "timepoint"

# DE with DESeq2
macro.de <- FindMarkers(pseudo_macro,
  ident.1 = "D02",
  ident.2 = "D05",
  test.use = "DESeq2",
  verbose = FALSE
)

head(macro.de, n = 10)
write.csv(macro.de, "DE_Macrophages_D02_vs_D05.csv")
macro.de$gene <- rownames(macro.de)

# Violin plot of C1qb and Isg15 genes

png("images/vlnplot_DE_genes.png", width = 800, height = 500)
VlnPlot(macro,
  features = c("C1qb", "Isg15"),
  group.by = "timepoint",
  pt.size = 0
)
dev.off()

# GSEA ----

# Convert gene symbols to Entrez IDs

gene_df <- bitr(rownames(macro.de),
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Mm.eg.db
)

macro.de <- left_join(macro.de, gene_df, by = c("gene" = "SYMBOL"))

# Build ranked gene list sorted by log2FC

gene_list <- macro.de$avg_log2FC
names(gene_list) <- macro.de$ENTREZID
gene_list <- sort(gene_list[!is.na(names(gene_list))], decreasing = TRUE)

gsea_result <- gseGO(
  geneList     = gene_list,
  OrgDb        = org.Mm.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",
  pvalueCutoff = 0.05
)

png("images/gsea_ridgeplot.png", width = 900, height = 1000)
ridgeplot(gsea_result) + ggtitle("GSEA: Macrophages D02 vs D05")
dev.off()

# ORA (Over-Representation Analysis) ----

# Upregulated at D02 (early infection)

sig_genes_D02 <- macro.de$gene[macro.de$p_val_adj < 0.05 & macro.de$avg_log2FC > 0]
sig_genes_D02 <- sig_genes_D02[!is.na(sig_genes_D02)]

ego_D02 <- enrichGO(
  gene         = sig_genes_D02,
  OrgDb        = org.Mm.eg.db,
  keyType      = "SYMBOL",
  ont          = "BP",
  pvalueCutoff = 0.05
)

# Upregulated at D05 (peak myeloid recruitment)

sig_genes_D05 <- macro.de$gene[macro.de$p_val_adj < 0.05 & macro.de$avg_log2FC < 0]
sig_genes_D05 <- sig_genes_D05[!is.na(sig_genes_D05)]

ego_D05 <- enrichGO(
  gene         = sig_genes_D05,
  OrgDb        = org.Mm.eg.db,
  keyType      = "SYMBOL",
  ont          = "BP",
  pvalueCutoff = 0.05
)

# Plot both

png("images/ora_dotplot_D02.png", width = 800, height = 600)
dotplot(ego_D02) + ggtitle("ORA: Genes Upregulated in Macrophages at D02")
dev.off()

png("images/ora_dotplot_D05.png", width = 800, height = 600)
dotplot(ego_D05) + ggtitle("ORA: Genes Upregulated in Macrophages at D05")
dev.off()

# ORA Comparison plot

compare_GO <- function(df) {
  up_genes <- df %>%
    filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
    pull(gene) %>%
    na.omit() %>%
    unique()

  down_genes <- df %>%
    filter(p_val_adj < 0.05 & avg_log2FC < 0) %>%
    pull(gene) %>%
    na.omit() %>%
    unique()

  compare_df <- compareCluster(
    geneCluster = list(
      Upregulated = up_genes,
      Downregulated = down_genes
    ),
    fun = "enrichGO",
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pvalueCutoff = 0.05
  )

  return(compare_df)
}


comp <- compare_GO(macro.de)

png("images/ora_compare_dotplot.png", width = 900, height = 700)
dotplot(comp) + ggtitle("GO Comparison: Up vs Downregulated Genes")
dev.off()
