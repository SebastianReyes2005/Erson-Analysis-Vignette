library(dplyr)
library(Seurat)
library(patchwork)
library(tibble)
library(SingleR)
library(SCINA)
library(celldex)
library(ggplot2)

# CreateSeuratObject
CreateSeuratObject <- function(fileDirectory, objectName, projectName) {
  rawCounts <- Read10X_h5(fileDirectory)
  seurat_obj <- Seurat::CreateSeuratObject(counts = rawCounts, 
                                           project = projectName, 
                                           min.cells = 3, 
                                           min.features = 200)
  assign(objectName, seurat_obj, envir = .GlobalEnv)
  return(objectName)
}

# AddMitoRna
AddMitoRna <- function(seuratObject){
  seuratObject[["percent.mt"]] <- Seurat::PercentageFeatureSet(seuratObject, pattern = "^MT-")
  return(seuratObject)
}

# MakeVlnPlot
MakeVlnPlot <- function(seuratObject) {
  print(VlnPlot(seuratObject, 
                features = c("nFeature_RNA", 
                             "nCount_RNA", 
                             "percent.mt"), 
                ncol = 3))
  return(seuratObject)
}

# QCSeurat
QCSeurat <- function(seuratObject, minFeature, maxFeature, percentMT) {
  seuratObject <- subset(seuratObject, subset = nFeature_RNA > minFeature & nFeature_RNA < maxFeature & percent.mt < percentMT)
  return(seuratObject)
}

# NormalizeAndVariable
NormalizeAndVariable <- function(seuratObject) {
  seuratObject <- NormalizeData(seuratObject)
  seuratObject <- FindVariableFeatures(seuratObject, selection.method = "vst", nfeatures = 2000)
  top10 <- head(VariableFeatures(seuratObject), 5)
  plot1 <- VariableFeaturePlot(seuratObject)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  print(plot2)
  return(seuratObject)
}

#Scale
Scale <- function(seuratObject) {
  seuratObject <- ScaleData(seuratObject, features = rownames(seuratObject))
  seuratObject <- ScaleData(seuratObject, vars.to.regress = "percent.mt",)
  return(seuratObject)
}

# PcaElbow
PcaElbow <- function(seuratObject) {
  seuratObject <- RunPCA(seuratObject, features = VariableFeatures(seuratObject))
  print(ElbowPlot(object = seuratObject))
  return(seuratObject)
}

# MakeHeatmap
MakeHeatmap <- function(seuratObject, PcNumberVector){
  DimHeatmap(seuratObject, dims = PcNumberVector, balanced = TRUE)
  return(seuratObject)
}

# MakeUMAP
MakeUMAP <- function(seuratObject1, PCs) {
  seuratObject2 <- RunUMAP(seuratObject1, dims = PCs)
  print(DimPlot(seuratObject2, reduction = "umap"))
  return(seuratObject2)
}

#FindPosMarkers
FindPosMarkers <- function(seuratObject){
  posMarkers <- FindAllMarkers(seuratObject, only.pos = TRUE)
  posMarkers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC >1)
  return(seuratObject)
}

# Top5Markers
Top5Markers <- function(seuratObject) {
  posMarkers <- FindAllMarkers(seuratObject, only.pos = TRUE)
  separatedMarkers <- posMarkers %>%
    dplyr::filter(avg_log2FC >1) %>%
    group_by(cluster) %>%
    slice_head(n = 5)
  print(separatedMarkers, n = Inf)
}
# SingleHPCA
SingleHPCA <- function(seuratObject){
  HpcaLabels <- SingleR(test = as.SingleCellExperiment(seuratObject), 
                        ref =    HumanPrimaryCellAtlasData(), 
                        labels = HumanPrimaryCellAtlasData()$label.main)
  seuratObject1 <- AddMetaData(seuratObject, 
                               metadata = HpcaLabels@listData[["labels"]], 
                               col.name = "HPCACellLabels")
  return(seuratObject1)
}

# SingleNeuro

NeuroMeta <- read.csv("Reference Datasets/NeuroMeta.csv.gz")

NeuroExpression <- read.csv("Reference Datasets/NeuroExpression.csv.gz")

columns_to_keep <- colnames(NeuroExpression) %in% NeuroMeta$X | colnames(NeuroExpression) == "X"

NeuroExpressionFiltered <- NeuroExpression[, columns_to_keep]

rm(NeuroExpression)

rownames(NeuroExpressionFiltered) <- NeuroExpressionFiltered[,1]
NeuroExpressionFiltered <- NeuroExpressionFiltered[,-1]

NeuroExpressionFiltered <- NeuroExpressionFiltered %>%
  mutate(across(everything(), as.numeric))

NeuroIDs <- as.character(NeuroMeta$cell_type)
names(NeuroIDs) <- colnames(NeuroExpressionFiltered)

rm(NeuroMeta)
NeuroExpressionFiltered <- as.matrix(NeuroExpressionFiltered)

SingleNeuro <- function(seuratObject) {
  # Function to run SingleR
  NeuroClusterLabels <- SingleR(test = as.SingleCellExperiment(seuratObject),
                                ref = NeuroExpressionFiltered,
                                labels = NeuroIDs)
  
  #To add SingleR output to Seurat Object
  seuratObject1 <- AddMetaData(seuratObject, metadata = NeuroClusterLabels@listData[["labels"]], col.name = "NeuroCellLabels")
  
  #To group everything by Neuro Cell Type
  Idents(seuratObject1) <- "NeuroCellLabels"
  return(seuratObject1)
}

# SingleImmune

ImmuneRef <- read.csv("Reference Datasets/ImmuneRef.csv")

ImmuneRef_list <- lapply(1:39, function(i) {
  genes <- ImmuneRef[, i]
  genes <- genes[!is.na(genes)]
  return(genes)
})

names(ImmuneRef_list) <- colnames(ImmuneRef)
ImmuneRef <- ImmuneRef_list
rm(ImmuneRef_list)
SingleImmune <- function(seuratObject) {
  
  SeuratObjectExpressionData <- as.matrix(GetAssayData(seuratObject, assay = "RNA", layer = "data"))
  
  results <- SCINA(SeuratObjectExpressionData, ImmuneRef, rm_overlap = FALSE)
  seuratObject1 <- AddMetaData(seuratObject, metadata = results$cell_labels, col.name = "ImmuneCellLabels")
  Idents(seuratObject1) <- "ImmuneCellLabels"
  return(seuratObject1)
}

# SingleUnbiasedNeuro
SingleUnbiasedNeuro <- function(seuratObject){
  
  cluster_label_counts <- table(seuratObject@meta.data$seurat_clusters, 
                                seuratObject@meta.data$NeuroCellLabels)
  
  most_abundant_labels <- apply(cluster_label_counts, 1, function(x) {
    names(x)[which.max(x)]
  })
  
  seuratObject@meta.data$UnbiasedNeuroCellLabels <- factor(seuratObject@meta.data$seurat_clusters, levels = names(most_abundant_labels), labels = most_abundant_labels)
  return(seuratObject)
}

# SingleUnbiasedImmune
SingleUnbiasedImmune <- function(seuratObject){
  cluster_label_counts <- table(seuratObject@meta.data$seurat_clusters, 
                                seuratObject@meta.data$ImmuneCellLabels)
  
  most_abundant_labels <- apply(cluster_label_counts, 1, function(x) {
    names(x)[which.max(x)] 
  })
  
  seuratObject@meta.data$UnbiasedImmuneCellLabels <- factor(seuratObject@meta.data$seurat_clusters, levels = names(most_abundant_labels), labels = most_abundant_labels)
  
  return(seuratObject)
}

# BarPlotIdPercent
BarPlotIdPercent <- function(seuratObject, Identity){
  
  TotalCells <- dim(seuratObject@meta.data)[1]
  
  cellIdentities <- seuratObject@meta.data[[Identity]]
  
  cellCounts <- table(cellIdentities)
  
  cellCountsDF <- as.data.frame(cellCounts)
  
  colnames(cellCountsDF) <- c("CellIdentity", "Count")
  
  ggplot(cellCountsDF, aes(x = CellIdentity, 
                           y = Count/TotalCells*100, 
                           fill = CellIdentity)) + 
    geom_bar(stat = "identity") +
    theme_minimal() + 
    labs(title = "Distribution of Cell Types", x = "Cell Identity", y = "Percentage of Cells") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(text = element_text("sans"))
}
