install.packages("bigmemory")
library(bigmemory)
install.packages("dplyr")
library(dplyr)
install.packages("Seurat")
library(Seurat)
install.packages("patchwork")
library(patchwork)
install.packages("tibble")
library(tibble)
install.packages("SingleR")
library(SingleR)
BiocManager::install("celldex")
library(celldex)
install.packages("ggplot2")
library(ggplot2)


#Creating functions
CreateSeuratObject <- function(fileDirectory, objectName, projectName) {
  rawCounts <- Read10X_h5(fileDirectory)
  seurat_obj <- Seurat::CreateSeuratObject(counts = rawCounts, project = projectName, min.cells = 3, min.features = 200)
  assign(objectName, seurat_obj, envir = .GlobalEnv)
  return(objectName)
}

MakeVlnPlot <- function(seuratObject) {
  print(VlnPlot(seuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  return(seuratObject)
}

AddMitoRna <- function(seuratObject){
  seuratObject[["percent.mt"]] <- Seurat::PercentageFeatureSet(seuratObject, pattern = "^MT-")
  return(seuratObject)
}

QCSeurat <- function(seuratObject, minFeature, maxFeature, percentMT) {
  seuratObject <- subset(seuratObject, subset = nFeature_RNA > minFeature & nFeature_RNA < maxFeature & percent.mt < percentMT)
  return(seuratObject)
}

MakeFeaturePlot <- function(seuratObject, feature1, feature2) {
  plot1 <- FeatureScatter(seuratObject, feature1 = feature1, feature2 = feature2)
  print(plot1)
}

NormalizeAndVariable <- function(seuratObject) {
  seuratObject <- NormalizeData(seuratObject)
  seuratObject <- FindVariableFeatures(seuratObject, selection.method = "vst", nfeatures = 2000)
  top10 <- head(VariableFeatures(seuratObject), 5)
  plot1 <- VariableFeaturePlot(seuratObject)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  print(plot2)
  return(seuratObject)
}

Scale <- function(seuratObject) {
  seuratObject <- ScaleData(seuratObject, features = rownames(seuratObject))
  seuratObject <- ScaleData(seuratObject, vars.to.regress = "percent.mt",)
  return(seuratObject)
}

PcaElbow <- function(seuratObject) {
  seuratObject <- RunPCA(seuratObject, features = VariableFeatures(seuratObject))
  print(ElbowPlot(object = seuratObject))
  return(seuratObject)
}

MakeHeatmap <- function(seuratObject, PcNumberVector){
  DimHeatmap(seuratObject, dims = PcNumberVector, balanced = TRUE)
  return(seuratObject)
}

FindCellClusters <- function(seuratObject, PCs, resolution) {
  seuratObject2 <- FindNeighbors(seuratObject, dims = PCs)
  seuratObject3 <- FindClusters(seuratObject2, resolution = resolution)
  return(seuratObject3)
}

MakeUMAP <- function(seuratObject1, PCs) {
  seuratObject2 <- RunUMAP(seuratObject1, dims = PCs)
  print(DimPlot(seuratObject2, reduction = "umap"))
  return(seuratObject2)
}

FindPosMarkers <- function(seuratObject){
  posMarkers <- FindAllMarkers(seuratObject, only.pos = TRUE)
  posMarkers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC >1)
  return(seuratObject)
}


Top5Markers <- function(seuratObject) {
  posMarkers <- FindAllMarkers(seuratObject, only.pos = TRUE)
  separatedMarkers <- posMarkers %>%
    dplyr::filter(avg_log2FC >1) %>%
    group_by(cluster) %>%
    slice_head(n = 5)
  print(separatedMarkers, n = Inf)
}

SingleHPCA <- function(seuratObject){
  HpcaLabels <- SingleR(test = as.SingleCellExperiment(seuratObject), ref = HumanPrimaryCellAtlasData(), labels = HumanPrimaryCellAtlasData()$label.main)
  #To add HPCA output to Seurat object
  seuratObject1 <- AddMetaData(seuratObject, metadata = HpcaLabels@listData[["labels"]], col.name = "HPCACellLabels")
  return(seuratObject1)
}


SingleImmune <- function(seuratObject) {
#Extraction of Seurat Expression Data
SeuratObjectExpressionData <- as.matrix(GetAssayData(seuratObject, assay = "RNA", layer = "data"))

#Assignment of CellIDs 
#NOTE: significant overlap of genetic signatures between cell types--had to set rm_overlap to FALSE
results <- SCINA(SeuratObjectExpressionData, ImmuneRef, rm_overlap = FALSE)
seuratObject1 <- AddMetaData(seuratObject, metadata = results$cell_labels, col.name = "ImmuneCellLabels")
Idents(seuratObject1) <- "ImmuneCellLabels"
return(seuratObject1)
}

SingleNeuro <- function(seuratObject) {
  # Function to run SingleR
  NeuroClusterLabels <- SingleR(test = as.SingleCellExperiment(seuratObject),
                                ref = expressionFiltered,
                                labels = CellIDs)
  
  #To add SingleR output to Seurat Object
  seuratObject1 <- AddMetaData(seuratObject, metadata = NeuroClusterLabels@listData[["labels"]], col.name = "NeuroCellLabels")
  
  #To group everything by Neuro Cell Type
  Idents(seuratObject1) <- "NeuroCellLabels"
  return(seuratObject1)
}

AssignCellIDs <- function(seuratObject, CellIdVector){
  names(CellIdVector) <- levels(seuratObject)
  seuratObject <- RenameIdents(seuratObject, CellIdVector)
  return(seuratObject)
}

BarPlot <- function(seuratObject, TotalCells){
  cellIdentities <- Idents(seuratObject)
  cellCounts <- table(cellIdentities)
  cellCountsDF <- as.data.frame(cellCounts)
  colnames(cellCountsDF) <- c("CellIdentity", "Count")
  ggplot(cellCountsDF, aes(x = CellIdentity, y = Count/TotalCells*100, fill = CellIdentity )) +
    geom_bar(stat = "identity") +
    theme_minimal() + 
    labs(title = "Distribution of Cell Types", x = "Cell Identity", y = "Percentage of Cells") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(text = element_text("sans"))
}

  








