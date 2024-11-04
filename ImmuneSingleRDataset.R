install.packages("SCINA")
library(SCINA)

#Reference dataset
.ImmuneRef <- read.csv("~/Documents/Erson_Lab/Datasets/CellIDReferences/Immune Cells/E-CURD-122.marker_genesymbols_inferred_cell_type_SCINA.csv")
.ImmuneRef_list <- lapply(1:ncol(.ImmuneRef), function(i) {
  genes <- .ImmuneRef[, i]
  genes <- genes[!is.na(genes)]
  return(genes)
})
names(.ImmuneRef_list) <- colnames(.ImmuneRef)
.ImmuneRef <- .ImmuneRef_list
rm(.ImmuneRef_list)





#BELOW CODE ALL INTEGRATED INTO CUSTOM FUNCTION (SingleImmune)

#Extraction of Seurat Expression Data
Seurat31ExpressionData <- as.matrix(GetAssayData(seurat31, assay = "RNA", layer = "data"))

#Assignment of CellIDs 
#NOTE: significant overlap of genetic signatures between cell types--had to set rm_overlap to FALSE
results <- SCINA(Seurat31ExpressionData, ImmuneRef, rm_overlap = FALSE)
seurat31 <- AddMetaData(seurat31, metadata = results$cell_labels, col.name = "ImmuneCellLabels")
Idents(seurat31) <- "ImmuneCellLabels"

#UMAP plot based off CellIDs
seurat31 <- RunUMAP(seurat31, dims = 1:9)
DimPlot(seurat31, reduction = "umap", group.by = "ImmuneCellLabels", label = FALSE, pt.size = 0.5) + 
  ggtitle("UMAP with Immune Cell Types")

