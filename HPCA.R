#To run SingleR using HPCA
HpcaLabels <- SingleR(test = as.SingleCellExperiment(seurat31), ref = HumanPrimaryCellAtlasData(), labels = HumanPrimaryCellAtlasData()$label.main)

#To add HPCA output to Seurat object
seurat31 <- AddMetaData(seurat31, metadata = HpcaLabels@listData[["labels"]], col.name = "HPCACellLabels")

#To create UMAP plot based off HPCA output
seurat31 <- RunUMAP(seurat31, dims = 1:9)
DimPlot(seurat31, reduction = "umap", group.by = "HPCACellLabels", label = FALSE, pt.size = 0.5) + 
  ggtitle("UMAP with HPCA Cell Types")   
