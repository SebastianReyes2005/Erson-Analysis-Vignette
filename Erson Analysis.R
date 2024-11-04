#For 31 : P1 and CD1 
CreateSeuratObject("~/Documents/Erson_Lab/Datasets/Asuzu Samples/Asuzu31.h5",
                   "seurat31", 
                   "P1&CD1" )
seurat31 <- AddMitoRna(seurat31)
MakeVlnPlot(seurat31)
dim(seurat31)
seurat31 <- QCSeurat(seurat31, 200, 3000, 20)
dim(seurat31)
MakeFeaturePlot(seurat31, "nCount_RNA", "percent.mt")
MakeFeaturePlot(seurat31, "nCount_RNA", "nFeature_RNA")
seurat31 <- NormalizeAndVariable(seurat31)
seurat31 <- Scale(seurat31)
seurat31 <- PcaElbow(seurat31)
MakeHeatmap(seurat31, 1:10)
seurat31 <- FindCellClusters(seurat31, 1:9, 0.5)
seurat31 <- RunUMAP(seurat31, dims = 1:9)
DimPlot(seurat31, reduction = "umap")
seurat31 <- FindPosMarkers(seurat31) 
seurat31 <- SingleNeuro(seurat31)
print(DimPlot(seurat31, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle("Neuro (P1/CD/Non-Tumor)"))
Idents(seurat31) <- seurat31@meta.data$NeuroCellLabels
BarPlot(seurat31, 2956) + ggtitle("Cell Type Distribution (P1/Cushing's Disease/Non-Tumor)")


#For immune characterization
SeuratFiltered <- subset(seurat32, subset = NeuroCellLabels == "Imm")
CategoryCounts <- SeuratFiltered@meta.data %>%
  group_by(ImmuneCellLabels) %>%
  tally()
ValidCategories <- CategoryCounts %>% 
  filter(n >= 5) %>%
  pull(ImmuneCellLabels)

print(DimPlot(subset(seurat33, subset = NeuroCellLabels == "Imm" & ImmuneCellLabels %in% ValidCategories),
              reduction = "umap",
              group.by = "ImmuneCellLabels",
              repel = FALSE,
              label = FALSE, 
              label.size = 6) + ggtitle("")) +
  theme(legend.text = element_text(size = 10))
FeaturePlot((subset(seurat38, subset = NeuroCellLabels == "Imm")),
            features = "FOS") + 
  ggtitle("FOS")





#For 32: P2 and CD2 Tumor
CreateSeuratObject("~/Documents/Erson_Lab/Datasets/Asuzu Samples/Asuzu32.h5",
                   "seurat32", 
                   "CD2tumor" )
seurat32 <- AddMitoRna(seurat32)
MakeVlnPlot(seurat32)
dim(seurat32)
seurat32 <- QCSeurat(seurat32, 200, 3000, 20)
dim(seurat32)
MakeFeaturePlot(seurat32, "nCount_RNA", "nFeature_RNA")
seurat32 <- NormalizeAndVariable(seurat32)
seurat32 <- Scale(seurat32)
seurat32 <- PcaElbow(seurat32)
seurat32 <- RunUMAP(seurat32, dims = 1:11)
seurat32 <- FindCellClusters(seurat32, 1:11, 0.5)
DimPlot(seurat32, reduction = "umap")
seurat32 <- FindPosMarkers(seurat32)
seurat32 <- SingleNeuro(seurat32)
print(DimPlot(seurat32, reduction = "umap", group.by = "NeuroCellLabels") + ggtitle("Neuro (P2/CD/Tumor)"))
Idents(seurat32) <- seurat32@meta.data$NeuroCellLabels
BarPlot(seurat32, 4559) + ggtitle("Cell Type Distribution (P2/Cushing's Disease/Tumor)")





#For 33: P2 and CD2 Normal
CreateSeuratObject("~/Documents/Erson_Lab/Datasets/Asuzu Samples/Asuzu33.h5",
                   "seurat33", 
                   "P2&CD2Normal" )
seurat33 <- AddMitoRna(seurat33)
MakeVlnPlot(seurat33)
dim(seurat33)
seurat33 <- QCSeurat(seurat33, 200, 3000, 20)
dim(seurat33)
MakeFeaturePlot(seurat33, "nCount_RNA", "percent.mt")
MakeFeaturePlot(seurat33, "nCount_RNA", "nFeature_RNA")
seurat33 <- NormalizeAndVariable(seurat33)
seurat33 <- Scale(seurat33)
seurat33 <- PcaElbow(seurat33)
MakeHeatmap(seurat33, 1:10)
seurat33 <- FindCellClusters(seurat33, 1:12, 0.5)
seurat33 <- RunUMAP(seurat33, dims = 1:12)
DimPlot(seurat33, reduction = "umap")
seurat33 <- FindPosMarkers(seurat33)
seurat33 <- SingleNeuro(seurat33)
print(DimPlot(seurat33, reduction = "umap", group.by = "NeuroCellLabels") + ggtitle("Neuro (P2/CD/Non-Tumor)"))
Idents(seurat33) <- seurat33@meta.data$seurat_clusters
BarPlot(seurat33, 4975) + ggtitle("Cell Type Distribution (P2/Cushing's Disease/Non-Tumor)")

VlnPlot(seurat33, features = "STAT6") +  NoLegend() + ggtitle("STAT6 (P2/CD/Non-Tumor)")

install.packages("clipr")
library(clipr)
ClusterDifferences <- FindMarkers(seurat33, ident.1 = 7, ident.2 = 11)
ClusterDifferences <- rownames_to_column(ClusterDifferences, var = "Genes")
write_clip(head(ClusterDifferences$Genes, n = 300))

#For 34: P3 and CD3 tumor
CreateSeuratObject("~/Documents/Erson_Lab/Datasets/Asuzu Samples/Asuzu34.h5",
                   "seurat34", 
                   "P3&CD3tumor" )
seurat34 <- AddMitoRna(seurat34)
MakeVlnPlot(seurat34)
dim(seurat34)
seurat34 <- QCSeurat(seurat34, 200, 3000, 20)
dim(seurat34)
MakeFeaturePlot(seurat34, "nCount_RNA", "percent.mt")
MakeFeaturePlot(seurat34, "nCount_RNA", "nFeature_RNA")
seurat34 <- NormalizeAndVariable(seurat34)
seurat34 <- Scale(seurat34)
seurat34 <- PcaElbow(seurat34)
MakeHeatmap(seurat34, 11:20)
seurat34 <- FindCellClusters(seurat34, 1:10, 0.5)
seurat34 <- FindPosMarkers(seurat34)
seurat34 <- SingleNeuro(seurat34)
seurat34 <- RunUMAP(seurat34, dims = 1:10)
print(DimPlot(seurat34, reduction = "umap", group.by = "NeuroCellLabels") + ggtitle("Neuro (P3/CD/Tumor)"))
Idents(seurat34) <- seurat34@meta.data$seurat_clusters
BarPlot(seurat34, 1958) + ggtitle("Cell Type Distribution (P3/Cushing's Disease/Tumor)")

VlnPlot(seurat34, features = "STAT6") +  NoLegend() + ggtitle("STAT6 (P3/CD/Tumor)")


#For 35: P3 and CD3 normal
CreateSeuratObject("~/Documents/Erson_Lab/Datasets/Asuzu Samples/Asuzu35.h5",
                   "seurat35", 
                   "P3&CD3Normal" )
seurat35 <- AddMitoRna(seurat35)
MakeVlnPlot(seurat35)
dim(seurat35)
seurat35 <- QCSeurat(seurat35, 200, 3000, 20)
dim(seurat35)
MakeFeaturePlot(seurat35, "nCount_RNA", "percent.mt")
MakeFeaturePlot(seurat35, "nCount_RNA", "nFeature_RNA")
seurat35 <- NormalizeAndVariable(seurat35)
seurat35 <- Scale(seurat35)
seurat35 <- PcaElbow(seurat35)
MakeHeatmap(seurat35, 1:15)
seurat35 <- FindCellClusters(seurat35, 1:11, 0.5)
seurat35 <- FindPosMarkers(seurat35)
seurat35 <- SingleNeuro(seurat35)
seurat35 <- RunUMAP(seurat35, dims = 1:11)
print(DimPlot(seurat35, reduction = "umap", group.by = "NeuroCellLabels") + ggtitle("Neuro (P3/CD/Non-Tumor)"))
Idents(seurat35) <- seurat35@meta.data$seurat_clusters
BarPlot(seurat35, 3763) + ggtitle("Cell Type Distribution (P3/Cushing's Disease/Non-Tumor)")

VlnPlot(seurat35, features = "STAT6") +  NoLegend() + ggtitle("STAT6 (P3/CD/Non-Tumor)")


#For 36: P4, and GH tumor
CreateSeuratObject("~/Documents/Erson_Lab/Datasets/Asuzu Samples/Asuzu36.h5",
                   "seurat36", 
                   "P4&GHtumor")
seurat36 <- AddMitoRna(seurat36)
MakeVlnPlot(seurat36)
dim(seurat36)
seurat36 <- QCSeurat(seurat36, 200, 3000, 20)
dim(seurat36)
MakeFeaturePlot(seurat36, "nCount_RNA", "percent.mt")
MakeFeaturePlot(seurat36, "nCount_RNA", "nFeature_RNA")
seurat36 <- NormalizeAndVariable(seurat36)
seurat36 <- Scale(seurat36)
seurat36 <- PcaElbow(seurat36)
MakeHeatmap(seurat36, 1:9)
seurat36 <- FindCellClusters(seurat36, 1:9, 0.5)
seurat36 <- FindPosMarkers(seurat36) 
seurat36 <- SingleNeuro(seurat36)
seurat36 <- RunUMAP(seurat36, dims = 1:9)
print(DimPlot(seurat36, reduction = "umap", group.by = "NeuroCellLabels") + ggtitle("Neuro (P4/GH/Tumor)"))
Idents(seurat36) <- seurat36@meta.data$NeuroCellLabels
BarPlot(seurat36, 3871) + ggtitle("Cell Type Distribution (P4/Growth Hormone/Tumor)")

VlnPlot(seurat36, features = "STAT6") +  NoLegend() + ggtitle("STAT6 (P4/GH/Tumor)")


#For 37: P4 and GH normal
CreateSeuratObject("~/Documents/Erson_Lab/Datasets/Asuzu Samples/Asuzu37.h5",
                   "seurat37", 
                   "P4&GHnormal")
seurat37 <- AddMitoRna(seurat37)
MakeVlnPlot(seurat37)
dim(seurat37)
seurat37 <- QCSeurat(seurat37, 200, 3000, 20)
dim(seurat37)
MakeFeaturePlot(seurat37, "nCount_RNA", "percent.mt")
MakeFeaturePlot(seurat37, "nCount_RNA", "nFeature_RNA")
seurat37 <- NormalizeAndVariable(seurat37)
seurat37 <- Scale(seurat37)
seurat37 <- PcaElbow(seurat37)
MakeHeatmap(seurat35, 1:13)
seurat37 <- FindCellClusters(seurat37, 1:11, 0.5)
seurat37 <- FindPosMarkers(seurat37) 
seurat37 <- SingleNeuro(seurat37)
seurat37 <- RunUMAP(seurat37, dims = 1:11)
print(DimPlot(seurat37, reduction = "umap", group.by = "NeuroCellLabels") + ggtitle("Neuro (P4/GH/Non-Tumor)"))
Idents(seurat37) <- seurat37@meta.data$seurat_clusters
BarPlot(seurat37, 2946) + ggtitle("Cell Type Distribution (P4/Growth Hormone/Non-Tumor)")
BarPlot(seurat37, 2946)

VlnPlot(seurat37, features = "STAT6") +  NoLegend() + ggtitle("STAT6 (P4/GH/Non-Tumor)")


#For 38: P5 and NFPA tumor
VlnPlot(seurat38, features = c("nFeature_RNA", "percent.mt"), ncol = 2)

CreateSeuratObject("~/Documents/Erson_Lab/Datasets/Asuzu Samples/Asuzu38.h5",
                   "seurat38", 
                   "P5&NFPAtumor" )
seurat38 <- AddMitoRna(seurat38)
MakeVlnPlot(seurat38)
dim(seurat38)
seurat38 <- QCSeurat(seurat38, 200, 3000, 20)
dim(seurat38)
MakeFeaturePlot(seurat38, "nCount_RNA", "percent.mt")
MakeFeaturePlot(seurat38, "nCount_RNA", "nFeature_RNA")
seurat38 <- NormalizeAndVariable(seurat38)
seurat38 <- Scale(seurat38)
seurat38 <- PcaElbow(seurat38)
MakeHeatmap(seurat38, 1:12)
seurat38 <- FindCellClusters(seurat38, 1:12, 0.5)
seurat38 <- FindPosMarkers(seurat38) 
seurat38 <- SingleNeuro(seurat38)
seurat38 <- RunUMAP(seurat38, dims = 1:12)
print(DimPlot(seurat38, reduction = "umap", group.by = "NeuroCellLabels") + ggtitle("Neuro (P5/NFPA/Tumor)"))
Idents(seurat38) <- seurat38@meta.data$NeuroCellLabels
BarPlot(seurat38, 8339) + ggtitle("Cell Type Distribution (P5/Non-Functional Pituitary Adenoma/Tumor Tissue)")
BarPlot(seurat38, 8339)


VlnPlot(seurat38, features = "STAT6") +  NoLegend() + ggtitle("STAT6 (P5/NFPA/Tumor)")

