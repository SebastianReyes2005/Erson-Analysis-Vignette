# Load necessary libraries
library(dplyr)
library(SingleR)
library(SingleCellExperiment)

# Load the metadata and expression data
.metadataOriginal <- read.csv("~/Documents/Erson_Lab/Datasets/CellIDReferences/Neuroendocrine Cells/GSE142653_pit_dev_CellInfo.csv.gz")
expressionOriginal <- read.csv("~/Documents/Erson_Lab/Datasets/CellIDReferences/Neuroendocrine Cells/GSE142653_pit_dev_5181_count.csv.gz")

# Ensure metadata in data frame
.metadataOriginal <- as.data.frame(.metadataOriginal)

# Filter expression data to match the metadata
.expressionFiltered <- expressionOriginal %>% select(X, one_of(.metadataOriginal$X))
rm(expressionOriginal) # Remove the original expression data from memory

# Set row names to the first column
rownames(.expressionFiltered) <- .expressionFiltered[,1]
.expressionFiltered <- .expressionFiltered[,-1]

# Convert expression data to numeric
.expressionFiltered <- .expressionFiltered %>%
  mutate(across(everything(), as.numeric))

# Ensure CellIDs are a character vector
.CellIDs <- as.character(.metadataOriginal$cell_type)
names(.CellIDs) <- colnames(.expressionFiltered)

# Convert metadata and expression data to matrices
.metadataOriginal <- as.matrix(.metadataOriginal)
.expressionFiltered <- as.matrix(.expressionFiltered)





#ALL BELOW FUNCTIONS IN CUSTOM FUNCTION (SingleNeuro)

# Function to run SingleR
  NeuroClusterLabels <- SingleR(test = as.SingleCellExperiment(seurat31),
                           ref = expressionFiltered,
                           labels = CellIDs)
  
#To add SingleR output to Seurat Object
seurat31 <- AddMetaData(seurat31, metadata = NeuroClusterLabels@listData[["labels"]], col.name = "NeuroCellLabels")

#Run UMAP based off SingleR output
seurat31 <- RunUMAP(seurat31, dims = 1:9)
DimPlot(seurat31, reduction = "umap", group.by = "NeuroCellLabels", label = FALSE, pt.size = 0.5) + 
  ggtitle("UMAP with Neuroendocrine Cell Types (P1/Cushing's Disease/Non-Tumor)")   
Idents(seurat31) <- "NeuroCellLabels"
VlnPlot(seurat31, features = "STAT6")
