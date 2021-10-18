## for Figure3-G-1

# import library

## for plot and color
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(scater)

## for sparse matrix
library(Matrix)

## for hash
library(hash)

## for Seurat
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(patchwork)

## for Monocle3
library(monocle3)

## for colors
library(viridis)

## load all barcodes coordinates
## PB10all_cluster_pos.txt from Xiaonan
All_barcode_coord = read.csv("/media/gulab/GUDR2/PB10_new/PB10all_cluster_pos.txt", header = FALSE, sep = "\t")

## load only neurons coordinates
## Seg_cluster_only_neuron_pos_ABC.txt from file_convert.ipynb
Coord_for_Neuron = read.csv("/media/gulab/GUDR2/PB10_new/Seg_cluster_only_neuron_pos_ABC.txt", header = FALSE, sep = "\t")

## plot Calca/Tac1/Nts cells distribution on PBN region (-5.20mm)
tiff("/media/gulab/GUDR2/PB10_new/Cell_type_PBN_A.tiff", height = 7, width = 7, units = 'in', res = 500)
par(mar = rep(0, 4))

## plot all reads as grey on PBN region (-5.20mm)
plot(All_barcode_coord[,1:2], 
     pch=15,cex=.5, col = 'grey',
     ylim = c(6200, 7500), xlim = c(21800, 23000), yaxt = "n", xaxt='n') 
par(new=TRUE)

## plot Calca cell as green on PBN region (-5.20mm)
plot(Coord_for_Neuron[Coord_for_Neuron$V5==5,1:2], 
     pch=15,cex=.5, col = '#33A02C',
     ylim = c(6200, 7500), xlim = c(21800, 23000), yaxt = "n", xaxt='n') 
par(new=TRUE)

## plot Tac1 Cell as brown on PBN region (-5.20mm)
plot(Coord_for_Neuron[Coord_for_Neuron$V5==9,1:2], 
     pch=15,cex=.5, col = '#B15928',
     ylim = c(6200, 7500), xlim = c(21800, 23000), yaxt = "n", xaxt='n') 
par(new=TRUE)
dev.off()

## plot Calca/Tac1/Nts cells distribution on PBN region (-5.35mm)
tiff("/media/gulab/GUDR2/PB10_new/Cell_type_PBN_B.tiff", height = 7, width = 7, units = 'in', res = 500)
par(mar = rep(0, 4))

## plot all reads as grey on PBN region (-5.35mm)
plot(All_barcode_coord[,1:2], 
     pch=15,cex=.5, col = 'grey',
     ylim = c(6900, 8300), xlim = c(36300, 37500), yaxt = "n", xaxt='n') 
par(new=TRUE)

## plot Calca cell as green on PBN region (-5.35mm)
plot(Coord_for_Neuron[Coord_for_Neuron$V5==5,1:2], 
     pch=15,cex=.5, col = '#33A02C',
     ylim = c(6900, 8300), xlim = c(36300, 37500), yaxt = "n", xaxt='n') 
par(new=TRUE)

## plot Tac1 cell as blue on PBN region (-5.35mm)
plot(Coord_for_Neuron[Coord_for_Neuron$V5==9,1:2], 
     pch=15,cex=.5, col = '#B15928',
     ylim = c(6900, 8300), xlim = c(36300, 37500), yaxt = "n", xaxt='n') 
par(new=TRUE)
dev.off()

## plot Calca/Tac1/Nts cells distribution on PBN region (-5.50mm)
tiff("/media/gulab/GUDR2/PB10_new/Cell_type_PBN_C.tiff", height = 7, width = 7, units = 'in', res = 500)
par(mar = rep(0, 4))

## plot all reads as grey on PBN region (-5.50mm)
plot(All_barcode_coord[,1:2], 
     pch=15,cex=.5, col = 'grey',
     ylim = c(6500, 7800), xlim = c(49200, 50500), yaxt = "n", xaxt='n') 
par(new=TRUE)

## plot Calca cell as green on PBN region (-5.50mm)
plot(Coord_for_Neuron[Coord_for_Neuron$V5==5,1:2], 
     pch=15,cex=.5, col = '#33A02C',
     ylim = c(6500, 7800), xlim = c(49200, 50500), yaxt = "n", xaxt='n') 
par(new=TRUE)

## plot Tac1 cell as blue on PBN region (-5.50mm)
plot(Coord_for_Neuron[Coord_for_Neuron$V5==9,1:2], 
     pch=15,cex=.5, col = '#B15928',
     ylim = c(6500, 7800), xlim = c(49200, 50500), yaxt = "n", xaxt='n') 
par(new=TRUE)

dev.off()

## load gene cell matrix from dge file
library(DropSeq.util)
## PB10_allABC_cellExp.txt from Xiaonan
PB10.all <- loadSparseDge("/media/gulab/GUDR2/PB10_new/PB10_allABC_cellExp.txt")

## create seurat object from dge matrix
PBN.dge <- CreateSeuratObject(
  counts = PB10.all,
  project = 'SlideSeq',
  assay = 'Spatial'
)

### load location information for dge file
## PB10_allABC_coord.txt from Xiaonan
position_all <- read.csv(
  file = "/media/gulab/GUDR2/PB10_new/PB10_allABC_coord.txt"
)

## process spatial coordinates for dge Seurat object
rownames(x = position_all) <- position_all$binID
position_all <- position_all[, 2:3]

## set spatial coordinates for each cell of dge Seurat object
PBN.dge[['images']] <-new(
  Class = 'SlideSeq',
  assay = "Spatial",
  coordinates = position_all
)

## filter cells if cell has counts less than 2^7 = 128
PBN.dge <- subset(PBN.dge, subset = nCount_Spatial > 512)

# Single Cell pre-processing
## SCT Normalization, select 3000 highly variable genes
PBN.dge <- SCTransform(PBN.dge, assay = "Spatial", verbose = TRUE, variable.features.n = 3000)
# 3000 feature, PCA default dims, UMAP 15 dims, Find Neighbor 15 dims, Cluster 0.5res

# Dimension Reduction
# PBN.dge <- ScaleData(object = PBN.dge)

## Run PCA
PBN.dge <- RunPCA(PBN.dge, assay = "SCT", verbose = TRUE, rev.pca = TRUE)
# VlnPlot(PBN.dge, features = c("Sncg"))

## Run UMAP
PBN.dge <- RunUMAP(PBN.dge, reduction = "pca", dims = 1:15)

## Find Neighbors
PBN.dge <- FindNeighbors(PBN.dge, reduction = "pca", dims = 1:15)
PBN.dge <- FindClusters(PBN.dge, resolution = 0.5, verbose = TRUE)

## isolate all neurons and create a new Seurat object
only_Neuron.matrix <- PBN.dge$Spatial@counts[, PBN.dge$SCT_snn_res.0.5 == 12 | PBN.dge$SCT_snn_res.0.5 == 2 | PBN.dge$SCT_snn_res.0.5 == 1 | PBN.dge$SCT_snn_res.0.5 == 9 | PBN.dge$SCT_snn_res.0.5 == 11 | PBN.dge$SCT_snn_res.0.5 == 13]
only_Neuron.PBN <- CreateSeuratObject(counts = only_Neuron.matrix,
                                      project = 'SlideSeq',
                                      assay = 'Spatial')
only_Neuron.positions <- PBN.dge[["images"]]@coordinates[PBN.dge$SCT_snn_res.0.5 == 12 | PBN.dge$SCT_snn_res.0.5 == 2 | PBN.dge$SCT_snn_res.0.5 == 1 | PBN.dge$SCT_snn_res.0.5 == 9 | PBN.dge$SCT_snn_res.0.5 == 11 | PBN.dge$SCT_snn_res.0.5 == 13, ]
only_Neuron.PBN[['images']] <-new(
  Class = 'SlideSeq',
  assay = "Spatial",
  coordinates = only_Neuron.positions
)

## SCT normalization for only_Neuron.PBN 
only_Neuron.PBN <- SCTransform(only_Neuron.PBN, assay = "Spatial", verbose = FALSE, variable.features.n = 3000)

## Run PCA for only_Neuron.PBN 
only_Neuron.PBN <- RunPCA(only_Neuron.PBN, assay = "SCT", verbose = FALSE, rev.pca = TRUE, npcs = 50)

## Run UMAP for only_Neuron.PBN 
only_Neuron.PBN <- RunUMAP(only_Neuron.PBN, reduction = "pca", dims = 1:15)

## Find neighbors for only_Neuron.PBN 
only_Neuron.PBN <- FindNeighbors(only_Neuron.PBN, reduction = "pca", dims = 1:15)

## Find clusters for only_Neuron.PBN 
only_Neuron.PBN <- FindClusters(only_Neuron.PBN, resolution = 0.5, verbose = FALSE)

## Calca cell counts on PBN region (-5.50mm)
dim(only_Neuron.PBN[ ,only_Neuron.PBN$seurat_clusters == 5 & (only_Neuron.PBN[["images"]]@coordinates$x >= 21800 & only_Neuron.PBN[["images"]]@coordinates$x <= 23000 & only_Neuron.PBN[["images"]]@coordinates$y >= 6200 & only_Neuron.PBN[["images"]]@coordinates$y <= 7500)])
## 83

## Tac1 cell counts on PBN region (-5.50mm)
dim(only_Neuron.PBN[ ,only_Neuron.PBN$seurat_clusters == 9 & (only_Neuron.PBN[["images"]]@coordinates$x >= 21800 & only_Neuron.PBN[["images"]]@coordinates$x <= 23000 & only_Neuron.PBN[["images"]]@coordinates$y >= 6200 & only_Neuron.PBN[["images"]]@coordinates$y <= 7500)])
## 21

## Calca cell counts on PBN region (-5.35mm)
dim(only_Neuron.PBN[ ,only_Neuron.PBN$seurat_clusters == 5 & (only_Neuron.PBN[["images"]]@coordinates$x >= 36300 & only_Neuron.PBN[["images"]]@coordinates$x <= 37500 & only_Neuron.PBN[["images"]]@coordinates$y >= 6900 & only_Neuron.PBN[["images"]]@coordinates$y <= 8300)])
## 141

## Tac1 cell counts on PBN region (-5.35mm)
dim(only_Neuron.PBN[ ,only_Neuron.PBN$seurat_clusters == 9 & (only_Neuron.PBN[["images"]]@coordinates$x >= 36300 & only_Neuron.PBN[["images"]]@coordinates$x <= 37500 & only_Neuron.PBN[["images"]]@coordinates$y >= 6900 & only_Neuron.PBN[["images"]]@coordinates$y <= 8300)])
## 44

## Calca cell counts on PBN region (-5.20mm)
dim(only_Neuron.PBN[ ,only_Neuron.PBN$seurat_clusters == 5 & (only_Neuron.PBN[["images"]]@coordinates$x >= 49200 & only_Neuron.PBN[["images"]]@coordinates$x <= 50500 & only_Neuron.PBN[["images"]]@coordinates$y >= 6500 & only_Neuron.PBN[["images"]]@coordinates$y <= 7800)])
## 17

## Tac1 cell counts on PBN region (-5.20mm)
dim(only_Neuron.PBN[ ,only_Neuron.PBN$seurat_clusters == 9 & (only_Neuron.PBN[["images"]]@coordinates$x >= 49200 & only_Neuron.PBN[["images"]]@coordinates$x <= 50500 & only_Neuron.PBN[["images"]]@coordinates$y >= 6500 & only_Neuron.PBN[["images"]]@coordinates$y <= 7800)])
## 65

## bar plot for cell counts
png("/media/gulab/GUDR2/PB10_new/barplot.png", height = 7, width = 3, units = 'in', res = 300)
par(mfrow=c(1,1),bg='white')
Values <- matrix(c(83, 21, 141, 44, 17, 65), nrow = 2, ncol = 3)
colors = c('#33A02C', '#B15928')
barplot(Values, col = colors)
dev.off()
