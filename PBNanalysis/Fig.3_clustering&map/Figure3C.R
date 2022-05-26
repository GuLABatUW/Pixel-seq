## for Figure3-C

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


## load PBN data from 4 section 30d1, 30d2, CN1, and CN2

## PBNct_S4.rds from Xiaonan
PBN_4 = readRDS("/media/gulab/GUDR2/PB10_new/PBNct_S4.rds")
PBN_4_df <- data.frame(PBN_4$integrated_snn_res.0.3)
write.csv(PBN_4_df, "/media/gulab/GUDR2/PB10_new/PBN_4_cluster.csv")

table(PBN_4[, PBN_4$orig.ident == '30d1']$seurat_clusters)
## set color for UMAP
mycol = c("#66C2A5", 
          "#FC8D62", 
          "#8DA0CB", 
          "#E78AC3", 
          "#FDBF6F", 
          "#FFD92F", 
          "#E5C494", 
          "#E31A1C", 
          "#A6CEE3", 
          "#B3B3B3", 
          "#B2DF8A", 
          "#33A02C", 
          "#FB9A99", 
          "#1F78B4", 
          "#A6D854", 
          "#FF7F00", 
          "#CAB2D6", 
          "#8DD3C7", 
          "#01665E", 
          "#B15928", 
          "#6A3D9A", 
          "#FFFFB3" 
          )

## generate UMAP plot

## dimension reduction and clustering 
PBN_4 <- RunPCA(PBN_4, verbose = TRUE, npcs = 50)
ElbowPlot(PBN_4)
PBN_4 <- RunUMAP(PBN_4, reduction = "pca", dims = 1:22)
PBN_4 <- FindNeighbors(PBN_4, reduction = "pca", dims = 1:22)
PBN_4 <- FindClusters(PBN_4, verbose = TRUE, resolution = 0.3)

## Marker list
PBN_4.markers <- FindAllMarkers(PBN_4, assay = "SCT", only.pos = TRUE)
View(PBN_4.markers)

## generate UMAP plot
png("/media/gulab/GUDR2/PB10_new/All4_UMAP.png", height = 7, width = 7, units = 'in', res = 300)
DimPlot(PBN_4, reduction = "umap", label = FALSE, cols = mycol, shuffle = TRUE, pt.size = .5)
dev.off()

png("/media/gulab/GUDR2/PB10_new/30d1_UMAP.png", height = 7, width = 7, units = 'in', res = 300)
DimPlot(PBN_4[, PBN_4$orig.ident=='30d1'], reduction = "umap", label = FALSE, cols = mycol, shuffle = TRUE, pt.size = .5)
dev.off()

png("/media/gulab/GUDR2/PB10_new/30d2_UMAP.png", height = 7, width = 7, units = 'in', res = 300)
DimPlot(PBN_4[, PBN_4$orig.ident=='30d2'], reduction = "umap", label = FALSE, cols = mycol, shuffle = TRUE, pt.size = .5)
dev.off()

png("/media/gulab/GUDR2/PB10_new/CNR1_UMAP.png", height = 7, width = 7, units = 'in', res = 300)
DimPlot(PBN_4[, PBN_4$orig.ident=='CNR1'], reduction = "umap", label = FALSE, cols = mycol, shuffle = TRUE, pt.size = .5)
dev.off()

png("/media/gulab/GUDR2/PB10_new/CNR2_UMAP.png", height = 7, width = 7, units = 'in', res = 300)
DimPlot(PBN_4[, PBN_4$orig.ident=='CNR2'], reduction = "umap", label = FALSE, cols = mycol, shuffle = TRUE, pt.size = .5)
dev.off()
