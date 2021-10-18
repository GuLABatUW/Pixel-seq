## for Figure 5-C

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
## neuronSeurat.rds from Xiaonan
## otherSeurat.rds from Xiaonan
Neuron = readRDS("/media/gulab/GUDR2/PB10_new/neuronSeurat.rds")
Other = readRDS("/media/gulab/GUDR2/PB10_new/otherSeurat.rds")
DefaultAssay(Neuron) <- "SCT"
DefaultAssay(Other) <- "SCT"
## check gene expression level in all clusters
## Apoe
index = 4
gene_name = 'Apoe'
(sum(Other[, Other$orig.ident=='CNR1' & Other$seurat_clusters == index & Other[, Other$orig.ident=='CNR1' & Other$seurat_clusters == index]@assays$SCT@data[gene_name, ] > 0]@assays$SCT@data[gene_name, ]) 
  + sum(Other[, Other$orig.ident=='30d1' & Other$seurat_clusters == index & Other[, Other$orig.ident=='30d1' & Other$seurat_clusters == index]@assays$SCT@data[gene_name, ] > 0]@assays$SCT@data[gene_name, ]))/
  ((1+sum(GetAssayData(object = Other[, Other$orig.ident=='CNR1' & Other$seurat_clusters == index], slot = "counts")[gene_name,]>0)) +
     (1+sum(GetAssayData(object = Other[, Other$orig.ident=='30d1' & Other$seurat_clusters == index], slot = "counts")[gene_name,]>0)))
## Astro1 1.256128
## Microglia 0.551664
## A2 0.806
## A3 0.559
## A4 0.431
## A5 0.491

index = 11
gene_name = 'Apoe'
(sum(Neuron[, Neuron$orig.ident=='CNR1' & Neuron$seurat_clusters == index & Neuron[, Neuron$orig.ident=='CNR1' & Neuron$seurat_clusters == index]@assays$SCT@data[gene_name, ] > 0]@assays$SCT@data[gene_name, ]) 
  + sum(Neuron[, Neuron$orig.ident=='30d1' & Neuron$seurat_clusters == index & Neuron[, Neuron$orig.ident=='30d1' & Neuron$seurat_clusters == index]@assays$SCT@data[gene_name, ] > 0]@assays$SCT@data[gene_name, ]))/
  ((1+sum(GetAssayData(object = Neuron[, Neuron$orig.ident=='CNR1' & Neuron$seurat_clusters == index], slot = "counts")[gene_name,]>0)) +
     (1+sum(GetAssayData(object = Neuron[, Neuron$orig.ident=='30d1' & Neuron$seurat_clusters == index], slot = "counts")[gene_name,]>0)))
## C1 0.517
## C2 0.395
## C3 0.366
## C4 0.371
## C6 0.598
## C7 0.676

## Mif
gene_name = 'Mif'
index = 11
(sum(Other[, Other$orig.ident=='CNR1' & Other$seurat_clusters == index & Other[, Other$orig.ident=='CNR1' & Other$seurat_clusters == index]@assays$SCT@data[gene_name, ] > 0]@assays$SCT@data[gene_name, ]) 
  + sum(Other[, Other$orig.ident=='30d1' & Other$seurat_clusters == index & Other[, Other$orig.ident=='30d1' & Other$seurat_clusters == index]@assays$SCT@data[gene_name, ] > 0]@assays$SCT@data[gene_name, ]))/
  ((1+sum(GetAssayData(object = Other[, Other$orig.ident=='CNR1' & Other$seurat_clusters == index], slot = "counts")[gene_name,]>0)) +
     (1+sum(GetAssayData(object = Other[, Other$orig.ident=='30d1' & Other$seurat_clusters == index], slot = "counts")[gene_name,]>0)))
## A1 0.751
## M 0.707
## A2 0.725
## A3 0.672
## A4 0.779
## A5 0.713

index = 11
(sum(Neuron[, Neuron$orig.ident=='CNR1' & Neuron$seurat_clusters == index & Neuron[, Neuron$orig.ident=='CNR1' & Neuron$seurat_clusters == index]@assays$SCT@data[gene_name, ] > 0]@assays$SCT@data[gene_name, ]) 
  + sum(Neuron[, Neuron$orig.ident=='30d1' & Neuron$seurat_clusters == index & Neuron[, Neuron$orig.ident=='30d1' & Neuron$seurat_clusters == index]@assays$SCT@data[gene_name, ] > 0]@assays$SCT@data[gene_name, ]))/
  ((1+sum(GetAssayData(object = Neuron[, Neuron$orig.ident=='CNR1' & Neuron$seurat_clusters == index], slot = "counts")[gene_name,]>0)) +
     (1+sum(GetAssayData(object = Neuron[, Neuron$orig.ident=='30d1' & Neuron$seurat_clusters == index], slot = "counts")[gene_name,]>0)))
## C1 1.088
## C2 1.156
## C3 1.045
## C4 1.118
## C6 1.124
## C7 1.137

## C1qb
index = 4
gene_name = 'C1qb'
(sum(Other[, Other$orig.ident=='CNR1' & Other$seurat_clusters == index & Other[, Other$orig.ident=='CNR1' & Other$seurat_clusters == index]@assays$SCT@data[gene_name, ] > 0]@assays$SCT@data[gene_name, ]) 
  + sum(Other[, Other$orig.ident=='30d1' & Other$seurat_clusters == index & Other[, Other$orig.ident=='30d1' & Other$seurat_clusters == index]@assays$SCT@data[gene_name, ] > 0]@assays$SCT@data[gene_name, ]))/
  ((1+sum(GetAssayData(object = Other[, Other$orig.ident=='CNR1' & Other$seurat_clusters == index], slot = "counts")[gene_name,]>0)) +
     (1+sum(GetAssayData(object = Other[, Other$orig.ident=='30d1' & Other$seurat_clusters == index], slot = "counts")[gene_name,]>0)))
## A1 0.029
## M 1.065
## A2 0
## A3 0
## A4 0
## A5 0.045

index = 11
gene_name = 'C1qb'
(sum(Neuron[, Neuron$orig.ident=='CNR1' & Neuron$seurat_clusters == index & Neuron[, Neuron$orig.ident=='CNR1' & Neuron$seurat_clusters == index]@assays$SCT@data[gene_name, ] > 0]@assays$SCT@data[gene_name, ]) 
  + sum(Neuron[, Neuron$orig.ident=='30d1' & Neuron$seurat_clusters == index & Neuron[, Neuron$orig.ident=='30d1' & Neuron$seurat_clusters == index]@assays$SCT@data[gene_name, ] > 0]@assays$SCT@data[gene_name, ]))/
  ((1+sum(GetAssayData(object = Neuron[, Neuron$orig.ident=='CNR1' & Neuron$seurat_clusters == index], slot = "counts")[gene_name,]>0)) +
     (1+sum(GetAssayData(object = Neuron[, Neuron$orig.ident=='30d1' & Neuron$seurat_clusters == index], slot = "counts")[gene_name,]>0)))
## C1 0.053
## C2 0.066
## C3 0
## C4 0
## C6 0.036
## C7 0.050

## Spp1
index = 4
gene_name = 'Spp1'
(sum(Other[, Other$orig.ident=='CNR1' & Other$seurat_clusters == index & Other[, Other$orig.ident=='CNR1' & Other$seurat_clusters == index]@assays$SCT@data[gene_name, ] > 0]@assays$SCT@data[gene_name, ]) 
  + sum(Other[, Other$orig.ident=='30d1' & Other$seurat_clusters == index & Other[, Other$orig.ident=='30d1' & Other$seurat_clusters == index]@assays$SCT@data[gene_name, ] > 0]@assays$SCT@data[gene_name, ]))/
  ((1+sum(GetAssayData(object = Other[, Other$orig.ident=='CNR1' & Other$seurat_clusters == index], slot = "counts")[gene_name,]>0)) +
     (1+sum(GetAssayData(object = Other[, Other$orig.ident=='30d1' & Other$seurat_clusters == index], slot = "counts")[gene_name,]>0)))
## A1 0.029
## M 1.065
## A2 
## A3 
## A4 
## A5 

index = 11
gene_name = 'Spp1'
(sum(Neuron[, Neuron$orig.ident=='CNR1' & Neuron$seurat_clusters == index & Neuron[, Neuron$orig.ident=='CNR1' & Neuron$seurat_clusters == index]@assays$SCT@data[gene_name, ] > 0]@assays$SCT@data[gene_name, ]) 
  + sum(Neuron[, Neuron$orig.ident=='30d1' & Neuron$seurat_clusters == index & Neuron[, Neuron$orig.ident=='30d1' & Neuron$seurat_clusters == index]@assays$SCT@data[gene_name, ] > 0]@assays$SCT@data[gene_name, ]))/
  ((1+sum(GetAssayData(object = Neuron[, Neuron$orig.ident=='CNR1' & Neuron$seurat_clusters == index], slot = "counts")[gene_name,]>0)) +
     (1+sum(GetAssayData(object = Neuron[, Neuron$orig.ident=='30d1' & Neuron$seurat_clusters == index], slot = "counts")[gene_name,]>0)))
## C1 0
## C2 0
## C3 0
## C4 0.150
## C6 0
## C7 0

## get cell cordinates and communication score
## Coord_30d1_clean.csv from Figure5A.ipynb
## Coord_CNR1_clean.csv from Figure5A.ipynb
CellCluster_30d1 = read.csv("/media/gulab/GUDR2/PB10_new/Coord_30d1_clean.csv", header = TRUE, sep = "\t")
CellCluster_CNR1 = read.csv("/media/gulab/GUDR2/PB10_new/Coord_CNR1_clean.csv", header = TRUE, sep = "\t")

## get average all communication score to do nomarlization
mean_Sall_CNR1 = mean(CellCluster_CNR1$Sall)
mean_Sall_30d1 = mean(CellCluster_30d1$Sall)
mean_Sall_CNR1
mean_Sall_30d1
## normalization
for (i in 1 : dim(CellCluster_30d1)[[1]]){
  CellCluster_30d1$Sall[i] = CellCluster_30d1$Sall[i] * mean_Sall_CNR1 / mean_Sall_30d1
}
mean(CellCluster_CNR1$Sall)
mean(CellCluster_30d1$Sall)
for (i in 1 : dim(CellCluster_30d1)[[1]]){
  CellCluster_30d1$SApoe[i] = CellCluster_30d1$SApoe[i] * mean_Sall_CNR1 / mean_Sall_30d1
}
for (i in 1 : dim(CellCluster_30d1)[[1]]){
  CellCluster_30d1$SMif[i] = CellCluster_30d1$SMif[i] * mean_Sall_CNR1 / mean_Sall_30d1
}
for (i in 1 : dim(CellCluster_30d1)[[1]]){
  CellCluster_30d1$SSpp1[i] = CellCluster_30d1$SSpp1[i] * mean_Sall_CNR1 / mean_Sall_30d1
}
for (i in 1 : dim(CellCluster_30d1)[[1]]){
  CellCluster_30d1$SC1qb[i] = CellCluster_30d1$SC1qb[i] * mean_Sall_CNR1 / mean_Sall_30d1
}

## generate communication score distribution plot for each gene
svg("/media/gulab/GUDR2/PB10_new/Apoe_density.svg", height = 10, width = 5)
plot(density(CellCluster_CNR1[CellCluster_CNR1$cluster %in% c(18), ]$SApoe), xlim = c(0, 10), ylim = c(0, 0.8), col = 'blue', lwd = 6, xaxt = "n", yaxt = "n")
par(new=TRUE)
plot(density(CellCluster_30d1[CellCluster_30d1$cluster %in% c(18), ]$SApoe), xlim = c(0, 10), ylim = c(0, 0.8), col = 'red', lwd = 6, xaxt = "n", yaxt = "n")
dev.off()

svg("/media/gulab/GUDR2/PB10_new/Mif_density.svg", height = 10, width = 5)
plot(density(CellCluster_CNR1[CellCluster_CNR1$cluster %in% c(1), ]$SMif), xlim = c(0, 5), ylim = c(0, 0.8), col = 'blue', lwd = 6, xaxt = "n", yaxt = "n")
par(new=TRUE)
plot(density(CellCluster_30d1[CellCluster_30d1$cluster %in% c(1), ]$SMif), xlim = c(0, 5), ylim = c(0, 0.8), col = 'red', lwd = 6, xaxt = "n", yaxt = "n")
dev.off()

svg("/media/gulab/GUDR2/PB10_new/Spp1_density.svg", height = 10, width = 5)
plot(density(CellCluster_CNR1[CellCluster_CNR1$cluster == 3, ]$SSpp1), xlim = c(0, 2), ylim = c(0, 3), col = 'blue', lwd = 6, xaxt = "n", yaxt = "n")
par(new=TRUE)
plot(density(CellCluster_30d1[CellCluster_30d1$cluster == 3, ]$SSpp1), xlim = c(0, 2), ylim = c(0, 3), col = 'red', lwd = 6, xaxt = "n", yaxt = "n")
dev.off()

svg("/media/gulab/GUDR2/PB10_new/C1qb_density.svg", height = 10, width = 5)
plot(density(CellCluster_CNR1[CellCluster_CNR1$cluster %in% c(26), ]$SC1qb), xlim = c(0, 2), ylim = c(0, 3), col = 'blue', lwd = 6, xaxt = "n", yaxt = "n")
par(new=TRUE)
plot(density(CellCluster_30d1[CellCluster_30d1$cluster %in% c(26), ]$SC1qb), xlim = c(0, 2), ylim = c(0, 3), col = 'red', lwd = 6, xaxt = "n", yaxt = "n")
dev.off()