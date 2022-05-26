## for Figure3-G-2

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

## load reads map and gene coord file
## PB10all_cluster_pos.txt from Xiaonan
## PB10_Calca_in_region.final from file_convert.ipynb
## PB10_Tac1_in_region.final from file_convert.ipynb
Pos_PB10_all = read.csv("/media/gulab/GUDR2/PB10_new/PB10all_cluster_pos.txt", header = FALSE, sep = "\t")
Calca_all = read.csv("/media/gulab/GUDR2/PB10_new/PB10_Calca_in_region.final", header = FALSE, sep = "\t")
Tac1_all = read.csv("/media/gulab/GUDR2/PB10_new/PB10_Tac1_in_region.final", header = FALSE, sep = "\t")

## figure 3F-2A for gene distribution on PBN region -5.2mm
tiff("/media/gulab/Hector/PB10_new/Gene_A.tiff", height = 7, width = 7, units = 'in', res = 500)
par(mar = rep(0, 4))

## display reads map as grey
plot(Pos_PB10_all[,1:2], 
     pch=15,cex=.5, col = 'grey',
     ylim = c(6200, 7500), xlim = c(21800, 23000), yaxt = "n", xaxt='n') 
par(new=TRUE)

## display Calca gene as green point
plot(Calca_all[,3:4], 
     pch=19,cex=1.5, col = '#33A02C',
     ylim = c(6200, 7500), xlim = c(21800, 23000), yaxt = "n", xaxt='n') 
par(new=TRUE)

## display Tac1 gene as blue point 
plot(Tac1_all[,3:4], 
     pch=19,cex=1.5, col = '#B15928',
     ylim = c(6200, 7500), xlim = c(21800, 23000), yaxt = "n", xaxt='n') 
dev.off()

## figure 3F-2B for gene distribution on PBN region -5.35mm
tiff("/media/gulab/Hector/PB10_new/Gene_B.tiff", height = 7, width = 7, units = 'in', res = 500)
par(mar = rep(0, 4))
plot(Pos_PB10_all[,1:2], 
     pch=15,cex=.5, col = 'grey',
     ylim = c(6900, 8300), xlim = c(36300, 37500), yaxt = "n", xaxt='n') 
par(new=TRUE)
plot(Calca_all[,3:4], 
     pch=19,cex=1.5, col = '#33A02C',
     ylim = c(6900, 8300), xlim = c(36300, 37500), yaxt = "n", xaxt='n') 
par(new=TRUE)
plot(Tac1_all[,3:4], 
     pch=19,cex=1.5, col = '#B15928',
     ylim = c(6900, 8300), xlim = c(36300, 37500), yaxt = "n", xaxt='n') 
dev.off()

## figure 3F-2C for gene distribution on PBN region -5.50mm
tiff("/media/gulab/Hector/PB10_new/Gene_C.tiff", height = 7, width = 7, units = 'in', res = 500)
par(mar = rep(0, 4))
plot(Pos_PB10_all[,1:2], 
     pch=15,cex=.5, col = 'grey',
     ylim = c(6500, 7800), xlim = c(49200, 50500), yaxt = "n", xaxt='n') 
par(new=TRUE)
plot(Calca_all[,3:4], 
     pch=19,cex=1.5, col = '#33A02C',
     ylim = c(6500, 7800), xlim = c(49200, 50500), yaxt = "n", xaxt='n') 
par(new=TRUE)
plot(Tac1_all[,3:4], 
     pch=19,cex=1.5, col = '#B15928',
     ylim = c(6500, 7800), xlim = c(49200, 50500), yaxt = "n", xaxt='n') 
dev.off()
