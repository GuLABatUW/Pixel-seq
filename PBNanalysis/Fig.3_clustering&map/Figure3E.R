## for Figure3-F

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

## load data and generate only neuron csv file
only_neuron_4 = readRDS("/media/GUDR2/Hector/PB10_new/only_neuron_4.rds")
only_neuron_4_df <- data.frame(only_neuron_4$integrated_snn_res.0.65)
write.csv(only_neuron_4_df, "/media/gulab/GUDR2/PB10_new/only_neuron_4.csv")

## Marker list
Only_Neuron_4.markers <- FindAllMarkers(only_neuron_4, assay = "SCT", only.pos = TRUE)
View(Only_Neuron_4.markers)

## load all barcodes coordinates
All_barcode_coord = read.csv("/media/gulab/GUDR2/PB10_new/Seg_cluster_prg2_pos.txt", header = FALSE, sep = "\t")

## load only neurons coordinates
## Seg_cluster_only_neuron_pos.txt from file_covert.ipynb
Coord_for_Neuron = read.csv("/media/gulab/GUDR2/PB10_new/Seg_cluster_only_neuron_pos.txt", header = FALSE, sep = "\t")

## plot Spatial cell positions 
tiff("/media/gulab/GUDR2/PB10_new/Cell_type_PBN.tiff", height = 15, width = 15, units = 'in', res = 500)
par(mar = rep(0, 4))

## plot all reads as grey
plot(All_barcode_coord[All_barcode_coord$V4 == '30d1_S2',1:2], 
     pch=15,cex=.4, col = "grey",
     ylim = c(4500, 8800), xlim = c(19500, 23800), yaxt = "n", xaxt='n') 
par(new=TRUE)

## plot Calca, Nts cell as green
plot(Coord_for_Neuron[Coord_for_Neuron$V4 == '30d1_S2' & Coord_for_Neuron$V5 == 14,1:2], 
     pch=15,cex=.4, col = "#33A02C",
     ylim = c(4500, 8800), xlim = c(19500, 23800), yaxt = "n", xaxt='n') 
par(new=TRUE)

## plot Resp18/Ctxn2, cell as orange
plot(Coord_for_Neuron[Coord_for_Neuron$V4 == '30d1_S2' & Coord_for_Neuron$V5 == 1,1:2], 
     pch=15,cex=.4, col = "#FC8D62",
     ylim = c(4500, 8800), xlim = c(19500, 23800), yaxt = "n", xaxt='n') 
par(new=TRUE)

## plot Tac1 Cell as brown
plot(Coord_for_Neuron[Coord_for_Neuron$V4 == '30d1_S2'& Coord_for_Neuron$V5 == 13,1:2], 
     pch=15,cex=.4, col = "#B15928",
     ylim = c(4500, 8800), xlim = c(19500, 23800), yaxt = "n", xaxt='n') 
par(new=TRUE)

## plot Calca, Sncg cell as blue
plot(Coord_for_Neuron[Coord_for_Neuron$V4 == '30d1_S2'& Coord_for_Neuron$V5 == 16,1:2], 
     pch=15,cex=.4, col = "#1F78B4",
     ylim = c(4500, 8800), xlim = c(19500, 23800), yaxt = "n", xaxt='n') 
par(new=TRUE)

## plot Sst, Resp18 cell as red
plot(Coord_for_Neuron[Coord_for_Neuron$V4 == '30d1_S2'& Coord_for_Neuron$V5 == 8,1:2], 
     pch=15,cex=.4, col = "#E31A1C",
     ylim = c(4500, 8800), xlim = c(19500, 23800), yaxt = "n", xaxt='n') 
par(new=TRUE)

## Plot Sncg, Uchl1 cell as yellow
plot(Coord_for_Neuron[Coord_for_Neuron$V4 == '30d1_S2'& Coord_for_Neuron$V5 == 5,1:2], 
     pch=15,cex=.4, col ="#FFD92F",
     ylim = c(4500, 8800), xlim = c(19500, 23800), yaxt = "n", xaxt='n') 
dev.off()

