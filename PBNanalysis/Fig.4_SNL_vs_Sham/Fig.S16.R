# For Sup Figure 13

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

library(viridis)

## load PBN data from 4 section 30d1, 30d2, CN1, and CN2
## PBNregion.rds from Xiaonan
PBN_all5 = readRDS("/media/gulab/GUDR2/PB10_new/PBNregion.rds")

## set color for Neurons UMAP
mycol = c("#A6D854", # cluster 0, light green for Fkbp3
          "#FC8D62", # cluster 1, orange for Resp18, Ctxn2
          "#8DA0CB", # cluster 2, deep blue for Gabrd, Nrep
          "#A6CEE3", # cluster 3, light blue for Pcp4, Pcp8
          "#E5C494", # cluster 4, wheat for Pcp2, Pvalb
          "#FFD92F", # cluster 5, yellow for Sncg, Uchl1, 
          "#CAB2D6", # cluster 6, light purple Timp4, Aldoc
          "#E78AC3", # cluster 7, pink for Fth1, Vps8
          "#E31A1C", # cluster 8, red for Sst, Resp18
          "#01665E", # cluster 9, deep green for Spp1
          "#8DD3C7", # cluster 10, light light cyan for Vps8, Mbp
          "#FF7F00", # cluster 11, pure orange for Uqcrb, Fcf1
          "#FB9A99", # cluster 12, rose for Fth1
          "#B15928", # cluster 13, blue for Tac1
          "#33A02C", # cluster 14, green for Calca, Nts
          "#66C2A5", # cluster 15, deep cyan for Bbip1
          "#1F78B4", # cluster 16, brown for Calca, Sncg
          "#6A0DAD" # cluster 17, deep purple Tmsb10
)


## load all barcodes coordinates
## Seg_cluster_prg2_pos_updated.txt from file_convert.ipynb
All_barcode_coord = read.csv("/media/gulab/GUDR2/PB10_new/Seg_cluster_prg2_pos_updated.txt", header = FALSE, sep = "\t")

## load only neurons coordinates
## Seg_cluster_only_neuron_pos.txt from file_convert.ipynb
Coord_for_Neuron = read.csv("/media/gulab/GUDR2/PB10_new/Seg_cluster_only_neuron_pos.txt", header = FALSE, sep = "\t")

## range of 30d1_S2
## 16500
## 26000
## 3000
## 11200

## plot Spatial cell positionsfor 30d1
tiff("/media/gulab/GUDR2/PB10_new/Cell_type_All_30d1.tiff", height = 1500, width = 1500)
par(mar = rep(0, 4))
plot(All_barcode_coord[All_barcode_coord$V4 == '30d1_S2',1:2], 
     pch=15,cex=.4, col = "grey",
     yaxt = "n", xaxt='n') 
par(new=TRUE)
for (i in 1:18){
  par(new=TRUE) 
  plot(Coord_for_Neuron[Coord_for_Neuron$V4 == '30d1_S2' & Coord_for_Neuron$V5 == i-1,1:2], 
       pch=15,cex=.25, col = mycol[i], 
       ylim = c(3000, 11200), xlim = c(16500, 26000), yaxt = "n", xaxt='n') 
}
dev.off()

## plot Spatial cell positionsfor CNR1
tiff("/media/gulab/GUDR2/PB10_new/Cell_type_All_CNR1.tiff", height = 1500, width = 1500)
par(mar = rep(0, 4))
plot(All_barcode_coord[All_barcode_coord$V4 == 'CNR1_S2',1:2], 
     pch=15,cex=.6, col = "grey",
     yaxt = "n", xaxt='n') 
par(new=TRUE)
for (i in 1:17){
  par(new=TRUE) 
  plot(Coord_for_Neuron[Coord_for_Neuron$V4 == 'CNR1_S2' & Coord_for_Neuron$V5 == i-1,1:2], 
       pch=15,cex=.25, col = mycol[i], 
       ylim = c(2000, 11200), xlim = c(29000, 38500), yaxt = "n", xaxt='n') 
}
dev.off()

## set color for UMAP for All Cells
mycol_all = c("#66C2A5", # cluster 0, deep cyan for Granule
              "#FC8D62", # cluster 1, orange for Neuron1
              "#8DA0CB", # cluster 2, deep blue for Oligo1
              "#E78AC3", # cluster 3, pink for Oligo2
              "#FDBF6F", # cluster 4, light orange for Astro
              "#FFD92F", # cluster 5, yellow for Neuron2
              "#E5C494", # cluster 6, light brown for RBCs
              "#E31A1C", # cluster 7, red for Neuron3
              "#A6CEE3", # cluster 8, light blue for Neuron4
              "#B3B3B3", # cluster 9, grey for VLMCs
              "#B2DF8A", # cluster 10, light light green for Neuron5
              "#33A02C", # cluster 11, green for Microglia
              "#FB9A99", # cluster 12, rose for EC
              "#1F78B4", # cluster 13, blue for Astro2
              "#A6D854", # cluster 14, light green for Neuron6
              "#FF7F00", # cluster 15, pure orange for Neuron7(CGRP)
              "#CAB2D6", # cluster 16, light purple for Neuron8
              "#8DD3C7", # cluster 17, light light cyan for Ependymal
              "#01665E", # cluster 18, deep green for Astro3
              "#B15928", # cluster 19, brown for Astro4
              "#6A3D9A", # cluster 20, light light cyan for unknown
              "#FFFFB3" # cluster 21, light light yellow for Neuron9
)

## save PBN_4 cluster infor for updating cluster info in All_barcode_coord
PBN_4_df <- data.frame(PBN_4$integrated_snn_res.0.3)
write.csv(PBN_4_df, "/media/gulab/Hector/PB10_new/PBN_4_cluster.csv")

## plot Spatial all cell positions for 30d1
tiff("/media/gulab/GUDR2/PB10_new/Cell_type_All_all_30d1.tiff", height = 1500, width = 1500)
par(mar = rep(0, 4))
for (i in 1:22){
  plot(All_barcode_coord[All_barcode_coord$V4 == '30d1_S2' & All_barcode_coord$V5 == i-1,1:2], 
       pch=15,cex=.25, col = mycol_all[i], 
       ylim = c(3000, 11200), xlim = c(16500, 26000), yaxt = "n", xaxt='n') 
  par(new=TRUE)
}
dev.off()

## plot Spatial all cell positions 
tiff("/media/gulab/GUDR2/PB10_new/Cell_type_All_all_CNR1.tiff", height = 1500, width = 1500)
par(mar = rep(0, 4))
for (i in 1:22){
  plot(All_barcode_coord[All_barcode_coord$V4 == 'CNR1_S2' & All_barcode_coord$V5 == i-1,1:2], 
       pch=15,cex=.25, col = mycol_all[i], 
       ylim = c(2000, 11200), xlim = c(29000, 38500), yaxt = "n", xaxt='n') 
  par(new=TRUE)
}
dev.off()
