## for Figure3-B

# import library

## for plot and color
library(ggplot2)
library(dplyr)
library(RColorBrewer)
# library(scater)

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

## load PBN data from 5 section 30d1, 30d2, CN1, and CN2
## PBNct_S4.rds from Xiaonan
PBN_4 = readRDS("/media/gulab/GUDR2/PB10_new/PBNct_S4.rds")

## set color for Violin Plot
mycol = c("#66C2A5", # cluster 0, deep cyan for Granule
          "#FC8D62", # cluster 1, orange for Neuron1
          "#8DA0CB", # cluster 2, deep blue for Oligo1
          "#E78AC3", # cluster 3, pink for Oligo2
          "#FDBF6F", # cluster 4, light orange for Astro1
          "#FFD92F", # cluster 5, yellow for Neuron2
          "#E5C494", # cluster 6, light brown for RBCs
          "#E31A1C", # cluster 7, red for Neuron3
          "#A6CEE3", # cluster 8, light blue for Neuron4
          "#B3B3B3", # cluster 9, grey for Neuron7
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
          "#6A3D9A", # cluster 20, light light cyan for Glia
          "#FFFFB3" # cluster 21, light light yellow for Neuron9
)

## Select clusters which have significant marker gene
selected_clusters = c(1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15, 17, 18, 19, 21)

## plot violin plot for each gene
DefaultAssay(PBN_4) <- "integrated"

### Pcp2
png("/media/gulab/GUDR2/PB10_new/Pcp2.png", height = 1, width = 15.78, units = 'in', res = 300)
VlnPlot(PBN_4, features = "Pcp2", cols = mycol[selected_clusters+1], pt.size = 0, idents = selected_clusters, same.y.lims = 1, flip = FALSE) + theme(legend.position = 'none', axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), plot.title = element_blank(), title = element_blank(), axis.line = element_line(size = 1.5), plot.margin = margin(0, 0, 0, 0)) + geom_line(size = 1.5)
dev.off()

### Cryab
png("/media/gulab/GUDR2/PB10_new/Cryab.png", height = 1, width = 15.78, units = 'in', res = 300)
VlnPlot(PBN_4, features = "Cryab", cols = mycol[selected_clusters+1], pt.size = 0, idents = selected_clusters, same.y.lims = 1, flip = FALSE) + theme(legend.position = 'none', axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), plot.title = element_blank(), title = element_blank(), axis.line = element_line(size = 1.5), plot.margin = margin(0, 0, 0, 0)) + geom_line(size = 1.5)
dev.off()

### Ptgds
png("/media/gulab/GUDR2/PB10_new/Ptgds.png", height = 1, width = 15.78, units = 'in', res = 300)
VlnPlot(PBN_4, features = "Ptgds", cols = mycol[selected_clusters+1], pt.size = 0, idents = selected_clusters, same.y.lims = 1, flip = FALSE) + theme(legend.position = 'none', axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), plot.title = element_blank(), title = element_blank(), axis.line = element_line(size = 1.5), plot.margin = margin(0, 0, 0, 0)) + geom_line(size = 1.5)
dev.off()

### Apoe
png("/media/gulab/GUDR2/PB10_new/Apoe.png", height = 1, width = 15.78, units = 'in', res = 300)
VlnPlot(PBN_4, features = "Apoe", cols = mycol[selected_clusters+1], pt.size = 0, idents = selected_clusters, same.y.lims = 1, flip = FALSE) + theme(legend.position = 'none', axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), plot.title = element_blank(), title = element_blank(), axis.line = element_line(size = 1.5), plot.margin = margin(0, 0, 0, 0)) + geom_line(size = 1.5)
dev.off()

### Resp18
png("/media/gulab/GUDR2/PB10_new/Resp18.png", height = 1, width = 15.78, units = 'in', res = 300)
VlnPlot(PBN_4, features = "Resp18", cols = mycol[selected_clusters+1], pt.size = 0, idents = selected_clusters, same.y.lims = 1, flip = FALSE) + theme(legend.position = 'none', axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), plot.title = element_blank(), title = element_blank(), axis.line = element_line(size = 1.5), plot.margin = margin(0, 0, 0, 0)) + geom_line(size = 1.5)
dev.off()

### Hba-a2
png("/media/gulab/GUDR2/PB10_new/Hba-a2.png", height = 1, width = 15.78, units = 'in', res = 300)
VlnPlot(PBN_4, features = "Hba-a2", cols = mycol[selected_clusters+1], pt.size = 0, idents = selected_clusters, same.y.lims = 1, flip = FALSE) + theme(legend.position = 'none', axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), plot.title = element_blank(), title = element_blank(), axis.line = element_line(size = 1.5), plot.margin = margin(0, 0, 0, 0)) + geom_line(size = 1.5)
dev.off()

### Sncg
png("/media/gulab/GUDR2/PB10_new/Sncg.png", height = 1, width = 15.78, units = 'in', res = 300)
VlnPlot(PBN_4, features = "Sncg", cols = mycol[selected_clusters+1], pt.size = 0, idents = selected_clusters, same.y.lims = 1, flip = FALSE) + theme(legend.position = 'none', axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), plot.title = element_blank(), title = element_blank(), axis.line = element_line(size = 1.5), plot.margin = margin(0, 0, 0, 0)) + geom_line(size = 1.5)
dev.off()

### Pcp4
png("/media/gulab/GUDR2/PB10_new/Pcp4.png", height = 1, width = 15.78, units = 'in', res = 300)
VlnPlot(PBN_4, features = "Pcp4", cols = mycol[selected_clusters+1], pt.size = 0, idents = selected_clusters, same.y.lims = 1, flip = FALSE) + theme(legend.position = 'none', axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), plot.title = element_blank(), title = element_blank(), axis.line = element_line(size = 1.5), plot.margin = margin(0, 0, 0, 0)) + geom_line(size = 1.5)
dev.off()

### B2m
png("/media/gulab/GUDR2/PB10_new/B2m.png", height = 1, width = 15.78, units = 'in', res = 300)
VlnPlot(PBN_4, features = "B2m", cols = mycol[selected_clusters+1], pt.size = 0, idents = selected_clusters, same.y.lims = 1, flip = FALSE) + theme(legend.position = 'none', axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), plot.title = element_blank(), title = element_blank(), axis.line = element_line(size = 1.5), plot.margin = margin(0, 0, 0, 0)) + geom_line(size = 1.5)
dev.off()

### Ly6c1
png("/media/gulab/GUDR2/PB10_new/Ly6c1.png", height = 1, width = 15.78, units = 'in', res = 300)
VlnPlot(PBN_4, features = "Ly6c1", cols = mycol[selected_clusters+1], pt.size = 0, idents = selected_clusters, same.y.lims = 1, flip = FALSE) + theme(legend.position = 'none', axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), plot.title = element_blank(), title = element_blank(), axis.line = element_line(size = 1.5), plot.margin = margin(0, 0, 0, 0)) + geom_line(size = 1.5)
dev.off()

### Ramp1
png("/media/gulab/GUDR2/PB10_new/Ramp1.png", height = 1, width = 15.78, units = 'in', res = 300)
VlnPlot(PBN_4, features = "Ramp1", cols = mycol[selected_clusters+1], pt.size = 0, idents = selected_clusters, same.y.lims = 1, flip = FALSE) + theme(legend.position = 'none', axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), plot.title = element_blank(), title = element_blank(), axis.line = element_line(size = 1.5), plot.margin = margin(0, 0, 0, 0)) + geom_line(size = 1.5)
dev.off()

### Sst
png("/media/gulab/GUDR2/PB10_new/Sst.png", height = 1, width = 15.78, units = 'in', res = 300)
VlnPlot(PBN_4, features = "Sst", cols = mycol[selected_clusters+1], pt.size = 0, idents = selected_clusters, same.y.lims = 1, flip = FALSE) + theme(legend.position = 'none', axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), plot.title = element_blank(), title = element_blank(), axis.line = element_line(size = 1.5), plot.margin = margin(0, 0, 0, 0)) + geom_line(size = 1.5)
dev.off()

### Calca
png("/media/gulab/GUDR2/PB10_new/Calca.png", height = 1, width = 15.78, units = 'in', res = 300)
VlnPlot(PBN_4, features = "Calca", cols = mycol[selected_clusters+1], pt.size = 0, idents = selected_clusters, same.y.lims = 1, flip = FALSE) + theme(legend.position = 'none', axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), plot.title = element_blank(), title = element_blank(), axis.line = element_line(size = 1.5), plot.margin = margin(0, 0, 0, 0)) + geom_line(size = 1.5)
dev.off()

### Tac1
png("/media/gulab/GUDR2/PB10_new/Tac1.png", height = 1, width = 15.78, units = 'in', res = 300)
VlnPlot(PBN_4, features = "Tac1", cols = mycol[selected_clusters+1], pt.size = 0, idents = selected_clusters, same.y.lims = 1, flip = FALSE) + theme(legend.position = 'none', axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), plot.title = element_blank(), title = element_blank(), axis.line = element_line(size = 1.5), plot.margin = margin(0, 0, 0, 0)) + geom_line(size = 1.5)
dev.off()

### Rarres2
png("/media/gulab/GUDR2/PB10_new/Rarres2.png", height = 1, width = 15.78, units = 'in', res = 300)
VlnPlot(PBN_4, features = "Rarres2", cols = mycol[selected_clusters+1], pt.size = 0, idents = selected_clusters, same.y.lims = 1, flip = FALSE) + theme(legend.position = 'none', axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), plot.title = element_blank(), title = element_blank(), axis.line = element_line(size = 1.5), plot.margin = margin(0, 0, 0, 0)) + geom_line(size = 1.5)
dev.off()

### Cyb5a
png("/media/gulab/GUDR2/PB10_new/Cyb5a.png", height = 1, width = 15.78, units = 'in', res = 300)
VlnPlot(PBN_4, features = "Cyb5a", cols = mycol[selected_clusters+1], pt.size = 0, idents = selected_clusters, same.y.lims = 1, flip = FALSE) + theme(legend.position = 'none', axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), plot.title = element_blank(), title = element_blank(), axis.line = element_line(size = 1.5), plot.margin = margin(0, 0, 0, 0)) + geom_line(size = 1.5)
dev.off()

### Tsc22d4
png("/media/gulab/GUDR2/PB10_new/Tsc22d4.png", height = 1, width = 15.78, units = 'in', res = 300)
VlnPlot(PBN_4, features = "Tsc22d4", cols = mycol[selected_clusters+1], pt.size = 0, idents = selected_clusters, same.y.lims = 1, flip = FALSE) + theme(legend.position = 'none', axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), plot.title = element_blank(), title = element_blank(), axis.line = element_line(size = 1.5), plot.margin = margin(0, 0, 0, 0)) + geom_line(size = 1.5)
dev.off()

### Prdx1
# png("/media/gulab/GUDR2/PB10_new/Prdx1.png", height = 1, width = 15.78, units = 'in', res = 300)
# VlnPlot(PBN_4, features = "Prdx1", cols = mycol[selected_clusters+1], pt.size = 0, idents = selected_clusters, same.y.lims = 1, flip = FALSE) + theme(legend.position = 'none', axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), plot.title = element_blank(), title = element_blank(), axis.line = element_line(size = 1.5), plot.margin = margin(0, 0, 0, 0)) + geom_line(size = 1.5)
# dev.off()

### Bbip1
png("/media/gulab/GUDR2/PB10_new/Bbip1.png", height = 1, width = 15.78, units = 'in', res = 300)
VlnPlot(PBN_4, features = "Bbip1", cols = mycol[selected_clusters+1], pt.size = 0, idents = selected_clusters, same.y.lims = 1, flip = FALSE) + theme(legend.position = 'none', axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), plot.title = element_blank(), title = element_blank(), axis.line = element_line(size = 1.5), plot.margin = margin(0, 0, 0, 0)) + geom_line(size = 1.5)
dev.off()

## get counts of all cell types
for (i in 1 : length(selected_clusters)){
  print(dim(PBN_4[,PBN_4$integrated_snn_res.0.3 == selected_clusters[[i]]]))
}

## 3000 6967
## 3000 6510
## 3000 4965
## 3000 4348
## 3000 3820
## 3000 3764
## 3000 2772
## 3000 2190
## 3000 2157
## 3000 2020
## 3000 2004
## 3000 1596
## 3000 1575
## 3000 1497
## 3000 1237
## 3000 1072
## 3000  314
