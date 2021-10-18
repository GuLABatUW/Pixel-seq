## for sup Figure12-B

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
PBN_all5 = readRDS("/home/gulab/Downloads/PBNregion.rds")

## set color for UMAP
mycol = c("#66C2A5", # cluster 0, deep cyan for Granule
          "#FC8D62", # cluster 1, orange for Neuron1
          "#8DA0CB", # cluster 2, deep blue for Oligo1
          "#E78AC3", # cluster 3, pink for RBCs
          "#A6D854", # cluster 4, light green for Neuron6
          "#FFD92F", # cluster 5, yellow for Oligo2
          "#E5C494", # cluster 6, wheat for Astro
          "#E31A1C", # cluster 7, red for Neuron3
          "#A6CEE3", # cluster 8, light blue for Ependymal
          "#B3B3B3", # cluster 9, grey for Neuron7
          "#B2DF8A", # cluster 10, light light green for Neuron9
          "#33A02C", # cluster 11, green for VEC
          "#FB9A99", # cluster 12, rose for Neuron8
          "#1F78B4", # cluster 13, blue for Neuron2
          "#FDBF6F", # cluster 14, earth yellow for Neuron5
          "#FF7F00", # cluster 15, pure orange for Neuron4(CGRP)
          "#CAB2D6", # cluster 16, light purple for Immune
          "#6A3D9A", # cluster 17, purple for Oligo3
          "#01665E", # cluster 18, deep green for Satellute glia
          "#B15928", # cluster 19, brown for VLMCS
          "#8DD3C7" # cluster 20, light light cyan for Microglia
)

## generate UMAP plot
png("/media/gulab/GUDR2/PB10_new/All5_UMAP_30d1.png", height = 7, width = 7, units = 'in', res = 300)
DimPlot(PBN_all5[,PBN_all5$orig.ident == '30d1'], reduction = "umap", label = FALSE, cols = mycol, shuffle = TRUE, pt.size = .6)
dev.off()

png("/media/gulab/GUDR2/PB10_new/All5_UMAP_30d2.png", height = 7, width = 7, units = 'in', res = 300)
DimPlot(PBN_all5[,PBN_all5$orig.ident == '30d2'], reduction = "umap", label = FALSE, cols = mycol, shuffle = TRUE, pt.size = .6)
dev.off()

png("/media/gulab/GUDR2/PB10_new/All5_UMAP_CNR1.png", height = 7, width = 7, units = 'in', res = 300)
DimPlot(PBN_all5[,PBN_all5$orig.ident == 'CNR1'], reduction = "umap", label = FALSE, cols = mycol, shuffle = TRUE, pt.size = .6)
dev.off()

png("/media/gulab/GUDR2/PB10_new/All5_UMAP_CNR2.png", height = 7, width = 7, units = 'in', res = 300)
DimPlot(PBN_all5[,PBN_all5$orig.ident == 'CNR2'], reduction = "umap", label = FALSE, cols = mycol, shuffle = TRUE, pt.size = .6)
dev.off()
