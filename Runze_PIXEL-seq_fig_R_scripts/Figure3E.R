## for Figure3-E

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

## load data from 4 section: 30d1, 30d2, CN1, and CN2
## PB_neuron4Sample.rds from Xiaonan
only_neuron_4 = readRDS("/media/gulab/GUDR2/PB10_new/PB_neuron4Sample.rds")

## set color for Neurons UMAP
mycol = c("#A6D854", # cluster 0, light green for Fkbp3
          "#FC8D62", # cluster 1, orange for Resp18, Ctxn2
          "#8DA0CB", # cluster 2, deep blue for Gabrd, Nrep
          "#A6CEE3", # cluster 3, light blue for Pcp4, Pcp2
          "#E5C494", # cluster 4, wheat for Pcp2, Pvalb
          "#FFD92F", # cluster 5, yellow for Sncg, Uchl1, 
          "#CAB2D6", # cluster 6, light purple Timp4, Aldoc
          "#E78AC3", # cluster 7, pink for Fth1, Vps8
          "#E31A1C", # cluster 8, red for Sst, Resp18
          "#B3B3B3", # cluster 9, grey for Pcp4, Spp1
          "#8DD3C7", # cluster 10, light light cyan for Vps8, Mbp
          "#FF7F00", # cluster 11, pure orange for Uqcrb, Fcf1
          "#FB9A99", # cluster 12, rose for Fth1
          "#B15928", # cluster 13, blue for Tac1
          "#33A02C", # cluster 14, green for Calca, Nts
          "#66C2A5", # cluster 15, deep cyan for Bbip1
          "#1F78B4", # cluster 16, brown for Calca, Sncg
          "#6A0DAD" # cluster 17, deep purple Vps8, Tmsb10
          )

## Run UMAP
only_neuron_4 <- RunUMAP(only_neuron_4, reduction = "pca", dims = 1:12)

## Find Neighbors
only_neuron_4 <- FindNeighbors(only_neuron_4, reduction = "pca", dims = 1:12)

## Clustering
only_neuron_4 <- FindClusters(only_neuron_4, resolution = 0.65, verbose = TRUE)

## Plot UMAP
png("/media/gulab/GUDR2/PB10_new/UMAP_neuron_4.png", height = 7, width = 7, units = 'in', res = 300)
DimPlot(only_neuron_4, reduction = "umap", label = FALSE, pt.size = .6, cols = mycol)
dev.off()

only_neuron_4_df <- data.frame(only_neuron_4$integrated_snn_res.0.3)
write.csv(only_neuron_4_df, "/media/gulab/GUDR2/PB10_new/neuron_4_cluster.csv")

# saveRDS(only_neuron_4, "/media/gulab/GUDR2/PB10_new/PB_neuron4Sample.rds")

FeaturePlot(only_neuron_4, features = c("Calca", "Sncg"))
DimPlot(only_neuron_4[,only_neuron_4$orig.ident=='30d1'], reduction = "umap", label = FALSE, pt.size = .6, cols = mycol)
DimPlot(only_neuron_4[,only_neuron_4$orig.ident=='30d2'], reduction = "umap", label = FALSE, pt.size = .6, cols = mycol)
DimPlot(only_neuron_4[,only_neuron_4$orig.ident=='CNR1'], reduction = "umap", label = FALSE, pt.size = .6, cols = mycol)
DimPlot(only_neuron_4[,only_neuron_4$orig.ident=='CNR2'], reduction = "umap", label = FALSE, pt.size = .6, cols = mycol)


## check feature plot
DefaultAssay(only_neuron_4) <- "SCT"
FeaturePlot(only_neuron_4, feature = c("Sst", "Calca", "Tac1", "Nts", "Sncg"))


## get marker gene lst
PBN_4_neuron.markers <- FindAllMarkers(only_neuron_4, assay = "SCT", only.pos = TRUE)
View(PBN_4_neuron.markers)

## spatially show cell distribution
only_neuron_4 <- only_neuron_4[,only_neuron_4$orig.ident=='30d1']
p2 <- SpatialDimPlot(only_neuron_4, label = FALSE, pt.size.factor = 2, cols = mycol, stroke = 1)
p2

