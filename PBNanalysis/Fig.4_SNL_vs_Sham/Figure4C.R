## for Figure4-C

# import library

## for plot and color
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(scater)
library(vioplot)

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
## PBS2_30d1uniq.coord from Xiaonan
## PBS2_CN1uniq.coord from Xiaonan
## Seg_cluster_prg2_pos.txt from Xiaonan
reads_table_30d1 = read.csv("/media/gulab/GUDR2/PB10_new/PBS2_30d1uniq.coord", sep = ',', header = TRUE)
reads_table_CN1 = read.csv("/media/gulab/GUDR2/PB10_new/PBS2_CN1uniq.coord", sep = ',', header = TRUE)
All_barcode_coord = read.csv("/media/gulab/GUDR2/PB10_new/Seg_cluster_prg2_pos.txt", header = FALSE, sep = "\t")

min(reads_table_30d1[,2])
# 16500.9
max(reads_table_30d1[,2])
# 25999.94
min(reads_table_30d1[,3])
# 3000.043
max(reads_table_30d1[,3])
# 11165.5

min(reads_table_CN1[,2])
# 29000.11
max(reads_table_CN1[,2])
# 38499.92
min(reads_table_CN1[,3])
# 2000.045
max(reads_table_CN1[,3])
# 11165.5

Calca_all = read.csv("/media/gulab/GUDR2/PB10_new/PBS2_calca.final", header = FALSE, sep = "\t")
Calca_all_30d1 = Calca_all[Calca_all$V3 >= 16500 & Calca_all$V3 <= 26000 & Calca_all$V4 >= 3000 & Calca_all$V4 <= 12000, ]
Calca_all_CNR1 = Calca_all[Calca_all$V3 >= 29000 & Calca_all$V3 <= 38500 & Calca_all$V4 >= 2000 & Calca_all$V4 <= 12000, ]

Penk_all = read.csv("/media/gulab/GUDR2/PB10_new/PBS2_penk.final", header = FALSE, sep = "\t")
Penk_all_30d1 = Penk_all[Penk_all$V3 >= 16500 & Penk_all$V3 <= 26000 & Penk_all$V4 >= 3000 & Penk_all$V4 <= 12000, ]
Penk_all_CNR1 = Penk_all[Penk_all$V3 >= 29000 & Penk_all$V3 <= 38500 & Penk_all$V4 >= 2000 & Penk_all$V4 <= 12000, ]

Scg2_all = read.csv("/media/gulab/GUDR2/PB10_new/PBS2_scg2.final", header = FALSE, sep = "\t")
Scg2_all_30d1 = Scg2_all[Scg2_all$V3 >= 16500 & Scg2_all$V3 <= 26000 & Scg2_all$V4 >= 3000 & Scg2_all$V4 <= 12000, ]
Scg2_all_CNR1 = Scg2_all[Scg2_all$V3 >= 29000 & Scg2_all$V3 <= 38500 & Scg2_all$V4 >= 2000 & Scg2_all$V4 <= 12000, ]

Cdc42_all = read.csv("/media/gulab/GUDR2/PB10_new/PBS2_Cdc42.final", header = FALSE, sep = "\t")
Cdc42_all_30d1 = Cdc42_all[Cdc42_all$V3 >= 16500 & Cdc42_all$V3 <= 26000 & Cdc42_all$V4 >= 3000 & Cdc42_all$V4 <= 12000, ]
Cdc42_all_CNR1 = Cdc42_all[Cdc42_all$V3 >= 29000 & Cdc42_all$V3 <= 38500 & Cdc42_all$V4 >= 2000 & Cdc42_all$V4 <= 12000, ]

Stmn2_all = read.csv("/media/gulab/GUDR2/PB10_new/PBS2_Stmn2.final", header = FALSE, sep = "\t")
Stmn2_all_30d1 = Stmn2_all[Stmn2_all$V3 >= 16500 & Stmn2_all$V3 <= 26000 & Stmn2_all$V4 >= 3000 & Stmn2_all$V4 <= 12000, ]
Stmn2_all_CNR1 = Stmn2_all[Stmn2_all$V3 >= 29000 & Stmn2_all$V3 <= 38500 & Stmn2_all$V4 >= 2000 & Stmn2_all$V4 <= 12000, ]

Cck_all = read.csv("/media/gulab/GUDR2/PB10_new/PBS2_Cck.final", header = FALSE, sep = "\t")
Cck_all_30d1 = Cck_all[Cck_all$V3 >= 16500 & Cck_all$V3 <= 26000 & Cck_all$V4 >= 3000 & Cck_all$V4 <= 12000, ]
Cck_all_CNR1 = Cck_all[Cck_all$V3 >= 29000 & Cck_all$V3 <= 38500 & Cck_all$V4 >= 2000 & Cck_all$V4 <= 12000, ]

Selenom_all = read.csv("/media/gulab/GUDR2/PB10_new/PBS2_Selenom.final", header = FALSE, sep = "\t")
Selenom_all_30d1 = Selenom_all[Selenom_all$V3 >= 16500 & Selenom_all$V3 <= 26000 & Selenom_all$V4 >= 3000 & Selenom_all$V4 <= 12000, ]
Selenom_all_CNR1 = Selenom_all[Selenom_all$V3 >= 29000 & Selenom_all$V3 <= 38500 & Selenom_all$V4 >= 2000 & Selenom_all$V4 <= 12000, ]

Ctxn2_all = read.csv("/media/gulab/GUDR2/PB10_new/PBS2_Ctxn2.final", header = FALSE, sep = "\t")
Ctxn2_all_30d1 = Ctxn2_all[Ctxn2_all$V3 >= 16500 & Ctxn2_all$V3 <= 26000 & Ctxn2_all$V4 >= 3000 & Ctxn2_all$V4 <= 12000, ]
Ctxn2_all_CNR1 = Ctxn2_all[Ctxn2_all$V3 >= 29000 & Ctxn2_all$V3 <= 38500 & Ctxn2_all$V4 >= 2000 & Ctxn2_all$V4 <= 12000, ]

Pin1_all = read.csv("/media/gulab/GUDR2/PB10_new/PBS2_Pin1.final", header = FALSE, sep = "\t")
Pin1_all_30d1 = Pin1_all[Pin1_all$V3 >= 16500 & Pin1_all$V3 <= 26000 & Pin1_all$V4 >= 3000 & Pin1_all$V4 <= 12000, ]
Pin1_all_CNR1 = Pin1_all[Pin1_all$V3 >= 29000 & Pin1_all$V3 <= 38500 & Pin1_all$V4 >= 2000 & Pin1_all$V4 <= 12000, ]

Crh_all = read.csv("/media/gulab/GUDR2/PB10_new/PBS2_Crh.final", header = FALSE, sep = "\t")
Crh_all_30d1 = Crh_all[Crh_all$V3 >= 16500 & Crh_all$V3 <= 26000 & Crh_all$V4 >= 3000 & Crh_all$V4 <= 12000, ]
Crh_all_CNR1 = Crh_all[Crh_all$V3 >= 29000 & Crh_all$V3 <= 38500 & Crh_all$V4 >= 2000 & Crh_all$V4 <= 12000, ]

Adcyap1_all = read.csv("/media/gulab/GUDR2/PB10_new/PBS2_Adcyap1.final", header = FALSE, sep = "\t")
Adcyap1_all_30d1 = Adcyap1_all[Adcyap1_all$V3 >= 16500 & Adcyap1_all$V3 <= 26000 & Adcyap1_all$V4 >= 3000 & Adcyap1_all$V4 <= 12000, ]
Adcyap1_all_CNR1 = Adcyap1_all[Adcyap1_all$V3 >= 29000 & Adcyap1_all$V3 <= 38500 & Adcyap1_all$V4 >= 2000 & Adcyap1_all$V4 <= 12000, ]



## plot 30d1 Calca
tiff("/media/gulab/GUDR2/PB10_new/30d1_Calca.tiff", height = 1500, width = 1500)
par(mar = rep(0, 4))
## plot all reads as grey
plot(reads_table_30d1[,2:3],pch=15,cex=.8,
     col=rgb(0.6,0.6,0.6,reads_table_30d1$UMIcounts/500),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
par(new=TRUE)
plot(Calca_all[,3:4], 
     pch=19,cex=3, col = rgb(0,1,0,0.4/1.0),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
dev.off()

## plot CNR1 Calca
tiff("/media/gulab/GUDR2/PB10_new/CNR1_Calca.tiff", height = 1500, width = 1500)
par(mar = rep(0, 4))
## plot all reads as grey
plot(reads_table_CN1[,2:3],pch=15,cex=.4,
     col=rgb(0.05,0.05,0.05,reads_table_CN1$UMIcounts/650),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n')
par(new=TRUE)
plot(Calca_all[,3:4], 
     pch=19,cex=3, col = rgb(0,1,0,0.4/1.4),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n') 
dev.off()

## plot 30d1 Penk
tiff("/media/gulab/GUDR2/PB10_new/30d1_Penk.tiff", height = 1500, width = 1500)
par(mar = rep(0, 4))
## plot all reads as grey
plot(reads_table_30d1[,2:3],pch=15,cex=.8,
     col=rgb(0.6,0.6,0.6,reads_table_30d1$UMIcounts/500),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
par(new=TRUE)
plot(Penk_all[,3:4], 
     pch=19,cex=3, col = rgb(1,0,0,0.4/1.0),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
dev.off()

## plot CNR1 Penk
tiff("/media/gulab/GUDR2/PB10_new/CNR1_Penk.tiff", height = 1500, width = 1500)
par(mar = rep(0, 4))
## plot all reads as grey
plot(reads_table_CN1[,2:3],pch=15,cex=.4,
     col=rgb(0.05,0.05,0.05,reads_table_CN1$UMIcounts/650),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n')
par(new=TRUE)
plot(Penk_all[,3:4], 
     pch=19,cex=3, col = rgb(1,0,0,0.4/1.4),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n') 
dev.off()

## plot 30d1 Scg2
tiff("/media/gulab/GUDR2/PB10_new/30d1_Scg2.tiff", height = 1500, width = 1500)
par(mar = rep(0, 4))
## plot all reads as grey
plot(reads_table_30d1[,2:3],pch=15,cex=.8,
     col=rgb(0.6,0.6,0.6,reads_table_30d1$UMIcounts/500),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
par(new=TRUE)
## plot Scg2 reads as blue
plot(Scg2_all[,3:4], 
     pch=19,cex=3, col = rgb(0,0,1,0.4/1.0),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
dev.off()

## plot CNR1 Scg2
tiff("/media/gulab/GUDR2/PB10_new/CNR1_Scg2.tiff", height = 1500, width = 1500)
par(mar = rep(0, 4))
## plot all reads as grey
plot(reads_table_CN1[,2:3],pch=15,cex=.4,
     col=rgb(0.05,0.05,0.05,reads_table_CN1$UMIcounts/650),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n')
par(new=TRUE)
## plot Scg2 reads as blue
plot(Scg2_all[,3:4], 
     pch=19,cex=3, col = rgb(0,0,1,0.4/1.4),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n') 
dev.off()

## plot 30d1 Cck
tiff("/media/gulab/GUDR2/PB10_new/30d1_Cck.tiff", height = 1500, width = 1500)
par(mar = rep(0, 4))
## plot all reads as grey
plot(reads_table_30d1[,2:3],pch=15,cex=.8,
     col=rgb(0.6,0.6,0.6,reads_table_30d1$UMIcounts/500),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
par(new=TRUE)
## plot Cck reads as magenta
plot(Cck_all[,3:4], 
     pch=19,cex=3, col = rgb(1,0,1,0.4/1.0),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
dev.off()

## plot CNR1 Cck
tiff("/media/gulab/GUDR2/PB10_new/CNR1_Cck.tiff", height = 1500, width = 1500)
par(mar = rep(0, 4))
## plot all reads as grey
plot(reads_table_CN1[,2:3],pch=15,cex=.4,
     col=rgb(0.05,0.05,0.05,reads_table_CN1$UMIcounts/650),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n')
par(new=TRUE)
## plot Cck reads as magenta
plot(Cck_all[,3:4], 
     pch=19,cex=3, col = rgb(1,0,1,0.4/1.4),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n') 
dev.off()

## plot 30d1 All
tiff("/media/gulab/GUDR2/PB10_new/30d1_All.tiff", height = 1500, width = 1500)
par(mar = rep(0, 4))
## plot all reads as grey
plot(reads_table_30d1[,2:3],pch=15,cex=.8,
     col=rgb(0.6,0.6,0.6,reads_table_30d1$UMIcounts/500),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
par(new=TRUE)
plot(Calca_all[,3:4], 
     pch=19,cex=3, col = rgb(0,1,0,0.4/1.4),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
par(new=TRUE)
plot(Penk_all[,3:4], 
     pch=19,cex=3, col = rgb(1,0,0,0.4/1.4),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
par(new=TRUE)
## plot Scg2 reads as blue
plot(Scg2_all[,3:4], 
     pch=19,cex=3, col = rgb(0,0,1,0.4/1.4),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
par(new=TRUE)
## plot Cck reads as magenta
plot(Cck_all[,3:4], 
     pch=19,cex=3, col = rgb(1,0,1,0.4/1.4),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
dev.off()

## plot CNR1 All
tiff("/media/gulab/GUDR2/PB10_new/CNR1_All.tiff", height = 1500, width = 1500)
par(mar = rep(0, 4))
## plot all reads as grey
plot(reads_table_CN1[,2:3],pch=15,cex=.4,
     col=rgb(0.05,0.05,0.05,reads_table_CN1$UMIcounts/650),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n')
par(new=TRUE)
plot(Calca_all[,3:4], 
     pch=19,cex=3, col = rgb(0,1,0,0.4/1.4),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n') 
par(new=TRUE)
plot(Penk_all[,3:4], 
     pch=19,cex=3, col = rgb(1,0,0,0.4/1.4),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n') 
par(new=TRUE)
## plot Scg2 reads as blue
plot(Scg2_all[,3:4], 
     pch=19,cex=3, col = rgb(0,0,1,0.4/1.4),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n') 
par(new=TRUE)
## plot Cck reads as magenta
plot(Cck_all[,3:4], 
     pch=19,cex=3, col = rgb(1,0,1,0.4/1.4),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n') 
dev.off()

## plot 30d1 Calca all
tiff("/media/gulab/GUDR2/PB10_new/30d1_Calca_all.tiff", height = 1200, width = 1600)
par(mar = rep(0, 4))
## plot all reads as grey
plot(reads_table_30d1[,2:3],pch=19,cex=.1,
     col=rgb(0.1,0.1,0.1,reads_table_30d1$UMIcounts/1000), 
     xlim = c(16500, 26000), ylim = c(3000, 11200), yaxt = "n", xaxt='n')
par(new=TRUE)
## plot Adcyap1 reads as blue
plot(Calca_all_30d1[,3:4], 
     pch=19,cex=1.5, col = rgb(0,1,0,0.2/2.0), 
     xlim = c(16500, 26000), ylim = c(3000, 11200), yaxt = "n", xaxt='n') 
dev.off()

## plot CNR1 Calca all
tiff("/media/gulab/GUDR2/PB10_new/CNR1_Calca_all.tiff", height = 1200, width = 1600)
par(mar = rep(0, 4))
## plot all reads as grey
plot(reads_table_CN1[,2:3],pch=19,cex=.1,
     col=rgb(0.1,0.1,0.1,reads_table_CN1$UMIcounts/1400), 
     xlim = c(29000, 38500), ylim = c(2000, 11200), yaxt = "n", xaxt='n') 
par(new=TRUE)
## plot Adcyap1 reads as blue
plot(Calca_all_CNR1[,3:4], 
     pch=19,cex=1.5, col = rgb(0,1,0,0.2/2.8), 
     xlim = c(29000, 38500), ylim = c(2000, 11200), yaxt = "n", xaxt='n') 
dev.off()

## plot 30d1 Cdc42
tiff("/media/gulab/GUDR2/PB10_new/30d1_Cdc42.tiff", height = 15, width = 15, units = 'in', res = 500)
par(mar = rep(0, 4))
## plot all reads as grey
plot(reads_table_30d1[,2:3],pch=15,cex=.8,
     col=rgb(0.8,0.8,0.8,reads_table_30d1$UMIcounts/500),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
par(new=TRUE)
plot(Cdc42_all[,3:4], 
     pch=19,cex=1.5, col = '#9C27B0',
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
dev.off()

## plot CNR1 Cdc42
tiff("/media/gulab/GUDR2/PB10_new/CNR1_Cdc42.tiff", height = 15, width = 15, units = 'in', res = 500)
par(mar = rep(0, 4))
## plot all reads as grey
plot(reads_table_CN1[,2:3],pch=15,cex=.4,
     col=rgb(0.6,0.6,0.6,reads_table_CN1$UMIcounts/650),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n')
par(new=TRUE)
plot(Cdc42_all[,3:4], 
     pch=19,cex=1.5, col = '#9C27B0',
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n') 
dev.off()

## plot 30d1 Stmn2
tiff("/media/gulab/GUDR2/PB10_new/30d1_Stmn2.tiff", height = 15, width = 15, units = 'in', res = 500)
par(mar = rep(0, 4))
## plot all reads as grey
plot(reads_table_30d1[,2:3],pch=15,cex=.8,
     col=rgb(0.8,0.8,0.8,reads_table_30d1$UMIcounts/500),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
par(new=TRUE)
plot(Stmn2_all[,3:4], 
     pch=19,cex=1.5, col = '#9C27B0',
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
dev.off()

## plot CNR1 Stmn2
tiff("/media/gulab/GUDR2/PB10_new/CNR1_Stmn2.tiff", height = 15, width = 15, units = 'in', res = 500)
par(mar = rep(0, 4))
## plot all reads as grey
plot(reads_table_CN1[,2:3],pch=15,cex=.4,
     col=rgb(0.6,0.6,0.6,reads_table_CN1$UMIcounts/650),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n')
par(new=TRUE)
plot(Stmn2_all[,3:4], 
     pch=19,cex=1.5, col = '#9C27B0',
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n') 
dev.off()

## plot 30d1 Cck
tiff("/media/gulab/GUDR2/PB10_new/30d1_Cck.tiff", height = 15, width = 15, units = 'in', res = 500)
par(mar = rep(0, 4))
## plot all reads as grey
plot(reads_table_30d1[,2:3],pch=15,cex=.8,
     col=rgb(0.8,0.8,0.8,reads_table_30d1$UMIcounts/500),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
par(new=TRUE)
plot(Cck_all[,3:4], 
     pch=19,cex=1.5, col = '#9C27B0',
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
dev.off()

## plot CNR1 Cck
tiff("/media/gulab/GUDR2/PB10_new/CNR1_Cck.tiff", height = 15, width = 15, units = 'in', res = 500)
par(mar = rep(0, 4))
## plot all reads as grey
plot(reads_table_CN1[,2:3],pch=15,cex=.4,
     col=rgb(0.6,0.6,0.6,reads_table_CN1$UMIcounts/650),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n')
par(new=TRUE)
plot(Cck_all[,3:4], 
     pch=19,cex=1.5, col = '#9C27B0',
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n') 
dev.off()

## plot 30d1 Selenom
tiff("/media/gulab/GUDR2/PB10_new/30d1_Selenom.tiff", height = 15, width = 15, units = 'in', res = 500)
par(mar = rep(0, 4))
## plot all reads as grey
plot(reads_table_30d1[,2:3],pch=15,cex=.8,
     col=rgb(0.8,0.8,0.8,reads_table_30d1$UMIcounts/500),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
par(new=TRUE)
plot(Selenom_all[,3:4], 
     pch=19,cex=1.5, col = '#9C27B0',
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
dev.off()

## plot CNR1 Selenom
tiff("/media/gulab/GUDR2/PB10_new/CNR1_Selenom.tiff", height = 15, width = 15, units = 'in', res = 500)
par(mar = rep(0, 4))
## plot all reads as grey
plot(reads_table_CN1[,2:3],pch=15,cex=.4,
     col=rgb(0.6,0.6,0.6,reads_table_CN1$UMIcounts/650),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n')
par(new=TRUE)
plot(Selenom_all[,3:4], 
     pch=19,cex=1.5, col = '#9C27B0',
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n') 
dev.off()

## plot 30d1 Ctxn2
tiff("/media/gulab/GUDR2/PB10_new/30d1_Ctxn2.tiff", height = 15, width = 15, units = 'in', res = 500)
par(mar = rep(0, 4))
## plot all reads as grey
plot(reads_table_30d1[,2:3],pch=15,cex=.8,
     col=rgb(0.8,0.8,0.8,reads_table_30d1$UMIcounts/500),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
par(new=TRUE)
plot(Ctxn2_all[,3:4], 
     pch=19,cex=1.5, col = '#9C27B0',
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
dev.off()

## plot CNR1 Ctxn2
tiff("/media/gulab/GUDR2/PB10_new/CNR1_Ctxn2.tiff", height = 15, width = 15, units = 'in', res = 500)
par(mar = rep(0, 4))
## plot all reads as grey
plot(reads_table_CN1[,2:3],pch=15,cex=.4,
     col=rgb(0.6,0.6,0.6,reads_table_CN1$UMIcounts/650),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n')
par(new=TRUE)
plot(Ctxn2_all[,3:4], 
     pch=19,cex=1.5, col = '#9C27B0',
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n') 
dev.off()

## plot 30d1 Pin1
tiff("/media/gulab/GUDR2/PB10_new/30d1_Pin1.tiff", height = 15, width = 15, units = 'in', res = 500)
par(mar = rep(0, 4))
## plot all reads as grey
plot(reads_table_30d1[,2:3],pch=15,cex=.8,
     col=rgb(0.8,0.8,0.8,reads_table_30d1$UMIcounts/500),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
par(new=TRUE)
plot(Pin1_all[,3:4], 
     pch=19,cex=1.5, col = '#9C27B0',
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
dev.off()

## plot CNR1 Pin1
tiff("/media/gulab/GUDR2/PB10_new/CNR1_Pin1.tiff", height = 15, width = 15, units = 'in', res = 500)
par(mar = rep(0, 4))
## plot all reads as grey
plot(reads_table_CN1[,2:3],pch=15,cex=.4,
     col=rgb(0.6,0.6,0.6,reads_table_CN1$UMIcounts/650),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n')
par(new=TRUE)
plot(Pin1_all[,3:4], 
     pch=19,cex=1.5, col = '#9C27B0',
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n') 
dev.off()

## plot 30d1 Crh
tiff("/media/gulab/GUDR2/PB10_new/30d1_Crh.tiff", height = 15, width = 15, units = 'in', res = 500)
par(mar = rep(0, 4))
## plot all reads as grey
plot(reads_table_30d1[,2:3],pch=15,cex=.8,
     col=rgb(0.8,0.8,0.8,reads_table_30d1$UMIcounts/500),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n')
par(new=TRUE)
## plot Crh reads as blue
plot(Crh_all[,3:4], 
     pch=19,cex=1.5, col = 'blue',
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
dev.off()

## plot CNR1 Crh
tiff("/media/gulab/GUDR2/PB10_new/CNR1_Crh.tiff", height = 15, width = 15, units = 'in', res = 500)
par(mar = rep(0, 4))
## plot all reads as grey
plot(reads_table_CN1[,2:3],pch=15,cex=.4,
     col=rgb(0.6,0.6,0.6,reads_table_CN1$UMIcounts/650),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n') 
par(new=TRUE)
## plot Crh reads as blue
plot(Crh_all[,3:4], 
     pch=19,cex=1.5, col = 'blue',
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n') 
dev.off()

## plot 30d1 Adcyap1
tiff("/media/gulab/GUDR2/PB10_new/30d1_Adcyap1.tiff", height = 1200, width = 1600)
par(mar = rep(0, 4))
## plot all reads as grey
plot(reads_table_30d1[,2:3],pch=19,cex=.1,
     col=rgb(0.1,0.1,0.1,reads_table_30d1$UMIcounts/2000), yaxt = "n", xaxt='n')
par(new=TRUE)
## plot Adcyap1 reads as blue
plot(Adcyap1_all[,3:4], 
     pch=19,cex=1.5, col = rgb(0,1,0,0.2/2.0), yaxt = "n", xaxt='n') 
dev.off()

## plot CNR1 Adcyap1
tiff("/media/gulab/GUDR2/PB10_new/CNR1_Adcyap1.tiff", height = 1200, width = 1600)
par(mar = rep(0, 4))
## plot all reads as grey
plot(reads_table_CN1[,2:3],pch=19,cex=.1,
     col=rgb(0.1,0.1,0.1,reads_table_CN1$UMIcounts/2800), yaxt = "n", xaxt='n') 
par(new=TRUE)
## plot Adcyap1 reads as blue
plot(Adcyap1_all[,3:4], 
     pch=19,cex=1.5, col = rgb(0,1,0,0.2/2.8), yaxt = "n", xaxt='n') 
dev.off()

## Seg_cluster_prg2_pos_updated.txt from file_convert
## neuronSeurat.rds from Xiaonan
Coord_for_All = read.csv("/media/gulab/GUDR2/PB10_new/Seg_cluster_prg2_pos_updated.txt", header = FALSE, sep = "\t")
only_neuron_4 = readRDS("/media/gulab/GUDR2/PB10_new/neuronSeurat.rds")
Calca_30d1_lst = colnames(PBN_4[,PBN_4@assays$RNA@counts['Calca',] > 3 & PBN_4$orig.ident=='30d1'])
Calca_30d1_lst_selected = unique(Coord_for_All[Coord_for_All$V6 %in% Calca_30d1_lst 
                                               & Coord_for_All$V1 >= 19300
                                               & Coord_for_All$V1 <= 23900
                                               & Coord_for_All$V2 >= 3800
                                               & Coord_for_All$V2 <= 8400,]$V6) 
Calca_CNR1_lst = colnames(PBN_4[,PBN_4@assays$RNA@counts['Calca',] > 3 & PBN_4$orig.ident=='CNR1'])
Calca_CNR1_lst_selected = unique(Coord_for_All[Coord_for_All$V6 %in% Calca_CNR1_lst 
                                               & Coord_for_All$V1 >= 30800
                                               & Coord_for_All$V1 <= 36000
                                               & Coord_for_All$V2 >= 3000
                                               & Coord_for_All$V2 <= 8200,]$V6) 
Penk_30d1_lst = colnames(PBN_4[,PBN_4@assays$RNA@counts['Penk',] > 3 & PBN_4$orig.ident=='30d1'])
Penk_30d1_lst_selected = unique(Coord_for_All[Coord_for_All$V6 %in% Penk_30d1_lst 
                                               & Coord_for_All$V1 >= 19300
                                               & Coord_for_All$V1 <= 23900
                                               & Coord_for_All$V2 >= 3800
                                               & Coord_for_All$V2 <= 8400,]$V6)
Penk_CNR1_lst = colnames(PBN_4[,PBN_4@assays$RNA@counts['Penk',] > 3 & PBN_4$orig.ident=='CNR1'])
Penk_CNR1_lst_selected = unique(Coord_for_All[Coord_for_All$V6 %in% Penk_CNR1_lst 
                                               & Coord_for_All$V1 >= 30800
                                               & Coord_for_All$V1 <= 36000
                                               & Coord_for_All$V2 >= 3000
                                               & Coord_for_All$V2 <= 8200,]$V6) 
Scg2_30d1_lst = colnames(PBN_4[,PBN_4@assays$RNA@counts['Scg2',] > 3 & PBN_4$orig.ident=='30d1'])
Scg2_30d1_lst_selected = unique(Coord_for_All[Coord_for_All$V6 %in% Scg2_30d1_lst 
                                              & Coord_for_All$V1 >= 19300
                                              & Coord_for_All$V1 <= 23900
                                              & Coord_for_All$V2 >= 3800
                                              & Coord_for_All$V2 <= 8400,]$V6)
Scg2_CNR1_lst = colnames(PBN_4[,PBN_4@assays$RNA@counts['Scg2',] > 3 & PBN_4$orig.ident=='CNR1'])
Scg2_CNR1_lst_selected = unique(Coord_for_All[Coord_for_All$V6 %in% Scg2_CNR1_lst 
                                               & Coord_for_All$V1 >= 30800
                                               & Coord_for_All$V1 <= 36000
                                               & Coord_for_All$V2 >= 3000
                                               & Coord_for_All$V2 <= 8200,]$V6) 
Cck_30d1_lst = colnames(PBN_4[,PBN_4@assays$RNA@counts['Cck',] > 3 & PBN_4$orig.ident=='30d1'])
Cck_30d1_lst_selected = unique(Coord_for_All[Coord_for_All$V6 %in% Cck_30d1_lst 
                                              & Coord_for_All$V1 >= 19300
                                              & Coord_for_All$V1 <= 23900
                                              & Coord_for_All$V2 >= 3800
                                              & Coord_for_All$V2 <= 8400,]$V6)
Cck_CNR1_lst = colnames(PBN_4[,PBN_4@assays$RNA@counts['Cck',] > 3 & PBN_4$orig.ident=='CNR1'])
Cck_CNR1_lst_selected = unique(Coord_for_All[Coord_for_All$V6 %in% Cck_CNR1_lst 
                                               & Coord_for_All$V1 >= 30800
                                               & Coord_for_All$V1 <= 36000
                                               & Coord_for_All$V2 >= 3000
                                               & Coord_for_All$V2 <= 8200,]$V6) 

tiff("/media/gulab/GUDR2/PB10_new/Calca_cell_30d1.tiff", height = 1500, width = 1500)
par(mar = rep(0, 4))
plot(reads_table_30d1[,2:3],pch=15,cex=.8,
     col=rgb(0.6,0.6,0.6,reads_table_30d1$UMIcounts/500),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
par(new=TRUE)
plot(Coord_for_All[Coord_for_All$V6 %in% Calca_30d1_lst,1:2], 
     pch=15,cex=.5, col = rgb(0,1,0),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
dev.off()

tiff("/media/gulab/GUDR2/PB10_new/Calca_cell_CNR1.tiff", height = 1500, width = 1500)
par(mar = rep(0, 4))
plot(reads_table_CN1[,2:3],pch=15,cex=.4,
     col=rgb(0.05,0.05,0.05,reads_table_CN1$UMIcounts/650),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n')
par(new=TRUE)
plot(Coord_for_All[Coord_for_All$V6 %in% Calca_CNR1_lst,1:2], 
     pch=15,cex=.5, col = rgb(0,1,0),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n') 
dev.off()

tiff("/media/gulab/GUDR2/PB10_new/Cck_cell_30d1.tiff", height = 1500, width = 1500)
par(mar = rep(0, 4))
plot(reads_table_30d1[,2:3],pch=15,cex=.8,
     col=rgb(0.6,0.6,0.6,reads_table_30d1$UMIcounts/500),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
par(new=TRUE)
plot(Coord_for_All[Coord_for_All$V6 %in% Cck_30d1_lst,1:2], 
     pch=15,cex=.5, col = rgb(1,0,1),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
dev.off()

tiff("/media/gulab/GUDR2/PB10_new/Cck_cell_CNR1.tiff", height = 1500, width = 1500)
par(mar = rep(0, 4))
plot(reads_table_CN1[,2:3],pch=15,cex=.4,
     col=rgb(0.05,0.05,0.05,reads_table_CN1$UMIcounts/650),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n')
par(new=TRUE)
plot(Coord_for_All[Coord_for_All$V6 %in% Cck_CNR1_lst,1:2], 
     pch=15,cex=.5, col = rgb(1,0,1),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n') 
dev.off()

tiff("/media/gulab/GUDR2/PB10_new/Penk_cell_30d1.tiff", height = 1500, width = 1500)
par(mar = rep(0, 4))
plot(reads_table_30d1[,2:3],pch=15,cex=.8,
     col=rgb(0.6,0.6,0.6,reads_table_30d1$UMIcounts/500),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
par(new=TRUE)
plot(Coord_for_All[Coord_for_All$V6 %in% Penk_30d1_lst,1:2], 
     pch=15,cex=.5, col = rgb(1,0,0),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
dev.off()

tiff("/media/gulab/GUDR2/PB10_new/Penk_cell_CNR1.tiff", height = 1500, width = 1500)
par(mar = rep(0, 4))
plot(reads_table_CN1[,2:3],pch=15,cex=.4,
     col=rgb(0.05,0.05,0.05,reads_table_CN1$UMIcounts/650),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n')
par(new=TRUE)
plot(Coord_for_All[Coord_for_All$V6 %in% Penk_CNR1_lst,1:2], 
     pch=15,cex=.5, col = rgb(1,0,0),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n') 
dev.off()

tiff("/media/gulab/GUDR2/PB10_new/Scg2_cell_30d1.tiff", height = 1500, width = 1500)
par(mar = rep(0, 4))
plot(reads_table_30d1[,2:3],pch=15,cex=.8,
     col=rgb(0.6,0.6,0.6,reads_table_30d1$UMIcounts/500),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
par(new=TRUE)
plot(Coord_for_All[Coord_for_All$V6 %in% Scg2_30d1_lst,1:2], 
     pch=15,cex=.5, col = rgb(0,0,1),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
dev.off()

tiff("/media/gulab/GUDR2/PB10_new/Scg2_cell_CNR1.tiff", height = 1500, width = 1500)
par(mar = rep(0, 4))
plot(reads_table_CN1[,2:3],pch=15,cex=.4,
     col=rgb(0.05,0.05,0.05,reads_table_CN1$UMIcounts/650),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n')
par(new=TRUE)
plot(Coord_for_All[Coord_for_All$V6 %in% Scg2_CNR1_lst,1:2], 
     pch=15,cex=.5, col = rgb(0,0,1),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n') 
dev.off()

tiff("/media/gulab/GUDR2/PB10_new/All_cell_30d1.tiff", height = 1500, width = 1500)
par(mar = rep(0, 4))
plot(reads_table_30d1[,2:3],pch=15,cex=.8,
     col=rgb(0.6,0.6,0.6,reads_table_30d1$UMIcounts/500),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
par(new=TRUE)
plot(Coord_for_All[Coord_for_All$V6 %in% Calca_30d1_lst,1:2], 
     pch=15,cex=.5, col = rgb(0,1,0),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
par(new=TRUE)
plot(Coord_for_All[Coord_for_All$V6 %in% Penk_30d1_lst,1:2], 
     pch=15,cex=.5, col = rgb(1,0,0),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n')
par(new=TRUE)
plot(Coord_for_All[Coord_for_All$V6 %in% Scg2_30d1_lst,1:2], 
     pch=15,cex=.5, col = rgb(0,0,1),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
par(new=TRUE)
plot(Coord_for_All[Coord_for_All$V6 %in% Cck_30d1_lst,1:2], 
     pch=15,cex=.5, col = rgb(1,0,1),
     ylim = c(3800, 8400), xlim = c(19300, 23900), yaxt = "n", xaxt='n') 
dev.off()

tiff("/media/gulab/GUDR2/PB10_new/All_cell_CNR1.tiff", height = 1500, width = 1500)
par(mar = rep(0, 4))
plot(reads_table_CN1[,2:3],pch=15,cex=.4,
     col=rgb(0.05,0.05,0.05,reads_table_CN1$UMIcounts/650),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n')
par(new=TRUE)
plot(Coord_for_All[Coord_for_All$V6 %in% Calca_CNR1_lst,1:2], 
     pch=15,cex=.5, col = rgb(0,1,0),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n') 
par(new=TRUE)
plot(Coord_for_All[Coord_for_All$V6 %in% Penk_CNR1_lst,1:2], 
     pch=15,cex=.5, col = rgb(1,0,0),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n') 
par(new=TRUE)
plot(Coord_for_All[Coord_for_All$V6 %in% Scg2_CNR1_lst,1:2], 
     pch=15,cex=.5, col = rgb(0,0,1),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n') 
par(new=TRUE)
plot(Coord_for_All[Coord_for_All$V6 %in% Cck_CNR1_lst,1:2], 
     pch=15,cex=.5, col = rgb(1,0,1),
     ylim = c(3000, 8200), xlim = c(30800, 36000), yaxt = "n", xaxt='n') 
dev.off()

DimPlot(only_neuron_4, label = TRUE)

table(only_neuron_4[,colnames(only_neuron_4) %in% Calca_CNR1_lst_selected]$seurat_clusters)
sum(table(only_neuron_4[,colnames(only_neuron_4) %in% Calca_CNR1_lst_selected]$seurat_clusters))
table(only_neuron_4[,colnames(only_neuron_4) %in% Calca_30d1_lst_selected]$seurat_clusters) 
sum(table(only_neuron_4[,colnames(only_neuron_4) %in% Calca_30d1_lst_selected]$seurat_clusters))

table(only_neuron_4[,colnames(only_neuron_4) %in% Penk_CNR1_lst_selected]$seurat_clusters)
sum(table(only_neuron_4[,colnames(only_neuron_4) %in% Penk_CNR1_lst_selected]$seurat_clusters))
table(only_neuron_4[,colnames(only_neuron_4) %in% Penk_30d1_lst_selected]$seurat_clusters)
sum(table(only_neuron_4[,colnames(only_neuron_4) %in% Penk_30d1_lst_selected]$seurat_clusters))

table(only_neuron_4[,colnames(only_neuron_4) %in% Scg2_CNR1_lst_selected]$seurat_clusters)
sum(table(only_neuron_4[,colnames(only_neuron_4) %in% Scg2_CNR1_lst_selected]$seurat_clusters))
table(only_neuron_4[,colnames(only_neuron_4) %in% Scg2_30d1_lst_selected]$seurat_clusters)
sum(table(only_neuron_4[,colnames(only_neuron_4) %in% Scg2_30d1_lst_selected]$seurat_clusters))

table(only_neuron_4[,colnames(only_neuron_4) %in% Cck_CNR1_lst_selected]$seurat_clusters)
sum(table(only_neuron_4[,colnames(only_neuron_4) %in% Cck_CNR1_lst_selected]$seurat_clusters))
table(only_neuron_4[,colnames(only_neuron_4) %in% Cck_30d1_lst_selected]$seurat_clusters)
sum(table(only_neuron_4[,colnames(only_neuron_4) %in% Cck_30d1_lst_selected]$seurat_clusters))
