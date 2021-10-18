# For Sup Figure 3-12A

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

## display UMI counts reads map
### section PBS2_CN1uniq
### for edge cut
## PBS2_CN1uniq.coord from Xiaonan
reads_table_CN1 = read.csv("/media/gulab/GUDR2/PB10_new/PBS2_CN1uniq.coord", sep = ',', header = TRUE)
plot(reads_table_CN1[,2:3],pch=15,cex=.2,col=rgb(1,0,0,reads_table_CN1$UMIcounts/2000))

## get range
min(reads_table_CN1$CoordX)
## 29000
max(reads_table_CN1$CoordX)
## 38500
min(reads_table_CN1$CoordY)
## 2000
max(reads_table_CN1$CoordY)
## 11200

### UMI counts plot
tiff(filename = "/media/gulab/GUDR2/PB10_new/PB10_CN1.tiff", width = 1050, height = 920, units = "px", pointsize = 12,compression="lzw")
par(bg='black')
p = ggplot(reads_table_CN1, aes(x=CoordX, y=CoordY,colour=log2(UMIcounts))) + geom_point(size=0.2)
p + scale_fill_gradientn(limits = c(0, 10), colours = rev(viridis(256, option = "A")),aesthetics = "colour")+ theme_void()
dev.off()

### section PBS2_30d1uniq
### for edge cut
## PBS2_30d1uniq.coord from Xiaonan
reads_table_30d1 = read.csv("/media/gulab/GUDR2/PB10_new/PBS2_30d1uniq.coord", sep = ',', header = TRUE)
plot(reads_table_30d1[,2:3],pch=15,cex=.2,col=rgb(1,0,0,reads_table_30d1$UMIcounts/2000))

## get range
min(reads_table_30d1$CoordX)
## 16500
max(reads_table_30d1$CoordX)
## 26000
min(reads_table_30d1$CoordY)
## 3000
max(reads_table_30d1$CoordY)
## 11200

### UMI counts plot
tiff(filename = "/media/gulab/GUDR2/PB10_new/PB10_30d1.tiff", width = 1050, height = 820, units = "px", pointsize = 12,compression="lzw")
par(bg='black')
p = ggplot(reads_table_30d1, aes(x=CoordX, y=CoordY,colour=log2(UMIcounts))) + geom_point(size=0.2)
p + scale_fill_gradientn(limits = c(0, 10), colours = rev(viridis(256, option = "A")),aesthetics = "colour")+ theme_void()
dev.off()

### section PBS3_30d2uniq
### for edge cut
## PBS3_30d2uniq.coord from Xiaonan
reads_table_30d2 = read.csv("/media/gulab/GUDR2/PB10_new/PBS3_30d2uniq.coord", sep = ',', header = TRUE)
plot(reads_table_30d2[,2:3],pch=15,cex=.2,col=rgb(1,0,0,reads_table_30d2$UMIcounts/2000))

## get range
min(reads_table_30d2$CoordX)
## 21000
max(reads_table_30d2$CoordX)
## 31000
max(reads_table_30d2$CoordY)
## 11200

reads_table_30d2.clean <- reads_table_30d2[reads_table_30d2$CoordY >= 2000,]

### UMI counts plot
tiff(filename = "/media/gulab/GUDR2/PB10_new/PB10_30d2.tiff", width = 1100, height = 920, units = "px", pointsize = 12,compression="lzw")
par(bg='black')
p = ggplot(reads_table_30d2.clean, aes(x=CoordX, y=CoordY,colour=log2(UMIcounts))) + geom_point(size=0.2)
p + scale_fill_gradientn(limits = c(0, 10), colours = rev(viridis(256, option = "A")),aesthetics = "colour")+ theme_void()
dev.off()

### section PBS3_CN2uniq
### for edge cut
## PBS3_CN2uniq.coord from Xiaonan
reads_table_CN2 = read.csv("/media/gulab/GUDR2/PB10_new/PBS3_CN2uniq.coord", sep = ',', header = TRUE)
plot(reads_table_CN2[,2:3],pch=15,cex=.2,col=rgb(1,0,0,reads_table_CN2$UMIcounts/2000))
min(reads_table_CN2$CoordX)
## 36000
max(reads_table_CN2$CoordX)
## 47000
max(reads_table_CN2$CoordY)
## 11200

## only keep clean region
reads_table_CN2.clean <- reads_table_CN2[reads_table_CN2$CoordY >= 2000,]

### UMI counts plot
tiff(filename = "/media/gulab/GUDR2/PB10_new/PB10_CN2.tiff", width = 1200, height = 920, units = "px", pointsize = 12,compression="lzw")
par(bg='black')
p = ggplot(reads_table_CN2.clean, aes(x=CoordX, y=CoordY,colour=log2(UMIcounts))) + geom_point(size=0.2)
p + scale_fill_gradientn(limits = c(0, 10), colours = rev(viridis(256, option = "A")),aesthetics = "colour")+ theme_void()
dev.off()
