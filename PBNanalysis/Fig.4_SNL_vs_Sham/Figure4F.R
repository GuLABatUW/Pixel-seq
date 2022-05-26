## for Figure4-B

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

finalPost = readRDS("/media/gulab/GUDR2/PB10_new/finalPost.rds")
PBNS2 = readRDS("/media/gulab/GUDR2/PB10_new/PBNS2.rds")
neuronSeurat = readRDS("/media/gulab/GUDR2/PB10_new/neuronSeurat.rds")
otherSeurat = readRDS("/media/gulab/GUDR2/PB10_new/otherSeurat.rds")

regionPBN_30d1 = rownames(finalPost[finalPost$condition=='30d1' & finalPost$region %in% c(1,2,3,4),])
regionPBN_CN = rownames(finalPost[finalPost$condition=='CNR1' & finalPost$region %in% c(1,2,3,4),])
regionV_30d1 = rownames(finalPost[finalPost$condition=='30d1' & finalPost$region %in% c(8),])
regionV_CN = rownames(finalPost[finalPost$condition=='CNR1' & finalPost$region %in% c(8),])

regionPBN_upper_30d1 = rownames(finalPost[finalPost$condition=='30d1' & finalPost$region %in% c(1,2),])
regionPBN_upper_CN = rownames(finalPost[finalPost$condition=='CNR1' & finalPost$region %in% c(1,2),])
regionPBN_lower_30d1 = rownames(finalPost[finalPost$condition=='30d1' & finalPost$region %in% c(3,4),])
regionPBN_lower_CN = rownames(finalPost[finalPost$condition=='CNR1' & finalPost$region %in% c(3,4),])

regionPBN_up1_30d1 = rownames(finalPost[finalPost$condition=='30d1' & finalPost$region %in% c(1),])
regionPBN_up1_CN = rownames(finalPost[finalPost$condition=='CNR1' & finalPost$region %in% c(1),])
regionPBN_up2_30d1 = rownames(finalPost[finalPost$condition=='30d1' & finalPost$region %in% c(2),])
regionPBN_up2_CN = rownames(finalPost[finalPost$condition=='CNR1' & finalPost$region %in% c(2),])

regionPBN_lw1_30d1 = rownames(finalPost[finalPost$condition=='30d1' & finalPost$region %in% c(3),])
regionPBN_lw1_CN = rownames(finalPost[finalPost$condition=='CNR1' & finalPost$region %in% c(3),])
regionPBN_lw2_30d1 = rownames(finalPost[finalPost$condition=='30d1' & finalPost$region %in% c(4),])
regionPBN_lw2_CN = rownames(finalPost[finalPost$condition=='CNR1' & finalPost$region %in% c(4),])


metaAgene = PBNS2@meta.data
SctValueAll = PBNS2@assays$SCT@data
SctValueNeu = neuronSeurat@assays$SCT@data
SctValueOth = otherSeurat@assays$SCT@data

geneAll = SctValueAll[rownames(SctValueAll)=='Calca',]
geneNeu = SctValueNeu[rownames(SctValueNeu)=='Calca',]
geneOth = SctValueOth[rownames(SctValueOth)=='Calca',]

#geneAll = geneAll[geneAll>0]
#geneNeu = geneNeu[geneNeu>0]
#geneOth = geneOth[geneOth>0]

geneAll = as.data.frame(geneAll)
geneNeu = as.data.frame(geneNeu)
geneOth = as.data.frame(geneOth)

mytypeAll = t(matrix(unlist(strsplit(rownames(geneAll),'_')),3,nrow(geneAll)))
geneAll$condition = mytypeAll[,1]

mytypeNeu = t(matrix(unlist(strsplit(rownames(geneNeu),'_')),3,nrow(geneNeu)))
geneNeu$condition = mytypeNeu[,1]
geneNeu$S2NEUcluster = 50

plot1neuronSeurat = DimPlot(neuronSeurat)
plot1otherSeurat = DimPlot(otherSeurat)

geneNeu[rownames(plot1neuronSeurat$data),]$S2NEUcluster = plot1neuronSeurat$data$ident

mytypeOth = t(matrix(unlist(strsplit(rownames(geneOth),'_')),3,nrow(geneOth)))
geneOth$condition = mytypeOth[,1]
geneOth$S2OTHcluster = 50

geneOth[rownames(plot1otherSeurat$data),]$S2OTHcluster = plot1otherSeurat$data$ident

S2_CNR1= geneAll[geneAll$condition=='CNR1','geneAll']
S2_30d1 = geneAll[geneAll$condition=='30d1','geneAll']

PBN_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_CN,'geneAll']
PBN_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_30d1,'geneAll']

V_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionV_CN,'geneAll']
V_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionV_30d1,'geneAll']

pUp1_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_up1_CN,'geneAll']
pUp1_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_up1_30d1,'geneAll']

pUp2_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_up2_CN,'geneAll']
pUp2_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_up2_30d1,'geneAll']

pDn1_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_lw1_CN,'geneAll']
pDn1_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_lw1_30d1,'geneAll']

pDn2_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_lw2_CN,'geneAll']   
pDn2_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_lw2_30d1,'geneAll']

V_CNR1_cNeu8 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_up2_CN & geneNeu$S2NEUcluster==8,'geneNeu']
V_30d1_cNeu8 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_up2_30d1 & geneNeu$S2NEUcluster==8,'geneNeu']

PBN_CNR1_cNeu1 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_up2_CN & geneNeu$S2NEUcluster==1,'geneNeu']
PBN_30d1_cNeu1 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_up2_30d1 & geneNeu$S2NEUcluster==1,'geneNeu']

PBN_CNR1_cNeu3 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==3,'geneNeu']
PBN_30d1_cNeu3 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==3,'geneNeu']

PBN_CNR1_cNeu7 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==7,'geneNeu']
PBN_30d1_cNeu7 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==7,'geneNeu']

V_CNR1_cNeu9 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionV_CN & geneNeu$S2NEUcluster==9,'geneNeu']
V_30d1_cNeu9 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionV_30d1 & geneNeu$S2NEUcluster==9,'geneNeu']

PBN_CNR1_cNeu11 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==11,'geneNeu']
PBN_30d1_cNeu11 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==11,'geneNeu']

PBN_CNR1_cNeu12 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==12,'geneNeu']
PBN_30d1_cNeu12 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==12,'geneNeu']

PBN_CNR1_cNeu16 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==16,'geneNeu']
PBN_30d1_cNeu16 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==16,'geneNeu']

PBN_CNR1_cOth1 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==1,'geneOth']
PBN_30d1_cOth1 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==1,'geneOth']

PBN_CNR1_cOth2 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==2,'geneOth']
PBN_30d1_cOth2 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==2,'geneOth']

PBN_CNR1_cOth10 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==10,'geneOth']
PBN_30d1_cOth10 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==10,'geneOth']
wilcox.test(S2_CNR1, S2_30d1)
wilcox.test(PBN_CNR1, PBN_30d1)
wilcox.test(V_CNR1, V_30d1)
mean(V_30d1)/mean(V_CNR1)
wilcox.test(pUp1_CNR1, pUp1_30d1)
mean(pUp1_30d1)/mean(pUp1_CNR1)
wilcox.test(pUp2_CNR1, pUp2_30d1)
wilcox.test(pDn1_CNR1, pDn1_30d1)
wilcox.test(pDn2_CNR1, pDn2_30d1)
wilcox.test(V_CNR1_cNeu8, V_30d1_cNeu8)
wilcox.test(PBN_CNR1_cNeu1, PBN_30d1_cNeu1)
wilcox.test(PBN_CNR1_cNeu3, PBN_30d1_cNeu3)
wilcox.test(PBN_CNR1_cNeu7, PBN_30d1_cNeu7)
wilcox.test(V_CNR1_cNeu9, V_30d1_cNeu9)
mean(V_30d1_cNeu9)/mean(V_CNR1_cNeu9)
wilcox.test(PBN_CNR1_cNeu11, PBN_30d1_cNeu11)
wilcox.test(PBN_CNR1_cNeu12, PBN_30d1_cNeu12)
wilcox.test(PBN_CNR1_cNeu16, PBN_30d1_cNeu16)
wilcox.test(PBN_CNR1_cOth1, PBN_30d1_cOth1)
wilcox.test(PBN_CNR1_cOth2, PBN_30d1_cOth2)
wilcox.test(PBN_CNR1_cOth10, PBN_30d1_cOth10)

mean(V_30d1)/mean(V_CNR1)
mean(V_30d1_cNeu9)/mean(V_CNR1_cNeu9)

geneValue = list(S2_CNR1= geneAll[geneAll$condition=='CNR1','geneAll'],
                  S2_30d1 = geneAll[geneAll$condition=='30d1','geneAll'],
                 
                 PBN_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_CN,'geneAll'],
                 PBN_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_30d1,'geneAll'],
                 
                 V_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionV_CN,'geneAll'],
                 V_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionV_30d1,'geneAll'],
                 
                 pUp1_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_up1_CN,'geneAll'],
                 pUp1_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_up1_30d1,'geneAll'],
                 
                 pUp2_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_up2_CN,'geneAll'],
                 pUp2_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_up2_30d1,'geneAll'],
                 
                 pDn1_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_lw1_CN,'geneAll'],
                 pDn1_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_lw1_30d1,'geneAll'],
                 
                 pDn2_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_lw2_CN,'geneAll'],    
                 pDn2_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_lw2_30d1,'geneAll'],
                 
                 V_CNR1_cNeu8 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_up2_CN & geneNeu$S2NEUcluster==8,'geneNeu'],
                 V_30d1_cNeu8 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_up2_30d1 & geneNeu$S2NEUcluster==8,'geneNeu'],
                 
                 PBN_CNR1_cNeu1 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_up2_CN & geneNeu$S2NEUcluster==1,'geneNeu'],
                 PBN_30d1_cNeu1 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_up2_30d1 & geneNeu$S2NEUcluster==1,'geneNeu'],
                 
                 PBN_CNR1_cNeu3 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==3,'geneNeu'],
                 PBN_30d1_cNeu3 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==3,'geneNeu'],
                 
                 PBN_CNR1_cNeu7 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==7,'geneNeu'],
                 PBN_30d1_cNeu7 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==7,'geneNeu'],
                 
                 V_CNR1_cNeu9 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionV_CN & geneNeu$S2NEUcluster==9,'geneNeu'],
                 V_30d1_cNeu9 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionV_30d1 & geneNeu$S2NEUcluster==9,'geneNeu'],
                 
                 PBN_CNR1_cNeu11 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==11,'geneNeu'],
                 PBN_30d1_cNeu11 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==11,'geneNeu'],
                 
                 PBN_CNR1_cNeu12 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==12,'geneNeu'],
                 PBN_30d1_cNeu12 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==12,'geneNeu'],
                 
                 PBN_CNR1_cNeu16 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==16,'geneNeu'],
                 PBN_30d1_cNeu16 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==16,'geneNeu']
                 
                 # PBN_CNR1_cOth1 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==1,'geneOth'],
                 # PBN_30d1_cOth1 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==1,'geneOth'],
                 # 
                 # PBN_CNR1_cOth2 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==2,'geneOth'],
                 # PBN_30d1_cOth2 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==2,'geneOth'],
                 # 
                 # PBN_CNR1_cOth10 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==10,'geneOth'],
                 # PBN_30d1_cOth10 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==10,'geneOth']
)

std <- function(x) sd(x)/sqrt(length(x))

png("/media/gulab/GUDR2/PB10_new/Calca.png", height = 2, width = 28, units = 'in', res = 300)
par(mar = rep(0, 4))
plot(c(0,31),c(0,3),type='n',xlab = 'F',las=1,axes = F)
axis(4,at=c(0, 1, 2, 3),las=1)
colVio = rep(c(rgb(0,0,0.8,0.3), rgb(0.8,0,0,0.3)),18)
for(i in 1:30){
  vioplot(geneValue[[i]],at=i,drawRect=F,add=T,col=colVio[i],lwd = 6)
  meanV =mean(geneValue[[i]])
  sde = std(geneValue[[i]])
  lines(c(i-0.2,i+0.2),c(meanV,meanV),lwd=5)
  lines(c(i,i),c(meanV-sde,meanV+sde),lwd=5)
  lines(c(i-0.1,i+0.1),c(meanV-sde,meanV-sde),lwd=3)
  lines(c(i-0.1,i+0.1),c(meanV+sde,meanV+sde),lwd=3)
  points(i,meanV,pch=15,cex=2)
}
dev.off()

geneAll = SctValueAll[rownames(SctValueAll)=='Penk',]
geneNeu = SctValueNeu[rownames(SctValueNeu)=='Penk',]
geneOth = SctValueOth[rownames(SctValueOth)=='Penk',]

# geneAll = geneAll[geneAll>0]
# geneNeu = geneNeu[geneNeu>0]
# geneOth = geneOth[geneOth>0]

geneAll = as.data.frame(geneAll)
geneNeu = as.data.frame(geneNeu)
geneOth = as.data.frame(geneOth)

mytypeAll = t(matrix(unlist(strsplit(rownames(geneAll),'_')),3,nrow(geneAll)))
geneAll$condition = mytypeAll[,1]

mytypeNeu = t(matrix(unlist(strsplit(rownames(geneNeu),'_')),3,nrow(geneNeu)))
geneNeu$condition = mytypeNeu[,1]
geneNeu$S2NEUcluster = 50

plot1neuronSeurat = DimPlot(neuronSeurat)
plot1otherSeurat = DimPlot(otherSeurat)

geneNeu[rownames(plot1neuronSeurat$data),]$S2NEUcluster = plot1neuronSeurat$data$ident

mytypeOth = t(matrix(unlist(strsplit(rownames(geneOth),'_')),3,nrow(geneOth)))
geneOth$condition = mytypeOth[,1]
geneOth$S2OTHcluster = 50

geneOth[rownames(plot1otherSeurat$data),]$S2OTHcluster = plot1otherSeurat$data$ident

S2_CNR1= geneAll[geneAll$condition=='CNR1','geneAll']
S2_30d1 = geneAll[geneAll$condition=='30d1','geneAll']

PBN_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_CN,'geneAll']
PBN_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_30d1,'geneAll']

V_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionV_CN,'geneAll']
V_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionV_30d1,'geneAll']

pUp1_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_up1_CN,'geneAll']
pUp1_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_up1_30d1,'geneAll']

pUp2_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_up2_CN,'geneAll']
pUp2_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_up2_30d1,'geneAll']

pDn1_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_lw1_CN,'geneAll']
pDn1_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_lw1_30d1,'geneAll']

pDn2_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_lw2_CN,'geneAll']   
pDn2_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_lw2_30d1,'geneAll']

V_CNR1_cNeu8 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_up2_CN & geneNeu$S2NEUcluster==8,'geneNeu']
V_30d1_cNeu8 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_up2_30d1 & geneNeu$S2NEUcluster==8,'geneNeu']

PBN_CNR1_cNeu1 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_up2_CN & geneNeu$S2NEUcluster==1,'geneNeu']
PBN_30d1_cNeu1 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_up2_30d1 & geneNeu$S2NEUcluster==1,'geneNeu']

PBN_CNR1_cNeu3 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==3,'geneNeu']
PBN_30d1_cNeu3 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==3,'geneNeu']

PBN_CNR1_cNeu7 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==7,'geneNeu']
PBN_30d1_cNeu7 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==7,'geneNeu']

V_CNR1_cNeu9 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionV_CN & geneNeu$S2NEUcluster==9,'geneNeu']
V_30d1_cNeu9 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionV_30d1 & geneNeu$S2NEUcluster==9,'geneNeu']

PBN_CNR1_cNeu11 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==11,'geneNeu']
PBN_30d1_cNeu11 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==11,'geneNeu']

PBN_CNR1_cNeu12 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==12,'geneNeu']
PBN_30d1_cNeu12 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==12,'geneNeu']

PBN_CNR1_cNeu16 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==16,'geneNeu']
PBN_30d1_cNeu16 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==16,'geneNeu']

PBN_CNR1_cOth1 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==1,'geneOth']
PBN_30d1_cOth1 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==1,'geneOth']

PBN_CNR1_cOth2 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==2,'geneOth']
PBN_30d1_cOth2 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==2,'geneOth']

PBN_CNR1_cOth10 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==10,'geneOth']
PBN_30d1_cOth10 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==10,'geneOth']

wilcox.test(S2_CNR1, S2_30d1)
mean(S2_30d1)/mean(S2_CNR1)

wilcox.test(PBN_CNR1, PBN_30d1)
mean(PBN_30d1)/mean(PBN_CNR1)

wilcox.test(V_CNR1, V_30d1)
wilcox.test(pUp1_CNR1, pUp1_30d1)
wilcox.test(pUp2_CNR1, pUp2_30d1)
mean(pUp2_30d1)/mean(pUp2_CNR1)

wilcox.test(pDn1_CNR1, pDn1_30d1)
mean(pDn1_30d1)/mean(pDn1_CNR1)

wilcox.test(pDn2_CNR1, pDn2_30d1)
mean(pDn2_30d1)/mean(pDn2_CNR1)

wilcox.test(V_CNR1_cNeu8, V_30d1_cNeu8)
mean(V_30d1_cNeu8)/mean(V_CNR1_cNeu8)

wilcox.test(PBN_CNR1_cNeu1, PBN_30d1_cNeu1)
mean(PBN_30d1_cNeu1)/mean(PBN_CNR1_cNeu1)

wilcox.test(PBN_CNR1_cNeu3, PBN_30d1_cNeu3)
wilcox.test(PBN_CNR1_cNeu7, PBN_30d1_cNeu7)
wilcox.test(V_CNR1_cNeu9, V_30d1_cNeu9)
wilcox.test(PBN_CNR1_cNeu11, PBN_30d1_cNeu11)
mean(PBN_30d1_cNeu11)/mean(PBN_CNR1_cNeu11)

wilcox.test(PBN_CNR1_cNeu12, PBN_30d1_cNeu12)
mean(PBN_30d1_cNeu12)/mean(PBN_CNR1_cNeu12)

wilcox.test(PBN_CNR1_cNeu16, PBN_30d1_cNeu16)
wilcox.test(PBN_CNR1_cOth1, PBN_30d1_cOth1)
wilcox.test(PBN_CNR1_cOth2, PBN_30d1_cOth2)
wilcox.test(PBN_CNR1_cOth10, PBN_30d1_cOth10)

mean(PBN_CNR1_cOth10_up2)
mean(PBN_30d1_cOth10_up2)

mean(V_30d1)/mean(V_CNR1)
mean(V_30d1_cNeu9)/mean(V_CNR1_cNeu9)

geneValue = list(S2_CNR1= geneAll[geneAll$condition=='CNR1','geneAll'],
                 S2_30d1 = geneAll[geneAll$condition=='30d1','geneAll'],
                 
                 PBN_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_CN,'geneAll'],
                 PBN_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_30d1,'geneAll'],
                 
                 V_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionV_CN,'geneAll'],
                 V_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionV_30d1,'geneAll'],
                 
                 pUp1_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_up1_CN,'geneAll'],
                 pUp1_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_up1_30d1,'geneAll'],
                 
                 pUp2_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_up2_CN,'geneAll'],
                 pUp2_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_up2_30d1,'geneAll'],
                 
                 pDn1_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_lw1_CN,'geneAll'],
                 pDn1_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_lw1_30d1,'geneAll'],
                 
                 pDn2_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_lw2_CN,'geneAll'],    
                 pDn2_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_lw2_30d1,'geneAll'],
                 
                 V_CNR1_cNeu8 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_up2_CN & geneNeu$S2NEUcluster==8,'geneNeu'],
                 V_30d1_cNeu8 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_up2_30d1 & geneNeu$S2NEUcluster==8,'geneNeu'],
                 
                 PBN_CNR1_cNeu1 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_up2_CN & geneNeu$S2NEUcluster==1,'geneNeu'],
                 PBN_30d1_cNeu1 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_up2_30d1 & geneNeu$S2NEUcluster==1,'geneNeu'],
                 
                 PBN_CNR1_cNeu3 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==3,'geneNeu'],
                 PBN_30d1_cNeu3 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==3,'geneNeu'],
                 
                 PBN_CNR1_cNeu7 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==7,'geneNeu'],
                 PBN_30d1_cNeu7 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==7,'geneNeu'],
                 
                 V_CNR1_cNeu9 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionV_CN & geneNeu$S2NEUcluster==9,'geneNeu'],
                 V_30d1_cNeu9 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionV_30d1 & geneNeu$S2NEUcluster==9,'geneNeu'],
                 
                 PBN_CNR1_cNeu11 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==11,'geneNeu'],
                 PBN_30d1_cNeu11 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==11,'geneNeu'],
                 
                 PBN_CNR1_cNeu12 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==12,'geneNeu'],
                 PBN_30d1_cNeu12 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==12,'geneNeu'],
                 
                 PBN_CNR1_cNeu16 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==16,'geneNeu'],
                 PBN_30d1_cNeu16 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==16,'geneNeu']
                 
                 # PBN_CNR1_cOth1 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==1,'geneOth'],
                 # PBN_30d1_cOth1 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==1,'geneOth'],
                 # 
                 # PBN_CNR1_cOth2 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==2,'geneOth'],
                 # PBN_30d1_cOth2 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==2,'geneOth'],
                 # 
                 # PBN_CNR1_cOth10 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==10,'geneOth'],
                 # PBN_30d1_cOth10 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==10,'geneOth']
)

library(vioplot)
std <- function(x) sd(x)/sqrt(length(x))

png("/media/gulab/GUDR2/PB10_new/Penk.png", height = 2, width = 30, units = 'in', res = 300)
par(mar = rep(0, 4))
plot(c(0,31),c(0,1.5),type='n',xlab = 'F',las=1,axes = F)
axis(2,at=c(0,0.75,1.5),las=1)
colVio = rep(c(rgb(0,0,0.8,0.3), rgb(0.8,0,0,0.3)),18)
for(i in 1:30){
  vioplot(geneValue[[i]],at=i,drawRect=F,add=T,col=colVio[i],lwd = 6)
  meanV =mean(geneValue[[i]])
  sde = std(geneValue[[i]])
  lines(c(i-0.2,i+0.2),c(meanV,meanV),lwd=5)
  lines(c(i,i),c(meanV-sde,meanV+sde),lwd=5)
  lines(c(i-0.1,i+0.1),c(meanV-sde,meanV-sde),lwd=3)
  lines(c(i-0.1,i+0.1),c(meanV+sde,meanV+sde),lwd=3)
  points(i,meanV,pch=15,cex=2)
}
dev.off()

geneAll = SctValueAll[rownames(SctValueAll)=='Scg2',]
geneNeu = SctValueNeu[rownames(SctValueNeu)=='Scg2',]
geneOth = SctValueOth[rownames(SctValueOth)=='Scg2',]

#geneAll = geneAll[geneAll>0]
#geneNeu = geneNeu[geneNeu>0]
#geneOth = geneOth[geneOth>0]

geneAll = as.data.frame(geneAll)
geneNeu = as.data.frame(geneNeu)
geneOth = as.data.frame(geneOth)

mytypeAll = t(matrix(unlist(strsplit(rownames(geneAll),'_')),3,nrow(geneAll)))
geneAll$condition = mytypeAll[,1]

mytypeNeu = t(matrix(unlist(strsplit(rownames(geneNeu),'_')),3,nrow(geneNeu)))
geneNeu$condition = mytypeNeu[,1]
geneNeu$S2NEUcluster = 50

plot1neuronSeurat = DimPlot(neuronSeurat)
plot1otherSeurat = DimPlot(otherSeurat)

geneNeu[rownames(plot1neuronSeurat$data),]$S2NEUcluster = plot1neuronSeurat$data$ident

mytypeOth = t(matrix(unlist(strsplit(rownames(geneOth),'_')),3,nrow(geneOth)))
geneOth$condition = mytypeOth[,1]
geneOth$S2OTHcluster = 50

geneOth[rownames(plot1otherSeurat$data),]$S2OTHcluster = plot1otherSeurat$data$ident

S2_CNR1= geneAll[geneAll$condition=='CNR1','geneAll']
S2_30d1 = geneAll[geneAll$condition=='30d1','geneAll']

PBN_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_CN,'geneAll']
PBN_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_30d1,'geneAll']

V_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionV_CN,'geneAll']
V_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionV_30d1,'geneAll']

pUp1_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_up1_CN,'geneAll']
pUp1_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_up1_30d1,'geneAll']

pUp2_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_up2_CN,'geneAll']
pUp2_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_up2_30d1,'geneAll']

pDn1_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_lw1_CN,'geneAll']
pDn1_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_lw1_30d1,'geneAll']

pDn2_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_lw2_CN,'geneAll']   
pDn2_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_lw2_30d1,'geneAll']

V_CNR1_cNeu8 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_up2_CN & geneNeu$S2NEUcluster==8,'geneNeu']
V_30d1_cNeu8 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_up2_30d1 & geneNeu$S2NEUcluster==8,'geneNeu']

PBN_CNR1_cNeu1 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_up2_CN & geneNeu$S2NEUcluster==1,'geneNeu']
PBN_30d1_cNeu1 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_up2_30d1 & geneNeu$S2NEUcluster==1,'geneNeu']

PBN_CNR1_cNeu3 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==3,'geneNeu']
PBN_30d1_cNeu3 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==3,'geneNeu']

PBN_CNR1_cNeu7 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==7,'geneNeu']
PBN_30d1_cNeu7 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==7,'geneNeu']

V_CNR1_cNeu9 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionV_CN & geneNeu$S2NEUcluster==9,'geneNeu']
V_30d1_cNeu9 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionV_30d1 & geneNeu$S2NEUcluster==9,'geneNeu']

PBN_CNR1_cNeu11 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==11,'geneNeu']
PBN_30d1_cNeu11 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==11,'geneNeu']

PBN_CNR1_cNeu12 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==12,'geneNeu']
PBN_30d1_cNeu12 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==12,'geneNeu']

PBN_CNR1_cNeu16 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==16,'geneNeu']
PBN_30d1_cNeu16 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==16,'geneNeu']

PBN_CNR1_cOth1 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==1,'geneOth']
PBN_30d1_cOth1 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==1,'geneOth']

PBN_CNR1_cOth2 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==2,'geneOth']
PBN_30d1_cOth2 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==2,'geneOth']

PBN_CNR1_cOth10 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==10,'geneOth']
PBN_30d1_cOth10 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==10,'geneOth']

wilcox.test(S2_CNR1, S2_30d1)
mean(S2_30d1)/mean(S2_CNR1)

wilcox.test(PBN_CNR1, PBN_30d1)
mean(PBN_30d1)/mean(PBN_CNR1)

wilcox.test(V_CNR1, V_30d1)
mean(V_30d1)/mean(V_CNR1)

wilcox.test(pUp1_CNR1, pUp1_30d1)
mean(pUp1_30d1)/mean(pUp1_CNR1)

wilcox.test(pUp2_CNR1, pUp2_30d1)
mean(pUp2_30d1)/mean(pUp2_CNR1)

wilcox.test(pDn1_CNR1, pDn1_30d1)
mean(pDn1_30d1)/mean(pDn1_CNR1)

wilcox.test(pDn2_CNR1, pDn2_30d1)
mean(pDn2_30d1)/mean(pDn2_CNR1)

wilcox.test(V_CNR1_cNeu8, V_30d1_cNeu8)
mean(V_30d1_cNeu8)/mean(V_CNR1_cNeu8)

wilcox.test(PBN_CNR1_cNeu1, PBN_30d1_cNeu1)
wilcox.test(PBN_CNR1_cNeu3, PBN_30d1_cNeu3)
wilcox.test(PBN_CNR1_cNeu7, PBN_30d1_cNeu7)
mean(PBN_30d1_cNeu7)/mean(PBN_CNR1_cNeu7)

wilcox.test(V_CNR1_cNeu9, V_30d1_cNeu9)
wilcox.test(PBN_CNR1_cNeu11, PBN_30d1_cNeu11)
mean(PBN_30d1_cNeu11)/mean(PBN_CNR1_cNeu11)

wilcox.test(PBN_CNR1_cNeu12, PBN_30d1_cNeu12)
mean(PBN_30d1_cNeu12)/mean(PBN_CNR1_cNeu12)

wilcox.test(PBN_CNR1_cNeu16, PBN_30d1_cNeu16)
wilcox.test(PBN_CNR1_cOth1, PBN_30d1_cOth1)
wilcox.test(PBN_CNR1_cOth2, PBN_30d1_cOth2)
mean(PBN_30d1_cOth2)/mean(PBN_CNR1_cOth2)

wilcox.test(PBN_CNR1_cOth10, PBN_30d1_cOth10)

mean(V_30d1)/mean(V_CNR1)
mean(V_30d1_cNeu9)/mean(V_CNR1_cNeu9)

geneValue = list(S2_CNR1= geneAll[geneAll$condition=='CNR1','geneAll'],
                 S2_30d1 = geneAll[geneAll$condition=='30d1','geneAll'],
                 
                 PBN_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_CN,'geneAll'],
                 PBN_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_30d1,'geneAll'],
                 
                 V_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionV_CN,'geneAll'],
                 V_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionV_30d1,'geneAll'],
                 
                 pUp1_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_up1_CN,'geneAll'],
                 pUp1_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_up1_30d1,'geneAll'],
                 
                 pUp2_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_up2_CN,'geneAll'],
                 pUp2_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_up2_30d1,'geneAll'],
                 
                 pDn1_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_lw1_CN,'geneAll'],
                 pDn1_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_lw1_30d1,'geneAll'],
                 
                 pDn2_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_lw2_CN,'geneAll'],    
                 pDn2_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_lw2_30d1,'geneAll'],
                 
                 V_CNR1_cNeu8 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_up2_CN & geneNeu$S2NEUcluster==8,'geneNeu'],
                 V_30d1_cNeu8 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_up2_30d1 & geneNeu$S2NEUcluster==8,'geneNeu'],
                 
                 PBN_CNR1_cNeu1 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_up2_CN & geneNeu$S2NEUcluster==1,'geneNeu'],
                 PBN_30d1_cNeu1 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_up2_30d1 & geneNeu$S2NEUcluster==1,'geneNeu'],
                 
                 PBN_CNR1_cNeu3 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==3,'geneNeu'],
                 PBN_30d1_cNeu3 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==3,'geneNeu'],
                 
                 PBN_CNR1_cNeu7 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==7,'geneNeu'],
                 PBN_30d1_cNeu7 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==7,'geneNeu'],
                 
                 V_CNR1_cNeu9 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionV_CN & geneNeu$S2NEUcluster==9,'geneNeu'],
                 V_30d1_cNeu9 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionV_30d1 & geneNeu$S2NEUcluster==9,'geneNeu'],
                 
                 PBN_CNR1_cNeu11 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==11,'geneNeu'],
                 PBN_30d1_cNeu11 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==11,'geneNeu'],
                 
                 PBN_CNR1_cNeu12 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==12,'geneNeu'],
                 PBN_30d1_cNeu12 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==12,'geneNeu'],
                 
                 PBN_CNR1_cNeu16 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==16,'geneNeu'],
                 PBN_30d1_cNeu16 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==16,'geneNeu']
                 
                 # PBN_CNR1_cOth1 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==1,'geneOth'],
                 # PBN_30d1_cOth1 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==1,'geneOth'],
                 # 
                 # PBN_CNR1_cOth2 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==2,'geneOth'],
                 # PBN_30d1_cOth2 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==2,'geneOth'],
                 # 
                 # PBN_CNR1_cOth10 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==10,'geneOth'],
                 # PBN_30d1_cOth10 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==10,'geneOth']
)

library(vioplot)
std <- function(x) sd(x)/sqrt(length(x))

png("/media/gulab/GUDR2/PB10_new/Scg.png", height = 2, width = 30, units = 'in', res = 300)
par(mar = rep(0, 4))
plot(c(0,31),c(0,1),type='n',xlab = 'F',las=1,axes = F)
axis(2,at=c(0,0.5,1),las=1)
colVio = rep(c(rgb(0,0,0.8,0.3), rgb(0.8,0,0,0.3)),18)
for(i in 1:30){
  vioplot(geneValue[[i]],at=i,drawRect=F,add=T,col=colVio[i],lwd = 6)
  meanV =mean(geneValue[[i]])
  sde = std(geneValue[[i]])
  lines(c(i-0.2,i+0.2),c(meanV,meanV),lwd=5)
  lines(c(i,i),c(meanV-sde,meanV+sde),lwd=5)
  lines(c(i-0.1,i+0.1),c(meanV-sde,meanV-sde),lwd=3)
  lines(c(i-0.1,i+0.1),c(meanV+sde,meanV+sde),lwd=3)
  points(i,meanV,pch=15,cex=2)
}
dev.off()

geneAll = SctValueAll[rownames(SctValueAll)=='Apoe',]
geneNeu = SctValueNeu[rownames(SctValueNeu)=='Apoe',]
geneOth = SctValueOth[rownames(SctValueOth)=='Apoe',]

#geneAll = geneAll[geneAll>0]
#geneNeu = geneNeu[geneNeu>0]
#geneOth = geneOth[geneOth>0]

geneAll = as.data.frame(geneAll)
geneNeu = as.data.frame(geneNeu)
geneOth = as.data.frame(geneOth)

mytypeAll = t(matrix(unlist(strsplit(rownames(geneAll),'_')),3,nrow(geneAll)))
geneAll$condition = mytypeAll[,1]

mytypeNeu = t(matrix(unlist(strsplit(rownames(geneNeu),'_')),3,nrow(geneNeu)))
geneNeu$condition = mytypeNeu[,1]
geneNeu$S2NEUcluster = 50

plot1neuronSeurat = DimPlot(neuronSeurat)
plot1otherSeurat = DimPlot(otherSeurat)

geneNeu[rownames(plot1neuronSeurat$data),]$S2NEUcluster = plot1neuronSeurat$data$ident

mytypeOth = t(matrix(unlist(strsplit(rownames(geneOth),'_')),3,nrow(geneOth)))
geneOth$condition = mytypeOth[,1]
geneOth$S2OTHcluster = 50

geneOth[rownames(plot1otherSeurat$data),]$S2OTHcluster = plot1otherSeurat$data$ident

S2_CNR1= geneAll[geneAll$condition=='CNR1','geneAll']
S2_30d1 = geneAll[geneAll$condition=='30d1','geneAll']

PBN_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_CN,'geneAll']
PBN_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_30d1,'geneAll']

V_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionV_CN,'geneAll']
V_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionV_30d1,'geneAll']

pUp1_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_up1_CN,'geneAll']
pUp1_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_up1_30d1,'geneAll']

pUp2_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_up2_CN,'geneAll']
pUp2_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_up2_30d1,'geneAll']

pDn1_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_lw1_CN,'geneAll']
pDn1_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_lw1_30d1,'geneAll']

pDn2_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_lw2_CN,'geneAll']   
pDn2_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_lw2_30d1,'geneAll']

V_CNR1_cNeu8 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_up2_CN & geneNeu$S2NEUcluster==8,'geneNeu']
V_30d1_cNeu8 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_up2_30d1 & geneNeu$S2NEUcluster==8,'geneNeu']

PBN_CNR1_cNeu1 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_up2_CN & geneNeu$S2NEUcluster==1,'geneNeu']
PBN_30d1_cNeu1 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_up2_30d1 & geneNeu$S2NEUcluster==1,'geneNeu']

PBN_CNR1_cNeu3 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==3,'geneNeu']
PBN_30d1_cNeu3 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==3,'geneNeu']

PBN_CNR1_cNeu7 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==7,'geneNeu']
PBN_30d1_cNeu7 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==7,'geneNeu']

V_CNR1_cNeu9 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionV_CN & geneNeu$S2NEUcluster==9,'geneNeu']
V_30d1_cNeu9 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionV_30d1 & geneNeu$S2NEUcluster==9,'geneNeu']

PBN_CNR1_cNeu11 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==11,'geneNeu']
PBN_30d1_cNeu11 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==11,'geneNeu']

PBN_CNR1_cNeu12 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==12,'geneNeu']
PBN_30d1_cNeu12 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==12,'geneNeu']

PBN_CNR1_cNeu16 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==16,'geneNeu']
PBN_30d1_cNeu16 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==16,'geneNeu']

PBN_CNR1_cOth1 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==1,'geneOth']
PBN_30d1_cOth1 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==1,'geneOth']

PBN_CNR1_cOth2 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==2,'geneOth']
PBN_30d1_cOth2 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==2,'geneOth']

PBN_CNR1_cOth10 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==10,'geneOth']
PBN_30d1_cOth10 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==10,'geneOth']

wilcox.test(S2_CNR1, S2_30d1)
mean(S2_30d1)/mean(S2_CNR1)

wilcox.test(PBN_CNR1, PBN_30d1)
mean(PBN_30d1)/mean(PBN_CNR1)

wilcox.test(V_CNR1, V_30d1)
mean(V_30d1)/mean(V_CNR1)

wilcox.test(pUp1_CNR1, pUp1_30d1)
mean(pUp1_30d1)/mean(pUp1_CNR1)

wilcox.test(pUp2_CNR1, pUp2_30d1)
mean(pUp2_30d1)/mean(pUp2_CNR1)

wilcox.test(pDn1_CNR1, pDn1_30d1)
mean(pDn1_30d1)/mean(pDn1_CNR1)

wilcox.test(pDn2_CNR1, pDn2_30d1)
mean(pDn2_30d1)/mean(pDn2_CNR1)

wilcox.test(V_CNR1_cNeu8, V_30d1_cNeu8)
wilcox.test(PBN_CNR1_cNeu1, PBN_30d1_cNeu1)
mean(PBN_30d1_cNeu1)/mean(PBN_CNR1_cNeu1)

wilcox.test(PBN_CNR1_cNeu3, PBN_30d1_cNeu3)
mean(PBN_30d1_cNeu3)/mean(PBN_CNR1_cNeu3)

wilcox.test(PBN_CNR1_cNeu7, PBN_30d1_cNeu7)
mean(PBN_30d1_cNeu7)/mean(PBN_CNR1_cNeu7)

wilcox.test(V_CNR1_cNeu9, V_30d1_cNeu9)
mean(V_30d1_cNeu9)/mean(V_CNR1_cNeu9)

wilcox.test(PBN_CNR1_cNeu11, PBN_30d1_cNeu11)
mean(PBN_30d1_cNeu11)/mean(PBN_CNR1_cNeu11)

wilcox.test(PBN_CNR1_cNeu12, PBN_30d1_cNeu12)
mean(PBN_30d1_cNeu12)/mean(PBN_CNR1_cNeu12)

wilcox.test(PBN_CNR1_cNeu16, PBN_30d1_cNeu16)
mean(PBN_30d1_cNeu16)/mean(PBN_CNR1_cNeu16)

wilcox.test(PBN_CNR1_cOth1, PBN_30d1_cOth1)
mean(PBN_30d1_cOth1)/mean(PBN_CNR1_cOth1)

wilcox.test(PBN_CNR1_cOth2, PBN_30d1_cOth2)
mean(PBN_30d1_cOth2)/mean(PBN_CNR1_cOth2)

wilcox.test(PBN_CNR1_cOth10, PBN_30d1_cOth10)
mean(PBN_30d1_cOth10)/mean(PBN_CNR1_cOth10)

mean(V_30d1)/mean(V_CNR1)
mean(V_30d1_cNeu9)/mean(V_CNR1_cNeu9)

geneValue = list(S2_CNR1= geneAll[geneAll$condition=='CNR1','geneAll'],
                 S2_30d1 = geneAll[geneAll$condition=='30d1','geneAll'],
                 
                 PBN_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_CN,'geneAll'],
                 PBN_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_30d1,'geneAll'],
                 
                 V_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionV_CN,'geneAll'],
                 V_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionV_30d1,'geneAll'],
                 
                 pUp1_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_up1_CN,'geneAll'],
                 pUp1_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_up1_30d1,'geneAll'],
                 
                 pUp2_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_up2_CN,'geneAll'],
                 pUp2_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_up2_30d1,'geneAll'],
                 
                 pDn1_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_lw1_CN,'geneAll'],
                 pDn1_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_lw1_30d1,'geneAll'],
                 
                 pDn2_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_lw2_CN,'geneAll'],    
                 pDn2_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_lw2_30d1,'geneAll'],
                 
                 V_CNR1_cNeu8 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_up2_CN & geneNeu$S2NEUcluster==8,'geneNeu'],
                 V_30d1_cNeu8 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_up2_30d1 & geneNeu$S2NEUcluster==8,'geneNeu'],
                 
                 PBN_CNR1_cNeu1 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_up2_CN & geneNeu$S2NEUcluster==1,'geneNeu'],
                 PBN_30d1_cNeu1 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_up2_30d1 & geneNeu$S2NEUcluster==1,'geneNeu'],
                 
                 PBN_CNR1_cNeu3 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==3,'geneNeu'],
                 PBN_30d1_cNeu3 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==3,'geneNeu'],
                 
                 PBN_CNR1_cNeu7 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==7,'geneNeu'],
                 PBN_30d1_cNeu7 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==7,'geneNeu'],
                 
                 V_CNR1_cNeu9 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionV_CN & geneNeu$S2NEUcluster==9,'geneNeu'],
                 V_30d1_cNeu9 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionV_30d1 & geneNeu$S2NEUcluster==9,'geneNeu'],
                 
                 PBN_CNR1_cNeu11 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==11,'geneNeu'],
                 PBN_30d1_cNeu11 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==11,'geneNeu'],
                 
                 PBN_CNR1_cNeu12 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==12,'geneNeu'],
                 PBN_30d1_cNeu12 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==12,'geneNeu'],
                 
                 PBN_CNR1_cNeu16 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==16,'geneNeu'],
                 PBN_30d1_cNeu16 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==16,'geneNeu'],
                 
                 PBN_CNR1_cOth1 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==1,'geneOth'],
                 PBN_30d1_cOth1 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==1,'geneOth'],
                 
                 PBN_CNR1_cOth2 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==2,'geneOth'],
                 PBN_30d1_cOth2 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==2,'geneOth'],
                 
                 PBN_CNR1_cOth10 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==10,'geneOth'],
                 PBN_30d1_cOth10 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==10,'geneOth']
)

std <- function(x) sd(x)/sqrt(length(x))

png("/media/gulab/Hector/PB10_new/Apoe.png", height = 2, width = 30, units = 'in', res = 300)
par(mar = rep(0, 4))
plot(c(0,37),c(0,2),type='n',xlab = 'F',las=1,axes = F)
axis(3,at=c(0, 1, 2),las=1)
colVio = rep(c(rgb(0,0,0.8,0.3), rgb(0.8,0,0,0.3)),18)
for(i in 1:36){
  vioplot(geneValue[[i]],at=i,drawRect=F,add=T,col=colVio[i],lwd = 6)
  meanV =mean(geneValue[[i]])
  sde = std(geneValue[[i]])
  lines(c(i-0.2,i+0.2),c(meanV,meanV),lwd=5)
  lines(c(i,i),c(meanV-sde,meanV+sde),lwd=5)
  lines(c(i-0.1,i+0.1),c(meanV-sde,meanV-sde),lwd=3)
  lines(c(i-0.1,i+0.1),c(meanV+sde,meanV+sde),lwd=3)
  points(i,meanV,pch=15,cex=2)
}
dev.off()

geneAll = SctValueAll[rownames(SctValueAll)=='Selenom',]
geneNeu = SctValueNeu[rownames(SctValueNeu)=='Selenom',]
geneOth = SctValueOth[rownames(SctValueOth)=='Selenom',]

#geneAll = geneAll[geneAll>0]
#geneNeu = geneNeu[geneNeu>0]
#geneOth = geneOth[geneOth>0]

geneAll = as.data.frame(geneAll)
geneNeu = as.data.frame(geneNeu)
geneOth = as.data.frame(geneOth)

mytypeAll = t(matrix(unlist(strsplit(rownames(geneAll),'_')),3,nrow(geneAll)))
geneAll$condition = mytypeAll[,1]

mytypeNeu = t(matrix(unlist(strsplit(rownames(geneNeu),'_')),3,nrow(geneNeu)))
geneNeu$condition = mytypeNeu[,1]
geneNeu$S2NEUcluster = 50

plot1neuronSeurat = DimPlot(neuronSeurat)
plot1otherSeurat = DimPlot(otherSeurat)

geneNeu[rownames(plot1neuronSeurat$data),]$S2NEUcluster = plot1neuronSeurat$data$ident

mytypeOth = t(matrix(unlist(strsplit(rownames(geneOth),'_')),3,nrow(geneOth)))
geneOth$condition = mytypeOth[,1]
geneOth$S2OTHcluster = 50

geneOth[rownames(plot1otherSeurat$data),]$S2OTHcluster = plot1otherSeurat$data$ident

S2_CNR1= geneAll[geneAll$condition=='CNR1','geneAll']
S2_30d1 = geneAll[geneAll$condition=='30d1','geneAll']

PBN_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_CN,'geneAll']
PBN_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_30d1,'geneAll']

V_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionV_CN,'geneAll']
V_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionV_30d1,'geneAll']

pUp1_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_up1_CN,'geneAll']
pUp1_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_up1_30d1,'geneAll']

pUp2_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_up2_CN,'geneAll']
pUp2_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_up2_30d1,'geneAll']

pDn1_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_lw1_CN,'geneAll']
pDn1_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_lw1_30d1,'geneAll']

pDn2_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_lw2_CN,'geneAll']   
pDn2_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_lw2_30d1,'geneAll']

V_CNR1_cNeu8 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_up2_CN & geneNeu$S2NEUcluster==8,'geneNeu']
V_30d1_cNeu8 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_up2_30d1 & geneNeu$S2NEUcluster==8,'geneNeu']

PBN_CNR1_cNeu1 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_up2_CN & geneNeu$S2NEUcluster==1,'geneNeu']
PBN_30d1_cNeu1 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_up2_30d1 & geneNeu$S2NEUcluster==1,'geneNeu']

PBN_CNR1_cNeu3 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==3,'geneNeu']
PBN_30d1_cNeu3 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==3,'geneNeu']

PBN_CNR1_cNeu7 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==7,'geneNeu']
PBN_30d1_cNeu7 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==7,'geneNeu']

V_CNR1_cNeu9 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionV_CN & geneNeu$S2NEUcluster==9,'geneNeu']
V_30d1_cNeu9 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionV_30d1 & geneNeu$S2NEUcluster==9,'geneNeu']

PBN_CNR1_cNeu11 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==11,'geneNeu']
PBN_30d1_cNeu11 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==11,'geneNeu']

PBN_CNR1_cNeu12 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==12,'geneNeu']
PBN_30d1_cNeu12 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==12,'geneNeu']

PBN_CNR1_cNeu16 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==16,'geneNeu']
PBN_30d1_cNeu16 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==16,'geneNeu']

PBN_CNR1_cOth1 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==1,'geneOth']
PBN_30d1_cOth1 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==1,'geneOth']

PBN_CNR1_cOth2 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==2,'geneOth']
PBN_30d1_cOth2 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==2,'geneOth']

PBN_CNR1_cOth10 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==10,'geneOth']
PBN_30d1_cOth10 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==10,'geneOth']

wilcox.test(S2_CNR1, S2_30d1)
mean(S2_30d1)/mean(S2_CNR1)

wilcox.test(PBN_CNR1, PBN_30d1)
mean(PBN_30d1)/mean(PBN_CNR1)

wilcox.test(V_CNR1, V_30d1)
mean(V_30d1)/mean(V_CNR1)

wilcox.test(pUp1_CNR1, pUp1_30d1)
mean(pUp1_30d1)/mean(pUp1_CNR1)

wilcox.test(pUp2_CNR1, pUp2_30d1)
wilcox.test(pDn1_CNR1, pDn1_30d1)
mean(pDn1_30d1)/mean(pDn1_CNR1)

wilcox.test(pDn2_CNR1, pDn2_30d1)
mean(pDn2_30d1)/mean(pDn2_CNR1)

wilcox.test(V_CNR1_cNeu8, V_30d1_cNeu8)
wilcox.test(PBN_CNR1_cNeu1, PBN_30d1_cNeu1)
wilcox.test(PBN_CNR1_cNeu3, PBN_30d1_cNeu3)
wilcox.test(PBN_CNR1_cNeu7, PBN_30d1_cNeu7)
wilcox.test(V_CNR1_cNeu9, V_30d1_cNeu9)
mean(V_30d1_cNeu9)/mean(V_CNR1_cNeu9)

wilcox.test(PBN_CNR1_cNeu11, PBN_30d1_cNeu11)
mean(PBN_30d1_cNeu11)/mean(PBN_CNR1_cNeu11)

wilcox.test(PBN_CNR1_cNeu12, PBN_30d1_cNeu12)
wilcox.test(PBN_CNR1_cNeu16, PBN_30d1_cNeu16)
mean(PBN_30d1_cNeu16)/mean(PBN_CNR1_cNeu16)

wilcox.test(PBN_CNR1_cOth1, PBN_30d1_cOth1)
mean(PBN_30d1_cOth1)/mean(PBN_CNR1_cOth1)

wilcox.test(PBN_CNR1_cOth2, PBN_30d1_cOth2)
wilcox.test(PBN_CNR1_cOth10, PBN_30d1_cOth10)

mean(V_30d1)/mean(V_CNR1)
mean(V_30d1_cNeu9)/mean(V_CNR1_cNeu9)

geneValue = list(S2_CNR1= geneAll[geneAll$condition=='CNR1','geneAll'],
                 S2_30d1 = geneAll[geneAll$condition=='30d1','geneAll'],
                 
                 PBN_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_CN,'geneAll'],
                 PBN_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_30d1,'geneAll'],
                 
                 V_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionV_CN,'geneAll'],
                 V_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionV_30d1,'geneAll'],
                 
                 pUp1_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_up1_CN,'geneAll'],
                 pUp1_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_up1_30d1,'geneAll'],
                 
                 pUp2_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_up2_CN,'geneAll'],
                 pUp2_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_up2_30d1,'geneAll'],
                 
                 pDn1_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_lw1_CN,'geneAll'],
                 pDn1_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_lw1_30d1,'geneAll'],
                 
                 pDn2_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_lw2_CN,'geneAll'],    
                 pDn2_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_lw2_30d1,'geneAll'],
                 
                 V_CNR1_cNeu8 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_up2_CN & geneNeu$S2NEUcluster==8,'geneNeu'],
                 V_30d1_cNeu8 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_up2_30d1 & geneNeu$S2NEUcluster==8,'geneNeu'],
                 
                 PBN_CNR1_cNeu1 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_up2_CN & geneNeu$S2NEUcluster==1,'geneNeu'],
                 PBN_30d1_cNeu1 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_up2_30d1 & geneNeu$S2NEUcluster==1,'geneNeu'],
                 
                 PBN_CNR1_cNeu3 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==3,'geneNeu'],
                 PBN_30d1_cNeu3 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==3,'geneNeu'],
                 
                 PBN_CNR1_cNeu7 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==7,'geneNeu'],
                 PBN_30d1_cNeu7 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==7,'geneNeu'],
                 
                 V_CNR1_cNeu9 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionV_CN & geneNeu$S2NEUcluster==9,'geneNeu'],
                 V_30d1_cNeu9 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionV_30d1 & geneNeu$S2NEUcluster==9,'geneNeu'],
                 
                 PBN_CNR1_cNeu11 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==11,'geneNeu'],
                 PBN_30d1_cNeu11 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==11,'geneNeu'],
                 
                 PBN_CNR1_cNeu12 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==12,'geneNeu'],
                 PBN_30d1_cNeu12 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==12,'geneNeu'],
                 
                 PBN_CNR1_cNeu16 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==16,'geneNeu'],
                 PBN_30d1_cNeu16 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==16,'geneNeu'],
                 
                 PBN_CNR1_cOth1 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==1,'geneOth'],
                 PBN_30d1_cOth1 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==1,'geneOth'],
                 
                 PBN_CNR1_cOth2 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==2,'geneOth'],
                 PBN_30d1_cOth2 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==2,'geneOth'],
                 
                 PBN_CNR1_cOth10 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==10,'geneOth'],
                 PBN_30d1_cOth10 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==10,'geneOth']
)

std <- function(x) sd(x)/sqrt(length(x))

png("/media/gulab/Hector/PB10_new/Selenom.png", height = 2, width = 30, units = 'in', res = 300)
par(mar = rep(0, 4))
plot(c(0,37),c(0,1.5),type='n',xlab = 'F',las=1,axes = F)
axis(4,at=c(0, 0.5, 1, 1.5),las=1)
colVio = rep(c(rgb(0,0,0.8,0.3), rgb(0.8,0,0,0.3)),18)
for(i in 1:36){
  vioplot(geneValue[[i]],at=i,drawRect=F,add=T,col=colVio[i],lwd = 6)
  meanV =mean(geneValue[[i]])
  sde = std(geneValue[[i]])
  lines(c(i-0.2,i+0.2),c(meanV,meanV),lwd=5)
  lines(c(i,i),c(meanV-sde,meanV+sde),lwd=5)
  lines(c(i-0.1,i+0.1),c(meanV-sde,meanV-sde),lwd=3)
  lines(c(i-0.1,i+0.1),c(meanV+sde,meanV+sde),lwd=3)
  points(i,meanV,pch=15,cex=2)
}
dev.off()

geneAll = SctValueAll[rownames(SctValueAll)=='Cck',]
geneNeu = SctValueNeu[rownames(SctValueNeu)=='Cck',]
geneOth = SctValueOth[rownames(SctValueOth)=='Cck',]

# geneAll = geneAll[geneAll>0]
# geneNeu = geneNeu[geneNeu>0]
# geneOth = geneOth[geneOth>0]

geneAll = as.data.frame(geneAll)
geneNeu = as.data.frame(geneNeu)
geneOth = as.data.frame(geneOth)

mytypeAll = t(matrix(unlist(strsplit(rownames(geneAll),'_')),3,nrow(geneAll)))
geneAll$condition = mytypeAll[,1]

mytypeNeu = t(matrix(unlist(strsplit(rownames(geneNeu),'_')),3,nrow(geneNeu)))
geneNeu$condition = mytypeNeu[,1]
geneNeu$S2NEUcluster = 50

plot1neuronSeurat = DimPlot(neuronSeurat)
plot1otherSeurat = DimPlot(otherSeurat)

geneNeu[rownames(plot1neuronSeurat$data),]$S2NEUcluster = plot1neuronSeurat$data$ident

mytypeOth = t(matrix(unlist(strsplit(rownames(geneOth),'_')),3,nrow(geneOth)))
geneOth$condition = mytypeOth[,1]
geneOth$S2OTHcluster = 50

geneOth[rownames(plot1otherSeurat$data),]$S2OTHcluster = plot1otherSeurat$data$ident

S2_CNR1= geneAll[geneAll$condition=='CNR1','geneAll']
S2_30d1 = geneAll[geneAll$condition=='30d1','geneAll']

PBN_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_CN,'geneAll']
PBN_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_30d1,'geneAll']

V_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionV_CN,'geneAll']
V_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionV_30d1,'geneAll']

pUp1_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_up1_CN,'geneAll']
pUp1_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_up1_30d1,'geneAll']

pUp2_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_up2_CN,'geneAll']
pUp2_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_up2_30d1,'geneAll']

pDn1_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_lw1_CN,'geneAll']
pDn1_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_lw1_30d1,'geneAll']

pDn2_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_lw2_CN,'geneAll']   
pDn2_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_lw2_30d1,'geneAll']

V_CNR1_cNeu8 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_up2_CN & geneNeu$S2NEUcluster==8,'geneNeu']
V_30d1_cNeu8 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_up2_30d1 & geneNeu$S2NEUcluster==8,'geneNeu']

PBN_CNR1_cNeu1 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_up2_CN & geneNeu$S2NEUcluster==1,'geneNeu']
PBN_30d1_cNeu1 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_up2_30d1 & geneNeu$S2NEUcluster==1,'geneNeu']

PBN_CNR1_cNeu3 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==3,'geneNeu']
PBN_30d1_cNeu3 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==3,'geneNeu']

PBN_CNR1_cNeu7 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==7,'geneNeu']
PBN_30d1_cNeu7 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==7,'geneNeu']

V_CNR1_cNeu9 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionV_CN & geneNeu$S2NEUcluster==9,'geneNeu']
V_30d1_cNeu9 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionV_30d1 & geneNeu$S2NEUcluster==9,'geneNeu']

PBN_CNR1_cNeu11 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==11,'geneNeu']
PBN_30d1_cNeu11 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==11,'geneNeu']

PBN_CNR1_cNeu12 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==12,'geneNeu']
PBN_30d1_cNeu12 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==12,'geneNeu']

PBN_CNR1_cNeu16 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==16,'geneNeu']
PBN_30d1_cNeu16 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==16,'geneNeu']

PBN_CNR1_cOth1 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==1,'geneOth']
PBN_30d1_cOth1 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==1,'geneOth']

PBN_CNR1_cOth2 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==2,'geneOth']
PBN_30d1_cOth2 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==2,'geneOth']

PBN_CNR1_cOth10 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==10,'geneOth']
PBN_30d1_cOth10 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==10,'geneOth']

wilcox.test(S2_CNR1, S2_30d1)
mean(S2_30d1)/mean(S2_CNR1)

wilcox.test(PBN_CNR1, PBN_30d1)
mean(PBN_30d1)/mean(PBN_CNR1)

wilcox.test(V_CNR1, V_30d1)
mean(V_30d1)/mean(V_CNR1)

wilcox.test(pUp1_CNR1, pUp1_30d1)
mean(pUp1_30d1)/mean(pUp1_CNR1)

wilcox.test(pUp2_CNR1, pUp2_30d1)
mean(pUp2_30d1)/mean(pUp2_CNR1)

wilcox.test(pDn1_CNR1, pDn1_30d1)
wilcox.test(pDn2_CNR1, pDn2_30d1)

wilcox.test(V_CNR1_cNeu8, V_30d1_cNeu8)
mean(V_30d1_cNeu8)/mean(V_CNR1_cNeu8)

wilcox.test(PBN_CNR1_cNeu1, PBN_30d1_cNeu1)
mean(PBN_30d1_cNeu1)/mean(PBN_CNR1_cNeu1)

wilcox.test(PBN_CNR1_cNeu3, PBN_30d1_cNeu3)
mean(PBN_30d1_cNeu3)/mean(PBN_CNR1_cNeu3)

wilcox.test(PBN_CNR1_cNeu7, PBN_30d1_cNeu7)

wilcox.test(V_CNR1_cNeu9, V_30d1_cNeu9)
mean(V_30d1_cNeu9)/mean(V_CNR1_cNeu9)

wilcox.test(PBN_CNR1_cNeu11, PBN_30d1_cNeu11)
wilcox.test(PBN_CNR1_cNeu12, PBN_30d1_cNeu12)
wilcox.test(PBN_CNR1_cNeu16, PBN_30d1_cNeu16)
wilcox.test(PBN_CNR1_cOth1, PBN_30d1_cOth1)
mean(PBN_30d1_cOth1)/mean(PBN_CNR1_cOth1)
wilcox.test(PBN_CNR1_cOth2, PBN_30d1_cOth2)
mean(PBN_30d1_cOth2)/mean(PBN_CNR1_cOth2)
wilcox.test(PBN_CNR1_cOth10, PBN_30d1_cOth10)
mean(PBN_30d1_cOth10)/mean(PBN_CNR1_cOth10)

mean(PBN_CNR1_cOth10_up2)
mean(PBN_30d1_cOth10_up2)

mean(V_30d1)/mean(V_CNR1)
mean(V_30d1_cNeu9)/mean(V_CNR1_cNeu9)

geneValue = list(S2_CNR1= geneAll[geneAll$condition=='CNR1','geneAll'],
                 S2_30d1 = geneAll[geneAll$condition=='30d1','geneAll'],
                 
                 PBN_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_CN,'geneAll'],
                 PBN_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_30d1,'geneAll'],
                 
                 V_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionV_CN,'geneAll'],
                 V_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionV_30d1,'geneAll'],
                 
                 pUp1_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_up1_CN,'geneAll'],
                 pUp1_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_up1_30d1,'geneAll'],
                 
                 pUp2_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_up2_CN,'geneAll'],
                 pUp2_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_up2_30d1,'geneAll'],
                 
                 pDn1_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_lw1_CN,'geneAll'],
                 pDn1_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_lw1_30d1,'geneAll'],
                 
                 pDn2_CNR1 = geneAll[geneAll$condition=='CNR1' & rownames(geneAll) %in% regionPBN_lw2_CN,'geneAll'],    
                 pDn2_30d1 = geneAll[geneAll$condition=='30d1' & rownames(geneAll) %in% regionPBN_lw2_30d1,'geneAll'],
                 
                 V_CNR1_cNeu8 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_up2_CN & geneNeu$S2NEUcluster==8,'geneNeu'],
                 V_30d1_cNeu8 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_up2_30d1 & geneNeu$S2NEUcluster==8,'geneNeu'],
                 
                 PBN_CNR1_cNeu1 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_up2_CN & geneNeu$S2NEUcluster==1,'geneNeu'],
                 PBN_30d1_cNeu1 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_up2_30d1 & geneNeu$S2NEUcluster==1,'geneNeu'],
                 
                 PBN_CNR1_cNeu3 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==3,'geneNeu'],
                 PBN_30d1_cNeu3 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==3,'geneNeu'],
                 
                 PBN_CNR1_cNeu7 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==7,'geneNeu'],
                 PBN_30d1_cNeu7 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==7,'geneNeu'],
                 
                 V_CNR1_cNeu9 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionV_CN & geneNeu$S2NEUcluster==9,'geneNeu'],
                 V_30d1_cNeu9 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionV_30d1 & geneNeu$S2NEUcluster==9,'geneNeu'],
                 
                 PBN_CNR1_cNeu11 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==11,'geneNeu'],
                 PBN_30d1_cNeu11 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==11,'geneNeu'],
                 
                 PBN_CNR1_cNeu12 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==12,'geneNeu'],
                 PBN_30d1_cNeu12 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==12,'geneNeu'],
                 
                 PBN_CNR1_cNeu16 = geneNeu[geneNeu$condition=='CNR1' & rownames(geneNeu) %in% regionPBN_CN & geneNeu$S2NEUcluster==16,'geneNeu'],
                 PBN_30d1_cNeu16 = geneNeu[geneNeu$condition=='30d1' & rownames(geneNeu) %in% regionPBN_30d1 & geneNeu$S2NEUcluster==16,'geneNeu']
                 
                 # PBN_CNR1_cOth1 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==1,'geneOth'],
                 # PBN_30d1_cOth1 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==1,'geneOth'],
                 # 
                 # PBN_CNR1_cOth2 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==2,'geneOth'],
                 # PBN_30d1_cOth2 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==2,'geneOth'],
                 # 
                 # PBN_CNR1_cOth10 = geneOth[geneOth$condition=='CNR1' & rownames(geneOth) %in% regionPBN_CN & geneOth$S2OTHcluster==10,'geneOth'],
                 # PBN_30d1_cOth10 = geneOth[geneOth$condition=='30d1' & rownames(geneOth) %in% regionPBN_30d1 & geneOth$S2OTHcluster==10,'geneOth']
)

library(vioplot)
std <- function(x) sd(x)/sqrt(length(x))

png("/media/gulab/GUDR2/PB10_new/Cck.png", height = 2, width = 30, units = 'in', res = 300)
par(mar = rep(0, 4))
plot(c(0,31),c(0,1),type='n',xlab = 'F',las=1,axes = F)
axis(2,at=c(0,0.5,1),las=1)
colVio = rep(c(rgb(0,0,0.8,0.3), rgb(0.8,0,0,0.3)),18)
for(i in 1:30){
  vioplot(geneValue[[i]],at=i,drawRect=F,add=T,col=colVio[i],lwd = 6)
  meanV =mean(geneValue[[i]])
  sde = std(geneValue[[i]])
  lines(c(i-0.2,i+0.2),c(meanV,meanV),lwd=5)
  lines(c(i,i),c(meanV-sde,meanV+sde),lwd=5)
  lines(c(i-0.1,i+0.1),c(meanV-sde,meanV-sde),lwd=3)
  lines(c(i-0.1,i+0.1),c(meanV+sde,meanV+sde),lwd=3)
  points(i,meanV,pch=15,cex=2)
}
dev.off()
