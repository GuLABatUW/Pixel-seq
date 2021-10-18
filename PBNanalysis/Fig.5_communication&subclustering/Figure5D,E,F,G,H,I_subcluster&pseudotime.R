library(DropSeq.util)
library('Seurat')
PBNall2 <- loadSparseDge("~/lab508/novaSeqPBN/PBNall/PBN_fullp2_cellExp.txt")
positions <- read.csv(
  file = "~/lab508/novaSeqPBN/PBNall/PBN_fullp2_coord.txt"
)
rownames(positions) = positions$binID
PBNct <- CreateSeuratObject(
  counts = PBNall2,
  project = 'SlideSeq',
  assay = 'Spatial'
)

PBNct[['image']] <-new(
  Class = 'SlideSeq',
  assay = "Spatial",
  coordinates = positions[,2:3]
)

PBNct = PBNct[,PBNct@meta.data$orig.ident %in% c('CNR1','30d1')]
PBNct = PBNct[,PBNct$nCount_Spatial>256]

otherSeurat = readRDS('~/lab508/novaSeqPBN/otherSeurat.rds')
neuroSeurat = readRDS('~/lab508/novaSeqPBN/neuronSeurat.rds')
DimPlot(neuroSeurat, reduction = "umap",label = T,split.by = 'orig.ident')
DimPlot(otherSeurat, reduction = "umap",label = T,split.by = 'orig.ident')

############################extract the astro and microglia cluster 
AstroN = otherSeurat@meta.data[otherSeurat@meta.data$seurat_clusters %in% c(1,6,10,11),]
neuroN = neuroSeurat@meta.data[neuroSeurat@meta.data$seurat_clusters %in% c(0,2,6,7,8,10,11,15),]
microgN = otherSeurat@meta.data[otherSeurat@meta.data$seurat_clusters %in% c(9),]

Astro = PBNct[,rownames(AstroN)]
neuro = neuroSeurat[,rownames(neuroN)]
microG = PBNct[,rownames(microgN)]

####################################################################
#######################################Astro and microglia Normalize
###################################################################
Astro = SCTransform(Astro, verbose = T, assay = 'Spatial',ncells = 3000)
neuro = SCTransform(neuro, verbose = T, assay = 'Spatial')
microG = SCTransform(microG, verbose = T, assay = 'Spatial')

##################################################################
######################################subCluster of microglia 
##################################################################

microG <- RunPCA(microG, verbose = FALSE,npcs = 50)
ElbowPlot(microG)
microG <- RunUMAP(microG, reduction = "pca", dims = 1:14)
microG <- FindNeighbors(microG, reduction = "pca", dims = 1:14)
microG <- FindClusters(microG, verbose = FALSE,resolution = 0.4)
#Astro@meta.data$s4cluster = AstroN$seurat_clusters
DimPlot(microG,reduction = 'umap',pt.size = 0.8,cols = mycol)


astroSpa = SpatialDimPlot(Astro)$data
mytype = t(matrix(unlist(strsplit(rownames(astroSpa),'_')),3,nrow(astroSpa)))
astroSpa$condition = mytype[,1]
astroSpa = astroSpa[,c(1,2,7,8)]

par(mar=c(2,2,1,1),mfrow=c(2,4))
for(i in 0:7){
  plot(astroSpa[astroSpa$condition=='30d1',1:2],pch=19,cex=.6,type='n',axes=F)
  points(astroSpa[astroSpa$condition=='30d1' & astroSpa$ident==i,1:2],pch=19,cex=.5,col=mycol[2])
}





par(mar=c(2,2,1,1))
plot(pos_30d1[,2:3],pch=19,cex=.1,col='gray80',main='',type='n')
points(testMicroG[testMicroG$seurat_clusters==0,1:2],pch=19,cex=.4,col=mycol[1])
points(testMicroG[testMicroG$seurat_clusters==1,1:2],pch=19,cex=.4,col=mycol[2])

plot(pos_cnr1[,2:3],pch=19,cex=.1,col='gray80',main='',type='n')
points(testMicroG[testMicroG$seurat_clusters==0,1:2],pch=19,cex=.4,col=mycol[1])
points(testMicroG[testMicroG$seurat_clusters==1,1:2],pch=19,cex=.4,col=mycol[2])


markersmicroG = FindAllMarkers(microG,only.pos = T,logfc.threshold = 0.15)
markersOther  = FindAllMarkers(otherSeurat,only.pos = T)

library(slingshot)
library(tradeSeq)
library(SingleCellExperiment)
library(SummarizedExperiment)
mgfb <- as.SingleCellExperiment(microG, assay = "Spatial")
mgfb <- slingshot(mgfb, reducedDim = 'UMAP')

set.seed(820)
library(condiments)
library(ggplot2)
library(dplyr)
dfm <- bind_cols(
  as.data.frame(reducedDims(mgfb)$UMAP),
  pseudotime = mgfb$slingPseudotime_1,
  treatment = mgfb$orig.ident)
dfm$treatment = factor(dfm$treatment)

curve <- slingCurves(mgfb)[[1]]
p4 <- ggplot(dfm, aes(x = UMAP_1, y = UMAP_2, col = pseudotime)) +
  geom_point(size = .7)+
  scale_color_viridis_c() +
  labs(col = "pseudotime") +
  geom_path(data = curve$s[curve$ord, ] %>% as.data.frame(),
            col = "black", size = 1.5)+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p4



scores <- condiments::imbalance_score(
  Object = dfm %>% select(UMAP_1, UMAP_2) %>% as.matrix(), 
  conditions = dfm$treatment,
  k = 20, smooth = 40)

dfm$scores <- scores$scaled_scores
p3 <- ggplot(dfm, aes(x = UMAP_1, y = UMAP_2, col = scores)) +
  geom_point(size = .7) +
  scale_color_viridis_c(option = "C") +
  labs(col = "Scores")+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p3


dfm$condition = 1
dfm[dfm$treatment=='30d1',]$condition = 'P30d'
dfm[dfm$treatment=='CNR1',]$condition = 'CNR1'

p5 <- ggplot(dfm, aes(x = pseudotime)) +
  geom_density(alpha = .8, aes(fill = treatment), col = "transparent") +
  geom_density(aes(col = treatment), fill = "transparent",
               guide = FALSE, size = 1.5) +
  labs(x = "Pseudotime", fill = "Treatment") +
  guides(col = FALSE, fill = guide_legend(
    override.aes = list(size = 1.5, col = c("#FF66B2", "#6666FF"))
  )) +
  scale_fill_brewer(palette = "Accent") +
  scale_color_brewer(palette = "Accent")+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

p5


ggplot(dfm, aes(x = pseudotime,fill=condition)) +geom_density(alpha=0.4)+
  geom_density(aes(col = condition), fill = "transparent", size = 1)+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


####compare the lineage
ks.test(dfm[dfm$treatment=='30d1',]$pseudotime,dfm[dfm$treatment=='CNR1',]$pseudotime)


pseuName = rownames(dfm[!(is.na(dfm$pseudotime)),])

mgfbMatrix = mgfb[,pseuName]

geneSum = rowSums(counts)
geneSum = geneSum[order(geneSum,decreasing = T)]
choose = names(geneSum[1:1000])
numOrder = which(rownames(counts) %in% choose)

mgfb <- fitGAM(counts = as.matrix(assays(mgfbMatrix)$counts),
               pseudotime =dfm[pseuName,]$pseudotime,
               cellWeights = mgfbMatrix$slingshot@assays@data$weights[,1],
               conditions = factor(dfm[pseuName,]$treatment),
               parallel = T,
               genes=numOrder,
               nknots = 5)

mean(rowData(mgfb)$tradeSeq$converged)



condRes <- conditionTest(mgfb, l2fc = 0.1)
condRes$padj <- p.adjust(condRes$pvalue, "fdr")
conditionGenes <- rownames(condRes)[condRes$padj <= 0.05]
conditionGenes <- conditionGenes[!is.na(conditionGenes)]
oo <- order(condRes$waldStat, decreasing = TRUE)
condRes[order(condRes$waldStat, decreasing = TRUE),]

p7 <- plotSmoothers(mgfb, assays(mgfb)$counts,
                    gene = 'C1qb',lwd=1,
                    alpha = 1, border = TRUE, curvesCols = c("#E69F00", "#999999")) +
  scale_color_manual(values = c("#E69F00", "#999999")) +
  ggtitle('C1qb')

p7





















##################################################################
######################################subCluster of microglia 
##################################################################
Astro <- RunPCA(Astro, verbose = FALSE,npcs = 50)
ElbowPlot(Astro)
Astro <- RunUMAP(Astro, reduction = "pca", dims = 1:10)
Astro <- FindNeighbors(Astro, reduction = "pca", dims = 1:10)
Astro <- FindClusters(Astro, verbose = FALSE,resolution = 0.4)
#Astro@meta.data$s4cluster = AstroN$seurat_clusters
DimPlot(Astro,reduction = 'umap',split.by = 'orig.ident',pt.size = 0.3,cols = mycol)
FeaturePlot(Astro,features = c('Aldoc','Cst3',
                               'Agt','Sparc',
                               'Apoe','Sparcl1',
                               'Apod','Ramp1',
                               'Nkain4','Cldn10','Tsc22d4','Ptgds'),pt.size = 0.3)

FeaturePlot(Astro,features = c('Ckb'),pt.size = 0.3)
FeaturePlot(Astro,features = c('Cst3'),pt.size = 0.3)
FeaturePlot(Astro,features = c('Agt'),pt.size = 0.3)
FeaturePlot(Astro,features = c('Sparc'),pt.size = 0.3)
FeaturePlot(Astro,features = c('Apoe'),pt.size = 0.3)
FeaturePlot(Astro,features = c('Sparcl1'),pt.size = 0.3)
FeaturePlot(Astro,features = c('Apod'),pt.size = 0.3)
FeaturePlot(Astro,features = c('Ramp1'),pt.size = 0.3)
FeaturePlot(Astro,features = c('Nkain4'),pt.size = 0.3)
FeaturePlot(Astro,features = c('Cldn10'),pt.size = 0.3)
FeaturePlot(Astro,features = c('Tsc22d4'),pt.size = 0.3)
FeaturePlot(Astro,features = c('Ptgds'),pt.size = 0.3)
FeaturePlot(Astro,features = c('Atp6v0e'),pt.size = 0.3)



DimPlot(otherSeurat,group.by = 's4cluster',split.by = 'orig.ident',pt.size = 0.5)
markersAstro = FindAllMarkers(Astro,only.pos = T)



Astro@meta.data$astroPseudotime = df$pseudotime
Astro@meta.data$phase = 0
Astro@meta.data[Astro@meta.data$astroPseudotime>=2 & Astro@meta.data$astroPseudotime<7,]$phase = 1
Astro@meta.data[Astro@meta.data$astroPseudotime>=7 & Astro@meta.data$astroPseudotime<11.5,]$phase = 2
Astro@meta.data[Astro@meta.data$astroPseudotime>=11.5 & Astro@meta.data$astroPseudotime<14.5,]$phase = 3
Astro@meta.data[Astro@meta.data$astroPseudotime>=14.5,]$phase = 4


markersPhaseAstro = FindAllMarkers(Astro,only.pos = T,group.by ='phase',logfc.threshold = 0.2)

spacExp = Astro@assays$SCT[c('Sparc','Sparcl1'),]
spacExp = t(as.matrix(spacExp))
spacExp = as.data.frame(spacExp)
spacExp$time = Astro@meta.data$astroPseudotime
spacExp$ratio = (spacExp$Sparcl1+1)/(spacExp$Sparc+1)

###########@@@@@@@@@@@@@@@@@@@@@@@@
###########@@@@@@@@@@@@@@@@@@@@@@@@
neuro <- RunPCA(neuro, verbose = FALSE,npcs = 50)
ElbowPlot(neuro)
neuro <- RunUMAP(neuro, reduction = "pca", dims = 1:16)
neuro <- FindNeighbors(neuro, reduction = "pca", dims = 1:16)
neuro <- FindClusters(neuro, verbose = FALSE,resolution = 0.6)
#neuro@meta.data$s4cluster = neuroN$seurat_clusters
DimPlot(neuro,split.by = 'orig.ident')
markersneuro = FindAllMarkers(neuro,only.pos = T)



library(slingshot)
library(tradeSeq)
library(SingleCellExperiment)
library(SummarizedExperiment)
tgfb <- as.SingleCellExperiment(Astro, assay = "Spatial")
tgfb <- slingshot(tgfb, reducedDim = 'UMAP')

set.seed(821)
library(condiments)
library(ggplot2)
library(dplyr)
df <- bind_cols(
  as.data.frame(reducedDims(tgfb)$UMAP),
  pseudotime = tgfb$slingPseudotime_1,
  treatment = tgfb$orig.ident)
df$treatment = factor(df$treatment)

curve <- slingCurves(tgfb)[[1]]
p4 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = pseudotime)) +
  geom_point(size = .7)+
  scale_color_viridis_c() +
  labs(col = "pseudotime") +
  geom_path(data = curve$s[curve$ord, ] %>% as.data.frame(),
            col = "black", size = 1.5)+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p4



scores <- condiments::imbalance_score(
  Object = df %>% select(UMAP_1, UMAP_2) %>% as.matrix(), 
  conditions = df$treatment,
  k = 20, smooth = 40)

df$scores <- scores$scaled_scores
p3 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = scores)) +
  geom_point(size = .7) +
  scale_color_viridis_c(option = "C") +
  labs(col = "Scores")+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p3


df$condition = 1
df[df$treatment=='30d1',]$condition = 'P30d'
df[df$treatment=='CNR1',]$condition = 'CNR1'

p5 <- ggplot(df, aes(x = pseudotime)) +
  geom_density(alpha = .8, aes(fill = treatment), col = "transparent") +
  geom_density(aes(col = treatment), fill = "transparent",
               guide = FALSE, size = 1.5) +
  labs(x = "Pseudotime", fill = "Treatment") +
  guides(col = FALSE, fill = guide_legend(
    override.aes = list(size = 1.5, col = c("#FF66B2", "#6666FF"))
  )) +
  scale_fill_brewer(palette = "Accent") +
  scale_color_brewer(palette = "Accent")+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

p5


ggplot(df, aes(x = pseudotime,fill=condition)) +geom_density(alpha=0.4)+
  geom_density(aes(col = condition), fill = "transparent", size = 1)+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
                                                                                                      

####compare the lineage
ks.test(df[df$treatment=='30d1',]$pseudotime,df[df$treatment=='CNR1',]$pseudotime)



library(tradeSeq)
set.seed(5)
icMat <- evaluateK(counts = as.matrix(assays(tgfbMatrix)$counts),
                   pseudotime =df[pseuName,]$pseudotime,
                   cellWeights = tgfbMatrix$slingshot@assays@data$weights[,1],
                   conditions = factor(df[pseuName,]$treatment),
                   nGenes = 200,parallel=T,
                   k = 3:15,plot = TRUE)




pseuName = rownames(df[!(is.na(df$pseudotime)),])

tgfbMatrix = tgfb[,pseuName]

geneSum = rowSums(counts)
geneSum = geneSum[order(geneSum,decreasing = T)]
choose = names(geneSum[1:3000])
numOrder = which(rownames(counts) %in% choose)

tgfb <- fitGAM(counts = as.matrix(assays(tgfbMatrix)$counts),
               pseudotime =df[pseuName,]$pseudotime,
               cellWeights = tgfbMatrix$slingshot@assays@data$weights[,1],
               conditions = factor(df[pseuName,]$treatment),
               parallel = T,
               genes=numOrder,
               nknots = 5)

mean(rowData(tgfb)$tradeSeq$converged)



condRes <- conditionTest(tgfb, l2fc = 0.1)
condRes$padj <- p.adjust(condRes$pvalue, "fdr")
conditionGenes <- rownames(condRes)[condRes$padj <= 0.05]
conditionGenes <- conditionGenes[!is.na(conditionGenes)]
oo <- order(condRes$waldStat, decreasing = TRUE)
condRes[order(condRes$waldStat, decreasing = TRUE),]

library(RColorBrewer)
scales <- brewer.pal(3, "Accent")[1:2]


for(i in 1:200){
  png(paste('~/Downloads/astro/knot8',i,'.png',sep=''),width = 600,height = 400)
  p6 <- plotSmoothers(tgfb, assays(tgfb)$counts,
                    gene = rownames(assays(tgfb)$counts)[oo[i]],
                    alpha = 1, border = TRUE, curvesCols = scales) +
  scale_color_manual(values = scales) +
  ggtitle(rownames(assays(tgfb)$counts)[oo[i]])
  
  print(p6)
  dev.off()
}

p7 <- plotSmoothers(tgfb, assays(tgfb)$counts,
                    gene = 'Id3',lwd=1,
                    alpha = 1, border = TRUE, curvesCols = c("#E69F00", "#999999")) +
  scale_color_manual(values = c("#E69F00", "#999999")) +
  ggtitle('Id3')

p7


p7 <- plotSmoothers(tgfb, assays(tgfb)$counts,
                    gene = 'Sparc',lwd=1,
                    alpha = 1, border = TRUE, curvesCols = c("#E69F00", "#999999")) +
  scale_color_manual(values = c("#E69F00", "#999999")) +
  ggtitle('Sparc')

p7

p7 <- plotSmoothers(tgfb, assays(tgfb)$counts,
                    gene = 'Apoe',lwd=1,
                    alpha = 1, border = TRUE, curvesCols = c("#E69F00", "#999999")) +
  scale_color_manual(values = c("#E69F00", "#999999")) +
  ggtitle('Apoe')

p7

p7 <- plotSmoothers(tgfb, assays(tgfb)$counts,
                    gene = 'Cldn10',lwd=1,
                    alpha = 1, border = TRUE, curvesCols = c("#E69F00", "#999999")) +
  scale_color_manual(values = c("#E69F00", "#999999")) +
  ggtitle('Cldn10')

p7

FeaturePlot(Astro,features = c('C1qb'),pt.size = 0.3)
p7 <- plotSmoothers(tgfb, assays(tgfb)$counts,
                    gene = 'C1qb',lwd=1,
                    alpha = 1, border = TRUE, curvesCols = c("#E69F00", "#999999")) +
  scale_color_manual(values = c("#E69F00", "#999999")) +
  ggtitle('C1qb')

p7

p7 <- plotSmoothers(tgfb, assays(tgfb)$counts,
                    gene = 'Pla2g7',lwd=1,
                    alpha = 1, border = TRUE, curvesCols = c("#E69F00", "#999999")) +
  scale_color_manual(values = c("#E69F00", "#999999")) +
  ggtitle('Pla2g7')

p7



write.csv(condRes,'~/lab508/novaSeqPBN/communication/astro_pseudotime.csv')


library(cowplot)
library(scales)
### based on mean smoother
yhatSmooth <- predictSmooth(tgfb, gene = conditionGenes, nPoints = 50, tidy = FALSE) %>%
  log1p()
yhatSmoothScaled <- t(apply(yhatSmooth,1, scales::rescale))
library(pheatmap)
heatSmooth_TGF <- pheatmap(yhatSmoothScaled[, 51:100],
                           cluster_cols = FALSE,
                           show_rownames = FALSE, show_colnames = FALSE, main = "CNR1", legend = FALSE,
                           silent = TRUE
)

matchingHeatmap_mock <- pheatmap(yhatSmoothScaled[heatSmooth_TGF$tree_row$order, 1:50],
                                 cluster_cols = FALSE, cluster_rows = FALSE,
                                 show_rownames = FALSE, show_colnames = FALSE, main = "30d1",
                                 legend = FALSE, silent = TRUE 
)

p9 <- plot_grid(heatSmooth_TGF[[4]], matchingHeatmap_mock[[4]], ncol = 2)
p9










#type spatial distribution
pos_cnr1$astroPseudotime = NA
nameList = intersect(rownames(df),rownames(pos_cnr1))
pos_cnr1[nameList,]$astroPseudotime = df[nameList,]$pseudotime

pos_30d1$astroPseudotime = NA
nameList = intersect(rownames(df),rownames(pos_30d1))
pos_30d1[nameList,]$astroPseudotime = df[nameList,]$pseudotime



ggplot(pos_cnr1, aes(x = xcoord, y = ycoord, col = astroPseudotime)) +
  geom_point(size = .2)+
  scale_color_viridis_c(na.value = 'gray80') +
  labs(col = "astroPseudotime")+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggplot(pos_30d1, aes(x = xcoord, y = ycoord, col = astroPseudotime)) +
  geom_point(size = 0.2)+
  scale_color_viridis_c(na.value = 'gray80') +
  labs(col = "astroPseudotime")+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

par(mfrow=c(3,2),bg='black')
for(i in 0:7){
  plot(pos_cnr1[,2:3],pch=19,cex=.2,col='gray80',main = paste('cluster',i))
  points(pos_cnr1[rownames(Astro@meta.data[Astro@meta.data$seurat_clusters==i,]),2:3],pch=19,cex=.4,col=mycol[i+1])
  
  plot(pos_30d1[,2:3],pch=19,cex=.2,col='gray80', main = paste('cluster',i) )
  points(pos_30d1[rownames(Astro@meta.data[Astro@meta.data$seurat_clusters==i,]),2:3],pch=19,cex=.4,col=mycol[i+1])
}

#three states
#2-7
#7-11.5
#11.5-14.5
colset <- brewer.pal(n = 8, name ="Set1")
par(mar=c(2,2,1,1))
plot(pos_cnr1[,2:3],pch=19,cex=.1,col='gray80',main='A')
points(pos_cnr1[pos_cnr1$astroPseudotime>=2 & pos_cnr1$astroPseudotime<7,2:3],pch=19,cex=.4,col=colset[2])

plot(pos_30d1[,2:3],pch=19,cex=.1,col='gray80',main='A',type='n')
points(pos_30d1[pos_30d1$astroPseudotime>=2 & pos_30d1$astroPseudotime<7,2:3],pch=19,cex=.4,col=colset[1])


plot(pos_cnr1[,2:3],pch=19,cex=.1,col='gray80',main='B',type='n')
points(pos_cnr1[pos_cnr1$astroPseudotime>=7 & pos_cnr1$astroPseudotime<11.5,2:3],pch=19,cex=.4,col=colset[2])

plot(pos_30d1[,2:3],pch=19,cex=.1,col='gray80',main='B',type='n')
points(pos_30d1[pos_30d1$astroPseudotime>=7 & pos_30d1$astroPseudotime<11.5,2:3],pch=19,cex=.4,col=colset[1])


plot(pos_cnr1[,2:3],pch=19,cex=.2,col='gray80',main='C',type='n')
points(pos_cnr1[pos_cnr1$astroPseudotime>=11.5 & pos_cnr1$astroPseudotime<14.5,2:3],pch=19,cex=.4,col=colset[2])

plot(pos_30d1[,2:3],pch=19,cex=.2,col='gray80',main='C',type='n')
points(pos_30d1[pos_30d1$astroPseudotime>=11.5 & pos_30d1$astroPseudotime<14.5,2:3],pch=19,cex=.4,col=colset[1])



cnr1_A = pos_cnr1[!(is.na(pos_cnr1$astroPseudotime)),]
d301_A = pos_30d1[!(is.na(pos_30d1$astroPseudotime)),]
selectA_C = rownames(cnr1_A[cnr1_A$astroPseudotime>=2 & cnr1_A$astroPseudotime<7,])
selectA_T = rownames(d301_A[d301_A$astroPseudotime>=2 & d301_A$astroPseudotime<7,])
selectB_C = rownames(cnr1_A[cnr1_A$astroPseudotime>=7 & cnr1_A$astroPseudotime<11.5,])
selectB_T = rownames(d301_A[d301_A$astroPseudotime>=7 & d301_A$astroPseudotime<11.5,])
selectC_C = rownames(cnr1_A[cnr1_A$astroPseudotime>=11.5 & cnr1_A$astroPseudotime<14.5,])
selectC_T = rownames(d301_A[d301_A$astroPseudotime>=11.5 & d301_A$astroPseudotime<14.5,])

length(intersect(selectA_C, regionPBN_CN))/nrow(cnr1_A)
length(intersect(selectA_T, regionPBN_30d1))/nrow(d301_A)

length(intersect(selectB_C, regionPBN_CN))/nrow(cnr1_A)
length(intersect(selectB_T, regionPBN_30d1))/nrow(d301_A)

length(intersect(selectC_C, regionPBN_CN))/nrow(cnr1_A)
length(intersect(selectC_T, regionPBN_30d1))/nrow(d301_A)




apoeComm_C = read.table('~/lab508/novaSeqPBN/communication/Coord_CNR1_clean.csv',header=T)
apoeComm_T = read.table('~/lab508/novaSeqPBN/communication/Coord_30d1_clean.csv',header=T)
rownames(apoeComm_C) = apoeComm_C$Barcode
rownames(apoeComm_T) = apoeComm_T$Barcode

nameC = intersect(rownames(cnr1_A), apoeComm_C$Barcode)
nameT = intersect(rownames(d301_A), apoeComm_T$Barcode)

selectTimeC = cnr1_A[nameC,]
selectTimeT = d301_A[nameT,]

selectApoeAC = apoeComm_C[nameC,]
selectApoeAT = apoeComm_T[nameT,]

selectApoeAC$pseudotime = selectTimeC$astroPseudotime
selectApoeAT$pseudotime = selectTimeT$astroPseudotime

par(mfrow=c(2,1),mar=c(4,4,2,1))
plot(selectApoeAC$pseudotime,selectApoeAC$SApoe,pch=19,cex=.8,col=rgb(1,0,0,0.3),ylab='communication',xlab='Astro_pseudotime')
abline(v=c(2,7,11.5,14.5),lty=2)
plot(selectApoeAT$pseudotime,selectApoeAT$SApoe,pch=19,cex=.8,col=rgb(0,0,1,0.3),ylab='communication',xlab='Astro_pseudotime')
abline(v=c(2,7,11.5,14.5),lty=2)


smoothScatter(selectApoeAC$pseudotime,selectApoeAC$SApoe,pch=19,cex=.8,col=rgb(1,0,0,0.3),ylab='communication',xlab='Astro_pseudotime')
abline(v=c(2,7,11.5,14.5),lty=2)
smoothScatter(selectApoeAT$pseudotime,selectApoeAT$SApoe,pch=19,cex=.8,col=rgb(0,0,1,0.3),ylab='communication',xlab='Astro_pseudotime')
abline(v=c(2,7,11.5,14.5),lty=2)

selectApoe_AC = apoeComm_C[apoeComm_C$Barcode %in% pseuName & apoeComm_C$SApoe>1.6 & apoeComm_C$cluster==5,]$Barcode
selectApoe_AT = apoeComm_T[apoeComm_T$Barcode %in% pseuName & apoeComm_T$SApoe*2.4>1.6 & apoeComm_T$cluster==5,]$Barcode

par(mfrow=c(3,1),mar=c(4,4,2,1))
plot(density(selectApoeAC[selectApoeAC$pseudotime>2 & selectApoeAC$pseudotime<7,]$SApoe))
lines(density(selectApoeAT[selectApoeAT$pseudotime>2 & selectApoeAT$pseudotime<7,]$SApoe*2.4),col=2)

plot(density(selectApoeAC[selectApoeAC$pseudotime>7 & selectApoeAC$pseudotime<11.5,]$SApoe))
lines(density(selectApoeAT[selectApoeAT$pseudotime>7 & selectApoeAT$pseudotime<11.5,]$SApoe*2.4),col=2)

plot(density(selectApoeAC[selectApoeAC$pseudotime>11.5 & selectApoeAC$pseudotime<14.5,]$SApoe))
lines(density(selectApoeAT[selectApoeAT$pseudotime>11.5 & selectApoeAT$pseudotime<14.5,]$SApoe*2.4),col=2)




length(intersect(selectA_C,selectApoe_AC))/length(selectA_C)
length(intersect(selectB_C,selectApoe_AC))/length(selectB_C)
length(intersect(selectC_C,selectApoe_AC))/length(selectC_C)

length(intersect(selectA_T,selectApoe_AT))/length(selectA_T)
length(intersect(selectB_T,selectApoe_AT))/length(selectB_T)
length(intersect(selectC_T,selectApoe_AT))/length(selectC_T)


plot(density(apoeComm_C[apoeComm_C$cluster==4,]$SApoe),col=1,lty=2)






#CNR2/30d2 microglia

micGName_section2 = rownames(PBNct.combined.sct@meta.data[PBNct.combined.sct@meta.data$seurat_clusters==11 & PBNct.combined.sct@meta.data$orig.ident %in% c('30d2','CNR2'),])
microG2 = PBNct[,micGName_section2]

########################
### Normalize
########################
microG2 = SCTransform(microG2, verbose = T, assay = 'Spatial')
##################
#microg
############################

microG2 <- RunPCA(microG2, verbose = FALSE,npcs = 50)
ElbowPlot(microG2)
microG2 <- RunUMAP(microG2, reduction = "pca", dims = 1:14)
microG2 <- FindNeighbors(microG2, reduction = "pca", dims = 1:14)
microG2 <- FindClusters(microG2, verbose = FALSE,resolution = 0.4)
#Astro@meta.data$s4cluster = AstroN$seurat_clusters
DimPlot(microG2,reduction = 'umap',pt.size = 0.8,cols = mycol)

markerG2all = FindAllMarkers(microG2)
