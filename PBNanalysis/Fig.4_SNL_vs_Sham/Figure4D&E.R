########################################
#################S3 umap
########################################
neuronCell = PBNct.combined.sct@meta.data[PBNct.combined.sct@meta.data$seurat_clusters %in% c(0,1,5,7,10,14,15,21,8,16)
                                          & PBNct.combined.sct@meta.data$orig.ident %in% c('CNR2','30d2'),]

otherCell = PBNct.combined.sct@meta.data[PBNct.combined.sct@meta.data$seurat_clusters %in% c(2,3,4,6,9,11,12,13,17,18,19,20)
                                         & PBNct.combined.sct@meta.data$orig.ident %in% c('CNR2','30d2'),]

neuronSeurat = PBNct[,rownames(neuronCell)]
neuronSeurat@meta.data$s4Cluster = neuronCell$seurat_clusters

otherSeurat = PBNct[,rownames(otherCell)]
otherSeurat@meta.data$s4Cluster = otherCell$seurat_clusters


neuronSeurat <- SCTransform(neuronSeurat, assay = "Spatial", ncells = 3000, verbose = TRUE)
neuronSeurat <- RunPCA(neuronSeurat, verbose = FALSE)
ElbowPlot(neuronSeurat,ndims = 50)
neuronSeurat <- FindNeighbors(neuronSeurat, reduction = "pca", dims = 1:16)
neuronSeurat <- FindClusters(neuronSeurat, verbose = FALSE,resolution = 0.3)
neuronSeurat <- RunUMAP(neuronSeurat, reduction = "pca", dims = 1:16)
DimPlot(neuronSeurat, reduction = "umap",cols=mycol)

otherSeurat <- SCTransform(otherSeurat, assay = "Spatial", ncells = 3000, verbose = TRUE)
otherSeurat <- RunPCA(otherSeurat, verbose = FALSE)
ElbowPlot(otherSeurat,ndims = 50)
otherSeurat <- FindNeighbors(otherSeurat, reduction = "pca", dims = 1:12)
otherSeurat <- FindClusters(otherSeurat, verbose = FALSE,resolution = 0.4)
otherSeurat <- RunUMAP(otherSeurat, reduction = "pca", dims = 1:12)
DimPlot(otherSeurat, reduction = "umap",cols=mycol)

markersNeuron = FindAllMarkers(neuronSeurat,only.pos = T)
markersOther = FindAllMarkers(otherSeurat,only.pos = T,logfc.threshold = 0.2)






###############################
#########Figure4 
###############################
#-----------------------------------------------------------------------------------------------neuron/other split clustering

neuronCell = PBNct.combined.sct@meta.data[PBNct.combined.sct@meta.data$seurat_clusters %in% c(0,1,5,7,10,14,15,21,8,16)
                                          & PBNct.combined.sct@meta.data$orig.ident %in% c('CNR1','30d1'),]

otherCell = PBNct.combined.sct@meta.data[PBNct.combined.sct@meta.data$seurat_clusters %in% c(2,3,4,6,9,11,12,13,17,18,19,20)
                                         & PBNct.combined.sct@meta.data$orig.ident %in% c('CNR1','30d1'),]

S2Cell = PBNct.combined.sct@meta.data[PBNct.combined.sct@meta.data$orig.ident %in% c('CNR1','30d1'),]

neuronSeurat = PBNct[,rownames(neuronCell)]
neuronSeurat@meta.data$s4Cluster = neuronCell$seurat_clusters

otherSeurat = PBNct[,rownames(otherCell)]
otherSeurat@meta.data$s4Cluster = otherCell$seurat_clusters

PBNS2 = PBNct[,rownames(S2Cell)]
PBNS2 <- SCTransform(PBNS2, assay = "Spatial", ncells = 3000, verbose = TRUE)

neuronSeurat <- SCTransform(neuronSeurat, assay = "Spatial", ncells = 3000, verbose = TRUE)
neuronSeurat <- RunPCA(neuronSeurat, verbose = FALSE)
ElbowPlot(neuronSeurat,ndims = 50)
neuronSeurat <- FindNeighbors(neuronSeurat, reduction = "pca", dims = 1:18)
neuronSeurat <- FindClusters(neuronSeurat, verbose = FALSE,resolution = 0.7)
neuronSeurat <- RunUMAP(neuronSeurat, reduction = "pca", dims = 1:18)
DimPlot(neuronSeurat, reduction = "umap",split.by = 'orig.ident',cols=mycol)

otherSeurat <- SCTransform(otherSeurat, assay = "Spatial", ncells = 3000, verbose = TRUE)
otherSeurat <- RunPCA(otherSeurat, verbose = FALSE)
ElbowPlot(otherSeurat,ndims = 50)
otherSeurat <- FindNeighbors(otherSeurat, reduction = "pca", dims = 1:22)
otherSeurat <- FindClusters(otherSeurat, verbose = FALSE,resolution = 0.4)
otherSeurat <- RunUMAP(otherSeurat, reduction = "pca", dims = 1:22)
DimPlot(otherSeurat, reduction = "umap",split.by = 'orig.ident',cols=mycol)


mypal <- brewer.pal(n = 12, name ="Paired")
mydal <- brewer.pal(n = 8, name ="Set2")
mycal <- brewer.pal(n = 8, name ="Set3")
mycol = c(mydal, mypal, mycal)

plot1otherSeurat = DimPlot(otherSeurat, reduction = "umap",cols = mycol,label = T)
plot1neuronSeurat = DimPlot(neuronSeurat, reduction = "umap",cols = mycol,label = T)

DefaultAssay(neuronSeurat)='SCT'
DefaultAssay(otherSeurat)='SCT'
FeaturePlot(neuronSeurat,features = c('Calca','Tac1','Nts','Sncg','Apoe','Penk'),min.cutoff = 0,ncol = 3)
FeaturePlot(otherSeurat,features = c('Calca','Tac1','Nts','Sncg','Apoe','Penk'),min.cutoff = 0,ncol = 3)

markersNeuron = FindAllMarkers(neuronSeurat,only.pos = T)
markersOther = FindAllMarkers(otherSeurat,only.pos = T)

otherSeurat = readRDS('~/lab508/novaSeqPBN/otherSeurat.rds')
neuronSeurat = readRDS('~/lab508/novaSeqPBN/neuronSeurat.rds')
write.csv(markersNeuron,'~/lab508/novaSeqPBN/PBNall/S2neuronMarkers.list')
write.csv(markersOther,'~/lab508/novaSeqPBN/PBNall/S2othersMarkers.list')

umappingOth = plot1otherSeurat$data
mytype = t(matrix(unlist(strsplit(rownames(umappingOth),'_')),3,nrow(umappingOth)))
umappingOth$condition = mytype[,1]
umappingOth$binid = mytype[,3]
rownames(positions) = positions$binID
posall = positions[rownames(umappingOth),]
umappingOth$x = posall$xcoord
umappingOth$y = posall$ycoord
umappingNeu = plot1neuronSeurat$data
mytype = t(matrix(unlist(strsplit(rownames(umappingNeu),'_')),3,nrow(umappingNeu)))
umappingNeu$condition = mytype[,1]
umappingNeu$binid = mytype[,3]
rownames(positions) = positions$binID
posall = positions[rownames(umappingNeu),]
umappingNeu$x = posall$xcoord
umappingNeu$y = posall$ycoord
umappingNeu$type = 'Neu'
umappingOth$type = 'Oth'
umapping = rbind(umappingNeu,umappingOth)

par(mfrow=c(1,2))
selcluster = umapping[umapping$condition=='30d1',]
plot(selcluster$x,selcluster$y,pch=19,cex=.2,col='grey80')
selcluster = umapping[umapping$condition=='30d1' & umapping$ident==1 & umapping$type=='Oth',]
points(selcluster$x,selcluster$y,pch=19,cex=.3,col=3)
selcluster = umapping[umapping$condition=='30d1' & umapping$ident==9 & umapping$type=='Oth',]
points(selcluster$x,selcluster$y,pch=19,cex=.3,col=2)


allcluster = umapping[umapping$condition=='30d1' & umapping$type=='Oth',]
neuorder = c(7,0,2,6,8,10,11,15,12,9,14,5,3,13,1,4)
othorder = c(0,3,9,1,10,11,6,4,2,5,8,7)
par(mfrow=c(7,4),mar=c(1,1,1,1))
for(i in othorder){
  selcluster = umapping[umapping$ident==i & umapping$condition=='30d1' & umapping$type=='Oth',]
  if(nrow(selcluster)>0){
    plot(allcluster$x,allcluster$y,main="",axes = F,type='n')
    points(selcluster$x,selcluster$y,pch=19,cex=.3,col=mydal[2])
  }
}








#------------------------
markersPerNeuCluster = FindMarkers(neuronSeurat, ident.1 = "30d1", ident.2 = "CNR1", group.by = 'orig.ident', subset.ident = 0,assay='SCT',logfc.threshold = 0)
markersPerNeuCluster$gene = rownames(markersPerNeuCluster)
markersPerNeuCluster$cluster = 0
for(i in 1:15){
  markers <- FindMarkers(neuronSeurat, ident.1 = "30d1", ident.2 = "CNR1", group.by = 'orig.ident', subset.ident = i,assay='SCT',logfc.threshold = 0)
  markers$gene = rownames(markers)
  markers$cluster = i
  markersPerNeuCluster = rbind(markersPerNeuCluster,markers)
}
write.csv(markersPerNeuCluster[markersPerNeuCluster$p_val<0.05,],'~/Dropbox/novaSeq/figures/PBN_pain/DiffList/S2_Neuron_diff.txt')

#-----------------------
markersPerOthCluster = FindMarkers(otherSeurat, ident.1 = "30d1", ident.2 = "CNR1", group.by = 'orig.ident', subset.ident = 0,assay='SCT',logfc.threshold = 0)
markersPerOthCluster$gene = rownames(markersPerOthCluster)
markersPerOthCluster$cluster = 0
for(i in 1:11){
  markers <- FindMarkers(otherSeurat, ident.1 = "30d1", ident.2 = "CNR1", group.by = 'orig.ident', subset.ident = i,assay='SCT',logfc.threshold = 0)
  markers$gene = rownames(markers)
  markers$cluster = i
  markersPerOthCluster = rbind(markersPerOthCluster,markers)
}

write.csv(markersPerOthCluster[markersPerOthCluster$p_val<0.05,],'~/Dropbox/novaSeq/figures/PBN_pain/DiffList/S2_Other_diff.txt')

