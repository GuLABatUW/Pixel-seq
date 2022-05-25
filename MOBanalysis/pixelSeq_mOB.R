###################################################################################################
###################################################################spatial heatmap of mouse OB data, Fig2A 
###################################################################################################
library(ggplot2)
library(viridis)
bgn36 = read.table('~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/OB36uniq.coord',header=T)
bgn36Clean = bgn36[bgn36$CoordX>24000 & bgn36$CoordX<40000,]
tiff(filename = paste("~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/","Fig2A",".tiff",sep=''),width = 3200, height = 2400, units = "px", pointsize = 12,compression="lzw");
par(bg='black')
p = ggplot(bgn36Clean, aes(x=CoordX, y=CoordY,colour=log2(UMIcounts))) + geom_point(size=0.2)
p + scale_fill_gradientn(colours = rev(viridis(256, option = "A")),aesthetics = "colour")+ theme_void()
dev.off()


####################################################################################################
#################################################################method comparison, FigS7B
####################################################################################################
hdst= read.table('~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/methodCompare/hdst/hdstraw.txt',header=F,skip = 1)

pixel2um_OB36 = read.csv('~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/methodCompare/pixelSeq/binStat_OB36_7pixelCoord.txt',header=T)
pixel2umClean = pixel2um_OB36[pixel2um_OB36$xcoord>24000 & pixel2um_OB36$xcoord<40000 & pixel2um_OB36$UMIcounts>4,]

summary(pixel2umClean$UMIcounts)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#5.00   16.00   33.00   46.81   63.00 1979.00


#####2um feature comparison " HDST"," Stereo-seq","PIXEL-seq"
mypal <- brewer.pal(n = 8, name ="Dark2")
df<-data.frame(lower=c(1,4.3,16),upper=c(4,10,63),middle=c(2,7.3,46.81),ymin=c(1,1,5),ymax=c(8,20,150)
               ,Category=as.factor(c(" HDST"," Stereo-seq","PIXEL-seq")))
ggplot(df,aes(x=Category,color=Category))+geom_boxplot(aes(lower=lower,upper=upper,middle=middle,ymin=ymin,ymax=ymax),stat="identity",width=0.5,alpha=0.7)+theme_bw() + scale_color_brewer(palette="Dark2")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")


#####10um feature comparison " HDST"," Stereo-seq","slide-seqv2", "PIXEL-seq"
#Stereo-seq no data, take it as x25 from 4um^2 data

#hdst
hdst5bin = read.csv('~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/methodCompare/hdst/HDSTBin5Coord.txt',header=T)
summary(hdst5bin$UMIcounts)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.00    4.00    9.00   12.07   16.00   97.00 

#slide-seqV2
slidev2 = read.csv('~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/methodCompare/slideseqV2/slideBin1Coord.txt',header=T)
plot(slidev2$xcoord,slidev2$ycoord,cex=.1,col=rgb(1,0,0,slidev2$UMIcounts/21000*0.99))
points(slidev2$xcoord,slidev2$ycoord,cex=.1,col=rgb(1,0,0,slidev2$UMIcounts/21000*0.99))
summary(slidev2$UMIcounts)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#10.0    67.0   239.0   465.2   625.0 18567.0 

#summary(slidev2$GENEcounts)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#6.0    58.0   194.0   344.3   489.0 10935.0 

p = ggplot(slidev2, aes(x=xcoord, y=ycoord,colour=UMIcounts/GENEcounts)) + geom_point(size=0.2)
p + scale_fill_gradientn(colours = rev(viridis(256, option = "A")),aesthetics = "colour")+ theme_void()


#pixel-seq
pixel10um_OB36 = read.csv('~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/methodCompare/pixelSeq/binStat_OB36_33pixelCoord.txt',header=T)
pixel10umClean = pixel10um_OB36[pixel10um_OB36$xcoord>24000 & pixel10um_OB36$xcoord<40000 & pixel10um_OB36$UMIcounts>256,]
summary(pixel10umClean$UMIcounts)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#257.0   551.0   831.0   976.8  1241.0  9595.0 


p = ggplot(feature10umClean, aes(x=xcoord, y=ycoord,colour=UMIcounts/GENEcounts)) + geom_point(size=0.2)
p + scale_fill_gradientn(colours = rev(viridis(256, option = "A")),aesthetics = "colour")+ theme_void()


df<-data.frame(lower=c(4,67,551),upper=c(16,625,1241),middle=c(12.07,465.2,976.8),ymin=c(1,10,257),ymax=c(97,2500,2500)
               ,Category=as.factor(c(" HDST"," slide-seqV2","PIXEL-seq")))
ggplot(df,aes(x=Category,color=Category))+geom_boxplot(aes(lower=lower,upper=upper,middle=middle,ymin=ymin,ymax=ymax),stat="identity",width=0.5,alpha=0.7)+theme_bw() + scale_color_brewer(palette="Dark2",direction = 2)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")


#####50um reature comparison " HDST"," Stereo-seq","slide-seqv2", "PIXEL-seq"
#hdst
hdst25bin = read.csv('~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/methodCompare/hdst/HDSTBin25Coord.txt',header=T)
summary(hdst25bin$UMIcounts)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1      98     182     204     277     821 

#visium
visium = read.csv('~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/methodCompare/visium/tissue_positions_list.csv',header=F)
h5 = Read10X_h5('~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/methodCompare/visium/filtered_feature_bc_matrix.h5', use.names = TRUE, unique.features = TRUE)
summary(h5@i)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0    7639   15392   15447   22927   31052 


#slide-seqv2
slidebin80v2 = read.csv('~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/methodCompare/slideseqV2/slideBin80Coord.txt',header=T)
slidebin80v2Clean = slidebin80v2[slidebin80v2$UMIcounts>256,]
summary(slidebin80v2Clean$UMIcounts)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#257    2067    3824    4526    6228   18740 

#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#206    1295    2186    2285    3149   10935 

#pixel-seq
pixel50um_OB36 = read.csv('~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/pixelSeq/binStat_OB36154Coord.txt',header=T)
pixel50umClean = pixel50um_OB36[pixel50um_OB36$xcoord>24000 & pixel50um_OB36$xcoord<40000 & pixel50um_OB36$UMIcounts>4096,]
summary(pixel50umClean$UMIcounts)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#4107   12506   18537   19975   25489   62336 

df<-data.frame(ymin=c(1,1,257,4107),lower=c(98,7639,2067,12506),middle=c(204,15447,4526,19975),upper=c(277,22927,6228,25489),ymax=c(821,31052,18740,62336)
               ,Category=as.factor(c(" HDST"," Visium"," Slide-seqV2","PIXEL-seq")))
ggplot(df,aes(x=Category,color=Category))+geom_boxplot(aes(lower=lower,upper=upper,middle=middle,ymin=ymin,ymax=ymax),stat="identity",width=0.5,alpha=0.7)+theme_bw() + scale_color_brewer(palette="Dark2")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position="none")




#################################################################################################################
#####################################################################spatial heatmap of region specific genes, Fig.S7C
#################################################################################################################
myorg <- brewer.pal(n = 9, name ="YlOrBr")
doc2g = read.table('~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/pattern/doc2g.final',header=F)
par(bg=myorg[1])#,mar=c(3,3,1,1))
plot(doc2g[,3:4],pch=19,cex=.1,col=rgb(1,0,0,0.1),xlim=c(24000,40000),axes=F)
box(col='black')

penk = read.table('~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/pattern/penk.final',header=F)
par(bg=myorg[1])
plot(penk[,3:4],pch=19,cex=.1,col=rgb(1,0,0,0.1),xlim=c(24000,40000),axes=F)
box(col='black')

Fabp7 = read.table('~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/pattern/fabp7.final',header=F)
par(bg=myorg[1])
plot(Fabp7[,3:4],pch=19,cex=.1,col=rgb(1,0,0,0.01),xlim=c(24000,40000),axes=F)
box(col='black')

ptgds = read.table('~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/pattern/ptgds.final',header=F)
par(bg=myorg[1])
plot(ptgds[,3:4],pch=19,cex=.1,col=rgb(1,0,0,0.01),xlim=c(24000,40000),axes=F)
box(col='black')


################################################################################################################
#####################################################################gene abundance correlation, Fig.S7A
################################################################################################################
olbExpress = read.table('~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/OBnovaSeq_geneAbundance.txt',header=T)
olbPexpress = olbExpress[apply(olbExpress[,2:4], 1, max)>=10,]

a = cbind(log2(olbPexpress$OB13+1), 25-log2(olbPexpress$OB25+1), rep(0,nrow(olbPexpress)))
b = cbind(log2(olbPexpress$OB13+1), rep(25,nrow(olbPexpress)), log2(olbPexpress$OB36+1))
c = cbind(rep(0,nrow(olbPexpress)), 25-log2(olbPexpress$OB25+1), log2(olbPexpress$OB36+1))
all = rbind(a,b,c)
library(scatterplot3d)
scatterplot3d(all,pch=19,cex.symbols=0.2,zlim=c(0,25), highlight.3d=TRUE,angle =30,label.tick.marks=F)





##################################################################################################################
#########################################################################Spatial pattern analysis at vSeg for OB36
##################################################################################################################
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

expr <- data.table::fread(
  input = "~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/segmentation/spt_3U_MaskCellGeneExp.txt",
  sep = ",",
  data.table = FALSE
)

rownames(x = expr) <- (x = expr[, 1])
colnames(x = expr) <- (x = colnames(x = expr))
expr <- t(x = expr[, -1])
ssmOLB <- CreateSeuratObject(
  counts = expr,
  project = 'SlideSeq',
  assay = 'Spatial'
)
positions <- read.csv(
  file = "~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/segmentation/spt_3U_MaskCoord.txt"
)

nn1 = positions
rownames(x = nn1) <- nn1$binID
positions <- nn1[, 2:3]

ssmOLB[['image']] <-new(
  Class = 'SlideSeq',
  assay = "Spatial",
  coordinates = positions
)

#output for scRNA-seq annotation
library(SeuratDisk)
SaveH5Seurat(ssmOLB, filename = "~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/scAnnotation/OB36_vSegN8003U.h5Seurat",overwrite =T)
Convert("~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/scAnnotation/OB36_vSegN8003U.h5Seurat", dest = "h5ad",overwrite =T)

ssmOLB$log_nCount_Spatial <- log2(ssmOLB$nCount_Spatial)
ssmOLB = ssmOLB[, ssmOLB$nCount_Spatial > 512]
ssmOLB <- SCTransform(ssmOLB, assay = "Spatial", ncells = 3000, verbose = TRUE)
ssmOLB <- RunPCA(ssmOLB,assay = "SCT", verbose = FALSE)

#check how many dimensions is better
ElbowPlot(ssmOLB)

ssmOLB <- FindNeighbors(ssmOLB, reduction = "pca", dims = 1:18)
ssmOLB <- FindClusters(ssmOLB, verbose = FALSE)
ssmOLB <- RunUMAP(ssmOLB, reduction = "pca", dims = 1:18)

plot1 <- DimPlot(ssmOLB, reduction = "umap", label = TRUE)
plot2 <- SpatialDimPlot(ssmOLB, label.size = 3)


mypal <- brewer.pal(n = 12, name ="Paired")
mydal <- brewer.pal(n = 8, name ="Set2")
mycal <- brewer.pal(n = 8, name ="Set3")
mycol = c(mydal, mypal, mycal)

#cluster plot and spatial cluster

DimPlot(ssmOLB, reduction = "umap",cols = mycol)

mycolset = c(colPurple[c(7,7,7)],colGreen[c(7,7,7,7)],colGreys[c(5,5,5,5)],colBlues[c(6,6,6,6,6,6)],colOrange[4],colPink[c(3,3,3)],colBrown[c(4,4)])
myOrder = c(0,8,12,1,10,16,22,2,7,11,13,3,6,20,5,21,17,9,15,14,19,18,4)
par(mfrow=c(4,6),mar=c(1,1,1,1))
for(i in 1:23){
  plot(devi[,29:30],type='n',axes = F)
  points(devi[devi$cluster==myOrder[i],29:30],col=mycolset[i],pch=19,cex=.1)
}


###################################################################################################################
############################################################################scRNA-seq guided prediction
###################################################################################################################

library(SeuratDisk)

OBscAll = CreateSeuratObject(counts = counts_full, project = "OBscAll", min.cells = 3, min.features = 100)
OBscAll@meta.data<-cbind(OBscAll@meta.data,meta_full,CellType)
SaveH5Seurat(OBscAll, filename = "~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/scAnnotation/OB_sc.h5Seurat",overwrite =T)
Convert("~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/scAnnotation/OB_sc.h5Seurat", dest = "h5ad",overwrite =T)


#-----------------prepare spatial data for 10um

expr <- data.table::fread(
  input = "~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/segmentation/OB36bin33dCellGeneExp.txt",
  sep = ",",
  data.table = FALSE
)
rownames(x = expr) <- (x = expr[, 1])
colnames(x = expr) <- (x = colnames(x = expr))
expr <- t(x = expr[, -1])
OB36_pixel <- CreateSeuratObject(
  counts = expr,
  project = 'SlideSeq',
  assay = 'Spatial'
)
positions <- read.csv(
  file = "~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/segmentation/OB36bin33Coordd.txt"
)


nn1 = positions[positions$xcoord>24000 & positions$xcoord <40000,]
nn1$xcoord = nn1$xcoord - 24000 
nn1$xcoord = nn1$xcoord/max(nn1$xcoord)
nn1$ycoord = nn1$ycoord/max(nn1$ycoord)
rownames(x = nn1) <- nn1$binID
positions <- nn1[, 2:3]

OB36_pixel[['image']] <-new(
  Class = 'SlideSeq',
  assay = "Spatial",
  coordinates = positions
)

#built-in UMI check plot
#25um cutoff 1024

OB36_pixel = OB36_pixel[, OB36_pixel$nFeature_Spatial >= 10]
OB36_pixel = OB36_pixel[, OB36_pixel$nCount_Spatial >= 512]
SaveH5Seurat(OB36_pixel, filename = "~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/scAnnotation/OB36_10um_512Cut.h5Seurat",overwrite =T)
Convert("~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/scAnnotation/OB36_10um_512Cut.h5Seurat", dest = "h5ad",overwrite =T)


#---------------Run devi for annotation
system(paste("~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/scAnnotation/devi.py ",'OB36_10um_512Cut.h5ad',' OB36_10um_512Cut_devi.csv',sep=''))
system(paste("~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/scAnnotation/devi.py ",'OB36_vSegN8003U.h5ad',' OB36_vSegN8003U_2_devi.csv',sep=''))

#---------------analyze the devi result
devi = read.csv('~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/scAnnotation/OB36_vSegN8003U_2_devi.csv',header = T,row.names = 1)
coord = read.csv('~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/segmentation/spt_3U_MaskCoord.txt',header=T)

subCoord = coord[coord$binID %in% rownames(devi),]
rownames(subCoord) = subCoord$binID
subCoord = subCoord[rownames(devi),]

subClustering = clustering[clustering$cells %in% rownames(devi), ]
subClustering = subClustering[rownames(devi),]

devi = as.data.frame(devi)
devi$assign = colnames(devi)[apply(devi[,1:26], 1,which.max)]
devi$score = apply(devi[,1:26], 1,max)
nrow(devi[devi$score>0.3,])/nrow(devi)

devi$xcoord = subCoord$xcoord
devi$ycoord = subCoord$ycoord
devi$total_counts = subCoord$total_counts
devi$B_counts = subCoord$BarcodeCount
devi$cluster = subClustering$ident

####################################################################################################################################################################
######################################################################################################cell size vs tUMI by cell type, Fig. 2F
####################################################################################################################################################################
scDensity = data.frame(bcount = (devi$B_counts/3.14)^0.5*2,tUMI = devi$total_counts,density = devi$total_counts/(devi$B_counts*1.44), celltype = devi$assign,score = devi$score,color = rep('white',nrow(devi)))
scDensity[scDensity$celltype=='n02.PGC.1',]$color='blue'
scDensity[scDensity$celltype=='Neuron.M.TC',]$color='green'


#range


par(mar=c(3,3,1,1))
plot(log10(scDensity[scDensity$color=='white',]$tUMI),scDensity[scDensity$color=='white',]$bcount,cex=.8,col=rgb(0.5,0.5,0.5,0.05),pch=19,axes=F,xlab='',ylab='',xlim=c(2,5),ylim=c(4,30))
points(log10(scDensity[scDensity$celltype=='Neuron.M.TC',]$tUMI),scDensity[scDensity$celltype=='Neuron.M.TC',]$bcount,cex=.05,col='#ff8080',pch=19)
points(log10(scDensity[scDensity$celltype=='n02.PGC.1',]$tUMI),scDensity[scDensity$celltype=='n02.PGC.1',]$bcount,cex=.1,col='#1aff66',pch=19)

axis(2,at=c(5,10,20,30),labels = c(5,10,20,30),lwd=0.2)
axis(1,at=c(2,3,4,5),labels = c(2,3,4,5),lwd=0.2)
box(lwd=0.2)


plot(density(log10(scDensity[scDensity$color=='white',]$tUMI)),xlim=c(2,5),main='',lwd=3,axes=F,ylim=c(0,1.2))
lines(density(log10(scDensity[scDensity$celltype=='Neuron.M.TC',]$tUMI)),col='#ff8080',lwd=3)
lines(density(log10(scDensity[scDensity$celltype=='n02.PGC.1',]$tUMI)),col='#1aff66',lwd=3)
axis(1,at=c(2,3,4,5),labels = c(2,3,4,5),lwd=0.2)
axis(2,at=c(0,0.4,0.8,1.2),labels = c(0,0.4,0.8,1.2),lwd=0.2,las=2)
box(lwd=0.2)




lines(spline(scDensity[scDensity$color=='white',]$bcount,log2(scDensity[scDensity$color=='white',]$density)),col='gray40')
smoothScatter(scDensity[scDensity$color=='white',]$bcount,log2(scDensity[scDensity$color=='white',]$tUMI),colramp = colorRampPalette(c("white", 'gray40')),axes=F)
points(scDensity[scDensity$celltype=='Neuron.M.TC',]$bcount,log2(scDensity[scDensity$celltype=='Neuron.M.TC',]$tUMI),cex=.2,col=colSet1[2],pch=19)
points(scDensity[scDensity$celltype=='n02.PGC.1',]$bcount,log2(scDensity[scDensity$celltype=='n02.PGC.1',]$tUMI),cex=.35,col=colSet1[5],pch=19)
lines(spline(scDensity[scDensity$color=='white',]$bcount,log2(scDensity[scDensity$color=='white',]$tUMI),n=2,method = "natural"),col='blue')
colSet1 <- brewer.pal(n = 9, name ="Set1")



typeUniq = unique(scDensity$celltype)
meanUniq = rep(0,length(typeUniq))
sdUniq = rep(0,length(typeUniq))
for(i in 1:length(typeUniq)){
  sdUniq[i] = sd(scDensity[scDensity$celltype==typeUniq[i],]$bcount)
  meanUniq[i] = mean(scDensity[scDensity$celltype==typeUniq[i],]$bcount)
}
cbind(typeUniq,meanUniq,sdUniq)

################################################################################################################################################
#####################################################################################################################cell density, Fig. 2G
################################################################################################################################################
par(mar=c(2,2,1,1),mfrow=c(1,1))
plot(c(0,27),c(0,60),type='n',axes=F)
cellTypeN = c("Mes","OEC","Neuron.OSNs","n03.GC.1","n02.PGC.1","n05.PGC.2","n08.PGC.3","n18.EPL.IN","Neuron.M.TC","n07.GC.2","n09.GC.3",
  "n10.GC.4","n11.GC.5","n12.GC.6","n14.GC.7","n06.Transition","Neuron.Immature","Astro","MicroG","OPC","MyOligo","Mural","Neuron.AstroLike","ImmunoCells","EC")


for(i in 1:length(cellTypeN)){
  selectType = scDensity[scDensity$celltype == cellTypeN[i],]
  vioplot(selectType$density,add=T,at=i,col='NA',rectCol='gray10',colMed='gray10',pchMed=15,lineCol = 'gray10',border='gray10',axes=F)
}
#axis(1,at=1:26,labels = cellTypeN,las=2)
axis(2,at=c(0,20,40,60),labels=c(0,20,40,60),las=2)



#--------------------------cell size
plot(c(0,27),c(0,40),type='n',axes=F)
cellTypeN = c("Mes","OEC","Neuron.OSNs","n03.GC.1","n02.PGC.1","n05.PGC.2","n08.PGC.3","n18.EPL.IN","Neuron.M.TC","n07.GC.2","n09.GC.3",
              "n10.GC.4","n11.GC.5","n12.GC.6","n14.GC.7","n06.Transition","Neuron.Immature","Astro","MicroG","OPC","MyOligo","Mural","Neuron.AstroLike","ImmunoCells","EC")

colSize = c(rep('gray60',9),rep(colBlues[5],16),'gray80')
for(i in 1:length(cellTypeN)){
  selectType = scDensity[scDensity$celltype == cellTypeN[i],]
  vioplot(selectType$bcount,add=T,at=i,col='NA',rectCol='gray10',colMed='gray10',pchMed=15,lineCol = 'gray10',border='gray10',axes=F)
}
axis(2,at=c(0,10,20,30,40,50),labels=c(0,10,20,30,40,50),las=2)



###-------------------cell tUMI


pAstro = scDensity[scDensity$celltype=='Astro',]
MicroG = scDensity[scDensity$celltype=='MicroG',]
MyOligo = scDensity[scDensity$celltype=='MyOligo',]
OPC = scDensity[scDensity$celltype=='OPC',]
EC = scDensity[scDensity$celltype=='EC',]
ImmunoCells  = scDensity[scDensity$celltype=='ImmunoCells',] 
Mural = scDensity[scDensity$celltype=='Mural',]
OEC = scDensity[scDensity$celltype=='OEC',]
Mes = scDensity[scDensity$celltype=='Mes',]
Neuron.OSNs = scDensity[scDensity$celltype=='Neuron.OSNs',]
Neuron.PGC = scDensity[scDensity$celltype %in% c('n08.PGC.3','n02.PGC.1','n05.PGC.2'),]
Neuron.EPL.IN = scDensity[scDensity$celltype=='n18.EPL.IN',]
Neuron.M.TC = scDensity[scDensity$celltype=='Neuron.M.TC',]
Neuron.GC = scDensity[scDensity$celltype %in% c('n03.GC.1','n07.GC.2','n09.GC.3','n10.GC.4','n11.GC.5','n12.GC.6','n14.GC.7'),]
Neuron.Transition = scDensity[scDensity$celltype=='n06.Transition',]
Neuron.Immature = scDensity[scDensity$celltype=='Neuron.Immature',]
Neuron.AstroLike = scDensity[scDensity$celltype=='Neuron.AstroLike',]
RBCs = scDensity[scDensity$celltype=='RBCs',]


par(mar=c(2,2,1,1))
plot(c(0,19),c(7,15),type='n',axes=F)
vioplot(log2(pAstro$tUMI),add=T,at=1,col = colSize[1])
vioplot(log2(MicroG$tUMI),add=T,at=2,col = colSize[2])
vioplot(log2(MyOligo$tUMI),add=T,at=3,col = colSize[3])
vioplot(log2(OPC$tUMI),add=T,at=4,col = colSize[4])
vioplot(log2(EC$tUMI),add=T,at=5,col = colSize[5])
vioplot(log2(ImmunoCells$tUMI),add=T,at=6,col = colSize[6])
vioplot(log2(Mural$tUMI),add=T,at=7,col = colSize[7])
vioplot(log2(OEC$tUMI),add=T,at=8,col = colSize[8])
vioplot(log2(Mes$tUMI),add=T,at=9,col = colSize[9])
vioplot(log2(Neuron.OSNs$tUMI),add=T,at=10,col = colSize[10])
vioplot(log2(Neuron.PGC$tUMI),add=T,at=11,col = colSize[11])
vioplot(log2(Neuron.EPL.IN$tUMI),add=T,at=12,col = colSize[12])
vioplot(log2(Neuron.M.TC$tUMI),add=T,at=13,col = colSize[13])
vioplot(log2(Neuron.GC$tUMI),add=T,at=14,col = colSize[14])
vioplot(log2(Neuron.Transition$tUMI),add=T,at=15,col = colSize[15])
vioplot(log2(Neuron.Immature$tUMI),add=T,at=16,col = colSize[16])
vioplot(log2(Neuron.AstroLike$tUMI),add=T,at=17,col = colSize[17])
vioplot(log2(RBCs$tUMI),add=T,at=18,col = colSize[18])

axis(2,at=c(5,10,20,30),labels = c(5,10,20,30))

###-------------------------------density

par(mar=c(2,2,1,1))
plot(c(0,19),c(0,60),type='n',axes=F)
vioplot(pAstro$density,add=T,at=1,col = colSize[1])
vioplot(MicroG$density,add=T,at=2,col = colSize[2])
vioplot(MyOligo$density,add=T,at=3,col = colSize[3])
vioplot(OPC$density,add=T,at=4,col = colSize[4])
vioplot(EC$density,add=T,at=5,col = colSize[5])
vioplot(ImmunoCells$density,add=T,at=6,col = colSize[6])
vioplot(Mural$density,add=T,at=7,col = colSize[7])
vioplot(OEC$density,add=T,at=8,col = colSize[8])
vioplot(Mes$density,add=T,at=9,col = colSize[9])
vioplot(Neuron.OSNs$density,add=T,at=10,col = colSize[10])
vioplot(Neuron.PGC$density,add=T,at=11,col = colSize[11])
vioplot(Neuron.EPL.IN$density,add=T,at=12,col = colSize[12])
vioplot(Neuron.M.TC$density,add=T,at=13,col = colSize[13])
vioplot(Neuron.GC$density,add=T,at=14,col = colSize[14])
vioplot(Neuron.Transition$density,add=T,at=15,col = colSize[15])
vioplot(Neuron.Immature$density,add=T,at=16,col = colSize[16])
vioplot(Neuron.AstroLike$density,add=T,at=17,col = colSize[17])
vioplot(RBCs$density,add=T,at=18,col = colSize[18])

axis(2,at=c(5,10,20,30),labels = c(5,10,20,30))


########################################################################################################################################
####################################################################################cell composition, Fig. 2H
########################################################################################################################################

scAstro	=3130+1780+336
scMicroG	=2260+790+951
scMyOligo=	858
scOPC	=397
scEC	=2524+1124
scImmuneCell=	623+572
scMural=	158+769
scOEC	=2252+883+3436+3240+1740
scMes	=230+653
scOSNs=	1200
scPGC	=1037+377+279
scEPL.IN	=161
scM.TC	=90+116+927
scGC	=1914+3608+273+234+1669+237+679
scTransition=	2085
scImmature=	3910
scAstroLike	=2950
scRBCs=	219

scOLBcomp = c(scAstro,scMicroG,scMyOligo,scOPC,scEC,scImmuneCell,scMural,scOEC,scMes,scOSNs,scPGC,scEPL.IN,scM.TC,scGC,scTransition,scImmature,scAstroLike)
pixelOLBcomp = c(nrow(pAstro),nrow(MicroG),nrow(MyOligo),nrow(OPC),nrow(EC),nrow(ImmunoCells),nrow(Mural),nrow(OEC),nrow(Mes),nrow(Neuron.OSNs),
                 nrow(Neuron.PGC),nrow(Neuron.EPL.IN),nrow(Neuron.M.TC),nrow(Neuron.GC),nrow(Neuron.Transition),nrow(Neuron.Immature),
                 nrow(Neuron.AstroLike))


cellCompostion = data.frame(scOLB =scOLBcomp/sum(scOLBcomp), colo = colSize)
rownames(cellCompostion) = c("Astro","MicroG","MyOligo","OPC","EC","ImmunoCells","Mural","OEC","Mes","Neuron.OSNs","Neuron.PGC","Neuron.EPL.IN","Neuron.M.TC","Neuron.GC","Neuron.Transition","Neuron.Immature",
"Neuron.AstroLike","RBCs")
cellCompostion = cellCompostion[order(cellCompostion$scOLB),]
cellCompostionR = cellCompostion[!(rownames(cellCompostion) %in% 'RBCs'),]
par(mar=c(5,5,1,1))
barplot(cellCompostionR$scOLB,horiz = T,col=cellCompostionR$colo,names.arg=rownames(cellCompostionR),las=2,xlim=c(0,0.25))


cellCompostionPixel = data.frame(pixel =pixelOLBcomp/sum(pixelOLBcomp), colo = colSize)
rownames(cellCompostionPixel) = c("Astro","MicroG","MyOligo","OPC","EC","ImmunoCells","Mural","OEC","Mes","Neuron.OSNs","Neuron.PGC","Neuron.EPL.IN","Neuron.M.TC","Neuron.GC","Neuron.Transition","Neuron.Immature",
                             "Neuron.AstroLike","RBCs")
cellCompostionPixel = cellCompostionPixel[order(cellCompostionPixel$pixel),]
cellCompostionPixelR = cellCompostionPixel[!(rownames(cellCompostionPixel) %in% 'RBCs'),]
par(mar=c(5,5,1,1))
barplot(cellCompostionPixelR$pixel,horiz = T,col=cellCompostionPixelR$colo,names.arg=rownames(cellCompostionPixelR),las=2,xlim=c(0,0.25))


########################################################################################################################
################################################################spatial pattern analysis of 10um segmentation, Fig. 2D
#######################################################################################################################
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

expr <- data.table::fread(
  input = "~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/segmentation/OB36bin33dCellGeneExp.txt",
  sep = ",",
  data.table = FALSE
)

rownames(x = expr) <- (x = expr[, 1])
colnames(x = expr) <- (x = colnames(x = expr))
expr <- t(x = expr[, -1])
ssmOLB10 <- CreateSeuratObject(
  counts = expr,
  project = 'SlideSeq',
  assay = 'Spatial'
)
positions <- read.csv(
  file = "~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/segmentation/OB36bin33Coordd.txt"
)

nn1 = positions[positions$xcoord>24000 & positions$xcoord <40000,]


rownames(x = nn1) <- nn1$binID
positions <- nn1[, 2:3]

ssmOLB10[['image']] <-new(
  Class = 'SlideSeq',
  assay = "Spatial",
  coordinates = positions
)

ssmOLB = ssmOLB[, ssmOLB$nCount_Spatial > 512]
#normalization and dimensions reduction
ssmOLB10 <- SCTransform(ssmOLB10, assay = "Spatial", ncells = 2000, verbose = TRUE)

ssmOLB10 <- RunPCA(ssmOLB10,assay = "SCT", verbose = FALSE)
ElbowPlot(ssmOLB10)
#spatial pattern analysis
ssmOLB10 <- FindNeighbors(ssmOLB10, reduction = "pca", dims = 1:15)
ssmOLB10 <- FindClusters(ssmOLB10, verbose = FALSE,resolution = 0.4)
ssmOLB10 <- RunUMAP(ssmOLB10, reduction = "pca", dims = 1:15)
#direct pattern Umap and plot
plot1 <- DimPlot(ssmOLB10, reduction = "umap", label = TRUE)
plot2 <- SpatialDimPlot(ssmOLB10, label.size = 1)
plot1 + plot2


#plot the spatial pattern
library(RColorBrewer)
clustering10 = plot2$data
clustering10$umi = nn1[rownames(clustering),]$total_counts

write.table(clustering10,'~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/clustering/OB36_10um_Vd_15PCA_re0.4.txt')


#----------------------------sc Annotation for 10um
devi = read.csv('~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/scAnnotation/OB36_10um_512Cut_devi.csv',header = T,row.names = 1)
coord = read.csv('~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/segmentation/OB36bin33Coordd.txt',header=T)
subCoord = coord[coord$binID %in% rownames(devi),]
rownames(subCoord) = subCoord$binID
subCoord = subCoord[rownames(devi),]

devi = as.data.frame(devi)
devi$assign = colnames(devi)[apply(devi[,1:26], 1,which.max)]
devi$score = apply(devi[,1:26], 1,max)
nrow(devi[devi$score>0.3,])/nrow(devi)

devi$xcoord = subCoord$xcoord
devi$ycoord = subCoord$ycoord

cellTypeNames = colnames(devi)
par(mfrow=c(5,6),mar=c(3,3,1,1))
for(i in 1:length(cellTypeNames)){
  subPosType = devi[devi$assign==cellTypeNames[i] & devi$score>0.3,]
  if(nrow(subPosType)>0){
    plot(subPosType$xcoord,subPosType$ycoord,col=2,pch=19,cex=.1,main=cellTypeNames[i])
  }
}

write.table(devi,'~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/scAnnotation/OB36_10um_Vd_scAnnotation_100_250.txt')






#######################################################################################################
#####################################################################cluster level _UMAP _scAnnotation, Fig. 2C
#######################################################################################################

#overview of the whole picture
devi = read.csv('~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/scAnnotation/OB36_vSegN8003U_2_devi.csv',header = T,row.names = 1)
coord = read.csv('~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/segmentation/spt_3U_MaskCoord.txt',header=T)

subCoord = coord[coord$binID %in% rownames(devi),]
rownames(subCoord) = subCoord$binID
subCoord = subCoord[rownames(devi),]
.
devi = as.data.frame(devi)
devi$assign = colnames(devi)[apply(devi[,1:26], 1,which.max)]
devi$score = apply(devi[,1:26], 1,max)
nrow(devi[devi$score>0.3,])/nrow(devi)

devi$xcoord = subCoord$xcoord
devi$ycoord = subCoord$ycoord
devi$total_counts = subCoord$total_counts
devi$B_counts = subCoord$B_counts
cellTypeNames = colnames(devi)
par(mfrow=c(5,6),mar=c(3,3,1,1))
for(i in 1:length(cellTypeNames)){
  subPosType = devi[devi$assign==cellTypeNames[i] & devi$score>0.3,]
  if(nrow(subPosType)>0){
    plot(subPosType$xcoord,subPosType$ycoord,col=2,pch=19,cex=.1,main=cellTypeNames[i])
  }
}

write.table(devi,"~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/scAnnotation/scAnnotation_cluster_vSegN8003U.txt")

#####ONL layer (types: OSN, Mes, OEC)
par(mfrow=c(1,6))
par(mar=c(2,2,1,1))
ONL_OSN = devi[devi$assign=='Neuron.OSNs' & devi$score>0.3,]
ONL_MES = devi[devi$assign=='Mes' & devi$score>0.3,]
ONL_OEC = devi[devi$assign=='OEC' & devi$score>0.3,]
col1 <- brewer.pal(n = 9, name ="PuBu")
plot(ONL_OSN$xcoord,ONL_OSN$ycoord,col=col1[4],pch=19,axes = F,xlab='',ylab='',cex=.02,xlim=c(24000,40000),ylim=c(0,12000))
points(ONL_OEC$xcoord,ONL_OEC$ycoord,col=col1[6],pch=19,cex=.01)
points(ONL_MES$xcoord,ONL_MES$ycoord,col=col1[8],pch=19,cex=.01)
#points(ONL_MES$xcoord,ONL_MES$ycoord,col=col1[6],pch=19,cex=.2)
box(col='gray80')


#####GL layer (type: PGC, n03.GC-1)
GL_n03.GC.1 = devi[devi$assign=='n03.GC.1' & devi$score>0.3,]
GL_n02.PGC.1 = devi[devi$assign=='n02.PGC.1' & devi$score>0.3,]
GL_n05.PGC.2 = devi[devi$assign=='n05.PGC.2' & devi$score>0.3,]
GL_n08.PGC.3 = devi[devi$assign=='n08.PGC.3' & devi$score>0.3,]
col2 <- brewer.pal(n = 11, name ="PiYG")
plot(GL_n03.GC.1$xcoord,GL_n03.GC.1$ycoord,col=col2[2],pch=19,axes = F,xlab='',ylab='',cex=.02,xlim=c(24000,40000),ylim=c(0,12000))
points(GL_n02.PGC.1$xcoord,GL_n02.PGC.1$ycoord,col=col2[3],pch=19,cex=.01)
points(GL_n05.PGC.2$xcoord,GL_n05.PGC.2$ycoord,col=col2[4],pch=19,cex=.01)
points(GL_n08.PGC.3$xcoord,GL_n08.PGC.3$ycoord,col=col2[5],pch=19,cex=.01)
box(col='gray80')

####EPL layer (type EPL-IN, Astro)
EPL_Astro = devi[devi$assign=='Astro' & devi$score>0.3,]
EPL_n18.EPL.IN = devi[devi$assign=='n18.EPL.IN' & devi$score>0.3,]
col3 = brewer.pal(n = 9, name ="Greys")
plot(EPL_n18.EPL.IN$xcoord,EPL_n18.EPL.IN$ycoord,col=col3[7],pch=19,axes = F,xlab='',ylab='',cex=.2,xlim=c(24000,40000),ylim=c(0,12000))
points(EPL_Astro$xcoord,EPL_Astro$ycoord,col=col3[5],pch=19,cex=.2)
box(col='gray80')

#####MCL layer (type: M/TC)
MCL_M.TC = devi[devi$assign=='Neuron.M.TC' & devi$score>0.3,]
col2 <- brewer.pal(n = 9, name ="Greens")
plot(MCL_M.TC$xcoord,MCL_M.TC$ycoord,col=col2[8],pch=19,xlab='',ylab='',cex=.2,axes = F,xlim=c(24000,40000),ylim=c(0,12000))
box(col='gray80')



#####GC layer (type: GC,transition)
GC_n11.GC.5 = devi[devi$assign=='n11.GC.5' & devi$score>0.3,]
GC_n12.GC.6 = devi[devi$assign=='n12.GC.6' & devi$score>0.3,]
GC_n06.Transition = devi[devi$assign=='n06.Transition' & devi$score>0.3,]
col5 <- brewer.pal(n = 9, name ="Purples")
plot(GC_n11.GC.5$xcoord,GC_n11.GC.5$ycoord,col=col5[6],pch=19,xlab='',ylab='',cex=.2,axes = F,xlim=c(24000,40000),ylim=c(0,12000))
points(GC_n12.GC.6$xcoord,GC_n12.GC.6$ycoord,col=col5[8],pch=19,cex=.2)
points(GC_n06.Transition$xcoord,GC_n06.Transition$ycoord,col=col5[4],pch=19,cex=.2)
box(col='gray80')

#RMS layer (type, Immumature)
RMS_Neuron.Immature = devi[devi$assign=='Neuron.Immature' & devi$score>0.3,]
col6 <- brewer.pal(n = 9, name ="Oranges")
plot(RMS_Neuron.Immature$xcoord,RMS_Neuron.Immature$ycoord,col=col6[6],pch=19,xlab='',ylab='',cex=.2,axes = F,xlim=c(24000,40000),ylim=c(0,12000))
box(col='gray80')



###################################################################################################
#########################################################zoom in comparison, Fig. 2D
###################################################################################################

clusterU_sc = read.table('~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/scAnnotation/OB36_clusterlevel_cluster_scA.txt',header=F)
colnames(clusterU_sc) = c('xcoord','ycoord','binseq','tUMI','bcount','clusterUmap','cellType','score')

fullPlot = clusterU_sc[clusterU_sc$xcoord>24000 & clusterU_sc$xcoord<40000 & clusterU_sc$ycoord>0 & clusterU_sc$ycoord<12000,]

subPlot = clusterU_sc[clusterU_sc$xcoord>24000 & clusterU_sc$xcoord<28500 & clusterU_sc$ycoord>5000 & clusterU_sc$ycoord<6000,]


#plot the check the pattern per cluster
cluster = read.table('~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/clustering/vSeg3UN800_PCA_UMAP_cluster.txt',header=T)
par(mfrow=c(5,6),mar=c(3,3,1,1))
for(i in 0:22){
  subPosType = cluster[cluster$ident==i,]
  if(nrow(subPosType)>0){
    plot(subPosType$x,subPosType$y,col=rgb(1,0,0,0.1),pch=19,cex=.2,main=i,xlim=c(24000,40000),ylim=c(0,12000))
  }
}


library(RColorBrewer)
par(mfrow=c(2,1),mar=c(3,3,1,1))
colSet <- brewer.pal(n = 9, name ="Set1")
subPlot$celltypeCol = rep('white',nrow(subPlot))
subPlot$ctName = rep('grey80',nrow(subPlot))

#cell type - color assignment
colBlues <- brewer.pal(n = 9, name ="Blues")
colGreen <- brewer.pal(n = 9, name ="Greens")
colPurple <- brewer.pal(n = 9, name ="Purples")
colOrange <- brewer.pal(n = 9, name ="Oranges")
colRed <- brewer.pal(n = 9, name ="Reds")
colBrown <- brewer.pal(n = 9, name ="BrBG")
colPink <- brewer.pal(n = 9, name ="PiYG")
colGreys <- brewer.pal(n = 9, name ="Greys")

#pattern ONL layer
subPlot[subPlot$cellType == 'OEC',9] = colBlues[6]
subPlot[subPlot$cellType == 'Mes',9] = colBlues[8]
subPlot[subPlot$cellType == 'Neuron.OSNs',9] = colBlues[4]
#---------------------
subPlot[subPlot$cellType == 'OEC',10] = 'cBlue4'
subPlot[subPlot$cellType == 'Mes',10] = 'cBlue6'
subPlot[subPlot$cellType == 'Neuron.OSNs',10] = 'cBlue8'

#pattern GL
subPlot[subPlot$cellType == 'n03.GC.1',9] = colPink[2]
subPlot[subPlot$cellType == 'n02.PGC.1',9] = colPink[3]
subPlot[subPlot$cellType == 'n05.PGC.2',9] = colPink[4]
subPlot[subPlot$cellType == 'n08.PGC.3',9] = colPink[5]
#pattern EPL
subPlot[subPlot$cellType == 'n18.EPL.IN',9] = colBrown[3]
#pattern MCL
subPlot[subPlot$cellType == 'Neuron.M.TC',9] = colGreen[8]
#---------------------
subPlot[subPlot$cellType == 'n03.GC.1',10] = 'colPink2'
subPlot[subPlot$cellType == 'n02.PGC.1',10] = 'colPink3'
subPlot[subPlot$cellType == 'n05.PGC.2',10] = 'colPink4'
subPlot[subPlot$cellType == 'n08.PGC.3',10] = 'colPink5'
subPlot[subPlot$cellType == 'n18.EPL.IN',10] = 'colBrown3'
subPlot[subPlot$cellType == 'Neuron.M.TC',10] = 'colGreen8'


subPlot[subPlot$cellType == 'n11.GC.5',9] = colPurple[8]
subPlot[subPlot$cellType == 'n12.GC.6',9] = colPurple[7]
subPlot[subPlot$cellType == 'n09.GC.3',9] = colPurple[6]
subPlot[subPlot$cellType == 'n10.GC.4',9] = colPurple[5]
subPlot[subPlot$cellType == 'n07.GC.2',9] = colPurple[4]
subPlot[subPlot$cellType == 'n14.GC.7',9] = colPurple[3]
subPlot[subPlot$cellType == 'n06.Transition',9] = colPurple[2]
#---------------------
subPlot[subPlot$cellType == 'n11.GC.5',10] = 'colPurple8'
subPlot[subPlot$cellType == 'n12.GC.6',10] = 'colPurple7'
subPlot[subPlot$cellType == 'n09.GC.3',10] = 'colPurple6'
subPlot[subPlot$cellType == 'n10.GC.4',10] = 'colPurple5'
subPlot[subPlot$cellType == 'n07.GC.2',10] = 'colPurple4'
subPlot[subPlot$cellType == 'n14.GC.7',10] = 'colPurple3'
subPlot[subPlot$cellType == 'n06.Transition',10] = 'colPurple2'


subPlot[subPlot$cellType == 'Neuron.Immature',9] = colOrange[6]
subPlot[subPlot$cellType == 'Neuron.Immature',10] = 'colOrange6'

subPlot[subPlot$cellType == 'Neuron.AstroLike',9] = colOrange[1]
subPlot[subPlot$cellType == 'Neuron.AstroLike',10] = 'colOrange1'


subPlot[subPlot$cellType == 'EC',9] = colRed[7]
subPlot[subPlot$cellType == 'ImmunoCells',9] = colRed[5]
subPlot[subPlot$cellType == 'Mural',9] = colRed[3]
#---------------------
subPlot[subPlot$cellType == 'EC',10] = 'colRed7'
subPlot[subPlot$cellType == 'ImmunoCells',10] = 'colRed5'
subPlot[subPlot$cellType == 'Mural',10] = 'colRed3'

subPlot[subPlot$cellType == 'Astro',9] = colGreys[8]
subPlot[subPlot$cellType == 'MicroG',9] = colGreys[7]
subPlot[subPlot$cellType == 'MyOligo',9] = colGreys[6]
subPlot[subPlot$cellType == 'OPC',9] = colGreys[5]
#----------------
subPlot[subPlot$cellType == 'Astro',10] = 'colGreys8'
subPlot[subPlot$cellType == 'MicroG',10] = 'colGreys7'
subPlot[subPlot$cellType == 'MyOligo',10] = 'colGreys6'
subPlot[subPlot$cellType == 'OPC',10] = 'colGreys5'

par(mar=c(2,2,1,1),mfrow=c(2,1))
plot(subPlot$xcoord,12000-subPlot$ycoord,col=subPlot$celltypeCol,pch=15,cex=.25,axes = F,xlab = '',ylab = '')
box(col='gray80')

#plot the legend for annotation
typeColorUniq = unique(subPlot[,c(7,9)])
plot(c(0,26),c(0,5),type='n')
text(1:nrow(typeColorUniq),seq(0,5,length.out = nrow(typeColorUniq)),labels = typeColorUniq$cellType,col = typeColorUniq$celltypeCol,cex=2)



###cluster color assignment
subPlot$clusterColor = rep(colSet[9],nrow(subPlot))
subPlot[subPlot$clusterUmap == 0,11] = colPurple[8] #GC
subPlot[subPlot$clusterUmap == 8,11] = colPurple[7] #GC
subPlot[subPlot$clusterUmap == 12,11] = colPurple[3]#GC

subPlot[subPlot$clusterUmap == 1,11] = colGreen[8]    #MT/C
subPlot[subPlot$clusterUmap == 10,11] = colGreen[5]   #MT/C
subPlot[subPlot$clusterUmap == 16,11] = colGreen[7]   #MT/C
subPlot[subPlot$clusterUmap == 22,11] = colGreen[4]   #M/TC

subPlot[subPlot$clusterUmap == 2,11] = colGreys[6]    # mainly EPL, astro
subPlot[subPlot$clusterUmap == 7,11] = colGreys[8]    # EPL
subPlot[subPlot$clusterUmap == 11,11] = colGreys[3]   # EPL
subPlot[subPlot$clusterUmap == 13,11] = colGreys[5]  # EPL

subPlot[subPlot$clusterUmap == 3,11] = colBlues[4]    #ONL
subPlot[subPlot$clusterUmap == 6,11] = colBlues[8]    #ONL
subPlot[subPlot$clusterUmap == 5,11] = colBlues[6]    #ONL
subPlot[subPlot$clusterUmap == 21,11] = colBlues[2]  #ONL
subPlot[subPlot$clusterUmap == 17,11] = colBlues[1] 

subPlot[subPlot$clusterUmap == 9,11] = colOrange[4]  #RMS
subPlot[subPlot$clusterUmap == 15,11] = colPink[3]  #GL
subPlot[subPlot$clusterUmap == 14,11] = colPink[4]   #GL
subPlot[subPlot$clusterUmap == 19,11] = colPink[5]   #GL

subPlot[subPlot$clusterUmap == 18,11] =  colBrown[4] #other
subPlot[subPlot$clusterUmap == 4,11] = colBrown[5] #other

plot(subPlot$xcoord,12000-subPlot$ycoord,col=subPlot$clusterColor,pch=15,cex=.25,axes = F,xlab = '',ylab='')
box(col='gray80')
typeColorUniq = unique(subPlot[,c(6,11)])
plot(c(0,26),c(0,5),type='n')
text(1:nrow(typeColorUniq),seq(0,5,length.out = nrow(typeColorUniq)),labels = typeColorUniq$clusterUmap,col = typeColorUniq$clusterColor,cex=2)



###cluster color assignment for fullplot figure
fullPlot$clusterColor = rep(colSet[9],nrow(fullPlot))
fullPlot[fullPlot$clusterUmap == 0,9] = colPurple[8] #GC
fullPlot[fullPlot$clusterUmap == 8,9] = colPurple[7] #GC
fullPlot[fullPlot$clusterUmap == 12,9] = colPurple[3]#GC

fullPlot[fullPlot$clusterUmap == 1,9] = colGreen[8]    #MT/C
fullPlot[fullPlot$clusterUmap == 10,9] = colGreen[5]   #MT/C
fullPlot[fullPlot$clusterUmap == 16,9] = colGreen[7]   #MT/C
fullPlot[fullPlot$clusterUmap == 22,9] = colGreen[4]   #M/TC

fullPlot[fullPlot$clusterUmap == 2,9] = colGreys[6]    # mainly EPL, astro
fullPlot[fullPlot$clusterUmap == 7,9] = colGreys[8]    # EPL
fullPlot[fullPlot$clusterUmap == 11,9] = colGreys[3]   # EPL
fullPlot[fullPlot$clusterUmap == 13,9] = colGreys[5]  # EPL

fullPlot[fullPlot$clusterUmap == 3,9] = colBlues[4]    #ONL
fullPlot[fullPlot$clusterUmap == 6,9] = colBlues[8]    #ONL
fullPlot[fullPlot$clusterUmap == 5,9] = colBlues[6]    #ONL
fullPlot[fullPlot$clusterUmap == 21,9] = colBlues[2]  #ONL
fullPlot[fullPlot$clusterUmap == 17,9] = colBlues[1] 

fullPlot[fullPlot$clusterUmap == 9,9] = colOrange[4]  #RMS
fullPlot[fullPlot$clusterUmap == 15,9] = colPink[3]  #GL
fullPlot[fullPlot$clusterUmap == 14,9] = colPink[4]   #GL
fullPlot[fullPlot$clusterUmap == 19,9] = colPink[5]   #GL

fullPlot[fullPlot$clusterUmap == 18,9] =  colBrown[4] #other
fullPlot[fullPlot$clusterUmap == 4,9] = colBrown[5] #other

tiff(filename = paste("~/Downloads/","OB36pixel_cluster",".tiff"),width = 4800, height = 3200, units = "px", pointsize = 12,compression="lzw");
par(bg='black',mar=c(0,0,0,0))
plot(fullPlot$xcoord,fullPlot$ycoord,col=fullPlot$clusterColor,pch=15,cex=.25,axes = F,xlab = '',ylab='')
dev.off()

box(col='gray80')
typeColorUniq = unique(subPlot[,c(6,11)])
plot(c(0,26),c(0,5),type='n')
text(1:nrow(typeColorUniq),seq(0,5,length.out = nrow(typeColorUniq)),labels = typeColorUniq$clusterUmap,col = typeColorUniq$clusterColor,cex=2)




##################################################################################################################zoom in 10um, Fig. 2D
bin10um_cluster = read.table('~/lab508/PIXEL-seq_Data&Code/mOB_dataAnalysis/clustering/OB36_bin10clusterlevel_UMAP.txt',header=F)
colnames(bin10um_cluster) = c('xcoord','ycoord','binseq','tUMI','clusterUmap','cellType','score')
sub10um = bin10um_cluster[bin10um_cluster$xcoord>24000 & bin10um_cluster$xcoord<28500 & bin10um_cluster$ycoord>5000 & bin10um_cluster$ycoord<6000,]
sub10um$celltypeCol = rep('white',nrow(sub10um))
sub10um$ctName = rep('white',nrow(sub10um))
#cell type - color assignment
colBlues <- brewer.pal(n = 9, name ="Blues")
colGreen <- brewer.pal(n = 9, name ="Greens")
colPurple <- brewer.pal(n = 9, name ="Purples")
colOrange <- brewer.pal(n = 9, name ="Oranges")
colRed <- brewer.pal(n = 9, name ="Reds")
colGreys <- brewer.pal(n = 9, name ="Greys")


sub10um[sub10um$cellType == 'OEC',8] = colBlues[6]
sub10um[sub10um$cellType == 'Mes',8] = colBlues[8]
sub10um[sub10um$cellType == 'Neuron.OSNs',8] = colBlues[4]
#---------------------
sub10um[sub10um$cellType == 'OEC',9] = 'cBlue4'
sub10um[sub10um$cellType == 'Mes',9] = 'cBlue6'
sub10um[sub10um$cellType == 'Neuron.OSNs',9] = 'cBlue8'

#pattern GL
sub10um[sub10um$cellType == 'n03.GC.1',8] = colPink[2]
sub10um[sub10um$cellType == 'n02.PGC.1',8] = colPink[3]
sub10um[sub10um$cellType == 'n05.PGC.2',8] = colPink[4]
sub10um[sub10um$cellType == 'n08.PGC.3',8] = colPink[5]
#pattern EPL
sub10um[sub10um$cellType == 'n18.EPL.IN',8] = colBrown[3]
#pattern MCL
sub10um[sub10um$cellType == 'Neuron.M.TC',8] = colGreen[8]
#---------------------

sub10um[sub10um$cellType == 'n03.GC.1',9] = 'colPink2'
sub10um[sub10um$cellType == 'n02.PGC.1',9] = 'colPink3'
sub10um[sub10um$cellType == 'n05.PGC.2',9] = 'colPink4'
sub10um[sub10um$cellType == 'n08.PGC.3',9] = 'colPink5'
sub10um[sub10um$cellType == 'n18.EPL.IN',9] = 'colBrown3'
sub10um[sub10um$cellType == 'Neuron.M.TC',9] = 'colGreen8'


sub10um[sub10um$cellType == 'n11.GC.5',8] = colPurple[8]
sub10um[sub10um$cellType == 'n12.GC.6',8] = colPurple[7]
sub10um[sub10um$cellType == 'n09.GC.3',8] = colPurple[6]
sub10um[sub10um$cellType == 'n10.GC.4',8] = colPurple[5]
sub10um[sub10um$cellType == 'n07.GC.2',8] = colPurple[4]
sub10um[sub10um$cellType == 'n14.GC.7',8] = colPurple[3]
sub10um[sub10um$cellType == 'n06.Transition',8] = colPurple[2]
#---------------------
sub10um[sub10um$cellType == 'n11.GC.5',9] = 'colPurple8'
sub10um[sub10um$cellType == 'n12.GC.6',9] = 'colPurple7'
sub10um[sub10um$cellType == 'n09.GC.3',9] = 'colPurple6'
sub10um[sub10um$cellType == 'n10.GC.4',9] = 'colPurple5'
sub10um[sub10um$cellType == 'n07.GC.2',9] = 'colPurple4'
sub10um[sub10um$cellType == 'n14.GC.7',9] = 'colPurple3'
sub10um[sub10um$cellType == 'n06.Transition',9] = 'colPurple2'


sub10um[sub10um$cellType == 'Neuron.Immature',8] = colOrange[6]
sub10um[sub10um$cellType == 'Neuron.Immature',9] = 'colOrange6'

sub10um[sub10um$cellType == 'Neuron.AstroLike',8] = colOrange[1]
sub10um[sub10um$cellType == 'Neuron.AstroLike',9] = 'colOrange1'


sub10um[sub10um$cellType == 'EC',8] = colRed[7]
sub10um[sub10um$cellType == 'ImmunoCells',8] = colRed[5]
sub10um[sub10um$cellType == 'Mural',8] = colRed[3]
#---------------------
sub10um[sub10um$cellType == 'EC',9] = 'colRed7'
sub10um[sub10um$cellType == 'ImmunoCells',9] = 'colRed5'
sub10um[sub10um$cellType == 'Mural',9] = 'colRed3'

sub10um[sub10um$cellType == 'Astro',8] = colGreys[8]
sub10um[sub10um$cellType == 'MicroG',8] = colGreys[7]
sub10um[sub10um$cellType == 'MyOligo',8] = colGreys[6]
sub10um[sub10um$cellType == 'OPC',8] = colGreys[5]
#----------------
sub10um[sub10um$cellType == 'Astro',9] = 'colGreys8'
sub10um[sub10um$cellType == 'MicroG',9] = 'colGreys7'
sub10um[sub10um$cellType == 'MyOligo',9] = 'colGreys6'
sub10um[sub10um$cellType == 'OPC',9] = 'colGreys5'

par(mfrow=c(1,1),mar=c(2,2,1,1))
plot(sub10um$xcoord,12000-sub10um$ycoord,col=sub10um$celltypeCol,pch=15,cex=.25,axes = F,xlab = '',ylab='')
box(col='gray80')



#------------------------------------10um cluster
sub10um$clusterColor = rep('white',nrow(sub10um))
sub10um[sub10um$clusterUmap == 2,10] = colBlues[4]    #ONL
sub10um[sub10um$clusterUmap == 7,10] = colBlues[8]    #ONL
sub10um[sub10um$clusterUmap == 4,10] = colBlues[6]    #ONL
sub10um[sub10um$clusterUmap == 6,10] = colBlues[3]  #ONL
sub10um[sub10um$clusterUmap == 8,10] = colBlues[2] 
sub10um[sub10um$clusterUmap == 11,10] = colBlues[1]  #ONL
sub10um[sub10um$clusterUmap == 13,10] = colBlues[3] 
sub10um[sub10um$clusterUmap == 16,10] = colBlues[2]  #ONL
sub10um[sub10um$clusterUmap == 17,10] = colBlues[1] 
sub10um[sub10um$clusterUmap == 18,10] = colBlues[3]  #ONL
sub10um[sub10um$clusterUmap == 20,10] = colBlues[2]  #ONL

sub10um[sub10um$clusterUmap == 0,10] = colPurple[8] #GC
sub10um[sub10um$clusterUmap == 5,10] = colPurple[6] #GC

sub10um[sub10um$clusterUmap == 1,10] = colPink[3]  #GL

sub10um[sub10um$clusterUmap == 3,10] = colGreys[8] #EPL

sub10um[sub10um$clusterUmap == 21,10] = colOrange[4]  #RMS

sub10um[sub10um$clusterUmap == 7,10] = colGreys[3]    # EPL
sub10um[sub10um$clusterUmap == 10,10] = colGreys[3]   # EPL
sub10um[sub10um$clusterUmap == 14,10] = colGreys[3]  # EPL
sub10um[sub10um$clusterUmap == 19,10] = colGreys[3]  # EPL
sub10um[sub10um$clusterUmap == 22,10] = colGreys[3]  # EPL

sub10um[sub10um$tUMI<512,10] = 'white'  # EPL

par(mfrow=c(1,1),mar=c(2,2,1,1))
plot(sub10um$xcoord,12000-sub10um$ycoord,col=sub10um$clusterColor,pch=15,cex=.25,axes = F,xlab = '',ylab='')
box(col='gray80')




##################################################################################################################################################################
##############################################################################################################################zoom in for segmentation
##################################################################################################################################################################

pixelrU_sc = read.table('~/lab508/novaSeqPBN/OB36/clusterPattern/OB36_subregion.txt',header=F)
colnames(pixelrU_sc) = c('xcoord','ycoord','px','py','nucleus','cellindex','binseq','tUMI','bcount','clusterUmap','cellType','score','clusterUMI')

#zoom in for pixel-level visualization
pixelrU_sc$celltypeCol = rep('white',nrow(pixelrU_sc))

pixelrU_sc[pixelrU_sc$cellType == 'OEC',14] = colBlues[6]
pixelrU_sc[pixelrU_sc$cellType == 'Mes',14] = colBlues[8]
pixelrU_sc[pixelrU_sc$cellType == 'Neuron.OSNs',14] = colBlues[4]
#---------------------
#pattern GL
pixelrU_sc[pixelrU_sc$cellType == 'n03.GC.1',14] = colPink[2]
pixelrU_sc[pixelrU_sc$cellType == 'n02.PGC.1',14] = colPink[3]
pixelrU_sc[pixelrU_sc$cellType == 'n05.PGC.2',14] = colPink[4]
pixelrU_sc[pixelrU_sc$cellType == 'n08.PGC.3',14] = colPink[5]
#pattern EPL
pixelrU_sc[pixelrU_sc$cellType == 'n18.EPL.IN',14] = colBrown[3]
#pattern MCL
pixelrU_sc[pixelrU_sc$cellType == 'Neuron.M.TC',14] = colGreen[8]
#---------------------

pixelrU_sc[pixelrU_sc$cellType == 'n11.GC.5',14] = colPurple[8]
pixelrU_sc[pixelrU_sc$cellType == 'n12.GC.6',14] = colPurple[7]
pixelrU_sc[pixelrU_sc$cellType == 'n09.GC.3',14] = colPurple[6]
pixelrU_sc[pixelrU_sc$cellType == 'n10.GC.4',14] = colPurple[5]
pixelrU_sc[pixelrU_sc$cellType == 'n07.GC.2',14] = colPurple[4]
pixelrU_sc[pixelrU_sc$cellType == 'n14.GC.7',14] = colPurple[3]
pixelrU_sc[pixelrU_sc$cellType == 'n06.Transition',14] = colPurple[2]
#---------------------
pixelrU_sc[pixelrU_sc$cellType == 'Neuron.Immature',14] = colOrange[6]
pixelrU_sc[pixelrU_sc$cellType == 'Neuron.AstroLike',14] = colOrange[1]
pixelrU_sc[pixelrU_sc$cellType == 'EC',14] = colRed[7]
pixelrU_sc[pixelrU_sc$cellType == 'ImmunoCells',14] = colRed[5]
pixelrU_sc[pixelrU_sc$cellType == 'Mural',14] = colRed[3]
#---------------------
pixelrU_sc[pixelrU_sc$cellType == 'Astro',14] = colGreys[8]
pixelrU_sc[pixelrU_sc$cellType == 'MicroG',14] = colGreys[7]
pixelrU_sc[pixelrU_sc$cellType == 'MyOligo',14] = colGreys[6]
pixelrU_sc[pixelrU_sc$cellType == 'OPC',14] = colGreys[5]
#----------------

###################################################################################################
#$$$$$$$$$$$$$$$$$$###############Example of cell boudary for the possible marker region
###################################################################################################
#marker gene
marker = read.table('~/lab508/novaSeqPBN/OB36/OB36_sampleRegion.final',header=F)
colnames(marker) = c('fovC','fovR','xcoord','ycoord','umi','spt','alignCode','chr','start','alignQ','uid','geneid','geneName','bioType','intronFrac')

plot(chooseRegion$pixelX,chooseRegion$pixelY,col=rgb(0,0,0,chooseRegion$tUMI/5000),cex=2,pch=15,axes=F,xlim=c(38500,38775),ylim=c(7765,7955))
points(marker[marker$geneName=='Apoe',]$xcoord,marker[marker$geneName=='Apoe',]$ycoord,col=rgb(1,0,0,0.3),cex=1.5,pch=15)
points(marker[marker$geneName=='Calb2',]$xcoord,marker[marker$geneName=='Calb2',]$ycoord,col=rgb(0,1,0,0.3),cex=1.5,pch=15)


################################################################################################################
#############################################################spatial pattern analysis of 10um segmentation
################################################################################################################
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

expr <- data.table::fread(
  input = "~/lab508/novaSeqPBN/OB36/OB36bin33dCellGeneExp.txt",
  sep = ",",
  data.table = FALSE
)

rownames(x = expr) <- (x = expr[, 1])
colnames(x = expr) <- (x = colnames(x = expr))
expr <- t(x = expr[, -1])
ssmOLB10 <- CreateSeuratObject(
  counts = expr,
  project = 'SlideSeq',
  assay = 'Spatial'
)
positions <- read.csv(
  file = "~/lab508/novaSeqPBN/OB36/OB36bin33Coordd.txt"
)

nn1 = positions[positions$xcoord>24000 & positions$xcoord <40000,]
#nn1$xcoord = nn1$xcoord - 24000 


rownames(x = nn1) <- nn1$binID
positions <- nn1[, 2:3]

ssmOLB10[['image']] <-new(
  Class = 'SlideSeq',
  assay = "Spatial",
  coordinates = positions
)

#ssmOLB = ssmOLB[, ssmOLB$nCount_Spatial > 512]
#normalization and dimensions reduction
ssmOLB10 <- SCTransform(ssmOLB10, assay = "Spatial", ncells = 2000, verbose = TRUE)

ssmOLB10 <- RunPCA(ssmOLB10,assay = "SCT", verbose = FALSE)
ElbowPlot(ssmOLB10)
#spatial pattern analysis
ssmOLB10 <- FindNeighbors(ssmOLB10, reduction = "pca", dims = 1:15)
ssmOLB10 <- FindClusters(ssmOLB10, verbose = FALSE,resolution = 0.4)
ssmOLB10 <- RunUMAP(ssmOLB10, reduction = "pca", dims = 1:15)
#direct pattern Umap and plot
plot1 <- DimPlot(ssmOLB10, reduction = "umap", label = TRUE)
plot2 <- SpatialDimPlot(ssmOLB10, label.size = 1)
plot1 + plot2

DimPlot(ssmOLB10, reduction = "umap", cols = mycol)

spaCluster = plot2$data
mycolset = c(colBlues[c(6,6,6,6,6,6,6,6,6,6,6)],colPurple[c(7,7)],colPink[c(3)],colGreys[c(5)],colOrange[4],colGreys[c(5,5,5,5,5)])
myOrder = c(2,7,4,6,8,11,13,16,17,18,20,0,5,1,3,21,7,10,14,19,22)
par(mfrow=c(4,6),mar=c(1,1,1,1))
for(i in 1:23){
  plot(spaCluster[,1:2],type='n',axes = F)
  points(spaCluster[spaCluster$ident==myOrder[i],1:2],col=mycolset[i],pch=19,cex=.1)
}
