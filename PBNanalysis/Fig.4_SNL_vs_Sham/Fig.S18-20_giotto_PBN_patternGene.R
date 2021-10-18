##################giotto
library(Giotto)
library(DropSeq.util)
#1. create & process Giotto object
PBS2 <- loadSparseDge("~/lab508/novaSeqPBN/PBNall/PBN_fullp2_cellExp.txt")
PBS3 <- loadSparseDge("~/lab508/novaSeqPBN/PBNall/PBN_fullp2_S3_cellExp.txt")
posS2 <- read.csv(
  file = "~/lab508/novaSeqPBN/PBNall/PBN_fullp2_S2_coord.txt"
)
posS3 <- read.csv(
  file = "~/lab508/novaSeqPBN/PBNall/PBN_fullp2_S3_coord.txt"
)

rownames(posS2) = posS2$binID
posS2 = posS2[,2:3]
pos30d1 = posS2[posS2$xcoord<26000,]
posCnr1 = posS2[posS2$xcoord>28000,]
exp30d1 = PBS2[,rownames(pos30d1)]
expCnr1 = PBS2[,rownames(posCnr1)]

rownames(posS3) = posS3$binID
posS3 = posS3[,2:3]
pos30d2 = posS3[posS3$xcoord<32000,]
posCnr2 = posS3[posS3$xcoord>32000,]
exp30d2 = PBS3[,rownames(pos30d2)]
expCnr2 = PBS3[,rownames(posCnr2)]

my_giotto_30d1 = createGiottoObject(raw_exprs = exp30d1,
                                      spatial_locs = pos30d1)
my_giotto_30d1 <- filterGiotto(gobject = my_giotto_30d1, 
                                 expression_threshold = 1, 
                                 gene_det_in_min_cells = 20, 
                                 min_det_genes_per_cell = 64)


my_giotto_30d1 <- normalizeGiotto(gobject = my_giotto_30d1, scalefactor = 6000, verbose = T)
my_giotto_30d1 <- addStatistics(gobject = my_giotto_30d1)
my_giotto_30d1 <- adjustGiottoMatrix(gobject = my_giotto_30d1, 
                                       expression_values = c('normalized'),
                                       covariate_columns = c('nr_genes', 'total_expr'))

#----3. dimension reduction
my_giotto_30d1 <- calculateHVG(gobject = my_giotto_30d1)


#---8. spatial network

#visualize information about the default Delaunay network
#create a spatial Delaunay network (default)
#create a spatial kNN network

plotStatDelaunayNetwork(gobject = my_giotto_30d1, maximum_distance = 100)
my_giotto_30d1 = createSpatialNetwork(gobject = my_giotto_30d1, minimum_k = 2, 
                                        maximum_distance_delaunay = 100)
my_giotto_30d1 = createSpatialNetwork(gobject = my_giotto_30d1, minimum_k = 2, 
                                        method = 'kNN', k = 10)
showNetworks(my_giotto_30d1)

#----9. spatial genes

#Visualize top 4 genes per method.

km_spatialgenes30d1 = binSpect(my_giotto_30d1)
rank_spatialgenes30d1 = binSpect(my_giotto_30d1, bin_method = 'rank')

#------12. cell neighborhood: cell-type/cell-type interactions

saveRDS(my_giotto_30d1,'~/lab508/novaSeqPBN/spatial/mygiotto30d1.rds')
#my_giotto_30d1 = readRDS('~/lab508/novaSeqPBN/spatial/mygiotto30d1.rds')

#######################################################################
#---------------------------------------------------------------------
#####################################################################
my_giotto_cnr1 = createGiottoObject(raw_exprs = expCnr1,
                                    spatial_locs = posCnr1)
my_giotto_cnr1 <- filterGiotto(gobject = my_giotto_cnr1, 
                               expression_threshold = 1, 
                               gene_det_in_min_cells = 20, 
                               min_det_genes_per_cell = 64)


my_giotto_cnr1 <- normalizeGiotto(gobject = my_giotto_cnr1, scalefactor = 6000, verbose = T)
my_giotto_cnr1 <- addStatistics(gobject = my_giotto_cnr1)
my_giotto_cnr1 <- adjustGiottoMatrix(gobject = my_giotto_cnr1, 
                                     expression_values = c('normalized'),
                                     covariate_columns = c('nr_genes', 'total_expr'))



#----3. dimension reduction
my_giotto_cnr1 <- calculateHVG(gobject = my_giotto_cnr1)


#---8. spatial network

#visualize information about the default Delaunay network
#create a spatial Delaunay network (default)
#create a spatial kNN network

plotStatDelaunayNetwork(gobject = my_giotto_cnr1, maximum_distance = 100)
my_giotto_cnr1 = createSpatialNetwork(gobject = my_giotto_cnr1, minimum_k = 2, 
                                      maximum_distance_delaunay = 100)
my_giotto_cnr1 = createSpatialNetwork(gobject = my_giotto_cnr1, minimum_k = 2, 
                                      method = 'kNN', k = 10)
showNetworks(my_giotto_cnr1)


#----9. spatial genes

#Visualize top 4 genes per method.

km_spatialgenesCNR1 = binSpect(my_giotto_cnr1)
rank_spatialgenesCNR1 = binSpect(my_giotto_cnr1, bin_method = 'rank')

#------12. cell neighborhood: cell-type/cell-type interactions

saveRDS(my_giotto_cnr1,'~/lab508/novaSeqPBN/spatial/mygiottoCNR1.rds')
#my_giotto_cnr1 = readRDS('~/lab508/novaSeqPBN/spatial/mygiottoCNR1.rds')





##################################################
#------------------------------------------------
################################################


my_giotto_30d2 = createGiottoObject(raw_exprs = exp30d2,
                                    spatial_locs = pos30d2)
my_giotto_30d2 <- filterGiotto(gobject = my_giotto_30d2, 
                               expression_threshold = 1, 
                               gene_det_in_min_cells = 20, 
                               min_det_genes_per_cell = 64)


my_giotto_30d2 <- normalizeGiotto(gobject = my_giotto_30d2, scalefactor = 6000, verbose = T)
my_giotto_30d2 <- addStatistics(gobject = my_giotto_30d2)
my_giotto_30d2 <- adjustGiottoMatrix(gobject = my_giotto_30d2, 
                                     expression_values = c('normalized'),
                                     covariate_columns = c('nr_genes', 'total_expr'))

#----3. dimension reduction
my_giotto_30d2 <- calculateHVG(gobject = my_giotto_30d2)


#---8. spatial network

#visualize information about the default Delaunay network
#create a spatial Delaunay network (default)
#create a spatial kNN network

plotStatDelaunayNetwork(gobject = my_giotto_30d2, maximum_distance = 100)
my_giotto_30d2 = createSpatialNetwork(gobject = my_giotto_30d2, minimum_k = 2, 
                                      maximum_distance_delaunay = 100)
my_giotto_30d2 = createSpatialNetwork(gobject = my_giotto_30d2, minimum_k = 2, 
                                      method = 'kNN', k = 10)
showNetworks(my_giotto_30d2)

#----9. spatial genes

#Visualize top 4 genes per method.

km_spatialgenes30d2 = binSpect(my_giotto_30d2)
rank_spatialgenes30d2 = binSpect(my_giotto_30d2, bin_method = 'rank')




#######################################################
#-----------------------------------------------------
#######################################################
my_giotto_cnr2 = createGiottoObject(raw_exprs = expCnr2,
                                    spatial_locs = posCnr2)
my_giotto_cnr2 <- filterGiotto(gobject = my_giotto_cnr2, 
                               expression_threshold = 1, 
                               gene_det_in_min_cells = 20, 
                               min_det_genes_per_cell = 64)


my_giotto_cnr2 <- normalizeGiotto(gobject = my_giotto_cnr2, scalefactor = 6000, verbose = T)
my_giotto_cnr2 <- addStatistics(gobject = my_giotto_cnr2)
my_giotto_cnr2 <- adjustGiottoMatrix(gobject = my_giotto_cnr2, 
                                     expression_values = c('normalized'),
                                     covariate_columns = c('nr_genes', 'total_expr'))



#----3. dimension reduction
my_giotto_cnr2 <- calculateHVG(gobject = my_giotto_cnr2)


#---8. spatial network

#visualize information about the default Delaunay network
#create a spatial Delaunay network (default)
#create a spatial kNN network

plotStatDelaunayNetwork(gobject = my_giotto_cnr2, maximum_distance = 100)
my_giotto_cnr2 = createSpatialNetwork(gobject = my_giotto_cnr2, minimum_k = 2, 
                                      maximum_distance_delaunay = 100)
my_giotto_cnr2 = createSpatialNetwork(gobject = my_giotto_cnr2, minimum_k = 2, 
                                      method = 'kNN', k = 10)
showNetworks(my_giotto_cnr2)


#----9. spatial genes

#Visualize top 4 genes per method.

km_spatialgenesCNR2 = binSpect(my_giotto_cnr2)

rank_spatialgenesCNR2 = binSpect(my_giotto_cnr2, bin_method = 'rank')

#------12. cell neighborhood: cell-type/cell-type interactions

saveRDS(my_giotto_cnr2,'~/lab508/novaSeqPBN/spatial/mygiottoCNR2.rds')
#my_giotto_cnr2 = readRDS('~/lab508/novaSeqPBN/spatial/mygiottoCNR2.rds')


library(ggvenn)
x <- list(
  A = rank_spatialgenesCNR1$genes[1:500], 
  B = rank_spatialgenes30d1$genes[1:500], 
  C = rank_spatialgenesCNR2$genes[1:500],
  D = rank_spatialgenes30d2$genes[1:500]
)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4,
  show_percentage = F 
)


commGene = intersect(intersect(rank_spatialgenesCNR1$genes[1:500],rank_spatialgenesCNR2$genes[1:500]),intersect(rank_spatialgenes30d1$genes[1:500],rank_spatialgenes30d2$genes[1:500]))
mCNR1 = as.data.frame(rank_spatialgenesCNR1[rank_spatialgenesCNR1$genes %in% commGene,])
mCNR2 = as.data.frame(rank_spatialgenesCNR2[rank_spatialgenesCNR2$genes %in% commGene,])
m30d1 = as.data.frame(rank_spatialgenes30d1[rank_spatialgenes30d1$genes %in% commGene,])
m30d2 = as.data.frame(rank_spatialgenes30d2[rank_spatialgenes30d2$genes %in% commGene,])

combineAllComm = cbind(mCNR1[order(mCNR1$genes),c(1,4:6)],
                       m30d1[order(m30d1$genes),4:6],
                       mCNR2[order(mCNR2$genes),4:6],
                       m30d2[order(m30d2$genes),4:6])


write.csv(combineAllComm,"~/Dropbox/novaSeq/figures/PBN_pain/giottoCommu/commPatternGene.txt")


######################################################################################################

###############spatial pattern gene
geneAbundance = rowSums(PBNct.combined.sct)


genelist = c('Mia','Uts2','Calca','Tac1','Uts2b','Gal','Nts','Sst','Grp','Ly86','Penk','Pltp','Psap','Calcb','Cd74','Crtam','Spp1','Dcn','Rarres2','Vtn',
             'Gabra6','Gabrd','Gng13','Prkcg','Cdhr4','Calb1','Calb2','Hpcal1','Mgp','Pvalb','Sncb','Vsnl1','Prkcg','Tuba1a','Fos','Slc10a4','Slc16a11',
             'Slc17a7','Slc18a3','Slc22a4')

PBS2_30d1 = PBS2[PBS2$CoordX>18200 & PBS2$CoordX<24200 & PBS2$CoordY>3000,]
PBS2_CN1 = PBS2[PBS2$CoordX>30000 & PBS2$CoordX<36000 & PBS2$CoordY>2000,]


for(i in 1:length(genelist)){
  tiff(paste('~/Dropbox/novaSeq/figures/Version3/patternGene/',genelist[i],'.tiff',sep=''),width = 1000,height = 400)
  par(mfrow=c(1,2),mar=c(1,1,1,1))
  system(paste("grep -P ",'\t',genelist[i],'\t ','~/lab508/novaSeqPBN/PBNall/PBS2.final','>~/lab508/novaSeqPBN/gene.final'))
  geneinfor = read.table('~/lab508/novaSeqPBN/gene.final',header=F)
  plot(PBS2_CN1[,1:2],cex=.1,col=rgb(0.1,0.1,0.1,PBS2_CN1$UMIcounts/6300),pch=19,axes=F)
  points(geneinfor[geneinfor$V3>30000 & geneinfor$V3<36000 & geneinfor$V4>2000,3:4],pch=15,cex=.8,col=rgb(0.5,0.1,0.8,0.2/1.4))
  plot(PBS2_30d1[,1:2],cex=.1,col=rgb(0.1,0.1,0.1,PBS2_30d1$UMIcounts/4500),pch=19,axes=F)
  points(geneinfor[geneinfor$V3>18200 & geneinfor$V3<24200 & geneinfor$V4>3000,3:4],pch=15,cex=.8,col=rgb(0.5,0.1,0.8,0.2/1.0))
  dev.off()
}
