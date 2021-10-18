
##############N-N distance calculation example
PBS3 = read.table('~/lab508/novaSeqPBN/PBNall/PBS3uniq.coord',header=T)
plot(PBS3$CoordX,PBS3$CoordY,col=rgb(1,0,0,PBS3$UMIcounts/2000),pch=15,cex=.3)
#PBS3_30d1 = PBS3[PBS3$CoordX>6000 & PBS3$CoordX<16000,]
PBS3_30d2 = PBS3[PBS3$CoordX>21000 & PBS3$CoordX<31000,]
PBS3_CN2 = PBS3[PBS3$CoordX>36000 & PBS3$CoordX<47000,]

nnDist_PBS3_30d2 = get.knnx(PBS3_30d2[,1:2],PBS3_30d2[,1:2],k=20)
options(scipen=999)
write.csv(PBS3_30d2,"~/lab508/novaSeqPBN/PBNall/PBS3_30d2uniq.coord")
write.csv(cbind(nnDist_PBS3_30d2$nn.index,round(nnDist_PBS3_30d2$nn.dist,2)),"~/lab508/novaSeqPBN/PBNall/PBS3_30d2_nnAll.txt")

nnDist_PBS3_CN2 = get.knnx(PBS3_CN2[,1:2],PBS3_CN2[,1:2],k=20)
options(scipen=999)
write.csv(PBS3_CN2,"~/lab508/novaSeqPBN/PBNall/PBS3_CN2uniq.coord")
write.csv(cbind(nnDist_PBS3_CN2$nn.index,round(nnDist_PBS3_CN2$nn.dist,2)),"~/lab508/novaSeqPBN/PBNall/PBS3_CN2_nnAll.txt")


