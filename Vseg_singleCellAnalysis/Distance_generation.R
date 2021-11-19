# import library
library(FNN)
library(igraph)
OB36 = read.table('OB36suniq.coord',header=T, sep = ',')
nnDist_OB36 = get.knnx(OB36[,2:3],OB36[,2:3],k=20) # use Knn to get 20 nearest barcode for each barcode
options(scipen=999)
write.csv(cbind(nnDist_OB36$nn.index,round(nnDist_OB36$nn.dist,2)),"OB36_nnAll.txt")
