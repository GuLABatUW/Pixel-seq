#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

keyset = args[1]
posfile = args[2]
palgorthm = args[3]

suppressPackageStartupMessages(library(igraph))
pbEdge = read.table(paste('~/lab508/novaSeqPBN/OB36/',posfile,sep=''),header=T)


outNet = graph.edgelist(as.matrix(pbEdge[,1:2]),directed=F)

#pEweight = pbEdge$V7/(log2(pbEdge$V10+1)+log2(pbEdge$V11+1))
#pEweight = 1/(log2(pbEdge$V8+1)+log2(pbEdge$V9+1))  #ws1
#pEweight = 1/log2(pbEdge$V8+pbEdge$V9+1)
#pEweight = 1/(pbEdge$V8+pbEdge$V9)

E(outNet)$weight=round(pbEdge$weight,2)

#E(outNet)$weight=as.numeric(pbEdge$V11)
#E(outNet)$weight=as.numeric(pbEdge$V8/(log2(pbEdge$V14+1)+pbEdge$V12))

if(palgorthm == 'edgeBetweenness'){
  wc <- edge.betweenness.community(outNet,weights = E(outNet)$weight, bridges=TRUE)
} else {
  wc = fastgreedy.community(outNet,weights = E(outNet)$weight)
}

tot = unique(rbind(as.matrix(pbEdge[,c(1,3,4)]),as.matrix(pbEdge[,c(2,5,6)])))

write.csv(cbind(tot[order(tot[,1]),],rep(keyset,length(wc$membership)),wc$membership),"~/lab508/novaSeqPBN/OB36/cellAsign3.csv")


