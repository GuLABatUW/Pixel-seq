#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

keyset = args[1]
posfile = args[2]
palgorthm = args[3]
outpath = args[4]

suppressPackageStartupMessages(library(igraph))
pbEdge = read.table(paste(outpath,posfile,sep=''),header=F)
outNet = graph.edgelist(as.matrix(pbEdge[,1:2]),directed=F)
E(outNet)$weight=as.numeric(pbEdge$V7)

if(palgorthm == 'edgeBetweenness'){
  wc <- edge.betweenness.community(outNet,weights = E(outNet)$value, bridges=TRUE)
} else {
  wc = fastgreedy.community(outNet,weights = E(outNet)$weight)
}

tot = unique(rbind(as.matrix(pbEdge[,c(1,3,4)]),as.matrix(pbEdge[,c(2,5,6)])))

write.csv(cbind(tot[order(tot[,1]),],rep(keyset,length(wc$membership)),wc$membership),"cellAsign2.csv")


