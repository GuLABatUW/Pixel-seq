
Region = allinfor[allinfor$xcoord>32000 & allinfor$xcoord<34000 & allinfor$ycoord>2000 & allinfor$ycoord<5000,]
Region = Region[order(Region$`Hba-a1`,decreasing = T),]
Region[1:10,]$X

gene2dis = read.table('~/lab508/novaSeqPBN/subcellular/genePosDistribution/gene2posDis.txt')
colnames(gene2dis) = c('gene','cellid','geneX','geneY','xcoord','ycoord','pos')
geneCount = table(gene2dis$gene)
geneCount = geneCount[geneCount>128]
geneCount = geneCount[order(geneCount,decreasing = T)]


#spatial plot of two chosed genes
plot(PBS2_sel[,1:2],col=rgb(0,0,0,PBS2_sel$UMIcounts/1200),cex=1,pch=15)
points(pbs2_calca[,3:4],col=rgb(0,0.8,0,0.15),cex=0.8,pch=15)
points(gene2dis[gene2dis$gene=='Cdc42',3:4],pch=15,cex=.8,col=rgb(1,0,0,0.3))


cellStat = table(gene2dis[gene2dis$gene=='Cdc42',]$cellid)
cellStat = cellStat[order(cellStat,decreasing = T)]
cellStat = cellStat[cellStat>10]
breaks <- seq(0,1,0.2)
# specify interval/bin labels
tags <- 1:5
# bucketing values into bins
for(i in 1:length(cellStat)){
  group_tags <- cut(gene2dis[gene2dis$gene=='Cdc42' & gene2dis$cellid == names(cellStat)[i],]$pos, 
                    breaks=breaks, 
                    include.lowest=TRUE, 
                    right=FALSE, 
                    labels=tags)
  if(i ==1){
    statHba = as.data.frame(table(group_tags))
  } else {
    statHba = cbind(statHba,as.data.frame(table(group_tags)))
  }
}
statHba = statHba[,seq(2,2*length(cellStat),2)]
colnames(statHba) = names(cellStat)[1:length(cellStat)]
statHba = apply(statHba, 2, function(x) x/sum(x))
meanHba = apply(statHba,1,mean)
sdHba = apply(statHba,1,std)

cellStat = table(gene2dis[gene2dis$gene=='Calca',]$cellid)
cellStat = cellStat[order(cellStat,decreasing = T)]
cellStat = cellStat[cellStat>10]
breaks <- seq(0,1,0.2)
# specify interval/bin labels
tags <- 1:5
# bucketing values into bins
for(i in 1:length(cellStat)){
  group_tags <- cut(gene2dis[gene2dis$gene=='Calca' & gene2dis$cellid == names(cellStat)[i],]$pos, 
                    breaks=breaks, 
                    include.lowest=TRUE, 
                    right=FALSE, 
                    labels=tags)
  if(i ==1){
    statCalca = as.data.frame(table(group_tags))
  } else {
    statCalca = cbind(statCalca,as.data.frame(table(group_tags)))
  }
}
statCalca = statCalca[,seq(2,2*length(cellStat),2)]
colnames(statCalca) = names(cellStat)[1:length(cellStat)]
statCalca = apply(statCalca, 2, function(x) x/sum(x))
meanCalca = apply(statCalca,1,mean)
sdCalca = apply(statCalca,1,std)


calcaDis = data.frame(pos = 1:5,meanSignal = meanCalca, stdSignal = sdCalca)
HbaDis = data.frame(pos = 1:5,meanSignal = meanHba, stdSignal = sdHba)

ggplot(calcaDis) +
  geom_bar( aes(x=pos, y=meanSignal), stat="identity", fill='green') +
  geom_errorbar( aes(x=pos, ymin=meanSignal-stdSignal, ymax=meanSignal+stdSignal), width=0.4, colour="black", alpha=0.9, size=0.8) +
  ggtitle("using standard error")+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ ylim(0, 0.6)

ggplot(HbaDis) +
  geom_bar( aes(x=pos, y=meanSignal), stat="identity", fill='red') +
  geom_errorbar( aes(x=pos, ymin=meanSignal-stdSignal, ymax=meanSignal+stdSignal), width=0.4, colour="black", alpha=0.9, size=0.8) +
  ggtitle("using standard error")+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ ylim(0, 0.6)

par(mfrow=c(1,2))
barplot(meanCalca,col = rgb(0,1,0,0.8),border = NA,ylim=c(0,0.5),las=1)
barplot(meanHba,col = rgb(1,0,0,0.8),border = NA,ylim=c(0,0.5),las=1)
axis(2,at=c(0,0.1,0.2,0.3,0.4,0.5),labels = c('','','','','',''))
