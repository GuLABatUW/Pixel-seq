#SPT analysis of RNAseqfish+
#seqFish+ source data downloaded from:
#https://github.com/CaiGroup/seqFISH-PLUS

RNAlocus = t(read.csv('~/Downloads/sourcedata/RNA_locations_cell_cortex/RNA_locations_cell_1.csv',header=F))
allRNAlocus = cbind(rep(1,nrow(RNAlocus)),RNAlocus)

for(i in 2:119){
  name = paste("RNA_locations_cell_",i,".csv",sep = "")
  rnalocus = t(read.csv(paste('~/Downloads/sourcedata/RNA_locations_cell_cortex/',name,sep=""),header=F))
  rnalocus = cbind(rep(i,nrow(rnalocus)),rnalocus)
  allRNAlocus = rbind(allRNAlocus,rnalocus)
}

library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colall = sample(rep(col_vector,2))
plot(allRNAlocus[allRNAlocus[,1]<119,3:4],col=colall[allRNAlocus[,1]],pch=15,cex=.5,xlim=c(0,2000),ylim=c(0,2000))
for(i in 1:118){
  select = allRNAlocus[allRNAlocus[,1]==i,]
  text(median(select[,3]),median(select[,4]),i,cex=1,col=1)
}


pixelDis = matrix(1:40000,200,200)
for(i in 1:200){
  a = (i-1)*10+1
  b = i*10
  for(j in 1:200){
    c = (j-1)*10+1
    d = j*10
    pixelDis[i,j] = nrow(as.matrix(allRNAlocus[allRNAlocus[,3]>a & allRNAlocus[,3]<b & allRNAlocus[,4]>c & allRNAlocus[,4]<d,]))
  }
}

plot(1:200,1:200,type="n")
for(i in 1:200){
  for(j in 1:200){
    points(i,j,cex=.5,pch=15,col=hsv(0,pixelDis[i,j]/max(pixelDis),pixelDis[i,j]/max(pixelDis)))
  }
}


library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colall = sample(rep(col_vector,2))
plot(allRNAlocus[,3:4],col=colall[allRNAlocus[,1]],pch=19,cex=.1)



plot(1:200,1:200,type="n")
for(i in 1:200){
  for(j in 1:200){
    points(i,j,cex=pixelDis[i,j]/max(pixelDis),pch=0)
  }
}


i =17
j = 36
a = (i-1)*10+1
b = i*10
c = (j-1)*10+1
d = j*10
subRNA = allRNAlocus[allRNAlocus[,3]>a & allRNAlocus[,3]<b & allRNAlocus[,4]>c & allRNAlocus[,4]<d,]
dim(subRNA)


cortox_SVZ = read.csv('~/Downloads/sourcedata/cortex_svz_counts.csv')




################check the ref RNA total sum dist distribution


library('FNN')

gene = 3444

coordA = selCell[selCell[,2]==gene,]
coordB = selCell[selCell[,2]!=gene,]

nn = get.knnx(selCell[,3:4],selCell[,3:4],2)
plot(density(nn$nn.dist[,2]),col=1,lwd=1.5,ylim=c(0,0.2))


##########################################circlar packing
#bead random distribution and coverage simulation

library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colall = rep(col_vector,73)

library(packcircles)
library(plotrix)
library('RANN')

#loading the data
RNAlocus = t(read.csv('~/Downloads/sourcedata/RNA_locations_cell_cortex/RNA_locations_cell_1.csv',header=F))
allRNAlocus = cbind(rep(1,nrow(RNAlocus)),RNAlocus)

for(i in 2:119){
  name = paste("RNA_locations_cell_",i,".csv",sep = "")
  rnalocus = t(read.csv(paste('~/Downloads/sourcedata/RNA_locations_cell_cortex/',name,sep=""),header=F))
  rnalocus = cbind(rep(i,nrow(rnalocus)),rnalocus)
  allRNAlocus = rbind(allRNAlocus,rnalocus)
}

beadssimulation = 10
beadDiameter=c(seq(0.2,1,0.1),seq(2,10,1))/2
beadDiameter=rev(beadDiameter)
beadNumPerFov = matrix(rep(0,length(beadDiameter)*beadssimulation),beadssimulation,length(beadDiameter))
beadAreaCovPerFov = matrix(rep(0,length(beadDiameter)*beadssimulation),beadssimulation,length(beadDiameter))
beadmRNACovPerFov = matrix(rep(0,length(beadDiameter)*beadssimulation),beadssimulation,length(beadDiameter))

for(j in 1:length(beadDiameter)){
  beadSize = beadDiameter[j]
  for(i in 1:beadssimulation){
    
    #bead packing part
    beadSd = beadSize*0.10
    beadNum = 200*200/pi/beadSize/beadSize*5
    radius <- rnorm(beadNum,beadSize,beadSd)
    packing <- circleProgressiveLayout(radius,sizetype = "radius")
    packing$x = packing$x-min(packing$x)/2
    packing$y = packing$y-min(packing$y)/2
    beads = packing[packing$x>0 & packing$x<200 & packing$y>0 & packing$y<200,]
    beadNumPerFov[i,j] = nrow(beads)
    beadAreaCovPerFov[i,j] = sum(pi*beads$radius*beads$radius)/200/200
    
    #mRNA assign part
    beadsRNAlocus = cbind(allRNAlocus,rep(0,nrow(allRNAlocus)))
    beadsRNAlocus[,3]=beadsRNAlocus[,3]*0.103/2
    beadsRNAlocus[,4]=beadsRNAlocus[,4]*0.103/2
    for (k in 1:nrow(beads)){
      #radical scan by radius
      a = nn2(beadsRNAlocus[,3:4],beads[k,1:2],radius = beads[2,3],k=nrow(beadsRNAlocus),searchtype="radius")
      b = a$nn.idx[a$nn.idx>0]
      beadsRNAlocus[b,5] = k
    }
    beadmRNACovPerFov[i,j] = nrow(beadsRNAlocus[beadsRNAlocus[,5]>0,])/nrow(beadsRNAlocus)
    #output simulated file
    mRNAFileName = paste("~/lab508/sptSimulation/simulatedBead_",beadSize,"um_rep",i,".txt",sep = "")
    beadFileName = paste("~/lab508/sptSimulation/simulatedBinfor_",beadSize,"um_rep",i,".txt",sep = "")
    write.table(beadsRNAlocus,mRNAFileName)
    write.table(beads,beadFileName)
  }
}









plot(1:200,1:200,type="n")
for(i in 1:nrow(beads)){
  draw.circle(beads[i,1],beads[i,2],radius = beads[i,3]*0.91,col="NA",border = "gray60")
}

for(i in 1:nrow(packing)){
  draw.circle(packing[i,1],packing[i,2],radius = packing[i,3]*0.91,col=colall[i],border = "NA")
}


beadDiameter_10um = beadDiameter
beadNumPerFov_10um = beadNumPerFov

beadDiameter_5um = beadDiameter
beadNumPerFov_5um = beadNumPerFov


beadDiameter_2.5um = beadDiameter
beadNumPerFov_2.5um = beadNumPerFov

beadDiameter_1um = beadDiameter
beadNumPerFov_1um = beadNumPerFov

#gene assigned to beads by locus
#convert to 1um per unit


library('FNN')
for (i in 1:nrow(beads)){
  #squreSelection
  a = nn2(beadsRNAlocus[,3:4],beads[i,1:2],radius = beads[2,3],k=nrow(beadsRNAlocus),searchtype="radius")
  b = a$nn.idx[a$nn.idx>0]
  beadsRNAlocus[b,5] = i
}




plot(1:200,1:200,type="n")
for(i in 1:nrow(beads)){
  draw.circle(beads[i,1],beads[i,2],radius = beads[i,3]*0.91,col=colall[i],border = "NA")
}



#stored information
beads1um
beads1um_RNAlocus



#overview of the beads-captured reads distribution
beadsinfor = read.table('~/lab508/sptSimulation/simulatedBinfor_0.5um_rep1.txt')
beadsRNA = read.table('~/lab508/sptSimulation/simulatedBead_0.5um_rep1.txt',skip = 1)

beadsReads = as.matrix(table(beadsRNA[,6]))
beadsReads = cbind(as.numeric(rownames(beadsReads)),beadsReads[,1])
beadsCol = rep(0,nrow(beads))

for(i in 1:nrow(beadsReads)){
  if(beadsReads[i,1] > 0){
    beadsCol[beadsReads[i,1]] = beadsReads[i,2]
  }
}
beadsCol = log2(beadsCol+1)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
beadsCol = range01(beadsCol)

par(bg = 'black')
plot(1:100,1:100,type="n")
for(i in 1:nrow(beadsinfor)){
  draw.circle(beadsinfor[i,1],beadsinfor[i,2],radius = beadsinfor[i,3]*0.9,col=rgb(1,0,0,beadsCol[i]))
}
box(col='white')



#########################beads purity

#calculate the reads number per cell

CellReads = as.matrix(table(allRNAlocus[,1]))
CellReads = cbind(as.numeric(rownames(CellReads)),CellReads[,1])


k = c(5,4.5,4,3.5,3,2.5,2,1.5,1,0.5)
meanCellCross = matrix(rep(0,40),10,4)
sdCellCross = matrix(rep(0,40),10,4)
for (m in 1:10){
  nn = m*2
  CrossCell = matrix(rep(0,40),10,4)
  for(j in 1:10){
    name = paste('~/lab508/sptSimulation/simulatedBead_',k[m],'um_rep',j,".txt",sep = "")
    #beadsinfor = read.table('~/lab508/sptSimulation/simulatedBinfor_5um_rep1.txt')
    beadsRNA = read.table(name,skip = 1)
    
    beadsToCell = as.matrix(t(table(beadsRNA[,c(2,6)])))
    
    beadsStat = matrix(rep(0,(nrow(beadsToCell)-1)*6),(nrow(beadsToCell)-1),6)
    #col1: beadsID, beadsTotReads, MaxCellReads, MaxCellToBeadPer, MaxCellToCellPer,TotCell
    
    for(i in 2:nrow(beadsToCell)){
      nc1 = i-1
      nc2 = sum(beadsToCell[i,])
      nc3 = max(beadsToCell[i,])
      nc4 = nc3/nc2
      mCellID = which.max(beadsToCell[i,])
      nc5 = nc3/CellReads[mCellID,2]
      nc6 = sum(beadsToCell[i,]>0)
      beadsStat[i-1,]=c(nc1,nc2,nc3,nc4,nc5,nc6)
    }
    
    CrossCell[j,1:4] = table(beadsStat[,6])[1:4]
  }
  CrossCell[is.na(CrossCell)] = 0
  meanCellCross[m,] = apply(CrossCell/apply(CrossCell,1,sum),2,mean)
  sdCellCross[m,] = apply(CrossCell/apply(CrossCell,1,sum),2,sd)
  #name2 = paste("~/Dropbox/crossCell",nn,"um.txt",sep = "")
  #write.table(CrossCell,name2)
}

meanCellCross = as.data.frame(meanCellCross)
sdCellCross = as.data.frame(sdCellCross)

library("RColorBrewer")
colp = brewer.pal(10,"Paired")
plot(c(0,10),c(0,1),type="n",axes = F,xlab='beadSize (um)',ylab='Percentage of Crossed Cells')
axis(side=1,at=1:10,labels=10:1)
axis(side=2,at=seq(0,1,0.2),labels=seq(0,1,0.2),las=1)
box()

polygon(c(1:10,10:1),c(meanCellCross$V1-sdCellCross$V1,rev(meanCellCross$V1+sdCellCross$V1)),col = colp[1],border = NA)
lines(c(1:10),meanCellCross$V1,col=colp[2],lwd=2)

polygon(c(1:10,10:1),c(meanCellCross$V2-sdCellCross$V2,rev(meanCellCross$V2+sdCellCross$V2)),col = colp[3],border = NA)
lines(c(1:10),meanCellCross$V2,col=colp[4],lwd=2)

polygon(c(1:10,10:1),c(meanCellCross$V3-sdCellCross$V3,rev(meanCellCross$V3+sdCellCross$V3)),col = colp[5],border = NA)
lines(c(1:10),meanCellCross$V3,col=colp[6],lwd=2)

polygon(c(1:10,10:1),c(meanCellCross$V4-sdCellCross$V4,rev(meanCellCross$V4+sdCellCross$V4)),col = colp[7],border = NA)
lines(c(1:10),meanCellCross$V4,col=colp[8],lwd=2)


#calculate the beads total reads count
CellReads = as.matrix(table(allRNAlocus[,1]))
CellReads = cbind(as.numeric(rownames(CellReads)),CellReads[,1])


library("RColorBrewer")
colp = brewer.pal(10,"RdYlBu")

k = c(5,4.5,4,3.5,3,2.5,2,1.5,1,0.5)
meanCellCross = matrix(rep(0,40),10,4)
sdCellCross = matrix(rep(0,40),10,4)
plot(c(0.5,10.5),c(0,5000),type="n",axes = F)
axis(side=2,at=c(0,1000,2000,3000,4000,5000),labels=c(0,1000,2000,3000,4000,5000))
axis(side=1,at=1:10,labels=c("10um","9um","8um","7um","6um","5um","4um","3um","2um","1um"))
box()
for (m in 1:10){
  nn = m*2
  CrossCell = matrix(rep(0,64),16,4)
  for(j in 1:16){
    name = paste('~/lab508/sptSimulation/simulatedBead_',k[m],'um_rep',j,".txt",sep = "")
    #beadsinfor = read.table('~/lab508/sptSimulation/simulatedBinfor_5um_rep1.txt')
    beadsRNA = read.table(name,skip = 1)
    
    beadsTotal = as.matrix(t(table(beadsRNA[,6])))
    if(j==1){
      readsTotal = beadsTotal[2:length(beadsTotal)]
    } else {
      readsTotal = c(readsTotal,beadsTotal[2:length(beadsTotal)])
    }
  }
  vioplot(readsTotal,col=colp[m],add = T,at=m)
}

##############beads ratio verus cell ratio
CellReads = as.matrix(table(allRNAlocus[,1]))
CellReads = cbind(as.numeric(rownames(CellReads)),CellReads[,1])


k = c(5,4.5,4,3.5,3,2.5,2,1.5,1,0.5)
meanCellCross = matrix(rep(0,40),10,4)
sdCellCross = matrix(rep(0,40),10,4)
for (m in 1:10){
  nn = m*2
  CrossCell = matrix(rep(0,64),16,4)
  for(j in 1:5){
    name = paste('~/lab508/sptSimulation/simulatedBead_',k[m],'um_rep',j,".txt",sep = "")
    #beadsinfor = read.table('~/lab508/sptSimulation/simulatedBinfor_5um_rep1.txt')
    beadsRNA = read.table(name,skip = 1)
    
    beadsToCell = as.matrix(t(table(beadsRNA[,c(2,6)])))
    
    beadsStat = matrix(rep(0,(nrow(beadsToCell)-1)*6),(nrow(beadsToCell)-1),6)
    #col1: beadsID, beadsTotReads, MaxCellReads, MaxCellToBeadPer, MaxCellToCellPer,TotCell
    
    for(i in 2:nrow(beadsToCell)){
      nc1 = i-1
      nc2 = sum(beadsToCell[i,])
      nc3 = max(beadsToCell[i,])
      nc4 = nc3/nc2
      mCellID = which.max(beadsToCell[i,])
      nc5 = nc3/CellReads[mCellID,2]
      nc6 = sum(beadsToCell[i,]>0)
      beadsStat[i-1,]=c(nc1,nc2,nc3,nc4,nc5,nc6)
    }
    
    if(j==1){
      ratio = beadsStat[,5]/beadsStat[,4]
    } else {
      ratio = c(ratio,beadsStat[,5]/beadsStat[,4])
    }
  }
  if(m==1){
    plot(ecdf(ratio),col=colp[m])
  } else {
    lines(ecdf(ratio),col=colp[m])
  }
}

legend("bottomright",c("10um","9um","8um","7um","6um","5um","4um","3um","2um","1um"),lwd=2,col=colp[1:10],bty="n")


###############percentage of cell recovered

CellReads = as.matrix(table(allRNAlocus[,1]))
CellReads = cbind(as.numeric(rownames(CellReads)),CellReads[,1])


k = c(5,4.5,4,3.5,3,2.5,2,1.5,1,0.5)
meanCellCross = matrix(rep(0,40),10,4)
sdCellCross = matrix(rep(0,40),10,4)
perCellCutoff = 0.66
CellRecovery = matrix(rep(0,119*10),119,10)
for (m in 1:10){
  nn = m*2
  for(j in 1:10){
    name = paste('~/lab508/sptSimulation/simulatedBead_',k[m],'um_rep',j,".txt",sep = "")
    #beadsinfor = read.table('~/lab508/sptSimulation/simulatedBinfor_5um_rep1.txt')
    beadsRNA = read.table(name,skip = 1)
    
    beadsToCell = as.matrix(t(table(beadsRNA[,c(2,6)])))
    
    beadsStat = matrix(rep(0,(nrow(beadsToCell)-1)*7),(nrow(beadsToCell)-1),7)
    #col1: beadsID, beadsTotReads, MaxCellReads, MaxCellToBeadPer, MaxCellToCellPer,TotCell, maxCellID
    
    for(i in 2:nrow(beadsToCell)){
      nc1 = i-1
      nc2 = sum(beadsToCell[i,])
      nc3 = max(beadsToCell[i,])
      nc4 = nc3/nc2
      mCellID = which.max(beadsToCell[i,])
      nc5 = nc3/CellReads[mCellID,2]
      nc6 = sum(beadsToCell[i,]>0)
      nc7 = mCellID
      beadsStat[i-1,]=c(nc1,nc2,nc3,nc4,nc5,nc6,nc7)
    }
    beadsStat=beadsStat[beadsStat[,4]>perCellCutoff,]
    detectCell = unique(beadsStat[,7])
    for(c in detectCell){
      CellRecovery[c,m] = CellRecovery[c,m]+sum(beadsStat[beadsStat[,7]==c,5])
    }
  }
}



beeswarm(list(x10um=CellReads[,2],x9um=CellReads[,2],x8um=CellReads[,2],x7um=CellReads[,2],x6um=CellReads[,2],x5um=CellReads[,2]
              ,x4um=CellReads[,2],x3um=CellReads[,2],x2um=CellReads[,2],x1um=CellReads[,2]),pwcol = hsv(range01(CellRecovery/10),range01(CellRecovery/10),range01(CellRecovery/10)),pch=19)



plot(ecdf(CellRecovery[,1]/10),col=colp[1],cex=.5,xlim=c(0,1.0))
for(i in 2:10){
  lines(ecdf(CellRecovery[,i]/10),col=colp[i],cex=.5)
}

plot(ecdf(CellRecovery[,1]/10*CellReads[,2]),col=colp[1],cex=.5)
for(i in 2:10){
  lines(ecdf(CellRecovery[,i]/10*CellReads[,2]),col=colp[i],cex=.5)
}

boxplot(CellRecovery/10,col=colp[1:10])

vioplot(CellRecovery/10,col=colp[1:10])

plot(c(0.5,10.5),c(0,9000),type="n",axes = F)
axis(side=2,at=c(0,3000,6000,9000),labels=c(0,3000,6000,9000))
axis(side=1,at=1:10,labels=c("10um","9um","8um","7um","6um","5um","4um","3um","2um","1um"))
box()

for(i in 1:10){
  vioplot(CellRecovery[,i]/10*CellReads[,2],col=colp[i],add=T,at=i)
}
