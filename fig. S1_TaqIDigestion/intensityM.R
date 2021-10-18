probeTaqI = read.table('~/lab508/pixelSeq_figures/TaqIDigestion/probeTaqI.csv',header=F)
probeCon = read.table('~/lab508/pixelSeq_figures/TaqIDigestion/probeCon.csv',header=F)

#background removal, defined by the minimum intensity

min1 = min(probeTaqI$V1)
min2 = min(probeCon$V1)

probeTaqI$V1 = probeTaqI$V1-min(c(min1,min2))
probeCon$V1 = probeCon$V1-min(c(min1,min2))

plot(ecdf(probeTaqI$V1/probeCon$V1),xlim=c(0,1),lwd=2,las=1)


intensityLeftPer = probeTaqI$V1/probeCon$V1

#0.2 as cutoff
length(intensityLeftPer[intensityLeftPer<=0.2])/length(intensityLeftPer)


#boxplot
boxplot(intensityLeftPer,ylim=c(0,0.2),pch=19,las=1,axes=F,outline = F)
axis(side=2,at=seq(0,0.2,0.05),labels=seq(0,0.2,0.05),las=1)


library(beeswarm)
vioplot(intensityLeftPer,ylim=c(0,1),col='white')
