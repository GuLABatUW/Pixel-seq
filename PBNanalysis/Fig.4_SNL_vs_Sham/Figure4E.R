library(viridis)
library(ggplot2)
goenrich = read.table('~/lab508/novaSeqPBN/GOBP/version3/enrichedBYenrichr.txt',header=F,sep="\t")
colnames(goenrich) = c('db','term','cluster','frac','count','pval')
goenrich$p.adjust = -log10(goenrich$pval)

goBP = as.data.frame(goenrich[goenrich$db=='BP',])

ggplot(goBP, aes(cluster,term, size = frac, colour = p.adjust,width=2.5)) +
  geom_point(shape = 19) +
    theme_linedraw(base_rect_size = 1.2) +scale_colour_gradientn(colours = rev(plasma(256))[80:256])+
  guides(x =  guide_axis(angle = 90))
