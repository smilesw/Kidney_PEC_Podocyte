library(Seurat)
library(fgsea)
library(ggplot2)
library(readxl)
#S1 is Seurat object comprising PEC and podocyte subpopulations
pathway=gmtPathways("D:/genesets_mouse.gmt")
for(i in 1:length(pathway)){
  geneset=pathway[[i]]
  geneset=intersect(geneset, rownames(S1))
  title=names(pathway)[[i]]
  S1 <- AddModuleScore(S1, features = list(geneset), name=title)
}

pathway_name=NULL
for(i in 1:length(pathway)){
  a=names(pathway)[[i]]
  a=paste(a, 1, sep = "")
  pathway_name=c(pathway_name,a)
}

VlnPlot(S1, features = pathway_name, stack=T, flip=T, fill.by="ident", 
        cols=c("#00b0f6","darkorange","darkorchid","navajowhite3", "tomato","forestgreen")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_blank(),
        legend.position = "none") + geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.5)
