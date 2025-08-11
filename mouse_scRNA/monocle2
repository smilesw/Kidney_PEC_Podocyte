library(dplyr)
library(Seurat)
library(patchwork)
library(DelayedArray)
library(monocle)
S1$ident=S1@active.ident
S1@meta.data=S1@meta.data[,1:4]
S1@meta.data=S1@meta.data[,-3]
x=GetAssayData(object = S1, slot = "counts")[rownames(GetAssayData(object = S1, slot = "counts")) %in% rownames(GetAssayData(object = S1)), colnames(GetAssayData(object = S1, slot = "counts")) %in% colnames(GetAssayData(object = S1))]
target=x[,colnames(x) %in% rownames(S1@meta.data)]
target=as.matrix(t(target))
name=merge(S1@meta.data, target, by="row.names")
rownames(name)=name$Row.names
name=name[,-1]
S1@meta.data=name[,1:4]
name=name[,-1:-3]
target=t(name)
g=as.data.frame(rownames(x))
rownames(g)=rownames(x)
colnames(g)=c("gene_short_name")

pd <- new("AnnotatedDataFrame", data = S1@meta.data)
fd <- new("AnnotatedDataFrame", data = g)
HSMM <- newCellDataSet(as.matrix(target), phenoData = pd, featureData = fd, expressionFamily=negbinomial.size())
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))

disp_table <- dispersionTable(HSMM)
ordering_genes <- subset(disp_table, mean_expression >= 0.05 & dispersion_empirical >= 2 * dispersion_fit)$gene_id
HSMM <- setOrderingFilter(HSMM, ordering_genes)

plot_ordering_genes(HSMM)
HSMM <- reduceDimension(HSMM, max_components=2, reduction_method = c("DDRTree"))
HSMM <- orderCells(HSMM, reverse=FALSE)
HSMM$ident=factor(HSMM$ident, levels = c("APEMP","TC","PCP","HP"))

#trajectory by subpopulation, group, and pseudotime
plot_cell_trajectory(HSMM, show_tree = FALSE ,show_backbone = FALSE, color_by="ident", show_branch_points = FALSE, cell_size=1.5) +
  scale_color_manual(values=c("#00b0f6","darkorange","darkorchid","navajowhite3","tomato","forestgreen")) +
  xlab("Component 1") + ylab("Component 2") +
  theme(text = element_text(size=30),
        axis.text.x = element_text(angle=0, size=20),
        axis.text.y = element_text(size=20),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text = element_text(size = 20, face = "italic"),
        legend.position = "right",
        legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 8)))

plot_cell_trajectory(HSMM, show_tree = FALSE ,show_backbone = FALSE, color_by="group", show_branch_points = FALSE, cell_size=1.5) +
  scale_color_manual(values=c("#1E90FF","#FF1493")) +
  xlab("Component 1") + ylab("Component 2") +
  theme(text = element_text(size=30),
        axis.text.x = element_text(angle=0, size=20),
        axis.text.y = element_text(size=20),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text = element_text(size = 20, face = "italic"),
        legend.position = "none",
        legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 8)))

plot_cell_trajectory(HSMM, show_tree = FALSE ,show_backbone = FALSE, color_by="Pseudotime", show_branch_points = FALSE, cell_size=1.5) +
  scale_color_gradient (low="navy", high="lightgoldenrod1") +
  xlab("Component 1") + ylab("Component 2") +
  theme(text = element_text(size=30),
        axis.text.x = element_text(angle=0, size=20),
        axis.text.y = element_text(size=20),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text = element_text(size = 20, face = "italic"),
        legend.position = "none",
        legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 8)))

#heatmap
genelist=c("Cers6","Sgms1","Ldhb","Myo1e","Tyro3","Loxl2","Bmp7","Miox","Axl","St6galnac3","Plod2","Actn4","Gpx3","Tspan2","Ucp2","Gstm1","Sod3","Itga3","Nr4a1","Mgat5","Tnfrsf12a","Tnfaip2","Nebl","Cyp1b1","Scin","Cldn1","Vcam1","Mafb","Dach1","Podxl","Aldob","Bcam","Dag1","Pros1")
pal <- colorRampPalette(c("deepskyblue","white", "red"))
plot_pseudotime_heatmap(HSMM_genuine[genelist,], num_clusters = 4, cores =1, cluster_rows = TRUE, show_rownames = T, hmcols=pal(1000)) 
