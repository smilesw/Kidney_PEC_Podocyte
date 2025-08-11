DefaultAssay(kidney_sub) <- 'peaks'
kidney_sub@assays$peaks@fragments[[1]]@path="D:/public_datasets/GSE195460_snATAC/cellranger_atac_aggr/outs/fragments.tsv.gz"
S1=subset(kidney_sub, idents=c("PEC","PODO"))

S1 <- RunTFIDF(S1)
S1 <- FindTopFeatures(S1, min.cutoff = 'q0')
S1 <- RunSVD(S1)
DepthCor(S1)
S1 <- RunUMAP(object = S1, reduction = 'lsi', dims = 2:30)
S1 <- FindNeighbors(object = S1, reduction = 'lsi', dims = 2:30)
S1 <- FindClusters(object = S1, verbose = FALSE, algorithm = 1, resolution = 1.0)
