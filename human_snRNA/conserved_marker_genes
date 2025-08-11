library(nichenetr)
library(dplyr)
library(Seurat)

# mouse Seurat object to convert human orthologs
S1_mouse$ident=S1_mouse@active.ident
exp_mtx <- as.matrix(S1_mouse@assays$RNA@counts)
con_df <- data.frame(mouse_orig = rownames(exp_mtx),
                     human = convert_mouse_to_human_symbols(rownames(exp_mtx)),
                     stringsAsFactors = F)

con_df <- con_df[!is.na(con_df$human),,F]
con_df=distinct(con_df, human, .keep_all = T)
exp_mtx <- exp_mtx[con_df$mouse_orig,]

# Now chnage the rownames of the matrix to the human gene names
rownames(exp_mtx) <- con_df$human

# Create the seurat object with human genes.
S1_mouse <- CreateSeuratObject(counts = exp_mtx, meta.data = S1_mouse@meta.data )
S1_mouse[["umap"]] <- S1_mouse[["umap"]]
S1_mouse <- NormalizeData(S1_mouse)

# S1 is a human Seurat object (snRNA-seq) from GSE183277, and S2 is a human Seurat object (snRNA-seq) from GSE195460 
S1 <- merge(S1_mouse, y = c(S1, S2))

# run standard analysis workflow
S1 <- NormalizeData(S1)
S1 <- FindVariableFeatures(S1)
S1 <- ScaleData(S1)
S1 <- RunPCA(S1)
S1 <- FindNeighbors(S1, dims = 1:30, reduction = "pca")
S1 <- FindClusters(S1, resolution = 2, cluster.name = "unintegrated_clusters")
S1 <- RunUMAP(S1, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
S1 <- IntegrateLayers(S1, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)

# re-join layers after integration
S1[["RNA"]] <- JoinLayers(S1[["RNA"]])
S1 <- FindNeighbors(S1, reduction = "integrated.cca", dims = 1:30)
S1 <- FindClusters(S1, resolution = 2.0)
S1 <- RunUMAP(S1, dims = 1:30, reduction = "integrated.cca")

marker_PEC1=FindConservedMarkers(S1, ident.1 = "PEC1", ident.2 = "Podocyte1", grouping.var = "sample", test.use = "MAST", logfc.threshold = 0.5, min.diff.pct = 0.1, min.pct = 0.1)
marker_PEC2=FindConservedMarkers(S1, ident.1 = "PEC2", ident.2 = "Podocyte1", grouping.var = "sample", test.use = "MAST", logfc.threshold = 0.5, min.diff.pct = 0.1, min.pct = 0.1)
marker_Prepodocyte=FindConservedMarkers(S1, ident.1 = "Prepodocyte", ident.2 = "Podocyte1", grouping.var = "sample", test.use = "MAST", logfc.threshold = 0.5, min.diff.pct = 0.1, min.pct = 0.1)
marker_Podocyte2=FindConservedMarkers(S1, ident.1 = "Podocyte2", ident.2 = "Podocyte1", grouping.var = "sample", test.use = "MAST", logfc.threshold = 0.5, min.diff.pct = 0.1, min.pct = 0.1)
marker_Podocyte3=FindConservedMarkers(S1, ident.1 = "Podocyte3", ident.2 = "Podocyte1", grouping.var = "sample", test.use = "MAST", logfc.threshold = 0.5, min.diff.pct = 0.1, min.pct = 0.1)
