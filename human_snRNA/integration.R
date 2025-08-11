library(Seurat)
library(SeuratWrappers)
library(SeuratData)

# KPMP data processing, dowloaded from https://cellxgene.cziscience.com/datasets (GSE183277)
kidney1@active.ident=kidney1$cell_type
rownames(kidney1@assays$RNA@counts)=kidney1@assays$RNA@meta.features$feature_name
rownames(kidney1@assays$RNA@data)=kidney1@assays$RNA@meta.features$feature_name
rownames(kidney1@assays$RNA@meta.features)=kidney1@assays$RNA@meta.features$feature_name

celltype=as.data.frame(names(table(kidney1@active.ident)))
colnames(celltype)="ident"
celltype=mutate(celltype, interest=ifelse(str_detect(celltype$ident, 'Podocyte')==TRUE, 'yes', 'no'))
a=celltype %>% filter(interest == "yes") 
celltype=mutate(celltype, interest=ifelse(str_detect(celltype$ident, 'Parietal Epithelial Cell')==TRUE, 'yes', 'no'))
b=celltype %>% filter(interest == "yes") 

S1=subset(kidney1, idents=c(a$ident, b$ident))
S1=subset(S1, cells=names(S1@active.ident[S1$condition.long %in% c("Normal Reference","Diabetic Kidney Disease (DKD)")]))

S1=SetIdent(S1, value = S1[["suspension_type"]])
S1=subset(S1, idents="nucleus")

S1[["RNA"]]@meta.features=data.frame(row.names= rownames(S1[["RNA"]]))
S1 = FindVariableFeatures(S1, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
S1 = ScaleData(S1, features=rownames(S1))
S1 = RunPCA(S1, npcs=30, verbose=T)
S1 = RunUMAP(S1,  dims=1:20)
S1 = FindNeighbors(S1,  dims=1:20)
S1 = FindClusters(S1, resolution = 0.5)

S1=UpdateSeuratObject(S1)
S1$sample="human"
S1@project.name="KPMP"
S1[["RNA"]] <- as(object = S1[["RNA"]], Class = "Assay5")
saveRDS(S1, file="./public_data/KPMP/Seurat_cell/Podo_PEC_normal_DKD.Rds")

# PNAS & Nature com data processing, dowloaded from https://cellxgene.cziscience.com/datasets (GSE195460)
kidney2@active.ident=kidney2$cell_type
rownames(kidney2@assays$RNA@counts)=kidney2@assays$RNA@meta.features$feature_name
rownames(kidney2@assays$RNA@data)=kidney2@assays$RNA@meta.features$feature_name
rownames(kidney2@assays$RNA@meta.features)=kidney2@assays$RNA@meta.features$feature_name

S2 = subset(kidney2, idents=c("parietal epithelial cell","podocyte"))
S2 = FindVariableFeatures(S2, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
S2 = ScaleData(S2, features=rownames(S2))
S2 = RunPCA(S2, npcs=30, verbose=T)
S2 = RunUMAP(S2,  dims=1:20)
S2 = FindNeighbors(S2,  dims=1:20)
S2 = FindClusters(S2, resolution = 1.0)

S2_final2=UpdateSeuratObject(S2_final2)
S2_final2$sample="human2"
S2_final2@project.name="Nature_com"
S2$group="dummy"
S2$group[which(S2$disease == "normal")] <- "CTL" 
S2$group[which(S2$disease == "type 2 diabetes mellitus")] <- "DKD" 
S2[["RNA"]] <- as(object = S2[["RNA"]], Class = "Assay5")

#integration
S1 <- merge(S1, y = S2)

# run standard anlaysis workflow
S1 <- NormalizeData(S1)
S1 <- FindVariableFeatures(S1)
S1 <- ScaleData(S1)
S1 <- RunPCA(S1)
S1 <- FindNeighbors(S1, dims = 1:30, reduction = "pca")
S1 <- FindClusters(S1, resolution = 2, cluster.name = "unintegrated_clusters")
S1 <- RunUMAP(S1, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(S1, reduction = "umap.unintegrated", group.by = c("sample", "seurat_clusters"))
S1 <- IntegrateLayers(S1, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)

# re-join layers after integration
S1[["RNA"]] <- JoinLayers(S1[["RNA"]])
S1 <- FindNeighbors(S1, reduction = "integrated.cca", dims = 1:30)
S1 <- FindClusters(S1, resolution = 2.0)
S1 <- RunUMAP(S1, dims = 1:30, reduction = "integrated.cca")
marker <- FindAllMarkers(S1, min.diff.pct = 0.3 , only.pos = FALSE)
saveRDS(S1, file="D:/Seurat_object/integration_larger_sample/integration_twohuman_mouse_cell.Rds")
