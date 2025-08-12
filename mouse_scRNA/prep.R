library(Seurat)
library(Matrix)
library(SeuratWrappers)
library(SeuratData)
library(data.table)
library(here)
sample_list <- list.files()
sampleNames = mapply(function(x) x[length(x)], strsplit(sample_list, split = '-'))
meta.info = data.table(file_Path = sample_list)
meta.info[, sampleID:=c("db/m_glom-1","db/m_glom-2","db/m_nonglom-1","db/m_nonglom-2","db/db_glom-1","db/db_glom-2","db/db_nonglom-1","db/db_nonglom-2")]
meta.info[, group:=c(rep("db/m",4),rep("db/db",4))]
meta.info[, sample:=mapply(function(x) x[1], strsplit(sampleNames, '_'))]
DT::datatable(meta.info)

norm_method = 'standard'
treat = FALSE
show = FALSE
batch_list = list()
for(i in 1:nrow(meta.info)){
  outs = meta.info$file_Path[i]
  sampleID = meta.info$sampleID[i]
  group = meta.info$group[i]
  sample = meta.info$sample[i]
  print(paste("Starting processing", sampleID, 'at', Sys.time()))
  mat <- Read10X_h5(here(outs,"raw_feature_bc_matrix_cellbender_filtered.h5"))
  
  kidney <- CreateSeuratObject(counts = mat, project = sampleID)
  
  kidney[["percent.mt"]] <- PercentageFeatureSet(kidney, pattern = "^mt-")
  kidney <- subset(kidney, percent.mt < 50 & nFeature_RNA <= 3000 & nFeature_RNA >= 300)
  
  if(norm_method == 'SCT'){
    kidney <- SCTransform(kidney, vars.to.regress = 'percent.mt')
  } else{
    kidney <- NormalizeData(kidney, verbose = F)
    kidney <- FindVariableFeatures(kidney, selection.method = "vst", 
                                   nfeatures = 2000, verbose = FALSE)
  }
  
  kidney <- RenameCells(object = kidney, add.cell.id = sampleID)  # add sample name as prefix
  if(treat == TRUE){
    kidney <- RunPCA(kidney, verbose = T)
    kidney <- RunUMAP(kidney, dims = 1:30, verbose = FALSE)
    kidney <- RunTSNE(kidney, dims = 1:30, verbose = FALSE)
    
    kidney <- FindNeighbors(kidney, dims = 1:30, verbose = FALSE)
    kidney <- FindClusters(kidney, verbose = FALSE)
    
    umap = DimPlot(kidney, label = TRUE) + NoLegend()
    tsne = TSNEPlot(kidney, label = TRUE) + NoLegend()
    
    if(show == TRUE){
      #pdf(paste0('~/Project/sc_ATC/figure/QC/',sampleID,'_Seurat_QC.pdf'),12,7.5)
      print(qc)
      plot_grid(plot1,plot2)
      dev.off()
      #pdf(paste0('~/Project/sc_ATC/figure/Clustering/',sampleID,'_Seurat_Dimension_reduction_and_Clustering.pdf'))
      print(umap)
      print(tsne)
    }
  }
  
  #dev.off()
  ### add information
  kidney@meta.data[,'group'] = group
  kidney@meta.data[,'sampleID'] = sampleID
  kidney@meta.data[,'sample'] = sample
  ### add batch data into a list to merge
  #assign(paste0(sampleID,'.Seurat'),kidney)
  batch_list <- append(batch_list, kidney)
  print(paste("Finishing processing", sampleID, 'at', Sys.time()))
}

kidney <- Reduce(function(x,y){merge(x,y,merge.data = TRUE)},batch_list)
kidney <- FindVariableFeatures(kidney, selection.method = "vst", 
                               nfeatures = 2000, verbose = FALSE)
all.genes=rownames(kidney)
kidney <- ScaleData(kidney, features=all.genes)
VlnPlot(kidney, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
kidney <- RunPCA(kidney, npcs=30, verbose=F)
kidney <- JackStraw(kidney, num.replicate = 100)
kidney <- ScoreJackStraw(kidney, dims = 1:20)
JackStrawPlot(kidney, dims = 1:20)
ElbowPlot(object = kidney)
num_PC = 30
kidney <- RunUMAP(kidney, reduction = 'pca', dims=1:num_PC)
kidney <- FindNeighbors(kidney, reduction = 'pca', dims=1:num_PC)
kidney <- FindClusters(kidney, resolution = 1)

#batch correction
kidney <- RunFastMNN(object.list = SplitObject(kidney, split.by = "SampleID"))
kidney <- RunUMAP(kidney, reduction = "mnn", dims = 1:30)
kidney <- FindNeighbors(kidney, reduction = "mnn", dims = 1:30)
kidney <- FindClusters(kidney, resolution = 1)
for(j in 0:31){
  
  n=table(kidney@active.ident)
  
  ident=j
  
  n[names(n) %in% c(0:ident)]=0
  
  n=n[-1]
  
  o=n[n>0]
  
  ident2=names(o[1])
  
  bim <- FindMarkers(kidney, ident.1 = ident,  ident.2 =ident2, only.pos = F, test.use = "bimod")
  
  deg=bim[bim$p_val_adj < 0.01,]
  
  deg1=deg[abs(deg$avg_logFC) >= 1,]
  
  degs=dim(deg)
  
  degs1=dim(deg1)
  
  
  for(i in c(names(o[-1]))){
    
    bim <- FindMarkers(kidney, ident.1 = ident,  ident.2 = i, only.pos = F, test.use = "bimod")
    
    deg=bim[bim$p_val_adj < 0.01,]
    
    deg1=deg[abs(deg$avg_logFC) >= 1,]
    
    degs=cbind(degs,dim(deg))
    
    degs1=cbind(degs1,dim(deg1))
    
  }
  
  print(degs1)
  
}
