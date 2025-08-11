library(Seurat)
library(Signac) 
library(EnsDb.Hsapiens.v86)
library(openxlsx)
library(here)
library(pheatmap)
library(BuenColors)
library(dplyr)
library(tibble)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(viridis)

DefaultAssay(S1) <- 'peaks'
S1@assays$peaks@fragments[[1]]@path="D:/public_datasets/GSE195460_snATAC/cellranger_atac_aggr/outs/fragments.tsv.gz"

# wrapper function for FindMarkers 
GetMarkers <- function(cluster, seurat) {
  print(paste0("Finding DAR for: ",cluster))
  dar <- FindMarkers(seurat, 
                     ident.1 = cluster,    
                     test.use = 'LR', 
                     latent.vars = "peak_region_fragments",
                     only.pos = TRUE,
                     min.pct = 0.1) # find all cluster-specific dars
  dar = dar %>% rownames_to_column("query_region")
  dar$fdr = p.adjust(dar$p_val, method='fdr')
  dar = dar %>% dplyr::filter(fdr < 0.05)
  dar=dplyr::mutate(dar, cellstate = cluster)
  cf <- ClosestFeature(seurat, regions=dar$query_region)
  return(merge(dar, cf, by="query_region"))
}
list.cluster.dar <- lapply(idents, function(x) {GetMarkers(x, seurat = S1)})
write.xlsx(list.cluster.dar, file = "D:/public_datasets/GSE195460_snATAC/DAR.xlsx", sheetName = idents, rowNames = T)

# identify all unique cell-type-specific peaks and filter for logfc > 0
all_dar <- bind_rows(list.cluster.dar) %>%
  dplyr::filter(avg_log2FC > 0) %>%
  dplyr::select("query_region") %>%
  dplyr::distinct()

dar_aver <- AverageExpression(S1, features = all_dar$coord, assays = "peaks")
pheatmap(dar_aver[["peaks"]], scale = 'row', cluster_rows = T, cluster_cols = F, show_rownames = FALSE,
         border_color = "white", legend = T,
         color = viridis(n = 256, alpha = 1, begin = 0, end = 1, option = "viridis"))
    
PEC1=list.cluster.dar[[1]] %>% select("query_region")
PEC1=data.frame(do.call('rbind', strsplit(as.character(PEC1$query_region), split='-', fixed=T))) 
write.table(APEMP, file = "./DAR_PEC1.txt", sep = "\t")

PEC2=list.cluster.dar[[2]] %>% select("query_region")
PEC2=data.frame(do.call('rbind', strsplit(as.character(PEC2$query_region), split='-', fixed=T))) 
write.table(TC, file = "./DAR_PEC2.txt", sep = "\t")

Prepodocyte=list.cluster.dar[[3]] %>% select("query_region")
Prepodocyte=data.frame(do.call('rbind', strsplit(as.character(Prepodocyte$query_region), split='-', fixed=T))) 
write.table(PCP, file = "./DAR_Prepodocyte.txt", sep = "\t")

Podocyte1=list.cluster.dar[[4]] %>% select("query_region")
Podocyte1=data.frame(do.call('rbind', strsplit(as.character(Podocyte1$query_region), split='-', fixed=T))) 
write.table(HP, file = "./DAR_Podocyte1.txt", sep = "\t")

Podocyte3=list.cluster.dar[[5]] %>% select("query_region")
Podocyte3=data.frame(do.call('rbind', strsplit(as.character(Podocyte3$query_region), split='-', fixed=T))) 
write.table(DP, file = "./DAR_Podocyte3.txt", sep = "\t")
