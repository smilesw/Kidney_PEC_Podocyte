library(Signac)
library(Seurat) 
library(JASPAR2020) 
library(TFBSTools) 
library(BSgenome.Hsapiens.UCSC.hg38) 
library(patchwork) 
library(motifmatchr) 
library(here)
library(chromVAR) 
library(future)
library(openxlsx) 
library(ggseqlogo)
library(tibble)
library(dplyr)
library(ggplot2)

DefaultAssay(S1) <- "peaks"
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

region=rownames(S1)
GL <- grep("^GL", rownames(S1), value = T)
KI <- grep("^KI", rownames(S1), value = T)
S1=subset(S1, features=region[!region %in% c(GL,KI)])

S1 <- AddMotifs(
  object = S1,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

S1 <- RegionStats(
  object = S1,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  sep = c("-", "-")
)

S1 <- RunChromVAR(
  object = S1,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

DefaultAssay(S1)="chromvar"
df <- FindMarkers(S1,
  ident.1 = 'PEC1',
  ident.2 = 'PEC2',
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

df$log10Padj=-log10(df$p_val_adj)
df <- df[order(df$log10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
df = df %>% rownames_to_column("motif")

ggplot(df, aes(rank, log10Padj, color = log10Padj)) + 
  geom_point(size = 2) +
  theme_classic() +
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") + 
  scale_color_gradientn(colors = rev(topo.colors(10)))
