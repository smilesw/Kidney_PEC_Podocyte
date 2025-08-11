library(Signac)
library(Seurat)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(patchwork)

# extract position frequency matrices for the motifs
pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# add motif information
S1 <- AddMotifs(S1, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm)

# gather the footprinting information for sets of motifs
S1@assays$peaks@fragments[[1]]@path="D:/public_datasets/GSE195460_snATAC/cellranger_atac_aggr/outs/fragments.tsv.gz"

S1 <- Footprint(
  object = S1,
  motif.name = c("HNF4A"),
  genome = BSgenome.Hsapiens.UCSC.hg38
)

PlotFootprint(S1, features = c("HNF4A"), label = F) & scale_colour_manual(values = c("#00b0f6","darkorange","darkorchid","navajowhite3","forestgreen")) 
