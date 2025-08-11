library(Signac) 
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork) 
library(cicero)
library(here) 
DefaultAssay(S1)="peaks"

# convert to cds monocle/cicero format
atac.cds <- as.cell_data_set(S1)
atac.cicero <- make_cicero_cds(atac.cds, reduced_coordinates=reducedDims(atac.cds)$UMAP)

# remove alt contigs
genomes <- seqlengths(S1)[1:25]

# convert chr sizes to df
genome.df <- data.frame("chr" = names(genomes), "length" = genomes)

# run cicero
conns <- run_cicero(atac.cicero, genomic_coords = genome.df, sample_num=100)

# generate ccan
ccans <- generate_ccans(conns)

# add to seurat obj
links <- ConnectionsToLinks(conns=conns, ccans=ccans)
Links(S1) <- links
