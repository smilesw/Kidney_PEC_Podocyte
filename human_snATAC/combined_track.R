library(bedr)
library(Signac)
library(Seurat)
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
library(openxlsx)
library(GenomicRanges)
library(readxl)

DefaultAssay(S1) <- "peaks"
# add ChIP-seq peak to Seurat object 
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
seqlevelsStyle(annotations) <- "UCSC"

chip_peak <- bed_to_granges("D:/public_datasets/ENCFF286UYA_HNF4A.bed")
chip_peak <- keepStanchip_peakdChromosomes(chip_peak, pruning.mode = "coarse")
chip_peak <- subsetByOverlaps(x = chip_peak, ranges = blacklist_hg38_unified, invert = TRUE)
chip_peak <- makeGRangesFromDataFrame(chip_peak)
peakwidths <- width(chip_peak)
chip_peak <- chip_peak[peakwidths  < 10000 & peakwidths > 20]

counts <- FeatureMatrix(
  fragments = Fragments(S1),
  features = chip_peak,
  cells = colnames(S1)
)
S1[["HNF4A_ChIP"]] <- CreateChromatinAssay(
  counts = counts,
  fragments = Fragments(S1),
  annotation = annotations
)

# add DAR to Seurat object 
ident <- getSheetNames("D:/public_datasets/GSE195460_snATAC/DAR.xlsx")
df <- read.xlsx("D:/public_datasets/GSE195460_snATAC/DAR.xlsx", sheet = ident[5])
peaks_df <- data.frame(peak = df$query_region)
peak_coords <- do.call(rbind, strsplit(peaks_df$peak, "-"))
colnames(peak_coords) <- c("chr", "start", "end")
peak_coords <- data.frame(chr = peak_coords[, 1],
                          start = as.numeric(peak_coords[, 2]),
                          end = as.numeric(peak_coords[, 3]))
peaks_gr <- GRanges(seqnames = peak_coords$chr,
                    ranges = IRanges(start = peak_coords$start,
                                     end = peak_coords$end))
peaks_gr <- keepStandardChromosomes(peaks_gr, pruning.mode = "coarse")
peaks_gr <- subsetByOverlaps(x = peaks_gr, ranges = blacklist_hg38_unified, invert = TRUE)
peaks_gr <- makeGRangesFromDataFrame(peaks_gr)
peakwidths <- width(peaks_gr)
peaks_gr <- peaks_gr[peakwidths  < 10000 & peakwidths > 20]

counts <- FeatureMatrix(
  fragments = Fragments(S1),
  features = peaks_gr,
  cells = colnames(S1)
)
S1[["DAR_podocyte3"]] <- CreateChromatinAssay(
  counts = counts,
  fragments = Fragments(S1),
  annotation = annotations
)

# coverage plot
region <- "chr1-155299000-155301462" #PKLR
region <- "chr9-101421388-101436287" #ALDOB
region <- "chr7-73591201-73626522" #MLXIPL

cov_plot <- CoveragePlot(object = S1, region = region, annotation = TRUE, peaks = TRUE, links = TRUE) & scale_fill_manual(values=c("#00b0f6","darkorange","darkorchid","navajowhite3","forestgreen"))
peak_plot2 <- PeakPlot(object = S1, assay = "DAR_DP", color = "red", region = region) + ylab("DAR") + theme(axis.title.y = element_text(size=10))
peak_plot3 <- PeakPlot(object = S1, assay = "HNF4A", color = "black", region = region) + ylab("HNF4A") + theme(axis.title.y = element_text(size=10))
coverage_plot <- BigwigTrack(region = region, bigwig = "/mnt/d/public_datasets/ENCFF531OYO.bigWig", smooth = 300, bigwig.scale = "separate")+
  scale_fill_manual(values = "black")+ theme(legend.position="none")+ theme(axis.title.y = element_text(size=0)) +theme(axis.title.y = element_blank())

CombineTracks(
  plotlist = list(cov_plot,
                  peak_plot2,
                  coverage_plot,
                  peak_plot3),
  heights = c(2,
              0.1,
              0.2,
              0.1)
) 
