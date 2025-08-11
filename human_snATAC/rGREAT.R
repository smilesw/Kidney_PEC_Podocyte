library(openxlsx)
library(GenomicRanges)
library(readxl)
library(rGREAT)

ident <- getSheetNames("D:/public_datasets/GSE195460_snATAC/DAR.xlsx")
output_dir <- "D:/public_datasets/GSE195460_snATAC/rGREAT"

for (cell_state in ident) {
  message("Processing sheet: ", cell_state)
  df <- read.xlsx("D:/public_datasets/GSE195460_snATAC/DAR.xlsx", sheet = cell_state)

  peaks_df <- data.frame(peak = df$query_region)
  peak_coords <- do.call(rbind, strsplit(peaks_df$peak, "-"))
  colnames(peak_coords) <- c("chr", "start", "end")
  peak_coords <- data.frame(chr = peak_coords[, 1],
                            start = as.numeric(peak_coords[, 2]),
                            end = as.numeric(peak_coords[, 3]))
  peaks_gr <- GRanges(seqnames = peak_coords$chr,
                      ranges = IRanges(start = peak_coords$start,
                                       end = peak_coords$end))
  job <- submitGreatJob(peaks_gr, species = "hg38")
  tb <- getEnrichmentTables(job)
  sheetname=c("GO_MF","GO_BP","GO_CC")
  
  out_file <- file.path(output_dir, paste0(cell_state, "_GREAT.xlsx"))
  write.xlsx(tb, file = out_file, sheetName = sheetname, rowNames = F)
}
