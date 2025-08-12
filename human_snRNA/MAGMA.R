library(MAGMA.Celltyping) 
library(dplyr)
library(data.table)
library(Seurat)
Sys.setenv(GITHUB_PAT = "token")

# GWAS raw data for eGFR
df <- fread("D:/bed/sumstats/eGFR_gwas.txt")  
colnames(df)=c("MarkerName","chr","position","Allele1","Allele2","Effect","P-value","n")
df_magma <- df[, .(SNP = MarkerName,
                   CHR = chr,
                   BP = position,
                   A1 = Allele1,
                   A2 = Allele2,
                   BETA = Effect,
                   P = `P-value`,
                   N = n)]

# GWAS raw data for microalbuminuria
df2 <- fread("D:/bed/sumstats/UACR_gwas_MAGMA.txt") 
df2=df2[,c(3,1,2,4,5,6,7,8)]
df_magma2 <- df2[, .(SNP = RSID,
                   CHR = Chr,
                   BP = Pos_b37,
                   A1 = Allele1,
                   A2 = Allele2,
                   BETA = Effect,
                   P = `P-value`,
                   N = n)]
df_magma2=na.omit(df_magma2)

# save
fwrite(df_magma, file = "D:/bed/sumstats/eGFR_MAGMA_input.tsv", sep = "\t")
fwrite(df_magma2, file = "D:/bed/sumstats/UACR_MAGMA_input.tsv", sep = "\t")

genesOutPath <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "D:/bed/sumstats/eGFR_MAGMA_input.tsv",
  genome_build = "GRCh37",
  upstream_kb = 10,
  downstream_kb = 1.5)
genesOutPath2 <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = "D:/bed/sumstats/UACR_MAGMA_input.tsv",
  genome_build = "GRCh37",
  upstream_kb = 10,
  downstream_kb = 1.5)

# make ctd
S1$ident=S1@active.ident
celltype_data <- EWCE::drop_uninformative_genes(
  exp = S1@assays$RNA$data, 
  input_species = "human",
  output_species = "human",
  level2annot = as.character(S1$ident)) 
annotLevels <- list(level1class=S1$ident)
ctd <- EWCE::generate_celltype_data(
  exp = celltype_data,
  annotLevels = annotLevels,
  groupName = "human_kidney")
ctd <- EWCE::load_rdata(ctd)
ctd <- MAGMA.Celltyping::prepare_quantile_groups(ctd = ctd, input_species = "human", output_species = "human")

# full_eGFR
ctAssocsTop <- MAGMA.Celltyping::calculate_celltype_associations(
  ctd = ctd,
  gwas_sumstats_path = "D:/bed/sumstats/eGFR_MAGMA_input.tsv",
  magma_dir = genesOutPath,
  EnrichmentMode = "Top 10%",
  upstream_kb = 10,
  downstream_kb = 1.5,
  force_new = TRUE)

FigsTopDecile <- MAGMA.Celltyping::plot_celltype_associations(
  ctAssocs = ctAssocsTop,
  ctd = ctd)

# full_UACR
ctAssocsTop <- MAGMA.Celltyping::calculate_celltype_associations(
  ctd = ctd,
  gwas_sumstats_path = "D:/bed/sumstats/UACR_MAGMA_input.tsv",
  magma_dir = genesOutPath2,
  EnrichmentMode = "Top 10%",
  upstream_kb = 10,
  downstream_kb = 1.5,
  force_new = TRUE)

FigsTopDecile <- MAGMA.Celltyping::plot_celltype_associations(
  ctAssocs = ctAssocsTop,
  ctd = ctd)
