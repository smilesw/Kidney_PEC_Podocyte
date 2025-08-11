library(data.table)
library(SCENIC)
library(AUCell)
library(RcisTarget)
library(Seurat)
library(NMF)
library(synchronicity)
library(dplyr)

exprMat <- as.matrix(GetAssayData(S1, slot = "data"))
ident = as.matrix(S1@active.ident)
cellInfo<-rbind(S1$nFeature_RNA, S1$nCount_RNA)
cellInfo<-t(as.matrix(cellInfo))
cellInfo<-cbind(ident,cellInfo)
colnames(cellInfo) <- c("CellType","nGene","nUMI")
cellInfo<-as.data.frame(cellInfo)
cellInfo$CellType=factor(cellInfo$CellType, levels=c("PEC1", "PEC2", "Prepodocyte", "Podocyte1", "Podocyte2", "Podocyte3"))
saveRDS(cellInfo, file="int/cellInfo.Rds")

# Color to assign to the variables (same format as for NMF::aheatmap)
colVars <- list(CellType=c("PEC1"="#00b0f6",
                           "PEC2"="darkorange",
                           "Prepodocyte"="darkorchid",
                           "Podocyte1"="navajowhite3",
                           "Podocyte2"="tomato", 
                           "Podocyte3"="forestgreen"))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
saveRDS(colVars, file="int/colVars.Rds")
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

org="mgi"
dbDir="/mnt/d/home/smilesw/program/cisTarget_databases/"
myDatasetTitle="SCENIC"
dbs <- defaultDbNames[[org]]
dbs[1] <-c("mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
dbs[2] <-c("mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10) 
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
interestingGenes <- c("Sox9", "Sox10", "Dlx5", "Pou3f3")
interestingGenes[which(!interestingGenes %in% genesKept)]
exprMat_filtered <- exprMat[genesKept, ]
ccc <- as.data.frame(exprMat_filtered)
ccc <- ccc[which(! row.names(ccc) %in% c("Hoxa11os","Ptges3l","Wbscr22","Ren1","Iqck","Gm12359","4430402I18Rik","9330151L19Rik", "Atr", "Pi4ka", "Sco2", "Med31", "Alkbh1", "Comtd1", "Scand1", "Ndufa13", "Dhfr","Pex1")),]
exprMat_filtered <- as.matrix(ccc)
rm(ccc)

set.seed(123)
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
exprMat_log <- log2(exprMat+1)
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat_filtered_log)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

#binary heatmap
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
data=read.csv("./regulon_list.csv", row.names = 1)
binary_heatmap=binaryRegulonActivity[rownames(data),]

state=as.data.frame(S1@active.ident)
group=as.data.frame(S1$group)

cellInfo<-cbind(state, group)
colnames(cellInfo) <- c("state","group")
cellInfo<-as.data.frame(cellInfo)

cellInfo$state=factor(cellInfo$state, levels = c("PEC1", "PEC2", "Prepodocyte", "Podocyte1", "Podocyte2", "Podocyte3"))
cellInfo$group=factor(cellInfo$group, levels = c("CTL","DKD"))

#cellInfo=cellInfo %>% select(CellType)
colVars <- list(state=c("PEC1"="#00b0f6", 
                        "PEC2"="darkorange",
                        "Prepodocyte"="darkorchid",
                        "Podocyte1"="navajowhite3",
                        "Podocyte2"="tomato",
                        "Podocyte3"="forestgreen"),
                group=c("CTL"="#1E90FF",
                        "DKD"="#FF1493"))

state1=cellInfo %>% filter(state == "PEC1" & group == "CTL")
state2=cellInfo %>% filter(state == "PEC1" & group == "DKD")
state3=cellInfo %>% filter(state == "PEC2" & group == "CTL")
state4=cellInfo %>% filter(state == "PEC2" & group == "DKD")
state5=cellInfo %>% filter(state == "Prepodocyte" & group == "CTL")
state6=cellInfo %>% filter(state == "Prepodocyte" & group == "DKD")
state7=cellInfo %>% filter(state == "Podocyte1" & group == "CTL")
state8=cellInfo %>% filter(state == "Podocyte1" & group == "DKD")
state9=cellInfo %>% filter(state == "Podocyte2" & group == "CTL")
state10=cellInfo %>% filter(state == "Podocyte2" & group == "DKD")
state11=cellInfo %>% filter(state == "Podocyte3" & group == "CTL")
state12=cellInfo %>% filter(state == "Podocyte3" & group == "DKD")
state_new=bind_rows(state1,state2,state3,state4,state5,state6,state7,state8,state9,state10,state11,state12)

m <- match(rownames(state_new), colnames(binary_heatmap))
binary_heatmap2=binary_heatmap[,m]
dev.off()
aheatmap(binary_heatmap2, scale="none", revC=TRUE, 
         annCol=state_new[,, drop=FALSE],
         annColor=colVars,
         labCol = NA,
         Rowv=NA,
         Colv=NA,
         color = c("white", "black"),
         legend=FALSE)
