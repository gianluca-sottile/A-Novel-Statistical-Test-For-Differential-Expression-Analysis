if(!requireNamespace("data.table", quietly = TRUE)){
  install.packages("data.table")
}
library(data.table)
if(!requireNamespace("preprocessCore", quietly = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("preprocessCore")
}
library(preprocessCore)
if(!requireNamespace("limma", quietly = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("limma")
}
library(limma)
if(!requireNamespace("siggenes", quietly = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("siggenes")
}
library("siggenes")
if(!requireNamespace("pbmcapply", quietly = TRUE)){
  install.packages("pbmcapply")
}
library(pbmcapply)
if(!requireNamespace("VennDiagram", quietly = TRUE)){
  install.packages("VennDiagram")
}
library("VennDiagram")
if(!requireNamespace("topGO", quietly = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("topGO")
}
library("topGO")
if(!requireNamespace("org.Hs.eg.db", quietly = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("org.Hs.eg.db")
}
library("org.Hs.eg.db")
if(!requireNamespace("GO.db", quietly = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("GO.db")
}
library("GO.db")
if(!requireNamespace("RISmed", quietly = TRUE)){
  install.packages("RISmed")
}
library("RISmed")

source("0.hy.test.R")

TISSUE <- "brca" # or kirc
if(TISSUE == "kirc") {
  files <- c("TCGA_KIRC67/geneExpression_inNormalSamples.txt",
             "TCGA_KIRC67/geneExpression_inTumorSamples.txt")  
} else {
  files <- c("TCGA_BRCA75/geneExpression_inNormalSamples.txt",
             "TCGA_BRCA75/geneExpression_inTumorSamples.txt")  
}

##### loading data #####
dati.geni <- fread(files[1], header=T)
dati.geni <- dati.geni[, 3:dim(dati.geni)[2]]
tmp <- fread(files[[2]], header=T, sep="\t")
dati.geni <- cbind(dati.geni, tmp[, 4:dim(tmp)[2]]) 
rm(tmp)

wdir <- paste0("Results and outputs/", TISSUE)
setwd(wdir)

dati.geni <- subset(dati.geni, BRCA_ensembl != "?") # delete unknown genes
global_geni <- data.frame(dati.geni)
rm(dati.geni)

func.double.rows <- function(data){
  
  nms <- as.character(data[,1])
  freq.geni <- table(data[,1])
  label <- names(which(freq.geni > 1))
  pos.doppi <- which(nms %in% label)
  label.doppi <- factor(nms[pos.doppi])
  nomi_unici <- setdiff(nms, label)
  
  if(length(label) != 0){
    A2 <- data[-pos.doppi,-1]; rownames(A2) <- nms[-pos.doppi]
    A2 <- rbind(apply(data[pos.doppi,-1], 2, function(.x) tapply(.x, label.doppi, mean)), A2)
  }else{
    A2 <- data[,-1]; rownames(A2) <- nms
  }
  
  return(A2)
}
global_geni <- func.double.rows(data = global_geni)
print(dim(global_geni))

nms <- dimnames(global_geni)
dati.geni <- normalize.quantiles(as.matrix(global_geni))
rm(global_geni)
dimnames(dati.geni) <- nms
print(dim(dati.geni))

screening <- function(data, perc = .95){
  
  sums <- apply(data, 1, sum)
  tokeep <- which(sums > 0)
  data <- data[tokeep, ]
  
  id0 <- apply(data > 0, 1, mean)
  tokeep <- which(id0 >= perc)
  data <- data[tokeep, ]
  
  # summary(sums[tokeep])
  
  return(data)
}
dati.geni <- screening(data = dati.geni, perc = 0)
dati.geni <- as.matrix(dati.geni)
print(dim(dati.geni))

minv <- min(dati.geni[dati.geni > 0]) / 2
dati.geni[dati.geni == 0] <- minv
dati.geni.log <- log2(dati.geni)
write.csv2(dati.geni, paste0(TISSUE, "_final_dataaset.csv"))
write.csv2(dati.geni.log, paste0(TISSUE, "_final_dataaset_log.csv"))

par(mfrow=c(2,3))
hist(dati.geni[1,]); hist(dati.geni[2,]); hist(dati.geni[3,])
hist(dati.geni.log[1,]); hist(dati.geni.log[2,]); hist(dati.geni.log[3,])
par(mfrow=c(1,1))


save.image(paste0(TISSUE, "_loading_and_preprocessing_", format(Sys.Date(), "%d%m%y"), ".RData"))