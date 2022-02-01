## Load requisite packages
library(lumi)
library(data.table)
library(ggplot2)
library(dplyr)
library(imputeTS)
library(glmnet)

## Read in CpG subsets (non-mQTLs) based on standard deviations 
top10 = read.csv("/Top10k_SD.csv")
top20 = read.csv("/Top20k_SD.csv")
top50 = read.csv("/Top50k_SD.csv")
all = read.csv("/All_nomqtls_SD.csv")


## Read in w3 methylation data 
meth = readRDS("/GS/GS_methylation/wave3_mvals.rds")

## Read in EPIC methylation annotation file 
anno = readRDS("/Daniel/EPIC_AnnotationObject.rds")
anno = as.data.frame(anno)

## Subset annotation file to 450k probes 
length(which(anno$Methyl450_Loci %in% "TRUE"))
anno1 = anno[which(anno$Methyl450_Loci %in% "TRUE"),]
meth1 = meth[which(row.names(meth) %in% anno1$Name),] 
meth1[which(is.nan(meth1))] <- NA
meth1[which(is.infinite(meth1))] <- NA


# Exclude cross-hybridising and polymorphic EPIC probes 
snps <- read.table("/mQTL_and_Illumina_annotation_files/epic_probes_to_exclude/snp_probes.txt", sep='\t', header=T)
snps <- snps[which(snps$EUR_AF >= 0.05), "IlmnID"] %>% as.character
ch1 <- read.table("/mQTL_and_Illumina_annotation_files/epic_probes_to_exclude/nonspecific_cph.txt", sep='\t', header=F)
ch1 <- as.character(ch1[,1])
ch2 <- read.table("/mQTL_and_Illumina_annotation_files/epic_probes_to_exclude/nonspecific_cpg.txt", sep='\t', header=F)
ch2 <- as.character(ch2[,1])
ychr <- anno[anno$chr=="chrY", "Name"]
exclude <- c(snps, ch1, ch2, ychr) %>% unique
meth1_resid = meth1[-which(row.names(meth1) %in% exclude),]

## Create subsets 

meth1_10k = meth1_resid[,which(colnames(meth1_resid) %in% top10$CpG)]
meth1_20k = meth1_resid[,which(colnames(meth1_resid) %in% top20$CpG)]
meth1_50k = meth1_resid[,which(colnames(meth1_resid) %in% top50$CpG)]
meth1_nomqtls = meth1_resid[,which(colnames(meth1_resid) %in% all$CpG)]

## Read in phenotypes 

phenos = readRDS("/GS_pheno_resids_for_OSCA_DNAm_VC_26May2021.rds")
names(phenos)[1] <- "Sample_Name"
id = read.csv("/GS/GS_methylation/wave3-final/samplesheet.final.csv")
id <- id[,c("Sample_Name", "Sample_Sentrix_ID")]
names(id) <- c("Sample_Name", "ID")
phenos <- merge(id, phenos, by = "Sample_Name")
phenos$Sample_Name <- NULL
phenos = phenos[which(phenos$ID %in% row.names(meth1_resid)),]

## make directories for biglasso 
setwd("/Rob/CpG_Array/Biglasso/")
dfs_names <- c("Top10k", "Top20k", "Top50k", "Meth_nomQTLs") 

dfs <- list(df1 = meth1_10k, df2 = meth1_20k, df3 = meth1_50k, df4 = meth1_nomqtls, df5 = meth1_resid) 
dfs_names <- c("Top10k", "Top20k", "Top50k", "Meth_nomQTLs", "Full_Meth") 

library(biglasso)

for(i in 1:length(dfs)){ 
  j = 2
  meth_tmp = dfs[[i]]
  
  tmp = phenos[,c(1,j)]
  tmp = tmp[!is.na(tmp[,2]),]
  
  a = which(row.names(meth_tmp) %in% tmp$ID)
  
  meth_tmp = meth_tmp[a,]
  
  ids = rownames(meth_tmp)
  tmp1 = tmp[match(ids, tmp$ID),]
  
  x1 = meth_tmp[,apply(meth_tmp, 2, function(tar) !any(is.na(tar)))]
  x1 = as.big.matrix(x1)
  
  #Cross-validation: find the best shrinkage value - the alpha=1 means it's a LASSO model
  
  y = as.numeric(tmp1[,2])
  cvfit <- tryCatch(
    {
      cv.biglasso(x1, y, seed = 1.234, alpha =1, family = "gaussian", nfolds = 10, ncores = 4)
    },
    error = function(cond) {
      cv.biglasso(x1, y, seed = 1.234, alpha = 1, family = "gaussian", nfolds = 10, ncores = 2)
    }
  )
  
  fit <- biglasso(x1, y, family = "gaussian", alpha = 1, lambda = cvfit$lambda.min)
  
  coefs = coef(cvfit)[which(coef(cvfit) != 0),]
  coefs = as.data.frame(coefs)
  saveRDS(coefs, paste0(dfs_names[[i]], "/", names(phenos)[[j]], "/", names(phenos)[[j]], "_biglasso.rds"))
  print(i)
  
} 

