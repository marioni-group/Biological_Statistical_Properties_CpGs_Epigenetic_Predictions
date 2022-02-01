## Load requisite packages 
library(lumi)
library(data.table)
library(ggplot2)
library(dplyr)
library(imputeTS)
library(glmnet)

## Read in CpG subsets (non-mQTLs) based on standard deviations 

top10 = read.table("/Probe_Lists_OSCA/Top10k_SD.list")
top20 = read.table("/Probe_Lists_OSCA/Top20k_SD.list")
top50 = read.table("/Probe_Lists_OSCA/Top50k_SD.list")
all = read.table("/Probe_Lists_OSCA/All_nomqtls_SD.list")

## Read in DNAm file and 
x = readRDS("w3_meth_agesexbatch_resid.rds")


## Create subsets
meth1_10k = meth1_resid[,which(colnames(meth1_resid) %in% top10$V1)]
meth1_20k = meth1_resid[,which(colnames(meth1_resid) %in% top20$V1)]
meth1_50k = meth1_resid[,which(colnames(meth1_resid) %in% top50$V1)]
meth1_nomqtls = meth1_resid[,which(colnames(meth1_resid) %in% all$V1)]


## Read in phenotypes 
phenos = readRDS("GS_pheno_resids_for_OSCA_DNAm_VC_26May2021.rds")
names(phenos)[1] <- "Sample_Name"
id = read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/wave3-final/samplesheet.final.csv")
id <- id[,c("Sample_Name", "Sample_Sentrix_ID")]
names(id) <- c("Sample_Name", "ID")
phenos <- merge(id, phenos, by = "Sample_Name")
phenos$Sample_Name <- NULL
phenos = phenos[which(phenos$ID %in% row.names(meth1_resid)),]


setwd("/Biglasso/")
dfs_names <- c("Top10k", "Top20k", "Top50k", "Meth_nomQTLs")
dfs <- list(df1 = meth1_10k, df2 = meth1_20k, df3 = meth1_50k, df4 = meth1_nomqtls) 

## Biglasso loop 

library(biglasso)
for(i in 1:length(dfs)){ 
  for(j in 3:ncol(phenos)){ 
print(j)
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
  

} 
} 