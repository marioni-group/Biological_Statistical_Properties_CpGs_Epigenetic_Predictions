## Load requisite libraries 
library(lumi)
library(data.table)
library(ggplot2)
library(dplyr)
library(imputeTS)
library(glmnet)


## Read in CpG subsets (non-mQTLs) based on standard deviations 
setwd("/Biglasso_Permutations/")
top10 = list.files("/Biglasso_Permutations/", ".")
top10 = top10[grep("_10k", top10)]
top20 = list.files("/Biglasso_Permutations/", ".")
top20 = top20[grep("_20k", top20)]
top50 = list.files("/Biglasso_Permutations/")
top50 = top50[grep("_50k", top50)]
top115 = list.files("/Biglasso_Permutations/")
top115 = top115[grep("_115k", top115)]

## Take 100 sub-samples for computation analysis 
top10 = top10[1:100]
top20 = top20[1:100]
top50 = top50[1:100]
top115 = top115[1:100]


## This script is for top 10k - replace with others for other CpG number 
cpgs <- lapply(top10,function(i){ read.table(i, header=F) })
names(cpgs) = gsub(".list", "",top10)

setwd("/Biglasso/Permutations/Top10k/")
for(i in 1:length(cpgs)){ 
  dir.create(names(cpgs)[i])
for(j in 2:ncol(phenos)){ 
  dir.create(paste0(names(cpgs)[i], "/", names(phenos)[j]))
   }     
}

## Read in methylation data 
x = readRDS("/w3_meth_agesexbatch_resid.rds")


## Read in phenotypes 
phenos = readRDS("/GS_pheno_resids_for_OSCA_DNAm_VC_26May2021.rds")
names(phenos)[1] <- "Sample_Name"
id = read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/wave3-final/samplesheet.final.csv")
id <- id[,c("Sample_Name", "Sample_Sentrix_ID")]
names(id) <- c("Sample_Name", "ID")


phenos <- merge(id, phenos, by = "Sample_Name")
phenos$Sample_Name <- NULL
phenos = phenos[which(phenos$ID %in% row.names(meth1_resid)),]


library(biglasso)
for(i in 1:length(cpgs)){ 
  cpg_tmp = as.data.frame(cpgs[i])
  print(paste("Iteration", i))
  for(j in 3:ncol(phenos)){ 
    print(j)
    
    meth_tmp = meth1_resid[,which(colnames(meth1_resid) %in% cpg_tmp$V1)]
    
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
    saveRDS(coefs, paste0(names(cpgs)[i], "/", names(phenos)[[j]], "/", names(phenos)[[j]], "_biglasso.rds"))
    
    
  } 
} 

