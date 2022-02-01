## Full Set ##
library(lumi)
library(data.table)
library(ggplot2)
library(dplyr)
library(imputeTS)


### Read in w3 methylation data 
meth = readRDS("wave3_mvals.rds")


## Read in EPIC methylation annotation file 
anno = readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject.rds")
anno = as.data.frame(anno)


## Subset annotation file to 450k probes 
length(which(anno$Methyl450_Loci %in% "TRUE"))
anno1 = anno[which(anno$Methyl450_Loci %in% "TRUE"),]
meth1 = meth[which(row.names(meth) %in% anno1$Name),] 
meth1[which(is.nan(meth1))] <- NA
meth1[which(is.infinite(meth1))] <- NA
meth1 = t(meth1)
meth1 = na_mean(meth1)
meth1 = as.data.frame(meth1)
osca_dat <- data.frame(FID=row.names(meth1), IID=row.names(meth1))
osca_dat <- cbind(osca_dat, meth1)
fwrite(osca_dat, "full_450k_w3.txt", row.names=F, sep=' ')


## In terminal, prepare BOD file 
cd /Cluster_Filespace/Marioni_Group/Rob/CpG_Array/
osca_Linux --efile full_450k_w3.txt --methylation-m --make-bod --out w3_profile_450k


## In R, correct .opi file from BOD analysis 
opi <- anno[colnames(osca_dat)[3:ncol(osca_dat)],c("chr","Name", "pos","UCSC_RefGene_Name", "strand")]
opi$chr <- gsub("chr", "", opi$chr)
opi$chr <- as.numeric(opi$chr)
opi$UCSC_RefGene_Name	<- as.factor(opi$UCSC_RefGene_Name) 
opi$strand <- as.factor(opi$strand)
opi[which(opi$UCSC_RefGene_Name	==""), "UCSC_RefGene_Name"] <- NA

 write.table(opi, file="w3_profile_450k.opi", 
	             col.names=F, 
 	             row.names=F, 
 	             quote=F, sep='\t')

 
## In terminal, run OSCA 

osca_Linux --befile w3_profile_450k --covar factors.cov --qcovar quant.qcov --adj-probe --make-bod --out w3_profile_450k_resid
osca_Linux --befile w3_profile_residualised  --make-orm --out w3_profile_residualised_orm

## REML 
stringList=bmi,pack_years,cholest,glucose,HDL,sodium,potassium,urea,creatinine,whr,fat,sBP,dBP,HR,FEV,FVC,Alc,PckYrs,g
for j in ${stringList//,/ }
do
osca_Linux --reml --orm w3_profile_residualised_orm --pheno Phenotypes/${j}.phen --out Full_Array/SD/${j}_full
done 

