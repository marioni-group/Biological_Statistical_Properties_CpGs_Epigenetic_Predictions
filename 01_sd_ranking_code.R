## Load packages
library(lumi)
library(data.table)
library(ggplot2)
library(dplyr)
library(imputeTS)

### Read in w3 methylation data 
meth = readRDS("wave3_mvals.rds")

## Read in EPIC methylation annotation file 
anno = readRDS("EPIC_AnnotationObject.rds")
anno = as.data.frame(anno)

## Subset annotation file to 450k probes 
length(which(anno$Methyl450_Loci %in% "TRUE"))
anno1 = anno[which(anno$Methyl450_Loci %in% "TRUE"),]
meth1 = meth[which(row.names(meth) %in% anno1$Name),] 
meth1[which(is.nan(meth1))] <- NA
meth1[which(is.infinite(meth1))] <- NA


# Exclude cross-hybridising and polymorphic EPIC probes 
snps <- read.table("mQTL_and_Illumina_annotation_files/epic_probes_to_exclude/snp_probes.txt", sep='\t', header=T)
snps <- snps[which(snps$EUR_AF >= 0.05), "IlmnID"] %>% as.character
ch1 <- read.table("mQTL_and_Illumina_annotation_files/epic_probes_to_exclude/nonspecific_cph.txt", sep='\t', header=F)
ch1 <- as.character(ch1[,1])
ch2 <- read.table("mQTL_and_Illumina_annotation_files/epic_probes_to_exclude/nonspecific_cpg.txt", sep='\t', header=F)
ch2 <- as.character(ch2[,1])
ychr <- anno[anno$chr=="chrY", "Name"]
exclude <- c(snps, ch1, ch2, ychr) %>% unique
meth1 = meth1[-which(row.names(meth1) %in% exclude),]

## Convert m values to beta values 
meth1 <- m2beta(meth1)

## Remove CpGs with mean methylation <10% or >90% in sample 
mean_mat = matrix(nrow = nrow(meth1), ncol = 2)
mean_mat = as.data.frame(mean_mat)

for(i in 1:nrow(meth1)){ 
  mean_mat[i,1] <- row.names(meth1)[i]
  mean_mat[i,2] <- mean(meth1[i,], na.rm = T)
  print(i)
}

names(mean_mat) <- c("CpG", "mean")
mean_mat1 = mean_mat[which(mean_mat$mean <= 0.1 | mean_mat$mean >= 0.9),]
meth2 = meth1[-which(row.names(meth1) %in% mean_mat1$CpG),]


## Read in mQTLs from GoDMC 
mqtls <- read.csv("GoDMC_mqtls_genomewide_CpGlist.csv")
## obtain non-mQTL CpGs 
meth2 = meth2[-which(row.names(meth2) %in% mqtls$GoDMC_GenomeWide_CpGs),]
mat = matrix(nrow = nrow(meth2), ncol = 2)
mat = as.data.frame(mat)


## Run loop to residualise data for age, sex and batch 
meth1_resid <- meth2

## Read in covariate information
id = read.csv("samplesheet.final.csv")
id <- id[,c("Sample_Name", "Sample_Sentrix_ID", "age", "sex", "Batch")]

## Merge in covariate info 
meth1_resid <- as.data.frame(meth1_resid)
id_meth = colnames(meth1_resid)
id1 <- id[match(id_meth, id$Sample_Sentrix_ID),]
table(id1$Sample_Sentrix_ID == colnames(meth1_resid))
meth1_resid = t(meth1_resid)

## Set up loop for residualisation 
for(i in 1:ncol(meth1_resid)){ 
  meth1_resid[,i] <- resid(lm(as.numeric(meth1_resid[,i]) ~ as.numeric(id1$age) + as.factor(id1$sex) + as.factor(id1$Batch), na.action = na.exclude))
  print(i)
} 
dat <- meth1_resid 
dat <- na_mean(dat)

## Run loop to get CpGs and statistics 
for(i in 1:nrow(meth2)){ 
  mat[i,1] <- row.names(meth2)[i]
  mat[i,2] <- sd(meth2[i,], na.rm = T)
  print(i)
}


## Tidy up dataframe
names(mat) <- c("CpG", "sd")
## Order the dataset by SD
mat1 <- mat[rev(order(mat$sd)),]

## Save out top 10,000, 20,000 and 50,000
write.csv(mat1[1:1e4,], "Top10k_SD.csv", row.names = F)
write.csv(mat1[1:2e4,], "Top20k_SD.csv", row.names = F)
write.csv(mat1[1:5e4,], "Top50k_SD.csv", row.names = F)
write.csv(mat1, "All_nomqtls_SD.csv", row.names = F)

