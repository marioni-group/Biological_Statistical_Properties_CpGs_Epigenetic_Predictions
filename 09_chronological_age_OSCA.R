#### Chronological age prediction analysis #####
library(lumi)
library(data.table)
library(ggplot2)
library(dplyr)
library(imputeTS)


### Read in w3 methylation data 
meth = readRDS("GS_methylation/wave3_mvals.rds")

## Read in EPIC methylation annotation file 
anno = readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject.rds")
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
meth1 = meth1[-which(row.names(meth1) %in% exclude),]
meth1_resid <- meth1


## prepare binary file 
dat <- meth1_resid 
dat <- na_mean(dat)
dat = as.data.frame(dat)
osca_dat <- data.frame(FID=row.names(meth1), IID=row.names(meth1))
osca_dat <- cbind(osca_dat, dat)
fwrite(osca_dat, "/DNAm_Age.txt", row.names=F, sep=' ')


## In terminal, prepare ORM  
cd /Cluster_Filespace/Marioni_Group/Rob/CpG_Array/
  osca_Linux --efile DNAm_Age.txt --methylation-m --make-bod --out w3_profile_450k


## in R, regenerate .opi file 
opi <- anno[colnames(osca_dat)[3:ncol(osca_dat)],c("chr","Name", "pos","UCSC_RefGene_Name", "strand")]
opi$chr <- gsub("chr", "", opi$chr)
opi$chr <- as.numeric(opi$chr)
opi$UCSC_RefGene_Name	<- as.factor(opi$UCSC_RefGene_Name) 
opi$strand <- as.factor(opi$strand)
opi[which(opi$UCSC_RefGene_Name	==""), "UCSC_RefGene_Name"] <- NA

write.table(opi, file="/w3_profile_450k.opi", 
            col.names=F, 
            row.names=F, 
            quote=F, sep='\t')



## In terminal, make ORM 
osca_Linux --befile w3_profile_450k  --make-orm --out w3_profile_450k_age_resid


## Make ORMs based on different probe sets 

cd /Probe_Lists_OSCA/
  
for i in *.list
do
A=$( echo $i | cut -d"/" -f5)
B=$( echo $A | cut -d. -f1)
C=$( echo $B | cut -d"_" -f -4)

osca_Linux --befile ../../w3_profile_450k --extract-probe $i --make-orm --out ../../Age_Lasso/ORM/${C}_orm

echo $i

done



## REML 
cd /Cluster_Filespace/Marioni_Group/Rob/CpG_Array/Age_Lasso/ORM/

for i in *.id
do

A=$( echo $i | cut -d"/" -f5)
B=$( echo $A | cut -d. -f1)
C=$( echo $B | cut -d"_" -f -4)

osca_Linux --reml --orm $B --pheno ../../age.phen --out ../REML/${C}_reml

echo $i

done
