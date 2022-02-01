## Make binary methylation file based on different CpG lists from CV analyses 
require(dplyr)
require(imputeTS)

## Read in methylation file 
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
snps <- read.table("/Cluster_Filespace/Marioni_Group/mQTL_and_Illumina_annotation_files/epic_probes_to_exclude/snp_probes.txt", sep='\t', header=T)
snps <- snps[which(snps$EUR_AF >= 0.05), "IlmnID"] %>% as.character
ch1 <- read.table("/Cluster_Filespace/Marioni_Group/mQTL_and_Illumina_annotation_files/epic_probes_to_exclude/nonspecific_cph.txt", sep='\t', header=F)
ch1 <- as.character(ch1[,1])
ch2 <- read.table("/Cluster_Filespace/Marioni_Group/mQTL_and_Illumina_annotation_files/epic_probes_to_exclude/nonspecific_cpg.txt", sep='\t', header=F)
ch2 <- as.character(ch2[,1])
ychr <- anno[anno$chr=="chrY", "Name"]
exclude <- c(snps, ch1, ch2, ychr) %>% unique
meth1 <- t(meth1)
meth1 = meth1[,-which(colnames(meth1) %in% exclude)]


## Read in covariate information
id = read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/wave3-final/samplesheet.final.csv")
id <- id[,c("Sample_Name", "Sample_Sentrix_ID", "age", "sex", "Batch")]
## Merge in covariate info 
meth1_resid <- as.data.frame(meth1_resid)
meth1_resid = meth1_resid[,which(!colnames(meth1_resid) %in% exclude)]
id_meth = row.names(meth1_resid)
id1 <- id[match(id_meth, id$Sample_Sentrix_ID),]
table(id1$Sample_Sentrix_ID == row.names(meth1_resid))

## Make qcovar and covar files 
 quant_cov <- data.frame(FID = id1$Sample_Sentrix_ID, 
                       IID = id1$Sample_Sentrix_ID,
                       age = id1$age)
write.table(quant_cov, file="quant.qcov", row.names=F, sep=' ', quote = F)


fact_cov <- data.frame(FID =id1$Sample_Sentrix_ID,
                      IID = id1$Sample_Sentrix_ID,
                      sex = id1$sex,
                      Batch = id1$Batch)
write.table(fact_cov, file="factors.cov", row.names=F, sep=' ', quote = F)


## Prepare data for binary methylation file 
osca_dat <- data.frame(IID=row.names(dat1), FID=row.names(dat1))
osca_dat <- cbind(osca_dat, dat1)
fwrite(osca_dat, file="w3_unresidualised_profile_for_bod.txt", row.names=F, sep=' ')

## Make phenotype files 
x = readRDS("/GS_pheno_resids_for_OSCA_DNAm_VC_26May2021.rds")
names(x)[1] <- "Sample_Name"
id = read.csv("samplesheet.final.csv")
id <- id[,c("Sample_Name", "Sample_Sentrix_ID")]
names(id)[2] <- "FID"
id$IID <- id$FID
x1 <- merge(x, id, by = "Sample_Name")
x1$Sample_Name <- NULL


## Make loop for phenotype file prep 

for(i in 2:19){ 
  tmp = x1[,c(20, 21, i)]
  if(length(which(is.na(x1[,i]))) > 0){ 
    tmp = tmp[-which(is.na(tmp[,3])),]
    tmp[,3] <- scale(tmp[,3])
    write.table(tmp, paste0("Phenotypes/", names(x1)[i], ".phen"), row.names=F, sep=' ', quote= F)
    
    } 
  tmp[,3] <- scale(tmp[,3])
  write.table(tmp, paste0("Phenotypes/",  names(x1)[i], ".phen"), row.names=F, sep=' ', quote = F)
  ids_phen = as.data.frame(tmp$IID)
  ids_phen$FID = ids_phen[,1]
  write.table(ids_phen, paste0("/Cluster_Filespace/Marioni_Group/Rob/CpG_Array/", names(x1)[i], ".list"), row.names = F, col.names = F, sep = "\t", quote =F)
  
}


### out of R in terminal - ideally open up new window, as remaking the .opi file will require the above files in R
### so its best to leave that window open as this runs for 20 minutes, makes a blank .opi file that needs overwriting
osca_Linux --efile w3_unresidualised_profile_for_bod.txt --methylation-m --covar factors.cov --qcovar quant.qcov --adj-probe --make-bod --out w3_profile


## open up R again 
## Manually make .opi file 

opi <- anno[colnames(osca_dat)[3:ncol(osca_dat)],c("chr","Name", "pos","UCSC_RefGene_Name", "strand")]
opi$chr <- gsub("chr", "", opi$chr)
opi$chr <- as.numeric(opi$chr)
opi$UCSC_RefGene_Name	 <- as.factor(opi$UCSC_RefGene_Name	)
opi$strand <- as.factor(opi$strand)
opi[which(opi$UCSC_RefGene_Name	==""), "UCSC_RefGene_Name"] <- NA

write.table(opi, file="w3_profile_unresidualised.opi", 
               col.names=F, 
               row.names=F, 
               quote=F, sep='\t')


## Make probe list files 
loop = list.files("CpG_Lists_SD/", ".csv")
for(i in 1:length(loop)){ 
## Read in probe list based on SD criteria 
tmp = read.csv(paste0("CpG_Lists_SD/", loop[i]))
## Take out CpGs
tmp_cpg = as.data.frame(tmp$CpG)
## Fix up file
names(tmp_cpg) <- "CV_cpgs"
## Take out file name
B = gsub(".csv", "", loop[i])
## save probe list file
write.table(tmp_cpg, paste0("Probe_Lists_OSCA/", B, ".list"), row.names = F, col.names = F, sep = "\t", quote =F)
print(i)
} 

## Make individual BODs for each probe list generated from SD analysis 
cd /Probe_Lists_OSCA/
resid=residualised
stringList=bmi,pack_years,cholest,glucose,HDL,sodium,potassium,urea,creatinine,whr,fat,sBP,dBP,HR,FEV,FVC,Alc,PckYrs,g


for r in ${resid//,/ }
do
for j in ${stringList//,/ }
do
   echo $j
for i in *.list
do
  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --befile ../../w3_profile --extract-probe $i --keep ../../${j}.list --make-bod --out ../../${j}/${r}/SD/BOD_Subsets/${B}_bod

echo $i

done

done

done 


## Make individual ORMs for each probe list generated from the SD analysis 

resid=residualised
stringList=bmi,pack_years,cholest,glucose,HDL,sodium,potassium,urea,creatinine,whr,fat,sBP,dBP,HR,FEV,FVC,Alc,PckYrs,g

for r in ${resid//,/ }
do

for j in ${stringList//,/ }
do
   echo $j

for i in ${j}/${r}/SD/BOD_Subsets/*.oii

do

  A=$( echo $i | cut -d"/" -f5)
  B=$( echo $A | cut -d. -f1)
  C=$( echo $B | cut -d"_" -f -4)

osca_Linux --befile ${j}/${r}/SD/BOD_Subsets/$B --make-orm --out ${j}/${r}/SD/ORM/${C}_orm

echo $i

done

done

done 


## REML analysis for each probe list generated from the SD analysis 

resid=residualised
stringList=bmi,pack_years,cholest,glucose,HDL,sodium,potassium,urea,creatinine,whr,fat,sBP,dBP,HR,FEV,FVC,Alc,PckYrs,g

for r in ${resid//,/ }
do

for j in ${stringList//,/ }
do
echo $j

for i in ${j}/${r}/SD/ORM/*.id

do

A=$( echo $i | cut -d"/" -f5)
B=$( echo $A | cut -d. -f1)
C=$( echo $B | cut -d"_" -f -4)

osca_Linux --reml --orm ${j}/${r}/SD/ORM/$B --pheno Phenotypes/${j}.phen --out ${j}/${r}/SD/REML/${C}_reml

echo $i

done

done

done


