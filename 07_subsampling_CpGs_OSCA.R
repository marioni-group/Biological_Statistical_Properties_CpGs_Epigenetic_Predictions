#### Sub-sampling analysis #### 
library(tidyverse)
meth = readRDS("wave3_mvals.rds")

## Read in EPIC methylation annotation file 
anno = readRDS("EPIC_AnnotationObject.rds")
anno = as.data.frame(anno)

## Subset annotation file to 450k probes 
length(which(anno$Methyl450_Loci %in% "TRUE"))
anno1 = anno[which(anno$Methyl450_Loci %in% "TRUE"),]
meth1 = meth[which(row.names(meth) %in% anno1$Name),] 


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


## Create function to randomly samples over 10,000, 20,000, 50,000 or ~115,000 rows or CpGs 
## 10k
get.samples.10k <- function(x) {  
  sample(row.names(x), 1e4)
}

## 20k
get.samples.20k <- function(x) {  
  sample(row.names(x), 2e4)
}

## 50k
get.samples.50k <- function(x) {  
  sample(row.names(x), 5e4)
}

## 115k
get.samples.115k <- function(x) {  
  sample(row.names(x), 115746)
}

## Set.seed for random sampling process 
set.seed(1.234)

## Generate a random sample 1000 times for 10k, 20k, 50k probes and 115k probes 

## 10k
df_10k <- replicate(1000,{
  get.samples.10k(meth1)
} 
) 

## 20k
df_20k <- replicate(1000,{
  get.samples.20k(meth1)
} 
) 

## 50k
df_50k <- replicate(1000,{
  get.samples.50k(meth1)
  } 
) 

## 115k
df_115k <- replicate(1000,{
  get.samples.115k(meth1)
} 
) 


## create dataframe for ease of saving out results 
df_10k = as.data.frame(df_10k)
df_20k = as.data.frame(df_20k)
df_50k = as.data.frame(df_50k)
df_115k = as.data.frame(df_115k)


## tidy up dataframe 
names(df_10k) <- paste("Perm", 1:ncol(df_10k), "10k", sep = "_")
names(df_20k) <- paste("Perm", 1:ncol(df_20k), "20k", sep = "_")
names(df_50k) <- paste("Perm", 1:ncol(df_50k), "50k", sep = "_")
names(df_115k) <- paste("Perm", 1:ncol(df_115k), "115k", sep = "_")


## Save out as individual probe sets to create ORMs and run REMLs based on them 

for(i in 1:ncol(df_10k)){ 
  tmp_cpg = as.data.frame(df_10k[,i])
  B = names(df_10k)[i]
  write.table(tmp_cpg, paste0("Permutations/Small_Sets/", B, ".list"), row.names = F, col.names = F, sep = "\t", quote =F)
}

for(i in 1:ncol(df_20k)){ 
  tmp_cpg = as.data.frame(df_20k[,i])
  B = names(df_20k)[i]
  write.table(tmp_cpg, paste0("/Permutations/Small_Sets/", B, ".list"), row.names = F, col.names = F, sep = "\t", quote =F)
}



for(i in 1:ncol(df_50k)){ 
tmp_cpg = as.data.frame(df_50k[,i])
B = names(df_50k)[i]
write.table(tmp_cpg, paste0("/Permutations/", B, ".list"), row.names = F, col.names = F, sep = "\t", quote =F)
}

for(i in 1:ncol(df_115k)){ 
  tmp_cpg = as.data.frame(df_115k[,i])
  B = names(df_115k)[i]
  write.table(tmp_cpg, paste0("/Permutations/", B, ".list"), row.names = F, col.names = F, sep = "\t", quote =F)
}




############# ORM ANALYSES ##################


## Out of R, set up ORMs based on permuted probe sets 

###### 50k and 115k #######################


## Group 1 - first 200 .... repeated up until group 10 
## replace x with number of group 
screen -S perm1
cd /Permutations/
  resid=residualised

for r in ${resid//,/ }
do

x=1
x1=`expr $x - 1`
y=`expr 200 \* $x1`

FILES=(*)
for i in "${FILES[@]:$y:200}"


do

A=$( echo $i | cut -d"/" -f3)
B=$( echo $A | cut -d. -f1)
C=$( echo $B | cut -d"_" -f -4)

osca_Linux --befile ../../w3_profile --extract-probe ../../CpG_Lists_SD/Permutations/$i --make-orm --out ../../Permutation_ORMs/${C}_orm

echo $i

done

done


###### 10k and 20k #######################


## Group 1 - first 200 .... repeated up until group 10 
## replace x with number of group 
screen -S perm1
cd /Permutations/Small_Sets/
  
  resid=residualised

for r in ${resid//,/ }
do

x=1
x1=`expr $x - 1`
y=`expr 200 \* $x1`

FILES=(*)
for i in "${FILES[@]:$y:200}"


do

A=$( echo $i | cut -d"/" -f3)
B=$( echo $A | cut -d. -f1)
C=$( echo $B | cut -d"_" -f -4)

osca_Linux --befile ../../../w3_profile --extract-probe ../../../CpG_Lists_SD/Permutations/Small_Sets/$i --make-orm --out ../../../Permutation_ORMs/Small_Sets/${C}_orm

echo $i

done

done


############# REML ANALYSES ##################

## 115k and 50k 
## two traits shown for reference on how to build loop 
cd /Permutation_ORMs/
  
  resid=residualised
stringList=bmi,cholest
for r in ${resid//,/ }
do

for j in ${stringList//,/ }
do
echo $j

for i in *.id
do

A=$( echo $i | cut -d"/" -f5)
B=$( echo $A | cut -d. -f1)
C=$( echo $B | cut -d"_" -f -3)
D=$( echo $C | cut -d"_" -f 3)


osca_Linux --reml --orm $B --pheno ../Phenotypes/${j}.phen --out ../${j}/${r}/Permutation_${D}/REML/${C}_reml

done

done

done



#### 10k and 20k ##########


screen -S perms1
cd /Permutation_ORMs/Small_Sets/
  
  resid=residualised
stringList=bmi,cholest
for r in ${resid//,/ }
do

for j in ${stringList//,/ }
do
echo $j

for i in *.id
do

A=$( echo $i | cut -d"/" -f5)
B=$( echo $A | cut -d. -f1)
C=$( echo $B | cut -d"_" -f -3)
D=$( echo $C | cut -d"_" -f 3)


osca_Linux --reml --orm $B --pheno ../../Phenotypes/${j}.phen --out ../../${j}/${r}/Permutation_${D}/REML/${C}_reml

done

done

done


########## COLLATE RESULTS ############

## Get phenotype names 
x = readRDS("/GS_pheno_resids_for_OSCA_DNAm_VC_26May2021.rds")
names(x)[1] <- "Sample_Name"

id = read.csv("/samplesheet.final.csv")
id <- id[,c("Sample_Name", "Sample_Sentrix_ID")]
names(id)[2] <- "FID"
id$IID <- id$FID
x1 <- merge(x, id, by = "Sample_Name")
x1$Sample_Name <- NULL
phenos = names(x1)[2:19]

## Obtain directories 
perm.loop = c("Permutation_10k", "Permutation_20k", "Permutation_50k", "Permutation_115k")
perm.loop.names = c("Top10k", "Top20k", "Top50k", "All_nomqtls")


final.list <- list()
for(i in perm.loop){ 
out_mat <- matrix(nrow = length(phenos), ncol = 6)
out_mat <- as.data.frame(out_mat)
names(out_mat) <- c("Trait", "Variance", "SE", "n", "Rank", "P")
for(j in 1:length(phenos)){ 

## Extract rsq files from relevant phenotype and permutation set   
tmp.dir = paste0(phenos[j],"/residualised/", i,"/REML/")
list.tmp = list.files(tmp.dir, ".rsq")

## Read all .rsq files 
tmp.dir = lapply(list.tmp, function(x) read.table(paste0(phenos[j],"/residualised/", i,"/REML/", x), header = T, fill = T))
full_array1 = lapply(tmp.dir, `[`,4,)
ns = lapply(tmp.dir, `[`,10,)

## Tidy up lists and combine
helper = gsub("Permutation_", "", i)
names_for_full = gsub(paste0("_",helper,".*", collapse = " "), "", list.tmp)
names_for_full = gsub("Perm_", "", names_for_full)
names(full_array1) <- names_for_full
names(ns) <- names_for_full
full_array2 <- do.call("rbind", full_array1)
ns2 <- do.call("rbind", ns)
names(ns2)[2] <- "n"
full_array2$Source = NULL 
ns2$Source = NULL 
ns2$SE = NULL 
full_array2$Permutation_Set <- row.names(full_array2)
ns2$Permutation_Set = row.names(full_array2)
full_df = merge(full_array2,ns2,by = "Permutation_Set")
full_df$Trait = as.character(phenos[j])
full_df1 <- full_df[,c(2,3,4,5,1)]

## Read in corresponding result from SD ranking (i.e. the non-permuted result) 
ind = which(perm.loop %in% i)
comp = read.table(paste0(phenos[j], "/residualised/SD/REML/", perm.loop.names[ind], "_SD_bod_orm_reml.rsq"), header = T, fill  = T)
comp1 = comp[4,c(2,3)]
comp1$n <- comp[10,2]
comp1$Trait <- phenos[j]
comp1$Permutation_Set <- "SD_ranked"
comp1 <- comp1[,c(5,1,2,3,4)]
total = rbind(comp1, full_df)

total$rank = round(rank(-total$Variance),0)

total$P = total$rank/nrow(total)

out_mat[j,1] <- phenos[j]
out_mat[j,c(2,3,4,5,6)] <- total[1,c(2,3,4,6,7)] 
} 
final.list[[i]] <- out_mat 
} 


perm.loop = c("Permutation_115k")
perm.loop.names = c("All_nomqtls")


final.list <- list()
for(i in perm.loop){ 
  out_mat <- matrix(nrow = length(phenos), ncol = 6)
  out_mat <- as.data.frame(out_mat)
  names(out_mat) <- c("Trait", "Variance", "SE", "n", "Rank", "P")
  for(j in 1:length(phenos)){ 
    
    ## Extract rsq files from relevant phenotype and permutation set   
    tmp.dir = paste0(phenos[j],"/residualised/", i,"/REML/")
    list.tmp = list.files(tmp.dir, ".rsq")
    
    ## Read all .rsq files 
    tmp.dir = lapply(list.tmp, function(x) read.table(paste0(phenos[j],"/residualised/", i,"/REML/", x), header = T, fill = T))
    full_array1 = lapply(tmp.dir, `[`,4,)
    ns = lapply(tmp.dir, `[`,10,)
    
    ## Tidy up lists and combine
    helper = gsub("Permutation_", "", i)
    names_for_full = gsub(paste0("_",helper,".*", collapse = " "), "", list.tmp)
    names_for_full = gsub("Perm_", "", names_for_full)
    names(full_array1) <- names_for_full
    names(ns) <- names_for_full
    full_array2 <- do.call("rbind", full_array1)
    ns2 <- do.call("rbind", ns)
    names(ns2)[2] <- "n"
    full_array2$Source = NULL 
    ns2$Source = NULL 
    ns2$SE = NULL 
    full_array2$Permutation_Set <- row.names(full_array2)
    ns2$Permutation_Set = row.names(full_array2)
    full_df = merge(full_array2,ns2,by = "Permutation_Set")
    full_df$Trait = as.character(phenos[j])
    full_df1 <- full_df[,c(2,3,4,5,1)]
    
    ## Read in corresponding result from SD ranking (i.e. the non-permuted result) 
    ind = which(perm.loop %in% i)
    comp = read.table(paste0(phenos[j], "/residualised/SD/REML/", perm.loop.names[ind], "_SD_bod_reml.rsq"), header = T, fill  = T)
    comp1 = comp[4,c(2,3)]
    comp1$n <- comp[10,2]
    comp1$Trait <- phenos[j]
    comp1$Permutation_Set <- "SD_ranked"
    comp1 <- comp1[,c(5,1,2,3,4)]
    total = rbind(comp1, full_df)
    
    total$rank = round(rank(-total$Variance),0)
    
    total$P = total$rank/nrow(total)
    
    out_mat[j,1] <- phenos[j]
    out_mat[j,c(2,3,4,5,6)] <- total[1,c(2,3,4,6,7)] 
  } 
  final.list[[i]] <- out_mat 
} 


### Age collation 

perm.loop = c("Permutation_10k", "Permutation_20k", "Permutation_50k", "Permutation_115k")
perm.loop.names = c("Top10k", "Top20k", "Top50k", "All_nomqtls")

out_mat <- matrix(nrow = length(perm.loop), ncol = 5)
out_mat <- as.data.frame(out_mat)
names(out_mat) <- c("Trait", "Variance", "SE", "Rank", "P")

for(i in perm.loop){ 
## Extract rsq files from relevant phenotype and permutation set   
tmp.dir = paste0("age","/residualised/", i,"/REML/")
list.tmp = list.files(tmp.dir, ".rsq")

## Read all .rsq files 
tmp.dir = lapply(list.tmp, function(x) read.table(paste0("age","/residualised/", i,"/REML/", x), header = T, fill = T))
full_array1 = lapply(tmp.dir, `[`,4,)

## Tidy up lists and combine
helper = gsub("Permutation_", "", i)
names_for_full = gsub(paste0("_",helper,".*", collapse = " "), "", list.tmp)
names_for_full = gsub("Perm_", "", names_for_full)
names(full_array1) <- names_for_full
full_array2 <- do.call("rbind", full_array1)
full_array2 = full_array2[which(full_array2$Variance > 0.010),]
full_array2$Source = NULL 
full_array2$Permutation_Set <- row.names(full_array2)
full_array2$Trait = as.character("age")
full_df <- full_array2

## Read in corresponding result from SD ranking (i.e. the non-permuted result) 
ind = which(perm.loop %in% i)
comp = read.table(paste0("Age/REML/", perm.loop.names[ind], "_SD_orm_reml.rsq"), header = T, fill  = T)
comp1 = comp[4,c(2,3)]
comp1$Trait <- "age"
comp1$Permutation_Set <- "SD_ranked"
comp1 <- comp1[,c(1,2,4,3)]
total = rbind(comp1, full_df)

total$rank = round(rank(-total$Variance),0)

total$P = total$rank/nrow(total)

out_mat[ind,1] <- as.character(i)
out_mat[ind,c(2,3,4,5)] <- total[1,c(3,4,5,6)] 
} 