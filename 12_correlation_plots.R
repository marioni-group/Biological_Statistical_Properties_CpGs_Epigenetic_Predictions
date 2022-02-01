## correlation matrices of tested variables 
library(corrplot)

## Read in phenotypes
phenos = readRDS("Phenos_Summary.rds")


## Set 2 
set2 = read.csv("/samplesheet.final.csv")
phenos.set2 = phenos[which(phenos$ID %in% set2$Sample_Name), ] 
row.names(phenos.set2) <- phenos.set2$ID 
phenos.set2$ID <- NULL 
cor.set2 = round(cor(phenos.set2, use = "complete.obs"),2)
ids = row.names(cor.set2)

pdf("/Corrplot_Set2_unadjusted.pdf", width = 8, height = 8)
corrplot(cor.set2, type = "upper", tl.col="black", tl.cex = 0.8, title = "Set 2 - Raw Phenotypes",addCoef.col = "black", number.cex = 0.72,mar=c(0,0,1,0))
dev.off() 


## Set 1 
unrel = read.csv("2578_unrelateds.csv")
set1 = readRDS("samples-5087.rds")
set1 = set1[which(set1$Sample_Name %in% unrel$Sample_Name), ]
phenos.set1 = phenos[which(phenos$ID %in% set1$Sample_Name), ] 
row.names(phenos.set1) <- phenos.set1$ID 
phenos.set1$ID <- NULL 
cor.set1 = round(cor(phenos.set1, use = "complete.obs"),2)

## prepare plots 
phenos.set1 = phenos.set1[match(ids, row.names(phenos.set1)),]
phenos.set1 = phenos.set1[,match(ids, colnames(phenos.set1))]


pdf("/Corrplot_Set1_unadjusted.pdf", width = 8, height = 8)
corrplot(cor.set1, type = "upper", tl.col="black", title = "Set 1 - Raw Phenotypes", tl.cex = 0.8, addCoef.col = "black", number.cex = 0.72,mar=c(0,0,1,0))
dev.off() 



## Adjusted Variables 

phenos = readRDS("GS_pheno_resids_for_OSCA_DNAm_VC_26May2021.rds") 


## Set 2 
set2 = read.csv("samplesheet.final.csv")

phenos.set2 = phenos[which(phenos$ID %in% set2$Sample_Name), ] 
row.names(phenos.set2) <- phenos.set2$ID 
phenos.set2$ID <- NULL 
cor.set2 = round(cor(phenos.set2, use = "complete.obs"),2)

## prepare plots 
phenos.set2 = phenos.set2[match(ids, row.names(phenos.set2)),]
phenos.set2 = phenos.set2[,match(ids, colnames(phenos.set2))]

pdf("/Corrplot_Set2_adjusted.pdf", width = 8, height = 8)
corrplot(cor.set2, type = "upper", tl.col="black", title = "Set 2 - Residualised Phenotypes", tl.cex = 0.8,addCoef.col = "black", number.cex = 0.72,mar=c(0,0,1,0))
dev.off() 


## Set 1 
unrel = read.csv("2578_unrelateds.csv")
set1 = readRDS("samples-5087.rds")
set1 = set1[which(set1$Sample_Name %in% unrel$Sample_Name), ]

phenos.set1 = phenos[which(phenos$ID %in% set1$Sample_Name), ] 
row.names(phenos.set1) <- phenos.set1$ID 
phenos.set1$ID <- NULL 
cor.set1 = round(cor(phenos.set1, use = "complete.obs"),2)

## prepare plots 
phenos.set1 = phenos.set1[match(ids, row.names(phenos.set1)),]
phenos.set1 = phenos.set1[,match(ids, colnames(phenos.set1))]


pdf("Corrplot_Set1_adjusted.pdf", width = 8, height = 8)
corrplot(cor.set1, type = "upper", tl.col="black", title = "Set 1 - Residualised Phenotypes", tl.cex = 0.8, addCoef.col = "black", number.cex = 0.72,mar=c(0,0,1,0))
dev.off() 
