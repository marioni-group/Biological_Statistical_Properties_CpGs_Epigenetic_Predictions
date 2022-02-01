## Load requisite packages 
library(lumi)
library(data.table)
library(ggplot2)
library(dplyr)
library(imputeTS)
library(wesanderson)


## Read in test set DNAm data 
meth = readRDS("GS_methylation/norm_mvals.rds")


## Subset to 5087 
subset = readRDS("/samples-5087.rds")
subset$ID = paste(subset$Sentrix_ID, subset$Sentrix_Position, sep = "_")
meth1 = meth[,which(colnames(meth) %in% subset$ID)] 


## Subset to unrelateds 
t1 = read.csv("2578_unrelateds.csv")
meth2 = meth1[,which(colnames(meth1) %in% t1$id)]

## Convert to beta values and tidy up dataframe 
meth2 <- m2beta(meth2)

## Read in phenotypes to prep for loop later 
phenos = readRDS("/Cluster_Filespace/Marioni_Group/Rob/CpG_Array/GS_pheno_resids_for_OSCA_DNAm_VC_26May2021.rds")
names(phenos)[1] <- "Sample_Name"
phenos$age <- NULL


## Read in CpGs and weights - again prep for loop 
setwd("Biglasso/")
dfs_names <- c("Full_Meth", "Meth_nomQTLs", "Top50k", "Top20k", "Top10k")

## Set up loop to calculate EpiScores for each trait and each methylation set 
## Use list to score outputs from different methylation sets 
mat <- matrix(nrow = ncol(meth2), ncol = ncol(phenos))
mat <- as.data.frame(mat)

for(i in 1:length(dfs_names)){ 
  mat[,1] <- colnames(meth2)
  for(j in 2:ncol(phenos)){ 
 ## Read in the scores from a given methylation set and trait   
 tmp = readRDS(paste0(dfs_names[[i]], "/", names(phenos)[[j]], "/", names(phenos)[[j]], "_biglasso.rds"))
 ## Several steps to tidy up dataframe just read in 
 tmp1 = as.data.frame(tmp)
 tmp1$CpG = row.names(tmp1)
 names(tmp1)[1] <- "beta"
 tmp1 = tmp1[-1,]
 tmp1 = tmp1[,c(2,1)]
 ## Subset methylation dataframe to those CpGs in the current predictor 
 tmp_meth = meth2[which(row.names(meth2) %in% tmp1$CpG),]
 ## Match order of CpGs in both datasets 
 ids = tmp1$CpG 
 tmp_meth = tmp_meth[match(ids, row.names(tmp_meth)), ]
 ## Calculation step - matrix multiplication 
 tmp_score = tmp_meth*tmp1$beta
 ## Get scores 
 scores = as.data.frame(colSums(tmp_score, na.rm = T))
 ## Tidy up scores dataframe 
 scores$ID <- row.names(scores)
 scores <- scores[,c(2,1)]
 names(scores)[2] <- names(phenos)[[j]]
 ## Store result in assigned results dataframe 
  mat[,j] <- scores[,2]
 names(mat)[[j]] <- names(scores)[2]
 names(mat)[1] <- "ID"
 write.csv(mat, paste0(dfs_names[[i]], "_scores_w1.csv"), row.names = F)
  } 
  
}

## R2 analyses 
meth_full = read.csv("/Full_Meth_scores_w1.csv")
meth_nomqtls = read.csv("/Meth_nomQTLs_scores_w1.csv")
meth_top50k = read.csv("/Top50k_scores_w1.csv")
meth_top20k = read.csv("/Top20k_scores_w1.csv")
meth_top10k = read.csv("/Top10k_scores_w1.csv")

## Prepare phenotypes for analysis 
phenos = read.csv("GS_phenos_unresidualised.csv")
names(phenos)[1] <- "Sample_Name"
phenos$age <- NULL 

## Subset phenotype file to 2578 unrelateds and extract demographic data 
id = readRDS("samples-5087.rds")
t1 = read.csv("2578_unrelateds.csv")
height = read.csv("GS_height.csv")
height = height[which(height$ID %in% t1$Sample_Name) ,]
id = id[which(id$Sample_Name %in% t1$Sample_Name),]
id$Basename = paste(id$Sentrix_ID, id$Sentrix_Position, sep = "_")
id <- id[,c("Sample_Name", "Basename", "age", "sex")]
phenos <- merge(id, phenos, by = "Sample_Name")
phenos$Sample_Name <- NULL
id2 = id$Sample_Name
height = height[match(id2, height$ID), ]


## Prepare loop for r squared analyses 
out_mat <- matrix(nrow = length(names(phenos)[4:21]),ncol = 4)
out_mat <- as.data.frame(out_mat)
names(out_mat) <- c("Predictor", "Null_rsq", "Predictor_rsq", "Incremental_rsq")
dfs_scores = list(df1 = meth_full, df2=meth_nomqtls, df3=meth_top50k, df4 = meth_top20k, df5 = meth_top10k)
  
for(i in 1:length(dfs_scores)){ 
    ## Match order of IDs in phenotype file and those in epigenetic predictor scores file 
    tmp = dfs_scores[[i]]
    ids = tmp$ID 
    phenos = phenos[match(ids, phenos$Basename), ]
    for(j in names(phenos)[4:21]){ 
      
      ## Obtain index of j amongst the predictors - used to store in output dataframe 
      ind = which(names(phenos)[4:21] %in% j)
      if(all(tmp[,j] == 0)){ 
        out_mat[ind, 1] <- j
        out_mat[ind, 2] <- 0
        out_mat[ind, 3] <- 0
        out_mat[ind, 4] <- 0
        } else { 
        if(j %in% c("FEV", "FVC")){ 
          pheno_r = summary(lm(scale(phenos[,j]) ~ phenos$age + factor(phenos$sex) + height$height, na.action = na.exclude))$r.squared
          epi_r = summary(lm(scale(phenos[,j]) ~ scale(tmp[,j]) + phenos$age + factor(phenos$sex) + height$height, na.action = na.exclude))$r.squared
          
          out_mat[ind, 1] <- j
          out_mat[ind, 2] <- pheno_r
          out_mat[ind, 3] <- epi_r
          out_mat[ind, 4] <- epi_r - pheno_r
          
          } else {   
          
        ## Calculate r squared from null model + model with epigenetic predictor 
      pheno_r = summary(lm(scale(phenos[,j]) ~ phenos$age + factor(phenos$sex), na.action = na.exclude))$r.squared
      epi_r = summary(lm(scale(phenos[,j]) ~ scale(tmp[,j]) + phenos$age + factor(phenos$sex), na.action = na.exclude))$r.squared
      
      out_mat[ind, 1] <- j
      out_mat[ind, 2] <- pheno_r
      out_mat[ind, 3] <- epi_r
      out_mat[ind, 4] <- epi_r - pheno_r
          }
    
  } 
      write.csv(out_mat, paste0(dfs_names[[i]], "_rsquared_w1.csv"), row.names = F)
      
    } 
  } 


## Plot results from incremental rsquared analyses 
meth_full = read.csv("/Full_Meth_rsquared_w1.csv")
meth_nomqtls = read.csv("/Meth_nomQTLs_rsquared_w1.csv")
meth_top50k = read.csv("/Top50k_rsquared_w1.csv")
meth_top20k = read.csv("/Top20k_rsquared_w1.csv")
meth_top10k = read.csv("/Top10k_rsquared_w1.csv")


## order the predictors by var explained 
meth_full = meth_full[order(meth_full$Incremental_rsq),]
ids = meth_full$Predictor 
meth_nomqtls = meth_nomqtls[match(ids, meth_nomqtls$Predictor),]
meth_top50k = meth_top50k[match(ids, meth_top50k$Predictor),]
meth_top20k = meth_top20k[match(ids, meth_top20k$Predictor),]
meth_top10k = meth_top10k[match(ids, meth_top10k$Predictor),]
total = cbind(cbind(cbind(cbind(meth_full, meth_nomqtls[,c(3,4)]), meth_top50k[,c(3,4)]), meth_top20k[,c(3,4)]), meth_top10k[,c(3,4)])
total[,2:ncol(total)] <- round((total[,2:ncol(total)]*100), 2)

## Get % of Full Array content 
total$Percentage_All_nomQTLs = round((total[,6]/total[,4])*100,2)
total$Percentage_Top50k = round((total[,8]/total[,4])*100,2)
total$Percentage_Top20k = round((total[,10]/total[,4])*100,2)
total$Percentage_Top10k = round((total[,12]/total[,4])*100,2)


write.csv(total, "rsquared_SuppTable.csv", row.names = F)


## Figure for main text 

pdf("Fig_R2.pdf", width = 11, height = 8.5)
var1 = data[which(data$trait %in% c("fat", "HDL", "bmi", "PckYrs")), ]
var1$trait = factor(var1$trait, levels = c("fat", "HDL", "bmi", "PckYrs"))
var1$condition = factor(var1$condition, levels = c("Top 10k non-mQTLs", "Top 20k non-mQTLs", "Top 50k non-mQTLs", "Variable non-mQTLs", "Full Array"))

x = ggplot(subset(var1, !is.na(var1$value)), aes(fill=condition, y=value, x=trait)) + 
  geom_bar(position="dodge", stat="identity") + scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100), limits = c(0,70))  + scale_x_discrete(labels=c("PckYrs" = "Smoking Pack Years", "bmi" = "Body Mass Index","HDL" = "HDL cholesterol", "fat" = "Body Fat %")) +
  xlab("Trait") +   ylab(bquote(paste("Incremental ", R^2))) + theme(legend.title.align = 0.5) + labs(fill = "CpG Set") + scale_fill_manual(labels = c("Top 10k non-mQTLs", "Top 20k non-mQTLs", "Top 50k non-mQTLs", "Variable non-mQTLs", "Full Array"), values = pal) + theme(legend.text=element_text(size=12)) + theme(axis.text=element_text(size=3),axis.title=element_text(size=14))
print(x)
dev.off()




