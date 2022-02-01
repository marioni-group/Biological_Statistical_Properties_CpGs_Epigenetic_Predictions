### Features of Sub-samples vs. Top 10k set #####
## For top 5k, simply subset 'top10k' file to first 5,000 ##

## Read in Top 10k non-mQTL CpGs  
top10k = read.table("Probe_Lists_OSCA/Top10k_SD.list")


## Read in EPIC array information 
anno = readRDS("EPIC_AnnotationObject.rds")
anno = as.data.frame(anno)
genes <- data.frame(ID=anno$Name, geneSymbol=anno$UCSC_RefGene_Name, CHR=anno$chr, MAPINFO=anno$pos, FEATURE=anno$UCSC_RefGene_Group, CpGISLAND=anno$Relation_to_Island)
genes1 <- data.frame(ID=anno1$Name, geneSymbol=anno1$UCSC_RefGene_Name, CHR=anno1$chr, MAPINFO=anno1$pos, FEATURE=anno1$UCSC_RefGene_Group, CpGISLAND=anno1$Relation_to_Island)
genes$FEATURE <- gsub(";.*", "" ,genes$FEATURE)
genes1$FEATURE <- gsub(";.*", "" ,genes1$FEATURE)


## Get regulatory features of top 10k probes 
anno_tmp = genes[which(genes$ID %in% top10k[,1]),]
## Get proportion of each feature in top 10k probes 
tmp_10k = as.data.frame(table(anno_tmp$FEATURE)/nrow(anno_tmp))
tmp_10k$Freq = tmp_10k$Freq*100

## Make dataframe of proportions of 1000 10k sub-samples  
setwd("/Permutations/")
## get 10k files 
loop = list.files(".", ".list")[grep("10k", list.files(".", ".list"))]

## make dataframe for output 
out_mat <- matrix(nrow = nrow(tmp_10k), ncol = length(loop))
out_mat = as.data.frame(out_mat)

for(i in loop){ 
  ind = which(loop %in% i)
  tmp = read.table(i, header = F)
  anno_perm = genes[which(genes$ID %in% tmp[,1]), ]
  ## Get proportion of each feature in top 10k probes 
  tmp_perm = as.data.frame(table(anno_perm$FEATURE)/nrow(anno_perm))
  tmp_perm$Freq = tmp_perm$Freq*100
  if(nrow(tmp_perm) == 7){ 
  ids = tmp_10k$Var1
  tmp_perm[match(ids, tmp_perm$Var1),]
  out_mat[,ind] <- tmp_perm[,2]
  print(ind)
  } else { 
    out_mat[,ind] <- NA
    print(ind)
    }
}

res <- out_mat 
## Rank observed proportions against sub-samples
res = as.data.frame(t(res))
names(res) <- tmp_10k$Var1
sds = as.data.frame(t(tmp_10k))[2,]
names(sds) <- tmp_10k$Var1
comb = rbind(sds, res)

