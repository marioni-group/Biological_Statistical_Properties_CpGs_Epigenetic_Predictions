## Load relevant libraries 
library(ggplot2)
library(wesanderson)

## Extract phenotype names and read in files with var. explained by full array
	full = list.files("Full_Array/SD/", ".")
	phenos = gsub("_full.rsq", "", full)

## Reading in and extract important information 
	full_array = lapply(full, function(x) read.table(paste0("Full_Array/SD/", x),header = T, fill = T))
	full_array1 = lapply(full_array, `[`,4,)
	ns = lapply(full_array, `[`,10,)

## Tidy up lists and combine
	names_for_full <- gsub("_full.*", "", full)
	names(full_array1) <- names_for_full
	names(ns) <- names_for_full
	full_array2 <- do.call("rbind", full_array1)
	ns2 <- do.call("rbind", ns)
	names(ns2)[2] <- "n"
	full_array2$Source = NULL 
	ns2$Source = NULL 
	ns2$SE = NULL 
	full_array2$Trait <- row.names(full_array2)
	ns2$Trait = row.names(full_array2)
	full_df = merge(full_array2,ns2,by = "Trait")
	full_df$Probe_Set = "450k"
	full_df1 <- full_df[,c(2,3,4,5,1)]

	
## Read in probe subset results 
## Create lists to store results 
	list_tmp <- list()
	subset_list <- list()

	
## Loop through phenotypes 
for(i in phenos){ 

## Directory to each phenotype and its results - here, residualised phenotypes 
	tmp = list.files(paste0(i, "/residualised/SD/REML/"), ".")
	tmp = tmp[-which(tmp %in% "Permutations")]
## Apply function to read files within each subdirectory (each phenotype)
	list_tmp = lapply(tmp, function(x) read.table(paste0(i, "/residualised/SD/REML/",x),header = T, fill = T))
## Extract n 
	n = list_tmp[[1]][10,2]
## Extract only variance/se for normalised variance explained 
	list1 = lapply(list_tmp, `[`,4, )
## Tidy up names of files within subdirectory
	tmp_trim = gsub("_SD.*", "", tmp)
	tmp_trim = gsub("_bod.*", "", tmp_trim)
	names(list1) <- tmp_trim
## Store variance results for each phenotype in a dataframe
	df = do.call("rbind", list1)
## Tidy up dataframe
	df$Source = NULL 
	df$n = n 
	df$Probe_Set = row.names(df)
	df$Trait = i 
## Store all phenotypes in a common list 
	subset_list[[i]] <- df 
} 

## Combine all the phenotypes together 
subset_list <- do.call("rbind", subset_list)
total_df = rbind(full_df1, subset_list)
total_df = total_df[order(total_df$Trait, total_df$Probe_Set),]
write.csv(total_df, "/REML_outputs.csv",row.names = F)


####### PLOTS ##########

total_df_plot = total_df[total_df$Probe_Set %in% c("450k", "All_nomqtls", "Top50k", "Top20k", "Top10k") ,]

var <- total_df_plot
var$Probe_Set <- factor(var$Probe_Set, levels = c("Top10k", "Top20k", "Top50k", "All_nomqtls", "450k"))
var$Variance <- var$Variance*100
var$SE = var$SE*100
var$lower = var$Variance-1.96*var$SE
var$upper = var$Variance+1.96*var$SE


## Figure for main text 
var1 = var[var$Trait %in% c("PckYrs", "bmi", "Alc", "fat"),]
var1$Trait = factor(var1$Trait, levels = c("fat", "Alc", "bmi", "PckYrs"))

pal = wes_palette(5, name = "Chevalier1", type = "continuous") 
pdf("Fig_REML.pdf", width = 11, height = 8.5)

x = ggplot(subset(var1, !is.na(var1$Variance)), aes(fill=Probe_Set, y=Variance, x=Trait)) + 
  geom_bar(position="dodge", stat="identity") + scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100), limits = c(0,100))  + scale_x_discrete(labels=c("PckYrs" = "Smoking Pack Years", "bmi" = "Body Mass Index","Alc" = "Alcohol consumption", "fat" = "Body Fat %")) +
  geom_errorbar(aes(ymin = lower, ymax = upper),position=position_dodge(.9), width = 0.25) + 
  xlab("Trait") + ylab("Proportion of Variance Explained %") + theme(legend.title.align = 0.5) + labs(fill = "CpG Set") + scale_fill_manual(labels = c("Top 10k non-mQTLs", "Top 20k non-mQTLs", "Top 50k non-mQTLs", "Variable non-mQTLs", "Full Array"), values = pal) + theme(text = element_text(size=14.5)) + theme(legend.text=element_text(size=12))
print(x)
dev.off()

