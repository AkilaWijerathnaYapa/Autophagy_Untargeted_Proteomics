# Load Libraries
library(tidyverse)
library(data.table)
library(DEP)
library(fitdistrplus)
library(SummarizedExperiment)
library(ggplot2)

##############################
##### Data preparation #######
##############################

## In data file Change Protein.ID column's AGI up to first (Keep only first AGI) !!! Do this manually
### collapsed the AGIs into genes in excel (removed the .1s)

# Load data
data <- fread('data/proteinGroups.txt')
data <- filter(data, Reverse != "+", `Potential contaminant` != "+")
colnames(data)

##############################
##### Data preparation #######
##############################

#extract gene names from fasta header and sort
data2 <- data %>%
  rowwise() %>% 
  mutate(Gene.names=paste(filter(as.tibble(unlist(strsplit(`Fasta headers`,'|',fixed = T))),str_detect(value,'AT.G')==F,str_detect(value,'Chr')==F,is.na(as.numeric(value)))$value, collapse = ';'))
data2 <- data2[,c(1,2,ncol(data2),3:(ncol(data2)-1))]

#check for duplicated gene names and remove isoform duplications in gene names          
data$Gene.names %>% duplicated() %>% any()

data2 <- data2 %>% 
  rowwise() %>% 
  mutate(Gene.names=paste(unique(unlist(strsplit(Gene.names, ';'))),collapse=';'))

#new gane names

new_gene.names <- fread('data/agi.gene.namesED.csv', header = T)


data2 <- separate_rows(data2,`Majority protein IDs`, convert = T) %>% 
  left_join(new_gene.names, by=c('Majority protein IDs'='AGI'))

## choose best name available (gene.names.y over gene.names.x)

data2 <- data2 %>% 
  mutate(Gene.names.x=ifelse(is.na(Gene.names.y)==T,Gene.names.x,Gene.names.y)) %>% 
  dplyr::rename(Gene.names=Gene.names.x) %>% 
  dplyr::select(-Gene.names.y)


# Make unique names using the annotation in the "Gene.names" column as primary names and the annotation in "Protein.IDs" as name for those that do not have an gene name.
data_unique <- make_unique(data2, "Majority protein IDs", "Gene.names", delim = ";")


####################################################
##### Generate a SummarizedExperiment object #######
####################################################


# Generate a SummarizedExperiment object by parsing condition information from the column names
LFQ_columns <- grep("LFQ.", colnames(data_unique)) # get LFQ column numbers
data_se <- make_se_parse(data_unique, LFQ_columns)
data_se


######################################
##### Filter on missing values #######
######################################

plot_frequency(data_se)

# Filter for proteins that are identified in all replicates of at least one condition
#data_filt <- filter_missval(data_se, thr = 0)

# Less stringent filtering:
# Filter for proteins that are identified in 2 out of 3 replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 1)


# Filter for proteins that are quantified in at least 2/3 of the samples.
#frac_filtered <- filter_proteins(data_se, "fraction", min = 0.66)
#data_filt <- frac_filtered 

plot_numbers(data_se)
# Plot a barplot of the number of identified proteins per samples
plot_numbers(data_filt)

# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_filt)

# Normalize the data
data_norm <- normalize_vsn(data_filt)

# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm)

#visualize normal distribution of data
#normal data
data_for_dist <- data_filt@assays$`.->data`@listData %>% as.data.frame() %>% na.omit() %>% gather(sample,value)
descdist(data_for_dist$value)

#normalised
data_for_dist_norm <- data_norm@assays$`.->data`@listData %>% as.data.frame() %>% na.omit() %>% gather(sample,value)
descdist(data_for_dist_norm$value)


############################################
##### Impute data for missing values #######
############################################

#Explore the pattern of missing values
##To explore the pattern of missing values in the data, a heatmap can be plotted indicating whether values are missing (0) or not (1). Only proteins with at least one missing value are visualized.

# Plot a heatmap of proteins with missing values
plot_missval(data_norm)

# Plot a heatmap of proteins with missing values
plot_missval(data_filt)

#The missing values seem to be randomly distributed across the samples (MAR). However, we do note a block of values that are missing in all control samples (bottom left side of the heatmap). These proteins might have missing values not at random (MNAR).

#To check whether missing values are biased to lower intense proteins, the densities and cumulative fractions are plotted for proteins with and without missing values.

# Plot intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(data_filt)
plot_detect(data_norm)
### Imputation options
#https://bioconductor.org/packages/release/bioc/vignettes/DEP/inst/doc/MissingValues.html
# Base Functions and Classes for Mass Spectrometry and Proteomics 
## https://bioconductor.org/packages/3.8/bioc/html/MSnbase.html

# All possible imputation methods are printed in an error, if an invalid function name is given.
impute(data_norm, fun = "")

## Error in match.arg(fun): 'arg' should be one of "bpca", "knn", "QRILC", "MLE", "MinDet", "MinProb", "man", "min", "zero", "mixed", "nbavg"
# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_filt, fun = "MinProb", q = 0.01)

# Impute missing data using random draws from a manually defined left-shifted Gaussian distribution (for MNAR)
data_imp_man <- impute(data_filt, fun = "man", shift = 1.8, scale = 0.3)

# Impute missing data using the k-nearest neighbour approach (for MAR)
data_imp_knn <- impute(data_filt, fun = "knn", rowmax = 0.9)

##??? Error in impute.knn(exprs(object), ...) : a column has more than 80 % missing values!


# Plot intensity distributions before and after imputation
plot_imputation(data_norm, data_imp)

#plot_imputation(data_norm, data_imp_knn) # imputation with method "k-nearest

plot_imputation(data_norm, data_imp_man) # imputation with method man

# check which imputation method is better
jpeg("imputation.jpeg", height = 6, width = 10, res = 300, units = "in") 
plot_imputation(data_se, data_imp, data_imp_man, data_imp_knn)
dev.off()
pdf("imputation.pdf", height = 6, width = 10)
plot_imputation(data_norm, data_imp, data_imp_man)#, data_imp_knn)
dev.off()


## According to the plots "data_imp" seems the best option for data imputation for this data set


# export imputated data and run statistic analysis
# extract imputated data and transformed to a dataframe

data_imp_190907 <- as.data.frame(assays(data_imp)) 

write.csv(data_imp_190907, "18APR_imp_190907.csv")

# DEP Package statistics
##############################################
##### Differential enrichment analysis #######
##############################################

## with "data_imp_knn" ##


# Differential enrichment analysis  based on linear models and empherical Bayes statistics

# Test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "wt_")

## Tested contrasts: Ubi4_vs_Ctrl, Ubi6_vs_Ctrl, Ubi1_vs_Ctrl
# Test all possible comparisons of samples
#data_diff_all_contrasts <- test_diff(data_imp, type = "all")

# Test manually defined comparisons
data_diff_manual <- test_diff(data_imp, type = "manual", 
                              test = c('atg2_0d__vs_WT_0d_', 
                                         'atg5_0d__vs_WT_0d_',
                                         'atg7_0d__vs_WT_0d_',
                                         'atg9_0d__vs_WT_0d_',
                                         'atg2_1d__vs_WT_1d_', 
                                         'atg5_1d__vs_WT_1d_',
                                         'atg7_1d__vs_WT_1d_',
                                         'atg9_1d__vs_WT_1d_',
                                         'atg2_3d__vs_WT_3d_', 
                                         'atg5_3d__vs_WT_3d_',
                                         'atg7_3d__vs_WT_3d_',
                                         'atg9_3d__vs_WT_3d_',
                                         'atg2_5d__vs_WT_5d_', 
                                         'atg5_5d__vs_WT_5d_',
                                         'atg7_5d__vs_WT_5d_',
                                         'atg9_5d__vs_WT_5d_'))

data_diff <- data_diff_manual

## Tested contrasts: e.g. Ubi4_vs_Ctrl, Ubi6_vs_Ctrl

# Denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))

######################
##### PCA PLOT #######
######################

# Plot the first and second principal components
plot_pca(dep, x = 1, n=nrow(dep),y = 2, point_size = 4)

plot_pca(dep, x = 1, n=nrow(dep),y = 2, indicate = "condition", point_size = 4)

plot_pca(dep, x = 1, n=nrow(dep),y = 2, indicate = "condition", label = TRUE, point_size = 4)

########################
##### Cor Matrix #######
########################

# Plot the Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

#####################
##### Heatmap #######
#####################

# Plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 2, show_row_names = T, row_font_size = 4)

# Plot a heatmap of all significant proteins (rows) and the tested contrasts (columns)
plot_heatmap(dep, type = "contrast", kmeans = TRUE, 
             k = 6, col_limit = 2, show_row_names = TRUE, row_font_size = 4)


##########################
##### Volcano Plot #######
##########################

# Plot a volcano plot for the contrast "Ubi6 vs Ctrl""
plot_volcano(dep, contrast = "atg2_0d__vs_WT_0d_", label_size = 3, add_names = TRUE)
plot_volcano(dep, contrast = "atg2_12h__vs_WT_12h_", label_size = 3, add_names = TRUE)
plot_volcano(dep, contrast = "atg2_24h__vs_WT_24h_", label_size = 3, add_names = TRUE)

plot_volcano(dep, contrast = "atg5_0h__vs_WT_0h_", label_size = 3, add_names = TRUE)
plot_volcano(dep, contrast = "atg5_12h__vs_WT_12h_", label_size = 3, add_names = TRUE)
plot_volcano(dep, contrast = "atg5_24h__vs_WT_24h_", label_size = 3, add_names = TRUE)

plot_volcano(dep, contrast = "atg7_0h__vs_WT_0h_", label_size = 3, add_names = TRUE)
plot_volcano(dep, contrast = "atg7_12h__vs_WT_12h_", label_size = 3, add_names = TRUE)
plot_volcano(dep, contrast = "atg7_24h__vs_WT_24h_", label_size = 3, add_names = TRUE)

plot_volcano(dep, contrast = "atg9_0h__vs_WT_0h_", label_size = 3, add_names = TRUE)
plot_volcano(dep, contrast = "atg9_12h__vs_WT_12h_", label_size = 3, add_names = TRUE)
plot_volcano(dep, contrast = "atg9_24h__vs_WT_24h_", label_size = 3, add_names = TRUE)

###############################################
##### Barplots of a protein of interest #######
###############################################

#change levels
dep_save <- dep
dep <- dep_save
#dep@colData$condition <- gsub('d_', 'd',dep@colData$condition)
dep@colData$condition <- gsub('h_', 'h',dep@colData$condition)

dep@colData$condition <- factor(dep@colData$condition, 
                                levels = c('WT_0d','WT_1d','WT_3d','WT_5d',
                                           'atg2_0d','atg2_1d','atg2_3d','atg2_5d',
                                           'atg5_0d','atg5_1d','atg5_3d','atg5_5d',
                                           'atg7_0d','atg7_1d','atg7_3d','atg7_5d',
                                           'atg9_0d','atg9_1d','atg9_3d','atg9_5d'))


# Plot a barplot for cruciferin 2 and cruciferin 3
plot_single(dep, proteins = 'AT1G03880.1')

plot_single(dep, proteins = 'ATCG00020.1', type = "contrast")

plot_single(dep, proteins = ' glutathione S-transferase 7 ', type = "contrast")

# Plot a barplot for the protein USP15 with the data centered
plot_single(dep, proteins = 'ATCG00020.1', type = "centered")

plot_single(dep, proteins = 'AT1G03880.1', type = "centered")

plot_single(dep, proteins = 'ATCG00020.1', type = "centered")

###########################
##### Results table #######
###########################


# Generate a results table
data_results <- get_results(dep)

# Number of significant proteins
data_results %>% filter(significant) %>% nrow()

# Column names of the results table
colnames(data_results)

get_results

write.csv(data_results, 'APR19_0.05_LFC1.5_190903.csv')


###########################
##### Merge Tables  #######
###########################


library(dplyr)
c <- read.csv("data/MapBIN_ED.csv", header = T, stringsAsFactors = F)
d <- read.csv("data/suba4_ED.csv", header = T, stringsAsFactors = F)
e <- read.csv("data/ED_NOV18_0.05_LFC1.5_190805.csv", header = T, stringsAsFactors = F)

MapBIN_SUBA_LFQ <- left_join(e, c, by = "AGI") %>% left_join(d, by = "AGI")

#MapBIN_SUBA <- dplyr::full_join(c,d, by = "AGI")

#MapBIN_SUBA_LFQ <- dplyr::full_join(MapBIN_SUBA,e, by = "AGI")

write.csv(MapBIN_SUBA_LFQ, "18NOVDEP_Results_190805.csv")


###################################
##custom volcano plot from jakob##
##################################

## custom mutant overlap volcano

#data

# Denote significant proteins based on user defined cutoffs
#dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))
#data_results <- get_results(dep)

#long format(gather)
data_long <- data_results %>% gather(sample,value,3:ncol(.)) %>% 
  mutate(type=sapply(strsplit(sample,'__'),'[',3),
         type=ifelse(is.na(type)==T,sapply(strsplit(sample,'__'),'[',2),type),
         sample=sapply(strsplit(sample,'__'),'[',1),
         genotype=sapply(strsplit(sample,'_'),'[',1),
         time=sapply(strsplit(sample,'_'),'[',2)) %>% 
  dplyr::select(-sample) %>%
  rowwise() %>% 
  filter(type != 'significant') %>% 
  spread(key = type,value = value)


#filter for both mutants under 0.1 pvalue, and use lowest FC and highest pvalue
data_volcano <-
  data_long %>%
  filter(genotype !='WT') %>%
  group_by(ID,time,genotype) %>%
  mutate(colour_p=ifelse(max(p.adj) <= 0.05,'red','blue')) %>%
  mutate(min_ratio=min(abs(ratio)),
         max_p=max(p.adj)) %>%
  filter(abs(ratio)==min_ratio) %>%
  dplyr::select(-min_ratio) %>%
  mutate(colour_r=ifelse(ratio <=-0.4 | ratio >= 0.4,'red','blue')) %>%
  mutate(sig=ifelse(colour_p=='blue'|colour_r=='blue','non_sig','sig')) %>%
  distinct(ID,fraction, .keep_all = T)



#levels
data_volcano$sig <- factor(data_volcano$sig, levels=c('sig','non_sig'))

library(ggrepel)

#volcano_atg9_24h
selected_gene_volcano <- filter(data_volcano,genotype=='atg9',time=='24h')
title=paste(unique(selected_gene_volcano$genotype),'vs WT at',(unique(selected_gene_volcano$time)),sep=' ')

p <- ggplot(selected_gene_volcano, aes(x=ratio,y=-log10(max_p),col=sig))+
  geom_point(pch=18,alpha=0.75)+
  geom_text_repel(data=filter(selected_gene_volcano,sig=='sig'),aes(label=ID),col='black',size=1.5, fontface='bold')+
  scale_colour_manual(values=c('#990000','#99ccff'))+
  theme(legend.position = 'none')+
  labs(title=title, x='log2 fold change', y='-log10 p-value')

ggsave(filename = paste0('volcano_atg9_24h','.png'),path = 'images',device = 'png',dpi=1080,plot = p)

#volcano_atg5_24h
selected_gene_volcano <- filter(data_volcano,genotype=='atg5',time=='24h')
title=paste(unique(selected_gene_volcano$genotype),'vs WT at',(unique(selected_gene_volcano$time)),sep=' ')

p <- ggplot(selected_gene_volcano, aes(x=ratio,y=-log10(max_p),col=sig))+
  geom_point(pch=18,alpha=0.75)+
  geom_text_repel(data=filter(selected_gene_volcano,sig=='sig'),aes(label=ID),col='black',size=1.5, fontface='bold')+
  scale_colour_manual(values=c('#990000','#99ccff'))+
  theme(legend.position = 'none')+
  labs(title=title, x='log2 fold change', y='-log10 p-value')

ggsave(filename = paste0('volcano_atg5_24h','.png'),path = 'images',device = 'png',dpi=1080,plot = p)


#hand curate some annotation

data_volcano <- data_volcano %>%
  ungroup() %>% 
  mutate(ID=ifelse(name=='AT4G16160.1.1','TIM17/TIM22/TIM23',ID),
         ifelse(ACC=='AT5G62270','mucin related AT5G62270',
                ifelse(ACC=='AT1G26460','PPR AT1G26460',
                       ifelse(ACC=='AT3G02650','PPR AT3G02650',
                              ifelse(ACC=='AT5G64670','PPR AT5G64670',desc)))))

