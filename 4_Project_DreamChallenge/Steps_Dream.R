#CHECK THAT ALL ROUTES ARE CORRECT

###################################################
# FUNCTION PCA
###################################################
PCAfunction <- function(all_genes, metadata, number){
  PCA  = prcomp(t(all_genes))
  PCA_df = as.data.frame(PCA$x)
  #adding a column with the names of samples because I want to merge the metadata
  PCA_df$SampleID = as.factor(rownames(PCA_df))
  PCA_df =  merge(PCA_df, metadata, by='SampleID')
  #some plots
  byTrain <- ggplot(PCA_df,aes(x=PC1, y = PC2, col=Train)) + geom_point() + ggtitle('Train/Test plot')
  variable <- paste("../../../../../home/esther/Desktop/PCA/byTrain", number, sep="_")
  ggsave(paste(variable, ".png"))
  byGroup <- ggplot(PCA_df,aes(x=PC1, y = PC2, col=Group)) + geom_point() + ggtitle('Group plot')
  variable <- paste("../../../../../home/esther/Desktop/PCA/byGroup", number, sep="_")
  ggsave(paste(variable, ".png"))
  byGA <- ggplot(PCA_df,aes(x=PC1, y = PC2, col=GA)) + geom_point() + ggtitle('Group plot')
  variable <- paste("../../../../../home/esther/Desktop/PCA/byGA", number, sep="_")
  ggsave(paste(variable, ".png"))
  outNULL <- ggplot(PCA_df[!is.na(PCA_df$Group),],aes(x=PC1, y = PC2, col=Group)) + geom_point() + ggtitle('Group plot')
  variable <- paste("../../../../../home/esther/Desktop/PCA/withoutNULL_Group", number, sep="_")
  ggsave(paste(variable, ".png"))
  
  #semiautomated plot : pay attention to (paste0 ,and aes_string in ggplot call) 
  # plot PCX and PC(X+1) . There are better ways in the above links  a completely different method
  
  for (pc_idx in 1:10){
    idx1 = paste0('PC',pc_idx)
    idx2 = paste0('PC',pc_idx+1)
    p =  ggplot(PCA_df[!is.na(PCA_df$Group),],aes_string(x=idx1, y = idx2, col='Group')) + geom_point() + 
      ggtitle(paste('Group plot',idx1,idx2))  + theme_bw()
    variable <- paste("../../../../../home/esther/Desktop/PCA/byGroupNoNULL", number, sep="_")
    ggsave(paste(variable, ".png"))
  }
}

###################################################
# FUNCTION RTSNE
###################################################
RTSNEfunction <- function(all_genes, metadata, number){
  tsne  = Rtsne(t(all_genes))
  tsne_df = as.data.frame(tsne$Y)
  #adding a column with the names of samples because I want to merge the metadata
  tsne_df$SampleID = as.factor(colnames(all_genes))
  tsne_df =  merge(tsne_df, metadata, by='SampleID')
  #some plots
  byTrain <- ggplot(tsne_df,aes(x=V1, y = V2, col=Train)) + geom_point() + ggtitle('Train/Test plot')
  variable <- paste("../../../../../home/esther/Desktop/RTSNE/tsne", number, sep="_")
  ggsave(paste(variable,  ".png"))
  byGroup <- ggplot(tsne_df,aes(x=V1, y = V2, col=Group)) + geom_point() + ggtitle('Group plot')
  variable <- paste("../../../../../home/esther/Desktop/RTSNE/tsne_group", number, sep="_")
  ggsave(paste(variable, ".png"))
  byGA <- ggplot(tsne_df,aes(x=V1, y = V2, col=GA)) + geom_point() + ggtitle('Group plot')
  variable <- paste("../../../../../home/esther/Desktop/RTSNE/tsne_byGA", number, sep="_")
  ggsave(paste(variable, ".png"))
  outNULL <- ggplot(tsne_df[!is.na(tsne_df$Group),],aes(x=V1, y = V2, col=Group)) + geom_point() + ggtitle('Group plot')
  variable <- paste("../../../../../home/esther/Desktop/RTSNE/tsne_withoutNULL_Group", number, sep="_")
  ggsave(paste(variable, ".png"))
  
  #semiautomated plot : pay attention to (paste0 ,and aes_string in ggplot call) 
  # plot PCX and PC(X+1) . There are better ways in the above links  a completely different method
  
  for (tsne_idx in 1:10){
    idx1 = paste0('TSNE',tsne_idx)
    idx2 = paste0('TSNE',tsne_idx+1)
    p =  ggplot(tsne_df[!is.na(tsne_df$Group),],aes_string(x=idx1, y = idx2, col='Group')) + geom_point() + 
      ggtitle(paste('Group plot',idx1,idx2))  + theme_bw()
    variable <- paste("../../../../../home/esther/Desktop/RTSNE/byGroupNoNULL", number, sep="_")
    ggsave(paste(variable, ".png"))
  }
}

###################################################
# FUNCTION FASTICA
###################################################
FASTICAfunction <- function(all_genes, metadata, number){
  fastica <- fastICA(all_genes, 2, alg.typ = "parallel", fun = "logcosh", alpha = 1, 
               method = "C", row.norm = FALSE, maxit = 200, 
               tol = 0.0001, verbose = TRUE)
  fastica_df = as.data.frame(fastica$X)
  #adding a column with the names of samples because I want to merge the metadata
  xx <- t(fastica_df)
  fastica_df$SampleID = as.factor(colnames(all_genes))
  fastica_df =  merge(fastica_df, metadata, by='SampleID')
  #some plots
  byTrain <- ggplot(tsne_df,aes(x=V1, y = V2, col=Train)) + geom_point() + ggtitle('Train/Test plot')
  variable <- paste("../../../../../home/esther/Desktop/FASTICA/fastica", number, sep="_")
  ggsave(paste(variable, ".png"))
  byGroup <- ggplot(tsne_df,aes(x=V1, y = V2, col=Group)) + geom_point() + ggtitle('Group plot')
  variable <- paste("../../../../../home/esther/Desktop/FASTICA/fastica_group", number, sep="_")
  ggsave(paste(variable, ".png"))
  byGA <- ggplot(tsne_df,aes(x=V1, y = V2, col=GA)) + geom_point() + ggtitle('Group plot')
  variable <- paste("../../../../../home/esther/Desktop/FASTICA/fastica_byGA", number, sep="_")
  ggsave(paste(variable, ".png"))
  outNULL <- ggplot(tsne_df[!is.na(tsne_df$Group),],aes(x=V1, y = V2, col=Group)) + geom_point() + ggtitle('Group plot')
  variable <- paste("../../../../../home/esther/Desktop/FASTICA/fastica_withoutNULL_Group", number, sep="_")
  ggsave(paste(variable, ".png"))
  
  #semiautomated plot : pay attention to (paste0 ,and aes_string in ggplot call) 
  # plot PCX and PC(X+1) . There are better ways in the above links  a completely different method
  
  for (fastica_idx in 1:10){
    idx1 = paste0('FASTICA',fastica_idx)
    idx2 = paste0('FASTICA',fastica_idx+1)
    p =  ggplot(fastica_df[!is.na(fastica_df$Group),],aes_string(x=idx1, y = idx2, col='Group')) + geom_point() + 
      ggtitle(paste('Group plot',idx1,idx2))  + theme_bw()
    variable <- paste("../../../../../home/esther/Desktop/FASTICA/byGroupNoNULL", number, sep="_")
    ggsave(paste(variable, ".png"))
  }
}
###################################################
# Remove all data loaded
###################################################
# > rm(list=ls())

###################################################
# How to install some libraries
###################################################
# > install.packages("tidyverse", dependencies = TRUE) #out this ggplot2 has an error

# > source("https://bioconductor.org/biocLite.R")
# > biocLite("preprocessCore")
# > biocLite("limma")

###################################################
# Libraries
###################################################
library(Rtsne)
library(stats)
library(ggplot2)
library(fastICA)
library(umap)
library(gridExtra)
library(corrplot)
library(dplyr)
library(tidyr)
library(preprocessCore)
library(limma)
library(fastICA)
library(Rtsne)

###################################################
# Load data 
###################################################
#all samples included in SC1 and SC2 involved in data preprocessing;
allSamplesSC1and2=read.csv("data/allSamplesSC1and2.csv",stringsAsFactors = FALSE)   

#defines what samples will be used in subchallenge 2 among those in allSamplesSC1and2
metadata = read.csv('data/anoSC2_v20_nokey.csv')

load('data/HTA20_RMA.RData') #this is the full HTA2.0 RMA preprocessed set  
#obtained from CELL files as described in preprocess_data_SC1.R (sub-challenge 1)  
HTA20 <- as.data.frame(eset_HTA20)

#see preprocess_SC2.r to see how this dataset was obtained from .CEL files in HuGene21ST folder
load("data/HuGene21ST_RMA.RData")
HuGene21ST <- as.data.frame(eset_HuGene21ST)

###################################################
# See data 
###################################################
head(metadata)
head(eset_HTA20)
head(eset_HuGene21ST)
summary(metadata)

###################################################
# Metadata plot
###################################################
# I'm setting Train as factor because otherwise it's treated as number
# and this is a problem if I want to use it to color stuff
metadata$Train = as.factor(metadata$Train)

# histogram - BY TRAIN
# > ggplot(data=metadata,aes(GA, fill=Train)) + geom_histogram(bins=40) + ggtitle('GA per train class') # GRAPHIC HISTOGRAM

#stacked histogram - BY GROUP
# > ggplot(metadata,aes(GA, fill=Group)) + geom_histogram(bins=40) + ggtitle('Group in train ') #  multiple population

#doing the same but separating the histogram
# facet_grid is used for that:   ' xxx ~ . ' = split the plot by xxx values and put them in rows
#                                ' - ~ xxx ' = split the plot by xxx values and put them in columns
# > ggplot(metadata,aes(GA, fill=Group)) + geom_histogram(bins=30) + facet_grid(Group ~ . ) + ggtitle('GA distribution') #GA
# > ggplot(metadata,aes(TTD, fill=Group)) + geom_histogram(bins=30) + facet_grid(Group ~ . ) + ggtitle('TTD distribution') #TTD

# this shows that the test is from 1 chip only while train is from 2 ()
# > table(metadata$Train,metadata$Platform)

### Back to metadata
# I sort (increasingly) by individual 1st and then by Gestational Age
metadata  = metadata[ with(metadata, order(IndividualID, GA)),]
# now I'll put a for loop here to retrieve the interval
cur_sample = metadata$IndividualID[1]
cur_GA     = metadata$GA[1]
metadata$Interval = 0
for (i in 1:nrow(metadata)){
  if (metadata$IndividualID[i] == cur_sample)  {
    #just to the difference 
    metadata$Interval[i] = metadata$GA[i] - cur_GA
  }
  cur_sample = metadata$IndividualID[i]
  cur_GA     = metadata$GA[i]
}
# > ggplot(metadata[metadata$Interval>0,], aes(Interval, fill=Group)) +geom_histogram(binwidth = 1)+ facet_grid(Group ~ . ) + ggtitle('Interval distribution')

###################################################
# Create a metadata out test (train = 0 = no GADel 'NA' value) -> metadata_two
###################################################
#Remove all columns  NA
delete.na <- function(DF, n=0) {
  DF[rowSums(is.na(DF)) <= n,]
}
#metadata_two is metadata out train = 0
metadata_two <- delete.na(metadata)

###################################################
# Could be omit it
###################################################
#See how many patients have been sampled 2,3,N times
# > number_of_samples_per_individual = as.data.frame(table(metadata_two$IndividualID)) #The times that appears the IndividualID
#change names so that Individual matches the metadata
# > names(number_of_samples_per_individual) = c('IndividualID','Nsamples')
#now add information just for the sake of plotting 
# > number_of_samples_per_individual = merge(number_of_samples_per_individual, metadata_two, by='IndividualID')
#plot the number of samples per group type (control - risk) 
# > ggplot(number_of_samples_per_individual, aes(Nsamples, fill=Group)) +geom_histogram()+ facet_grid(Group ~ . ) + ggtitle('Nsample distribution')

#Do a table  sampleid and individualid
# > genes__individual = as.data.frame(table(metadata_two$IndividualID, metadata_two$SampleID)) #The times that appears the IndividualID
# > genes__individual <- genes_with_individual[!(genes_with_individual$Freq==0),] #remove freq=0 as not exist because no appear
# > col_to_remove = c('Freq')
# > genes_with_individual = genes_with_individual[, !(colnames(genes_with_individual) %in% col_to_remove)] #remove Freq column
# > names(genes_with_individual) = c('IndividualID','SampleID')

###################################################
# Creating metadata_HTA20 and metadata_HuGene21ST
###################################################
#metadata_HTA20 es el metadata pero sin nulls y solo los que encontrare en el dataframe de HTA20
headerHTA20 <- as.data.frame(colnames(HTA20))
names(headerHTA20) = c('SampleID') #rename column to merge it
metadata_HTA20 <- merge(headerHTA20, metadata_two, by='SampleID')

#metadata_HuGene21ST es el metadata pero sin nulls y solo los que encontrare en el dataframe de HuGene21ST
headerHuGene21ST <- as.data.frame(colnames(HuGene21ST))
names(headerHuGene21ST) = c('SampleID')
metadata_HuGene21ST <- merge(headerHuGene21ST, metadata_two, by='SampleID')

remove(headerHTA20, headerHuGene21ST)

#see plots and save it 
# > plot(metadata_two)
# > plot(metadata)
# > plot(metadata_HTA20)
# > plot(metadata_HuGene21ST)

###################################################
# Remove interval = 0 and do the plot
###################################################
#remove interval = 0
# > metadata_two_0 <- metadata_two[!(metadata_two$Interval==0),]
# > metadata_0 <- metadata[!(metadata$Interval==0),]
# > metadata_HTA20_0 <- metadata_HTA20[!(metadata_ge_df$Interval==0),]
# > metadata_HuGene21ST_0 <- metadata_HuGene21ST[!(metadata_ge_df2$Interval==0),]

#see plots without interval and save it
# > plot(metadata_two_0)
# > plot(metadata_0)
# > plot(metadata_HTA20_0)
# > plot(metadata_HuGene21ST_0)

###################################################
# Create all_genes
# See HTA20 and HuGene21ST and how many equal Genes they have + merge by,...
###################################################
HTA20$Genes <- rownames_HTA20 <- rownames(HTA20) #Insert in HTA20 a column call Genes.

HuGene21ST$Genes <-rownames_HuGene21ST <- rownames(HuGene21ST)

all_genes = merge(HTA20, HuGene21ST, by='Genes') #Put together the two tables only those genes who are in both tables.

rownames(all_genes)<- all_genes$Genes #The column Genes will be the rownames

#Remove Genes column
#In all_genes we have all the genes that are in both files(HTA20, HuGene21ST) -> Y and the Sample ID that they have -> X
all_genes$Genes = NULL 
HTA20$Genes = NULL
HuGene21ST$Genes = NULL

###################################################
# Remove innecesary things
###################################################
remove(cur_GA, cur_sample, i, rownames_HTA20, rownames_HuGene21ST, delete.na, HTA20, HuGene21ST)

###################################################
# Convert into dataframe to order
###################################################
eset_HTA20 <- as.data.frame(eset_HTA20)
eset_HuGene21ST <- as.data.frame(eset_HuGene21ST)

###################################################
# Order dataframes
###################################################
all_genes <- all_genes[ , order(names(all_genes))]
eset_HTA20 <- eset_HTA20[ , order(names(eset_HTA20))]
eset_HuGene21ST <- eset_HuGene21ST[ , order(names(eset_HuGene21ST))]

###################################################
# Adapt the data and do PCA and do it
###################################################
# DO PCAfunction -> PCA
PCAfunction(all_genes, metadata, 1)
RTSNEfunction(all_genes, metadata, 1)
FASTICAfunction(all_genes, metadata, 1)

#quantile normalize over the two sets
#take rownames in a vector/table
cgenes = c(rownames(eset_HuGene21ST),rownames(eset_HTA20))
cgenes=names(table(cgenes)[table(cgenes)==2])
all_genes = data.matrix(all_genes) #convert dataframe into a matrix
all_genes = normalize.quantiles(all_genes) #obtain quantiles
colnames <- c(colnames(eset_HuGene21ST),colnames(eset_HTA20)) #put all colnames in a vector
colnames <- sort(colnames) #order the vector like all_genes was order
colnames(all_genes) = colnames #substitute colnames from all_genes by the others
rownames(all_genes) = cgenes #substitute rownames from all_genes by the others
#extract gestatioan age value for all samples used in analysis
ga=allSamplesSC1and2$GA[match(colnames(all_genes),allSamplesSC1and2$SampleID)]

# DO PCAfunction -> PCA
PCAfunction(all_genes, metadata, 2)
RTSNEfunction(all_genes, metadata, 2)
FASTICAfunction(all_genes, metadata, 2)

#remove platform effect and GA effect 
all_genes = removeBatchEffect(all_genes,batch=ifelse(substr(colnames(all_genes),1,3)=="GSM",0,1),covariates=cbind(ga,ga^2))
# DO PCAfunction -> PCA
PCAfunction(all_genes, metadata, 3)
RTSNEfunction(all_genes, metadata, 3)
FASTICAfunction(all_genes, metadata, 3)

# > remove(allSamplesSC1and2, eset_HTA20, eset_HuGene21ST, metadata_HTA20, metadata_HuGene21ST, cgenes, colnames, ga, PCAfunction)
