#Dream SC2 data load and description
install.packages("corrplot", dependencies = TRUE)
# libraries 
library(Rtsne)
library(stats)
library(ggplot2)
library(fastICA)
library(umap)
library(gridExtra)
library(corrplot)


###################################################
# load data 
###################################################

metadata = read.csv('Desktop/DUAL/4_Project_DreamChallenge/data/anoSC2_v20_nokey.csv')
load('Desktop/DUAL/4_Project_DreamChallenge/data/HTA20_RMA.RData')
ge_df = as.data.frame(eset_HTA20) #eset_HTA20 is the name that have the result of load(...). Convert large matrix into data frame

#See data
head(metadata)
head(ge_df)
  
###################################################
# Metadata plot
###################################################
# I'm setting Train as factor because otherwise it's treated as number
# and this is a problem if I want to use it to color stuff
metadata$Train = as.factor(metadata$Train)

# histogram - BY TRAIN
ggplot(data=metadata,aes(GA, fill=Train)) + geom_histogram() + ggtitle('GA per train class') # GRAPHIC HISTOGRAM

#stacked histogram - BY GROUP
ggplot(metadata,aes(GA, fill=Group)) + geom_histogram(bins=40) + ggtitle('Group in train ') # With multiple population

#doing the same but separating the histogram
# facet_grid is used for that:   ' xxx ~ . ' = split the plot by xxx values and put them in rows
#                                ' - ~ xxx ' = split the plot by xxx values and put them in columns
ggplot(metadata,aes(GA, fill=Group)) + geom_histogram(bins=30) + facet_grid(Group ~ . ) + ggtitle('GA distribution')
ggplot(metadata,aes(TTD, fill=Group)) + geom_histogram(bins=30) + facet_grid(Group ~ . ) + ggtitle('TTD distribution')

# show the chip used in train and test
# this shows that the test is from 1 chip only while train is from 2 ()
table(metadata$Train,metadata$Platform)

EnsurePackage<-function(x){
  x<-as.character(x)
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x, repos="http://cran.r-project.org")
    require(x, character.only=TRUE)
  }
}
EnsurePackage("ggmap")

#try to get aggregated values for each sample 
# meaning see how many patients have been sampled 2,3,N times

genes_with_individual = as.data.frame(table(metadata_two$IndividualID, metadata_two$SampleID)) #The times that appears the IndividualID
col_to_remove = c('Freq')
genes_with_individual = genes_with_individual[, !(colnames(genes_with_individual) %in% col_to_remove)]
names(genes_with_individual) = c('IndividualID','SampleID')

number_of_samples_per_individual = as.data.frame(table(metadata$IndividualID)) #The times that appears the IndividualID
#change names so that Individual matches the metadata -> cambia el numero de las columnas de la tabla.
names(number_of_samples_per_individual) = c('IndividualID','Nsamples', 'remove')
#now add information just for the sake of plotting 
# some of the new columns doesn't make sense because there are many to one associations so I'll remove them
number_of_samples_per_individual  =merge(number_of_samples_per_individual,metadata,by='IndividualID')
cols_to_remove = c('GA','GADel','TTD')
number_of_samples_per_individual = number_of_samples_per_individual[, !(colnames(number_of_samples_per_individual) %in% cols_to_remove)]
names(number_of_samples_per_individual)

#plot the number of samples per group type (control - risk) 
ggplot(number_of_samples_per_individual, aes(Nsamples, fill=Group)) +geom_histogram()+ facet_grid(Group ~ . ) + ggtitle('Nsample distribution')
# do you think this is informative or has some trends ? 


### Back to metadata ## Add a column with the time interval between two subsequent samples
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
  ggplot(metadata[metadata$Interval>0,], aes(Interval, fill=Group)) +geom_histogram(binwidth = 1)+ facet_grid(Group ~ . ) + ggtitle('Interval distribution')

  
  
  
###################################################################
#  DIMENSIONALITY REDUCTION PLOTS
###################################################################
  
  # 1) PCA
   # https://www.datacamp.com/community/tutorials/pca-analysis-r
   # https://github.com/kevinblighe/PCAtools <<-- For a much better analysis 
   #    This guy has several tutorias for RNA data analysis that are good material to work with
    PCA  = prcomp(t(ge_df))
    PCA_df = as.data.frame(PCA$x)
    #adding a column with the names of samples because I want to merge the metadata
    PCA_df$SampleID = as.factor(rownames(PCA_df))
    PCA_df =  merge(PCA_df, metadata, by='SampleID')
    #some plots
    ggplot(PCA_df,aes(x=PC1, y = PC2, col=Train)) + geom_point() + ggtitle('Train/Test plot')
    ggplot(PCA_df,aes(x=PC1, y = PC2, col=Group)) + geom_point() + ggtitle('Group plot')
    ggplot(PCA_df[!is.na(PCA_df$Group),],aes(x=PC1, y = PC2, col=Group)) + geom_point() + ggtitle('Group plot')
  
    #semiautomated plot : pay attention to (paste0 ,and aes_string within ggplot call) 
    # plot PCX and PC(X+1) . There are better ways in the above links with a completely different method
    
    for (pc_idx in 1:10){
      idx1 = paste0('PC',pc_idx)
      idx2 = paste0('PC',pc_idx+1)
      p =  ggplot(PCA_df[!is.na(PCA_df$Group),],aes_string(x=idx1, y = idx2, col='Group')) + geom_point() + 
        ggtitle(paste('Group plot',idx1,idx2))  + theme_bw()
      print(p)
    }
    
    
  # 2) ICA Independet component analysis [fastICA but feel free to search for others]
  
  
  # 3) tSNE (Rtsne package, used in bioinfo. There are better python implementation like tsneOPEN that also allow 
             # plotting the test set in the generated embedding . If we see trends in the training set this feature could
             # be really useful)
  tsne_plot = Rtsne(t(ge_df),dims = 2, perplexity = 15)
  #adding a column with the names of samples because I want to merge the metadata
  tsne_df = as.data.frame(tsne_plot$Y)
  tsne_df$SampleID = names(ge_df)
  tsne_df =  merge(tsne_df, metadata, by='SampleID')
    
  ggplot(tsne_df,aes(x=V1, y = V2, col=Train)) + geom_point() + ggtitle('Train/Test plot')
  ggplot(tsne_df[!is.na(tsne_df$Group),],aes(x=V1, y = V2, col=Group)) + geom_point() + ggtitle('Group plot')
  
  # 4) UMAP  (similar to Rtsne but faster and better. More used in ML than bioinfo)
  #    Try to figure out how to generate more than 2 columns from the documentation (if it's possible at all)
  umap_plot = umap(t(ge_df))
  umap_df = as.data.frame(umap_plot$layout)
  umap_df$SampleID = rownames(umap_df)
  umap_df =  merge(umap_df, metadata, by='SampleID')
  
  ggplot(umap_df,aes(x=V1, y = V2, col=Train)) + geom_point() + ggtitle('Train/Test plot')
  ggplot(umap_df[!is.na(umap_df$Group),],aes(x=V1, y = V2, col=Group)) + geom_point() + ggtitle('Group plot')
  
  # 5) MDS algorithm (check for it on the net as MDS dimesionality reduction R bioinfo)
  
  
  #6) Filter the input matrix and repeat the plots 
      # usual filters are by low-expression, but for Blood I'd add very high average expression too
      # in our matrix this can be done by row-wise sum/average . Check online a fast way to do it
  
  
#corplot 
M <- cor(mtcars)  
plot(M, method = "pie")
corrplot.mixed(M, lower.col = "black", number.cex = .7)

plot(mtcars)
