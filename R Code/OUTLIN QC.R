#############################################################
# PART 2: Quality filtering: Doublet removal and mitochondrial threshold filtering.
#############################################################

## Install and load packages
### Seurat v5.1.0 does not recognize the counts layer in the SCT assay. Install v5.0.3 
remotes::install_version("SeuratObject", "5.0.1", repos = c("https://satijalab.r-universe.dev", getOption("repos")), force = TRUE)
remotes::install_version("Seurat", "5.0.3", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
BiocManager::install('glmGamPoi')

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scDblFinder")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BiocParallel")

library(glmGamPoi)
library(Seurat)
library(SeuratObject)
library(cowplot)
library(ggplot2)
library(BiocManager)
library(tidyverse)
library(dplyr)
library(scDblFinder)
library(BiocParallel)
options(MulticoreParam=MulticoreParam(workers=4))


## Load in filtered feature counts for two OTULIN patients and three healthy controls. 

setwd("/Volumes/Markle Lab/Data/Christine Mariskanish/10604-JM experiment/")
samples <- c('SoupX_0012', 'SoupX_0013', 'SoupX_0014','SoupX_0015','SoupX_0016')
path <- 'OTULIN/SoupX Output'

d10x <- sapply(samples, function(i){
  d10x <- ReadMtx(mtx = file.path(path,i, 'matrix.mtx'), features =  file.path(path,i, "genes.tsv"), cells =  file.path(path,i, "barcodes.tsv"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x), split = "-"), '[[',1L),i,sep = "-")
  d10x})


OTULIN<- do.call("cbind", d10x)

merged <- CreateSeuratObject(
  OTULIN,
  project = "count_10604-JM",
  min.cells = 3,
  min.features = 200,
  names.field = 2
)

table(merged$orig.ident, merged$orig.ident)
    #0012  0013  0014  0015  0016
#0012 20246     0     0     0     0
#0013     0 25622     0     0     0
#0014     0     0 70609     0     0
#0015     0     0     0 57155     0
#0016     0     0     0     0 21666
## Add mitochondrial and ribosomal percentages. Create scatter plots to access quality of cells by sample.

merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")
merged[["percent.ribo"]] <- PercentageFeatureSet(merged, pattern = "^RP[SL]")

FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  NoLegend() + ggtitle('Number of UMIs per cell vs. Number of Genes per cell') + 
  scale_y_log10() + scale_x_log10() + facet_wrap(~merged$orig.ident)


FeatureScatter(merged, feature1 = "percent.mt", feature2 = "percent.ribo") +
  NoLegend() + ggtitle('Percent mito per cell vs. Percent ribo per cell') + 
  scale_y_log10() + scale_x_log10() + facet_wrap(~merged$orig.ident)
# Mitochondrial threshold at 10% would remove a small percentage of each sample.

##Split the object for sample-specific quality filtering
#Create a seurat object list
splitobj <- SplitObject(merged, split.by = "orig.ident")


#Create a list for faster processing with a for loop
ids <- c('patient12', 'patient13', 'control14', 'control15','control16')

#Create new ids
splitobj$patient12 <- splitobj$'0012'
splitobj$patient13 <- splitobj$'0013'
splitobj$control14 <- splitobj$'0014'
splitobj$control15 <- splitobj$'0015'
splitobj$control16 <- splitobj$'0016'

#Remove old IDs
splitobj$'0012' <- NULL
splitobj$'0013' <- NULL
splitobj$'0014' <- NULL
splitobj$'0015' <- NULL
splitobj$'0016' <- NULL
## Remove low quality cells with a threshold cutoff of the lowest 5% of features 
VlnPlot(splitobj$patient12, features = c('nFeature_RNA'))+
  geom_hline(yintercept = 700, linetype = 2)
VlnPlot(splitobj$patient13, features = c('nFeature_RNA'))+
  geom_hline(yintercept = 700, linetype = 2)
VlnPlot(splitobj$control14, features = c('nFeature_RNA', 'nCount_RNA'))
VlnPlot(splitobj$control15, features = c('nFeature_RNA', 'nCount_RNA'))
VlnPlot(splitobj$control16, features = c('nFeature_RNA')) +
  geom_hline(yintercept = 800, linetype = 2)

## 
minCov=1000 #If the minimum nCount is less than 1000, make the cutoff threshold greater than 1% of the quantile distribution. 
#Make a countHIGH variable that has a cutoff lower than 99% of the quantile distribution to remove empty drops, dead cells, or multiplets. 
#Remove the lower quantile of nFeature_RNA which indicate dead or dying cells. 
splitobj_test<- list()
for (i in ids){
  if(min(splitobj[[i]]$nCount_RNA)>=minCov){
    countLOW=min(splitobj[[i]]$nCount_RNA)
  }else{
    countLOW=quantile(splitobj[[i]]$nCount_RNA, prob=c(0.05))  
  }
  countHIGH=quantile(splitobj[[i]]$nCount_RNA, prob=0.95)
  featureLOW=quantile(splitobj[[i]]$nFeature_RNA, prob=0.05)
  
  ##subset
  splitobj_test[[i]] <- subset(splitobj[[i]], subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt <= 10)
  
}


## Remove doublets using scDblFinder

for (i in ids){
  sce <- as.SingleCellExperiment(splitobj_test[[i]])
  
  sce <- scDblFinder(sce, clusters = TRUE, BPPARAM = MulticoreParam(4), samples = "orig.ident")
  
  print(table(sce$orig.ident, sce$scDblFinder.class))
  
  all.equal(Cells(splitobj_test[[i]]), colnames(sce))
  
  splitobj_test[[i]] <- AddMetaData(splitobj_test[[i]], sce$scDblFinder.class, col.name= 'scDblFinder.class')
  
  rm(sce) ; gc()
  
  splitobj_test[[i]] <- subset(splitobj_test[[i]], subset = scDblFinder.class == 'singlet')
}

# Make certain that the doublets identified are removed from the object.
table(splitobj_test$patient12$scDblFinder.class)
table(splitobj_test$patient13$scDblFinder.class)
table(splitobj_test$control14$scDblFinder.class)
table(splitobj_test$control15$scDblFinder.class)
table(splitobj_test$control16$scDblFinder.class)

VlnPlot(splitobj_test$control16, features = c('nFeature_RNA')) +
  geom_hline(yintercept = 700, linetype = 2)

splitobj_test$control16 <- subset(splitobj_test$control16, subset = nFeature_RNA > 700)
## Save RDS object for later analyses and then clear your global environment to free up RAM.
getwd()
setwd('OTULIN/')
saveRDS(splitobj_test,'split_OTULIN.rds')

