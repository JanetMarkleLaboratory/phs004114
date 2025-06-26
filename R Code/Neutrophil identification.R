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
library(WriteXLS)
options(MulticoreParam=MulticoreParam(workers=4))
mem.maxNSize(1e10)

## Load in filtered feature counts for two OTULIN patients and three healthy controls. 

setwd("/Volumes/Markle Lab/Data/Christine Mariskanish/10604-JM experiment/OTULIN")
samples <- c('count_10406-JM-0012', 'count_10406-JM-0013', 'count_10406-JM-0014','count_10406-JM-0015','count_10406-JM-0016')


d10x <- sapply(samples, function(i){
  d10x <- Read10X_h5(file.path(i, 'filtered_feature_bc_matrix.h5'))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x), split = "-"), '[[',1L),i,sep = "-")
  d10x})


OTULIN <- do.call("cbind", d10x)

merged <- CreateSeuratObject(
  OTULIN,
  project = "count_10604-JM",
  min.cells = 3,
  min.features = 200,
  names.field = 2,
  names.delim = "\\-JM-"
)

table(merged$orig.ident, merged$orig.ident)

###### 0012  0013  0014  0015  0016
#0012 20247     0     0     0     0
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

##Split the object for sample-specific quality filtering

#Create a seurat object list
splitobj <- SplitObject(merged, split.by = "orig.ident")


#Create a list for faster processing with a for loop
ids <- c('patient12', 'patient13', 'control14', 'control15', 'control16')

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
VlnPlot(splitobj$patient12, features = c('nFeature_RNA', 'nCount_RNA'))
VlnPlot(splitobj$patient13, features = c('nFeature_RNA', 'nCount_RNA'))
VlnPlot(splitobj$control14, features = c('nFeature_RNA', 'nCount_RNA'))
VlnPlot(splitobj$control15, features = c('nFeature_RNA', 'nCount_RNA'))
VlnPlot(splitobj$control16, features = c('nFeature_RNA', 'nCount_RNA'))

for (i in ids){
  
  ## filter out low quality cells similar to methods used by Zhang et al. 2023.
splitobj[[i]] <- subset(splitobj[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 10000  & percent.mt <= 20)}

VlnPlot(splitobj$patient12, features = c('nFeature_RNA', 'nCount_RNA'))
VlnPlot(splitobj$patient13, features = c('nFeature_RNA', 'nCount_RNA'))
VlnPlot(splitobj$control14, features = c('nFeature_RNA', 'nCount_RNA'))
VlnPlot(splitobj$control15, features = c('nFeature_RNA', 'nCount_RNA'))
VlnPlot(splitobj$control16, features = c('nFeature_RNA', 'nCount_RNA'))

## Save RDS object for later analyses and then clear your global environment to free up RAM.
getwd()

saveRDS(splitobj,'split_OTULIN_neutrophils.rds')
splitobj<- readRDS('split_OTULIN_neutrophils.rds')
################################
## Differential Expression Analysis: Aggregate Expression to identify whether neutrophils are present in the patient PBMC population. 
## We know that neutrophils are not present in the healthy donor population.
## If we identify differentially expressed neutrophil markers that are higher in patients, then neutrophils are likely to be present.
################################
## min.cells.group = 2
merge_object <- merge(splitobj[[1]], splitobj[2:length(splitobj)])
merge_object$new.ident <- NA
merge_object@meta.data$new.ident[which(merge_object@meta.data$orig.ident == "0012")] <- "R57C"
merge_object@meta.data$new.ident[which(merge_object@meta.data$orig.ident == "0013")] <- "R57C"
merge_object@meta.data$new.ident[which(merge_object@meta.data$orig.ident == "0014")] <- "HC"
merge_object@meta.data$new.ident[which(merge_object@meta.data$orig.ident == "0015")] <- "HC"
merge_object@meta.data$new.ident[which(merge_object@meta.data$orig.ident == "0016")] <- "HC"

pseudo_otulin <- AggregateExpression(merge_object, assays = "RNA", return.seurat = T, group.by = c("new.ident", "orig.ident"))
Cells(pseudo_otulin)


Idents(pseudo_otulin) <- "new.ident"

bulk.de.selected <- FindMarkers(object = pseudo_otulin, min.cells.group = 2, ident.1 = "R57C",ident.2 = "HC", 
                            test.use = "DESeq2")
bulk.de.selected %>%
  dplyr::filter(p_val_adj < 0.05)%>%
  ungroup() -> bulk.de.sig

#### Results from this analysis suggest that neutrophils are not significantly present in patients. 
### Marker genes CSF3R, S100A9 and FCGR3B were used to detect neutrophils.
### According to the Human Protein Atlas (HPA), CSF3R is the most reliable marker of the three genes for neutrophils.
write.csv(bulk.de.selected, 'bulk.de.csv')


################################
## Aggregate Expression to compare DEGs found from NanoString. 
## We want values 1) by genotype and 2) by sample.
################################
list_of_genes <- c('BST2','HLA-A','HLA-B',
                   
                   'HLA-C','IFI35','IFIT2',
                   
                   'IFITM1',#'IFNA1/13','IFNA2',
                   
                   'IFNAR1','IFNAR2','IFNB1',
                   
                   'IRF1','IRF3','IRF4',
                   
                   'IRF5','IRF7','IRF8',
                   
                   'JAK1','MX1','PSMB8',
                   
                   'PTPN6','SOCS1','SOCS3',
                   
                   'STAT1','STAT2','TYK2')

#Check to see if any synonyms to the genes unidentifiedc through scRNA-seq exist within this dataset.
rownames(merge_object)[grepl('IFNA',rownames(merge_object))]
# "IFNA5"  "IFNAR2" "IFNAR1"

pseudo_otulin <- AggregateExpression(merge_object, assays = "RNA", return.seurat = T, group.by = c("new.ident", "orig.ident"), features = list_of_genes)
#The following 2 features were not found in the RNA assay: IFNA1/13, IFNA2
# I attempted to identify these genes, by placing different aliases (IFNA1, IFNA2B, IFNA13) but these transcripts were not amplified within this assay.

Cells(pseudo_otulin)
#"HC_0014"   "HC_0015"   "HC_0016"   "R57C_0012" "R57C_0013"


data <- pseudo_otulin@assays$RNA$scale.data
data <- as.data.frame(data)
data$genes <- rownames(data)

# Order genes for ease of data visualization
data <- data[order(data$genes, list_of_genes),]

## By sample
setwd("/Volumes/pmi/Data/Research/Markle Lab/Data/Christine Mariskanish/10604-JM experiment/OTULIN/Supplementary")
WriteXLS(data, 'scRNAseq z-scaled data of gene list by sample.xlsx')






