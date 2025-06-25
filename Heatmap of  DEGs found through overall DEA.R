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
library(ComplexHeatmap)
library(colorRamp2)
options(MulticoreParam=MulticoreParam(workers=4))
mem.maxVSize(vsize = 100000)
library(readxl)

###########
# Determine contribution of cells expressing the top 15 DEGs
##########
setwd('/Volumes/Markle Lab/Data/Christine Mariskanish/10604-JM experiment/OTULIN/')

object <- readRDS('object_postfiltering.rds')
object <- JoinLayers(object, assay = "RNA")
DEA_across_all_cells_ranked_by_FC <- read.csv("/Volumes/Markle Lab/Data/Christine Mariskanish/10604-JM experiment/OTULIN/Supplementary/bulk.de.ALL.csv")


pseudo_object <- AggregateExpression(object, assays = "RNA", return.seurat = T, group.by = c(""))
Cells(pseudo_object)

DEA_across_all_cells_ranked_by_FC %>%
  dplyr::filter(p_val_adj < 0.05) -> Sig_DEGs #985 DEGs


data <- pseudo_object@assays$RNA$data
data <- data[rownames(data) %in% Sig_DEGs$X,]
scaled <- t(scale(t(data)))
scaled


BlWhRd <-viridis(3)


col_ann <- data.frame(Genotype = c(rep('HD',3), rep('R57C/R57C', 2)),
                      Sample = c('HD 1','HD 2','HD 3','R57C/R57C 1', 'R57C/R57C 2'))
rownames(col_ann) <- colnames(scaled)

ann_color = list(
  'Genotype' = c("HD"="blue", "R57C/R57C"="red"),
  'Sample' = c('HD 1'="blue",'HD 2'="blue",'HD 3'="blue",'R57C/R57C 1'="red", 'R57C/R57C 2'="red"))



library(viridis)

# Create text annotation for samples
text_annot <- HeatmapAnnotation(Sample = anno_block(gp = gpar(fill = NA),
                                                    labels = c("HD 1", "HD 2","HD 3", "R57C/R57C 1","R57C/R57C 2"), show_name = T,
                                                    labels_gp = gpar(col = "black", fontsize = 10, family = 'arial')))



# Create a block annotation for genotypes
genotype_annotation <- data.frame(Genotype = c(rep('HD',3), rep("R57C/R57C",2)))
rownames(genotype_annotation) <- colnames(scaled)
genotype_colors <- c("HD" = "blue","HD" = "blue","HD" = "blue",
                     "R57C/R57C" = "red", "R57C/R57C"= "red")

block_annot <- HeatmapAnnotation(df = genotype_annotation,
                                 col = list(Genotype = genotype_colors),
                                 which = "col",
                                 annotation_name_side = "right")
col_fun = colorRamp2(c(-1,0, 1),  viridis(3))

# Create the heatmap with combined annotations
Heatmap(scaled,show_column_names = F, column_title = 'Total DEGs between Patient and HD cells', 
        show_row_dend = F, 
        col = col_fun,
        name= 'z-score', cluster_rows = T,
        top_annotation = c(block_annot, text_annot),  # Combine block and text annotations
        cluster_columns = FALSE,  # Maintain specified order
        column_split = colnames(scaled),  # Split columns by sample
        show_row_names = F)

