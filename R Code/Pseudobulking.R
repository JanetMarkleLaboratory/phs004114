#############################################################
# PART 4: DEA: pseudobulking
#############################################################
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

# We are interested in discovering DEGs between patients vs. HC monocytes. 
setwd('/Volumes/Markle Lab/Data/Christine Mariskanish/10604-JM experiment/OTULIN/')

object <- readRDS('object_postfiltering.rds')
################################
## Differential Expression Analysis: Aggregate Expression
################################
## min.cells.group = 2

pseudo_otulin <- AggregateExpression(object, assays = "RNA", return.seurat = T, group.by = c("new.ident", "newer.ident","orig.ident"))
Cells(pseudo_otulin)


Idents(pseudo_otulin) <- "newer.ident"

bulk.de <- FindMarkers(object = pseudo_otulin, min.cells.group = 2, ident.1 = "R57C",ident.2 = "HC", 
                                test.use = "DESeq2")
bulk.de %>%
  dplyr::filter(p_val_adj < 0.05)%>%
  ungroup() -> bulk.sig

write.csv(bulk.de, 'bulk.de.ALL.csv')

Idents(object) <- 'predicted.celltype.l1'
monocytes <- subset(object, idents ='Mono')
pseudo_otulin <- AggregateExpression(monocytes, assays = "RNA", return.seurat = T, group.by = c("newer.ident","new.ident", "predicted.celltype.l1"))
Cells(pseudo_otulin)

Idents(pseudo_otulin) <- "new.ident"
bulk.de.selected <- FindMarkers(object = pseudo_otulin, min.cells.group = 2, ident.1 = "R57C", ident.2 = "HC",
                                   test.use = "DESeq2")
#Top 15 genes
bulk.de.selected %>%
  dplyr::filter(p_val_adj < 0.05)%>%
  dplyr::filter(avg_log2FC>0)%>%
  arrange(-avg_log2FC)%>%
  slice_head(n = 15)-> bulk.de.sig
#285 significantly DEGs found from DEA

#### CXCL8 was found upregulated in patient monocytes. Create a heatmap of significant DEGs for monocytes.
write.csv(bulk.de.selected, 'bulk.de.monocytes.csv')
DefaultAssay(monocytes) <- "RNA"

monocytes <- NormalizeData(monocytes, assay = "RNA")
monocytes <- FindVariableFeatures(monocytes, assay = "RNA")
monocytes <- ScaleData(monocytes,  assay = "RNA")
monocytes <- RunPCA(monocytes, verbose = FALSE, assay = "RNA")

DimPlot(monocytes, reduction = "pca", split.by = 'new.ident')

FeaturePlot(monocytes, features = 'percent.mt')

ElbowPlot(monocytes, ndims = 50)
monocytes <- RunUMAP(monocytes, dims = 1:30, verbose = FALSE, assay = "RNA")
monocytes <- FindNeighbors(monocytes, dims = 1:30, verbose = FALSE, assay = "RNA")
monocytes <- FindClusters(monocytes, verbose = FALSE, resolution =seq(0, 2, by=0.05))

# Create a heatmap
meta_data_monocytes <- monocytes@meta.data
meta_data_monocytes <- meta_data_monocytes[,c("predicted.celltype.l1","new.ident","newer.ident")]
meta_data_monocytes$Sample<- meta_data_monocytes$new.ident
meta_data_monocytes$new.ident<- NULL
meta_data_monocytes$Celltype<- meta_data_monocytes$predicted.celltype.l1
meta_data_monocytes$predicted.celltype.l1<- NULL
meta_data_monocytes$Genotype <- meta_data_monocytes$newer.ident
meta_data_monocytes$newer.ident <-NULL
# order of annotations/colors are defined here
ordered_meta_data_monocytes <- meta_data_monocytes[order(meta_data_monocytes$Sample), ]
ordered_meta_data_monocytes$celltype.genotype.l1<- NULL
ordered_meta_data_monocytes$Celltype <- NULL

monocytes <- JoinLayers(monocytes, "RNA")
my_data_monocytes <- monocytes@assays$RNA$data
my_data_monocytes <- my_data_monocytes[rownames(my_data_monocytes)%in% rownames(bulk.de.sig),]
my_data_monocytes <- t(scale(t(my_data_monocytes)))
# Define sample groups for column splitting
sample_groups <- meta_data_monocytes$Sample
sample_order <- factor(sample_groups, levels = c("HC 1", "HC 2","HC 3", "R57C 1","R57C 2"))
col_fun = colorRamp2(c(-2,0, 2),  viridis(3))

genotypes <- meta_data_monocytes$Genotype
genotype_colors <- c("HC" = "blue","HC" = "blue","HC" = "blue",
              "R57C" = "red", "R57C"= "red")
# Reorder columns based on sample order
ordered_indices <- order(factor(sample_groups, levels = c("HC 1", "HC 2","HC 3", "R57C 1","R57C 2")))
ordered_matrix <- my_data_monocytes[, ordered_indices]

# Create text annotation for samples
text_annot <- HeatmapAnnotation(Sample = anno_block(gp = gpar(fill = NA),
                                                    labels = c("HC 1", "HC 2","HC 3", "R57C","R57C"), show_name = T,
                                                    labels_gp = gpar(col = "black", fontsize = 10, family = 'arial')))

# Create a block annotation for genotypes
genotype_annotation <- data.frame(Genotype = genotypes[ordered_indices])
rownames(genotype_annotation) <- colnames(ordered_matrix)

block_annot <- HeatmapAnnotation(df = genotype_annotation,
                                 col = list(Genotype = genotype_colors),
                                 which = "col",
                                 annotation_name_side = "right")

png("15_DEG_Monocytes xcell.png",width=7,height=7,units="in",res=300)
# Create the heatmap with combined annotations
Heatmap(ordered_matrix,show_column_names = F, column_title = 'Top 15 DEGs between patient and HC monocytes', show_row_dend = F, col = col_fun,
        name= 'Expression',row_order = rownames(bulk.de.sig),
        top_annotation = c(block_annot, text_annot),  # Combine block and text annotations
        cluster_columns = FALSE,  # Maintain specified order
        column_split = sample_groups[ordered_indices],  # Split columns by sample
        show_row_names = T)
dev.off()
########################################
### Create a heatmap for pseudobulked monocytes
########################################
# Create a heatmap
meta_data_monocytes <- pseudo_otulin@meta.data
meta_data_monocytes <- meta_data_monocytes[,c("new.ident","newer.ident")]
meta_data_monocytes$Sample<- meta_data_monocytes$newer.ident
meta_data_monocytes$newer.ident<- NULL
meta_data_monocytes$Genotype <- meta_data_monocytes$new.ident
meta_data_monocytes$new.ident <-NULL
# order of annotations/colors are defined here
ordered_meta_data_monocytes <- meta_data_monocytes[order(meta_data_monocytes$Sample), ]

my_data_monocytes <- pseudo_otulin@assays$RNA$data
my_data_monocytes <- my_data_monocytes[rownames(my_data_monocytes)%in% rownames(bulk.de.sig),]
my_data_monocytes <- t(scale(t(my_data_monocytes)))
# Define sample groups for column splitting
sample_groups <- meta_data_monocytes$Sample
sample_order <- factor(sample_groups, levels = c("HC 1", "HC 2","HC 3", "R57C 1","R57C 2"))


genotypes <- meta_data_monocytes$Genotype
genotype_colors <- c("HC" = "blue","HC" = "blue","HC" = "blue",
                     "R57C" = "red", "R57C"= "red")
# Reorder columns based on sample order
ordered_indices <- order(factor(sample_groups, levels = c("HC 1", "HC 2","HC 3", "R57C 1","R57C 2")))
ordered_matrix <- my_data_monocytes[, ordered_indices]
colnames(ordered_matrix) <- c("HC 1", "HC 2","HC 3", "R57C 1","R57C 2")
# Create text annotation for samples
text_annot <- HeatmapAnnotation(Sample = anno_block(gp = gpar(fill = NA),
                                                    labels = c("HC 1", "HC 2","HC 3", "R57C 1","R57C 2"), show_name = T,
                                                    labels_gp = gpar(col = "black", fontsize = 10, family = 'arial')))
col_fun = colorRamp2(c(-1.5,0, 1.5),  viridis(3))
# Create a block annotation for genotypes
genotype_annotation <- data.frame(Genotype = genotypes[ordered_indices])
rownames(genotype_annotation) <- colnames(ordered_matrix)

block_annot <- HeatmapAnnotation(df = genotype_annotation,
                                 col = list(Genotype = genotype_colors),
                                 which = "col",
                                 annotation_name_side = "right")

png("15_DEG_Monocytes.png",width=7,height=7,units="in",res=300)
# Create the heatmap with combined annotations
Heatmap(ordered_matrix,show_column_names = F, column_title = 'Top 15 DEGs between patient and HC monocytes', show_row_dend = F, col = col_fun,
        name= 'Expression', clustering_method_rows = "ward.D2", row_order = rownames(bulk.de.sig),
        top_annotation = c(block_annot, text_annot),  # Combine block and text annotations
        cluster_columns = FALSE,  # Maintain specified order
        column_split = sample_groups[ordered_indices],  # Split columns by sample
        show_row_names = T)
dev.off()



