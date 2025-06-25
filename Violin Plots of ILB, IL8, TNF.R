#############################################################
# Create violin plots for 3 genes
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
mem.maxVSize(vsize = 100000)
library(readxl)

# We are interested in discovering DEGs between patients vs. HC monocytes. Plot the top 15 ranked genes from pseudobulked data
## to see where (in level 1 identities) the up-regulation originates

DEA_across_all_cells_ranked_by_FC <- read_excel("/Volumes/Markle Lab/Data/Christine Mariskanish/10604-JM experiment/OTULIN/Supplementary/DEA across all cells ranked by FC.xlsx")
View(DEA_across_all_cells_ranked_by_FC)

setwd('/Volumes/Markle Lab/Data/Christine Mariskanish/10604-JM experiment/OTULIN/')

object <- readRDS('object_postfiltering.rds')

colors <- c(
  "HC" = "blue",
  "R57C" = "red"
)
object <- JoinLayers(object,assay = "RNA")

log2_expression_data <- object@assays$RNA$data
meta.data <- object@meta.data

rm(object)
gc() 
genes_fit <- DEA_across_all_cells_ranked_by_FC %>%
  dplyr::slice_head(n =15)

subset_log2_expression_data <- log2_expression_data[rownames(log2_expression_data) %in% genes_fit$Gene,]

subset_log2_expression_data <- as.data.frame(as.table(as.matrix(subset_log2_expression_data)))
colnames(subset_log2_expression_data) <- c("gene", "cell", "Expression")


meta.data$cell <- rownames(meta.data)
length(unique(subset_log2_expression_data$cell))
length(unique(rownames(meta.data)))

merge_log2_expression_data<- merge(subset_log2_expression_data, meta.data, by = "cell")


merge_log2_expression_data$predicted.celltype <- NA
merge_log2_expression_data$predicted.celltype[which(merge_log2_expression_data$predicted.celltype.l1 == "B")] <- "B cells"
merge_log2_expression_data$predicted.celltype[which(merge_log2_expression_data$predicted.celltype.l1 == "NK")] <- "NK cells"
merge_log2_expression_data$predicted.celltype[which(merge_log2_expression_data$predicted.celltype.l1 == "CD4 T")] <- "CD4+ T cells"
merge_log2_expression_data$predicted.celltype[which(merge_log2_expression_data$predicted.celltype.l1 == "CD8 T")] <- "CD8+ T cells"
merge_log2_expression_data$predicted.celltype[which(merge_log2_expression_data$predicted.celltype.l1 == "Mono")] <- "Monocytes"
merge_log2_expression_data$predicted.celltype[which(merge_log2_expression_data$predicted.celltype.l1 == "DC")] <- "Dendritic cells"
merge_log2_expression_data$predicted.celltype[which(merge_log2_expression_data$predicted.celltype.l1 == "other")] <- "Others"
merge_log2_expression_data$predicted.celltype[which(merge_log2_expression_data$predicted.celltype.l1 == "other T")] <- "Other T cells"
table(merge_log2_expression_data$predicted.celltype, merge_log2_expression_data$predicted.celltype.l1)


for (i in genes_fit$Gene){
  name <- paste("merge_log2_expression_data", i, sep = '')
  assign(x =name, value = merge_log2_expression_data[merge_log2_expression_data$gene == i,])
}

scaleFUN <- function(x) sprintf("%.2f", x)

setwd('Figures')
png("IL1B expression by L1 violin plots.png", width = 11, height = 7, units = 'in', res = 300)
ggplot(merge_log2_expression_dataIL1B, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = 'Expression')+ xlab('')+
  facet_wrap(~predicted.celltype, ncol = 4) +
  coord_cartesian(ylim = c(0, 6))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.x  = element_text(size = 12, color = "black", family = 'arial', angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 10, color = "black", family = 'arial'),
        axis.title.y = element_text(size = 15, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
dev.off()

png("CXCL8 expression by L1 violin plots.png", width = 11, height = 7, units = 'in', res = 300)
ggplot(merge_log2_expression_dataCXCL8, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+facet_wrap(~predicted.celltype, ncol = 4) +
  scale_fill_manual(values = colors)+ labs(x = ' ',y = 'Expression')+ xlab('')+
  facet_wrap(~predicted.celltype, ncol = 4) +
  coord_cartesian(ylim = c(0, 6))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.x  = element_text(size = 12, color = "black", family = 'arial', angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 10, color = "black", family = 'arial'),
        axis.title.y = element_text(size = 15, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
dev.off()

png("TNF expression by L1 violin plots.png", width = 11, height = 7, units = 'in', res = 300)
ggplot(merge_log2_expression_dataTNF, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+facet_wrap(~predicted.celltype, ncol = 4) +
  scale_fill_manual(values = colors)+ labs(x = ' ',y = 'Expression')+ xlab('')+
  facet_wrap(~predicted.celltype, ncol = 4) +
  coord_cartesian(ylim = c(0, 6))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.x  = element_text(size = 12, color = "black", family = 'arial', angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 10, color = "black", family = 'arial'),
        axis.title.y = element_text(size = 15, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
dev.off()


png("IL1B expression by L1 violin plots.png", width = 11, height = 7, units = 'in', res = 300)
ggplot(merge_log2_expression_dataIL1B, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = 'Expression')+ xlab('')+
  facet_wrap(~predicted.celltype, ncol = 4) +
  coord_cartesian(ylim = c(0, 6))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.x  = element_text(size = 12, color = "black", family = 'arial', angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 10, color = "black", family = 'arial'),
        axis.title.y = element_text(size = 15, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
dev.off()