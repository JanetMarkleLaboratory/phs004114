#############################################################
# Create violin plots for 23 genes that overlap NanoString
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

setwd('/Volumes/Markle Lab/Data/Christine Mariskanish/10604-JM experiment/OTULIN/Supplementary')
data <- readxl::read_excel('scRNAseq z-scaled data of gene list by sample.xlsx')

object <- readRDS('/Volumes/Markle Lab/Data/Christine Mariskanish/10604-JM experiment/OTULIN/object_postfiltering.rds')

colors <- c(
  "HC" = "blue",
  "R57C" = "red"
)
object <- JoinLayers(object,assay = "RNA")

log2_expression_data <- object@assays$RNA$data
meta.data <- object@meta.data

rm(object)
gc() 

genes <- data$genes

subset_log2_expression_data <- log2_expression_data[rownames(log2_expression_data) %in% genes,]

subset_log2_expression_data <- as.data.frame(as.table(as.matrix(subset_log2_expression_data)))
colnames(subset_log2_expression_data) <- c("gene", "cell", "Expression")


meta.data$cell <- rownames(meta.data)
length(unique(subset_log2_expression_data$cell))
length(unique(rownames(meta.data)))

merge_log2_expression_data<- merge(subset_log2_expression_data, meta.data, by = "cell")



for (i in genes){
  name <- paste("merge_log2_expression_data", i, sep = '')
  assign(x =name, value = merge_log2_expression_data[merge_log2_expression_data$gene == i,])
}

scaleFUN <- function(x) sprintf("%.2f", x)

setwd('/Volumes/Markle Lab/Data/Christine Mariskanish/10604-JM experiment/OTULIN/Figures/25 genes overlapping with NanoString')
png("BST2 expression by sample violin plots.png", width = 11, height = 7, units = 'in', res = 300)
ggplot(merge_log2_expression_dataBST2, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = 'Expression')+ xlab('')+
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

png("HLA-A expression by sample violin plots.png", width = 11, height = 7, units = 'in', res = 300)
ggplot(`merge_log2_expression_dataHLA-A`, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = 'Expression')+ xlab('')+
  
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

png("HLA-B expression by sample violin plots.png", width = 11, height = 7, units = 'in', res = 300)
ggplot(`merge_log2_expression_dataHLA-B`, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = 'Expression')+ xlab('')+
  
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


png("HLA-C expression by sample violin plots.png", width = 11, height = 7, units = 'in', res = 300)
ggplot(`merge_log2_expression_dataHLA-C`, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = 'Expression')+ xlab('')+
  
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

png("IFI35 expression by sample violin plots.png", width = 11, height = 7, units = 'in', res = 300)
ggplot(merge_log2_expression_dataIFI35, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = 'Expression')+ xlab('')+
  
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

png("IFIT2 expression by sample violin plots.png", width = 11, height = 7, units = 'in', res = 300)
ggplot(merge_log2_expression_dataIFIT2, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = 'Expression')+ xlab('')+
  
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

png("IFITM1 expression by sample violin plots.png", width = 11, height = 7, units = 'in', res = 300)
ggplot(merge_log2_expression_dataIFITM1, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = 'Expression')+ xlab('')+
  
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

png("IFNAR1 expression by sample violin plots.png", width = 11, height = 7, units = 'in', res = 300)
ggplot(merge_log2_expression_dataIFNAR1, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = 'Expression')+ 
  xlab('')+
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


png("IFNAR2 expression by sample violin plots.png", width = 11, height = 7, units = 'in', res = 300)
ggplot(merge_log2_expression_dataIFNAR2, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = 'Expression')+ xlab('')+
  
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

png("IFNB1 expression by sample violin plots.png", width = 11, height = 7, units = 'in', res = 300)
ggplot(merge_log2_expression_dataIFNB1, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = 'Expression')+ xlab('')+
  
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

png("IRF1 expression by sample violin plots.png", width = 11, height = 7, units = 'in', res = 300)
ggplot(merge_log2_expression_dataIRF1, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = 'Expression')+ xlab('')+
  
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


png("IRF3 expression by sample violin plots.png", width = 11, height = 7, units = 'in', res = 300)
ggplot(merge_log2_expression_dataIRF3, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = 'Expression')+ xlab('')+
  
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

png("IRF4 expression by sample violin plots.png", width = 11, height = 7, units = 'in', res = 300)
ggplot(merge_log2_expression_dataIRF4, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = 'Expression')+ xlab('')+
  
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


png("IRF5 expression by sample violin plots.png", width = 11, height = 7, units = 'in', res = 300)
ggplot(merge_log2_expression_dataIRF5, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = 'Expression')+ xlab('')+
  
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

png("IRF7 expression by sample violin plots.png", width = 11, height = 7, units = 'in', res = 300)
ggplot(merge_log2_expression_dataIRF7, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = 'Expression')+ xlab('')+
  
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

png("IRF8 expression by sample violin plots.png", width = 11, height = 7, units = 'in', res = 300)
ggplot(merge_log2_expression_dataIRF8, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = 'Expression')+ xlab('')+
  
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


png("TYK2 expression by sample violin plots.png", width = 11, height = 7, units = 'in', res = 300)
ggplot(merge_log2_expression_dataTYK2, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = 'Expression')+ xlab('')+
  
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


