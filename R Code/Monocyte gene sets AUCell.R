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
library(viridis)
library(AUCell)
###########
# Determine contribution of cells expressing the top 15 DEGs
##########
setwd('/Volumes/Markle Lab/Data/Christine Mariskanish/10604-JM experiment/OTULIN/')

object <- readRDS('object_postfiltering.rds')
object <- JoinLayers(object, assay = "RNA")

classical <- c('CD14','S100A8','S100A9','S100A12','CD99')
nonclassical <- c('FCGR3A', 'CXCR3','CDKN1C','HES4','CD79B','RHOC')
IFNresponse <- c('IFI44L','IFIT1','IFIT2','MX1','ISG15')
HLAhigh <- c('HLA-DPB1','HLA-DQA1','HLA-DPA1','HLA-DMA')


#Check to see that these gene names are within/detected in the seurat object
intersect(HLAhigh, rownames(object))
intersect(classical, rownames(object))
intersect(nonclassical, rownames(object))
intersect(IFNresponse, rownames(object))

#Extract counts matrix
expreMatrix <- object@assays$RNA$counts

## Run AUCell on the classical, non-classical, HLA-high and IFN-response gene sets
AUcell_output<- AUCell_run(expreMatrix, geneSets = classical)
head(AUcell_output[])
names(AUcell_output@assays@data@listData) <- "Classical monocyte gene set"

object <- AddMetaData(object = object, metadata = AUcell_output@assays@data@listData, col.name = "Classical monocyte gene set")
# non-classical
AUcell_output<- AUCell_run(expreMatrix, geneSets = nonclassical)
head(AUcell_output[])
names(AUcell_output@assays@data@listData) <- "Non-classical monocyte gene set"

object <- AddMetaData(object = object, metadata = AUcell_output@assays@data@listData, col.name = "Non-classical monocyte gene set")

#HLA-high
AUcell_output<- AUCell_run(expreMatrix, geneSets = HLAhigh)
head(AUcell_output[])
names(AUcell_output@assays@data@listData) <- "HLA-high monocyte gene set"

object <- AddMetaData(object = object, metadata = AUcell_output@assays@data@listData, col.name = "HLA-high monocyte gene set")

#IFN-responsive
AUcell_output<- AUCell_run(expreMatrix, geneSets = IFNresponse)
head(AUcell_output[])
names(AUcell_output@assays@data@listData) <- "IFN-responsive monocyte gene set"

object <- AddMetaData(object = object, metadata = AUcell_output@assays@data@listData, col.name = "IFN-responsive monocyte gene set")



metadata <- object@meta.data
setwd("/Volumes/Markle Lab/Data/Christine Mariskanish/10604-JM experiment/OTULIN/Supplementary/AUCell Monocyte Markers 06172025")

write.csv(metadata, 'AUCell monocyte markers metadata.csv')

# Plot density UMAPs
Mono <- subset(object, subset = predicted.celltype.l1 =='Mono')
axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(3, "cm")
)
setwd('/Volumes/Markle Lab/Data/Christine Mariskanish/10604-JM experiment/OTULIN/Figures/AUCell Monocyte Markers Density Feature UMAPs/Total L1 cells')
png("Classical Monocytes Feature Density.png", width = 11, height = 7, units = 'in', res = 300)
Nebulosa::plot_density(object, features = c("Classical monocyte gene set"),size = 2) + 
  scale_colour_gradientn(colours = c("grey",viridis(1),"#FDE725FF"),limits = c(0.0,0.2), breaks = seq(0.0, 0.2, 0.05), na.value = "grey")+
  ggtitle("Classical monocyte gene set")+guides(x = axis, y = axis)+ xlab('UMAP 1')+ylab('UMAP 2')+
  #coord_equal(ylim =c(-3.936416,  6.919397), xlim=c(-11.501336,   -2))+
  coord_equal(ylim =c(-16.61637,15.09251), xlim=c(-15.39274,14.75713))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        axis.line = element_line(arrow = arrow(type='closed', length = unit(3,'pt'))),
        axis.title = element_text(hjust = 0, size = 12, family = 'arial'),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks = element_blank())+   theme(aspect.ratio = 1)
dev.off()

# non-classical
png("Non-classical Monocytes Feature Density.png", width = 11, height = 7, units = 'in', res = 300)
Nebulosa::plot_density(object, features = c("Non-classical monocyte gene set"),size = 2) + 
  scale_colour_gradientn(colours = c("grey",viridis(1),"#FDE725FF"),limits = c(0.0,0.2), breaks = seq(0.0, 0.2, 0.05), na.value = "grey")+
  ggtitle("Non-classical monocyte gene set")+guides(x = axis, y = axis)+ xlab('UMAP 1')+ylab('UMAP 2')+
  #coord_equal(ylim =c(-3.936416,  6.919397), xlim=c(-11.501336,   -2))+
  coord_equal(ylim =c(-16.61637,15.09251), xlim=c(-15.39274,14.75713))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        axis.line = element_line(arrow = arrow(type='closed', length = unit(3,'pt'))),
        axis.title = element_text(hjust = 0, size = 12, family = 'arial'),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks = element_blank())+   theme(aspect.ratio = 1)
dev.off()


# IFN response
png("IFN-responsive Monocytes Feature Density .png", width = 11, height = 7, units = 'in', res = 300)
Nebulosa::plot_density(object, features = c("IFN-responsive monocyte gene set"),size = 2) + 
  scale_colour_gradientn(colours = c("grey",viridis(1),"#FDE725FF"),limits = c(0.0,0.5), breaks = seq(0.0, 0.5, 0.1), na.value = "grey")+
  ggtitle("IFN-responsive monocyte gene set")+guides(x = axis, y = axis)+ xlab('UMAP 1')+ylab('UMAP 2')+
  #coord_equal(ylim =c(-3.936416,  6.919397), xlim=c(-11.501336,   -2))+
  coord_equal(ylim =c(-16.61637,15.09251), xlim=c(-15.39274,14.75713))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        axis.line = element_line(arrow = arrow(type='closed', length = unit(3,'pt'))),
        axis.title = element_text(hjust = 0, size = 12, family = 'arial'),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks = element_blank())+   theme(aspect.ratio = 1)
dev.off()

# HLA-high
png("HLA-high Monocytes Feature Density .png", width = 11, height = 7, units = 'in', res = 300)
Nebulosa::plot_density(object, features = c("HLA-high monocyte gene set"),size = 2) + 
  scale_colour_gradientn(colours = c("grey",viridis(1),"#FDE725FF"),limits = c(0.0,0.2), breaks = seq(0.0, 0.2, 0.05), na.value = "grey")+
  ggtitle("HLA-high monocyte gene set")+guides(x = axis, y = axis)+ xlab('UMAP 1')+ylab('UMAP 2')+
  #coord_equal(ylim =c(-3.936416,  6.919397), xlim=c(-11.501336,   -2))+
  coord_equal(ylim =c(-16.61637,15.09251), xlim=c(-15.39274,14.75713))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        axis.line = element_line(arrow = arrow(type='closed', length = unit(3,'pt'))),
        axis.title = element_text(hjust = 0, size = 12, family = 'arial'),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks = element_blank())+   theme(aspect.ratio = 1)
dev.off()

saveRDS(object, '/Volumes/Markle Lab/Data/Christine Mariskanish/10604-JM experiment/OTULIN/AUCell object monocyte markers.rds')