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
DEA_across_all_cells_ranked_by_FC <- read_excel("/Volumes/Markle Lab/Data/Christine Mariskanish/10604-JM experiment/OTULIN/Supplementary/DEA across all cells ranked by FC.xlsx")

colors <- c(
  "HC" = "blue",
  "R57C" = "red"
)

genes_fit <- DEA_across_all_cells_ranked_by_FC %>%
  dplyr::slice_head(n =370)

object_small<- AddModuleScore(object, features = list(genes_fit$Gene), name = "module.score")
head(object_small[])

axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(3, "cm")
)

setwd('/Volumes/Markle Lab/Data/Christine Mariskanish/10604-JM experiment/OTULIN/Unofficial figures')
png("Gene Program Module Score 50 Genes viridis.png", width = 11, height = 7, units = 'in', res = 300)
FeaturePlot(object_small, features = "module.score1", label = F, repel = F,pt.size = 1, raster = T) +
  scale_colour_gradientn(colours = viridis(3), limits = c(-0.2,0.4), breaks = seq(-0.20, 0.4, 0.2), na.value = "grey", name = "Module Score")+
  guides(x = axis, y = axis)+ xlab('UMAP 1')+ylab('UMAP 2')+
  coord_equal(ylim =c(-16.37798,14.92378), xlim=c(-14.94919,14.49963))+
  theme(plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        axis.line = element_line(arrow = arrow(type='closed', length = unit(3,'pt'))),
        axis.title = element_text(hjust = 0, size = 12, family = 'arial'),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks = element_blank())
dev.off()

Patients<- subset(object_small, subset= newer.ident == "R57C")

png("R57C Gene Program Module Score 50 DEGs viridis.png", width = 11, height = 7, units = 'in', res = 300)
FeaturePlot(Patients, features = "module.score1", label = F, repel = F,pt.size = 1, raster = T) +
  scale_colour_gradientn(colours = viridis(3), limits = c(-0.2,0.4), breaks = seq(-0.20, 0.4, 0.2), na.value = "grey", name = "Module Score")+
  ggtitle("R57C Module Scores")+
  guides(x = axis, y = axis)+ xlab('UMAP 1')+ylab('UMAP 2')+
  coord_equal(ylim =c(-16.37798,14.92378), xlim=c(-14.94919,14.49963))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        axis.line = element_line(arrow = arrow(type='closed', length = unit(3,'pt'))),
        axis.title = element_text(hjust = 0, size = 12, family = 'arial'),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks = element_blank())
dev.off()

HCs <- subset(object_small, subset= newer.ident == "HC")
png("HC Gene Program Module Score 50 DEGs viridis.png", width = 11, height = 7, units = 'in', res = 300)
FeaturePlot(HCs, features = "module.score1", label = F, repel = F,pt.size = 1, raster = T) +
  scale_colour_gradientn(colours = viridis(3), limits = c(-0.2,0.4), breaks = seq(-0.20, 0.4, 0.2), na.value = "grey", name = "Module Score")+
  ggtitle("HC Module Scores")+
  guides(x = axis, y = axis)+ xlab('UMAP 1')+ylab('UMAP 2')+
  coord_equal(ylim =c(-16.37798,14.92378), xlim=c(-14.94919,14.49963))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        axis.line = element_line(arrow = arrow(type='closed', length = unit(3,'pt'))),
        axis.title = element_text(hjust = 0, size = 12, family = 'arial'),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks = element_blank())
dev.off()

meta.data370 <- object_small@meta.data
meta.data50 <- object_small@meta.data
meta.data15 <- object_small@meta.data
png("Cell type level 1 50 DEGs.png", width = 11, height = 7, units = 'in', res = 300)
ggplot(meta.data50, aes(y = module.score1, x = newer.ident, fill = factor(newer.ident)))+
  geom_violin()+ facet_grid(~predicted.celltype.l1)+ xlab('')+
  ylab('Module Score')+scale_fill_manual(values = colors)+
  labs(fill='Genotype') +
  theme(
    strip.text = element_text(size = 12, color = "black", family = 'arial'),axis.text.x  = element_blank(),
axis.text.y = element_text(size = 10, color = "black", family = 'arial'),
axis.title.y = element_text(size = 15, color = "black", family = 'arial'),
axis.ticks.x = element_blank(),
panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
geom = "crossbar",
width = 0.5,
#fatten = ,
color = "black")+ggtitle('50 DEGs AddModuleScore')
dev.off()
png("Cell type level 1 15 DEGs.png", width = 11, height = 7, units = 'in', res = 300)
ggplot(meta.data15, aes(y = module.score1, x = newer.ident, fill = factor(newer.ident)))+
  geom_violin()+ facet_grid(~predicted.celltype.l1)+ xlab('')+
  ylab('Module Score')+scale_fill_manual(values = colors)+
  labs(fill='Genotype') +
  theme(
    strip.text = element_text(size = 12, color = "black", family = 'arial'),axis.text.x  = element_blank(),
    axis.text.y = element_text(size = 10, color = "black", family = 'arial'),
    axis.title.y = element_text(size = 15, color = "black", family = 'arial'),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar",
               width = 0.5,
               #fatten = ,
               color = "black")+ggtitle('15 DEGs AddModuleScore')
dev.off()


# Function arguments
object = object.sct.filtered
features = list(genes_fit$Gene)
pool = rownames(object)
nbin = 24
ctrl = 100
k = FALSE
name = "module.score"
seed = 1

# Find how many gene lists were provided. In this case just one.
cluster.length <- length(x = features)

# Pull the expression data from the provided Seurat object
assay.data <- GetAssayData(object = object)
# For all genes, get the average expression across all cells (named vector)
data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
# Order genes from lowest average expression to highest average expression
data.avg <- data.avg[order(data.avg)]

# Use ggplot2's cut_number function to make n groups with (approximately) equal numbers of observations. 
# The 'rnorm(n = length(data.avg))/1e+30' part adds a tiny bit of noise to the data, presumably to break ties.
data.cut <- ggplot2::cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30,
                                n = nbin,
                                labels = FALSE,
                                right = FALSE)

# Set the names of the cuts as the gene names
names(x = data.cut) <- names(x = data.avg)

# Create an empty list the same length as the number of input gene sets. This will contain the names of the control genes
ctrl.use <- vector(mode = "list", length = cluster.length)

# For each of the input gene lists:
for (i in 1:cluster.length) {
  # Get the gene names from the input gene set as a character vector  
  features.use <- features[[i]]
  
  # Loop through the provided genes (1:num_genes) and for each gene, find ctrl (default=100) genes from the same expression bin (by looking in data.cut):
  for (j in 1:length(x = features.use)) {
    # Within this loop, 'data.cut[features.use[j]]' gives us the expression bin number. 
    # We then sample `ctrl` genes from that bin without replacement and add the gene names to ctrl.use.
    ctrl.use[[i]] <- c(ctrl.use[[i]],
                       names(x = sample(x = data.cut[which(x = data.cut == data.cut[features.use[j]])],
                                        size = ctrl,
                                        replace = FALSE)))
  }
}

# Have a quick look at what's in ctrl.use:
class(ctrl.use)
## [1] "list"
length(ctrl.use)
## [1] 1
class(ctrl.use[[1]])
## [1] "character"
# There should be length(features.use)*ctrl genes (i.e. 15*100):
length(ctrl.use[[1]])
## [1] 1500

# Plot the bins that have been created to split genes based on their average expression
plot(data.avg, pch=16, ylab="Average expression across all cells", xlab="All genes, ranked")
for(i in unique(data.cut)){
  cut_pos <- which(data.cut==i)
  if(i%%2==0){
    rect(xleft = cut_pos[1], ybottom = min(data.avg), xright = cut_pos[length(cut_pos)], ytop = max(data.avg), col=scales::alpha("grey", 0.3))
  } else {
    rect(xleft = cut_pos[1], ybottom = min(data.avg), xright = cut_pos[length(cut_pos)], ytop = max(data.avg), col=scales::alpha("white", 0.3))
  }
}

# Add red points for selected control genes
points(which(names(data.avg)%in%ctrl.use[[1]]), data.avg[which(names(data.avg)%in%ctrl.use[[1]])], pch=16, col="red")

# Add blue points for genes in the input gene list
points(which(names(data.avg)%in%features[[1]]), data.avg[which(names(data.avg)%in%features[[1]])], pch=16, col="blue")

# Add a legend
legend(x = "topleft",
       legend = c("gene", "selected control gene", "gene in geneset"),
       col = c("black", "red", "blue"),
       pch = 16)

# Remove any repeated gene names - even though we set replace=FALSE when we sampled genes from the same expression bin, there may be more than two genes in our input gene list that fall in the same expression bin, so we can end up sampling the same gene more than once.
ctrl.use <- lapply(X = ctrl.use, FUN = unique)


## Get control gene scores

# Create an empty matrix with dimensions;
# number of rows equal to the number of gene sets (just one here)
# number of columns equal to number of cells in input Seurat object
ctrl.scores <- matrix(data = numeric(length = 1L),
                      nrow = length(x = ctrl.use),
                      ncol = ncol(x = object))

# Loop through each provided gene set and add to the empty matrix the mean expression of the control genes in each cell
for (i in 1:length(ctrl.use)) {
  # Get control gene names as a vector  
  features.use <- ctrl.use[[i]]
  # For each cell, calculate the mean expression of *all* of the control genes 
  ctrl.scores[i, ] <- Matrix::colMeans(x = assay.data[features.use,])
}


## Get scores for input gene sets

# Similar to the above, create an empty matrix
features.scores <- matrix(data = numeric(length = 1L),
                          nrow = cluster.length,
                          ncol = ncol(x = object))

# Loop through input gene sets and calculate the mean expression of these genes for each cell
for (i in 1:cluster.length) {
  features.use <- features[[i]]
  data.use <- assay.data[features.use, , drop = FALSE]
  features.scores[i, ] <- Matrix::colMeans(x = data.use)
}

# Subtract the control scores from the feature scores - the idea is that if there is no enrichment of the genes in the geneset in a cell, 
# then the result of this subtraction should be ~ 0
features.scores.use <- features.scores - ctrl.scores

# Name the result the "name" variable + whatever the position the geneset was in the input list, e.g. "Cluster1"
rownames(x = features.scores.use) <- paste0(name, 1:cluster.length)

# Change the matrix from wide to long
features.scores.use <- as.data.frame(x = t(x = features.scores.use))

# Give the rows of the matrix, the names of the cells
rownames(x = features.scores.use) <- colnames(x = object)

# Add the result as a metadata column to the input Seurat object 
object[[colnames(x = features.scores.use)]] <- features.scores.use

# Voila!
VlnPlot(object,
        features = "module.score1", group.by = 'newer.ident', split.by = 'predicted.celltype.l1') +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

object <- AddModuleScore(object,features = list(genes_fit$Gene), name = "module.score")
meta.data <- object@meta.data
x <- table(object$newer.ident, object$predicted.celltype.l1)
x

meta.data$predicted.celltype <- NA
meta.data$predicted.celltype[which(meta.data$predicted.celltype.l1 == "B")] <- "B cells"
meta.data$predicted.celltype[which(meta.data$predicted.celltype.l1 == "NK")] <- "NK cells"
meta.data$predicted.celltype[which(meta.data$predicted.celltype.l1 == "CD4 T")] <- "CD4+ T cells"
meta.data$predicted.celltype[which(meta.data$predicted.celltype.l1 == "CD8 T")] <- "CD8+ T cells"
meta.data$predicted.celltype[which(meta.data$predicted.celltype.l1 == "Mono")] <- "Monocytes"
meta.data$predicted.celltype[which(meta.data$predicted.celltype.l1 == "DC")] <- "Dendritic cells"
meta.data$predicted.celltype[which(meta.data$predicted.celltype.l1 == "other")] <- "Others"
meta.data$predicted.celltype[which(meta.data$predicted.celltype.l1 == "other T")] <- "Other T cells"

png('module score violin plots by level 1 identity and genotype.png', width = 11, height = 7, units = 'in', res = 300)
ggplot(meta.data, aes(y = module.score1, x = newer.ident, fill = factor(newer.ident)))+
  geom_violin()+ facet_wrap(~predicted.celltype, nrow = 2)+ xlab('')+ylab('Module Score')+scale_fill_manual(values = colors)+
  labs(fill='Genotype') +
  theme(
    strip.text = element_text(size = 12, color = "black", family = 'arial'),
    axis.text.x  = element_blank(),
    axis.text.y = element_text(size = 10, color = "black", family = 'arial'),
    axis.title.y = element_text(size = 15, color = "black", family = 'arial'),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")
dev.off()

###########
# Determine contribution of cells expressing the top X DEGs using AUCell
##########

BiocManager::install("AUCell", force = T)
library(AUCell)
genes_fit <- DEA_across_all_cells_ranked_by_FC %>%
  dplyr::slice_head(n =15)

geneSets <- list(DEGs_15=genes_fit$Gene)
t <- merge(object, assay = "RNA", merge.data = TRUE)
expreMatrix <- t@assays$RNA$counts

AUcell_output <- AUCell_run(expreMatrix, geneSets = geneSets)
set.seed(333)
par(mfrow=c(3,3)) 
cells_assignment <- AUCell_exploreThresholds(AUcell_output, plotHist=TRUE, assign=TRUE) 
cells_assignment$DEGs_15$aucThr$thresholds
## 15 genes
#threshold nCells
#Global_k1 0.1215571  24412
#L_k2      0.2116392   2345
#R_k3      0.0340422 100688

##50 genes 
#threshold nCells
#Global_k1 0.12161346  24485
#L_k2      0.21191785   2342
#R_k3      0.03218605 101943
cells_assignment$DEGs_15$aucThr$selected
## 15 genes
#L_k2 
#0.2116392 

## 50 genes
#L_k2 
#0.2119178
CellsAssigned <- cells_assignment$DEGs_15$assignment
length(CellsAssigned)
#15
#2345 cells

#50
#2342


t<- AddMetaData(object = object, metadata = AUcell_output@assays@data@listData, col.name = "AUC")
VlnPlot(object = t, features = "AUC", group.by = 'predicted.celltype.l1',split.by = 'newer.ident')

meta.dataAUC15 <- t@meta.data
Patients_AUC15 <-meta.dataAUC15[meta.dataAUC15$newer.ident =='R57C', ]

meta.dataAUC50 <- t@meta.data


png("Cell type level 1 15 DEGs AUCell.png", width = 11, height = 7, units = 'in', res = 300)
ggplot(meta.dataAUC15, aes(y = AUC, x = newer.ident, fill = factor(newer.ident)))+
  geom_violin()+ facet_grid(~predicted.celltype.l1)+ xlab('')+
  ylab('Module Score')+scale_fill_manual(values = colors)+
  labs(fill='Genotype') +
  theme(
    strip.text = element_text(size = 12, color = "black", family = 'arial'),axis.text.x  = element_blank(),
    axis.text.y = element_text(size = 10, color = "black", family = 'arial'),
    axis.title.y = element_text(size = 15, color = "black", family = 'arial'),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar",
               width = 0.5,
               #fatten = ,
               color = "black")+ggtitle('15 DEGs AUCell')+stat_compare_means(comparisons = my_comparisons, method = 't.test')
dev.off()

png("Cell type level 1 50 DEGs AUCell.png", width = 11, height = 7, units = 'in', res = 300)
ggplot(meta.dataAUC50, aes(y = AUC, x = newer.ident, fill = factor(newer.ident)))+
  geom_violin()+ facet_grid(~predicted.celltype.l1)+ xlab('')+
  ylab('Module Score')+scale_fill_manual(values = colors)+
  labs(fill='Genotype') +
  theme(
    strip.text = element_text(size = 12, color = "black", family = 'arial'),axis.text.x  = element_blank(),
    axis.text.y = element_text(size = 10, color = "black", family = 'arial'),
    axis.title.y = element_text(size = 15, color = "black", family = 'arial'),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar",
               width = 0.5,
               #fatten = ,
               color = "black")+ggtitle('50 DEGs AUCell')
dev.off()

#PateintsAUCell50 <- subset(t, subset = newer.ident == 'R57C')
#HCsAUCell50 <- subset(t, subset = newer.ident == 'HC')

PateintsAUCell15 <- subset(t, subset = newer.ident == 'R57C')
HCsAUCell15 <- subset(t, subset = newer.ident == 'HC')


FeaturePlot(PateintsAUCell50, features = "AUC", label = F, repel = F,pt.size = 1, raster = T) +
  scale_colour_gradient(low = 'grey', high = 'blue', limits = c(-0.2,0.4), breaks = seq(-0.20, 0.4, 0.2), na.value = "grey", name = "AUC")+
  ggtitle("AUC")+
  guides(x = axis, y = axis)+ xlab('UMAP 1')+ylab('UMAP 2')+
  coord_equal(ylim =c(-16.37798,14.92378), xlim=c(-14.94919,14.49963))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        axis.line = element_line(arrow = arrow(type='closed', length = unit(3,'pt'))),
        axis.title = element_text(hjust = 0, size = 12, family = 'arial'),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks = element_blank())

png("UMAP Cell type level 1 50 DEGs R57C AUCell.png", width = 11, height = 7, units = 'in', res = 300)
FeaturePlot(PateintsAUCell50, features = "AUC", label = F, repel = F,pt.size = 1, raster = T) +
  scale_colour_gradientn(colours = viridis(3), limits = c(0,0.3), breaks = seq(0, 0.3, 0.1), na.value = viridis(1), name = "AUC")+
  ggtitle("R57C AUC")+
  guides(x = axis, y = axis)+ xlab('UMAP 1')+ylab('UMAP 2')+
  coord_equal(ylim =c(-16.37798,14.92378), xlim=c(-14.94919,14.49963))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        axis.line = element_line(arrow = arrow(type='closed', length = unit(3,'pt'))),
        axis.title = element_text(hjust = 0, size = 12, family = 'arial'),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks = element_blank())
dev.off()
png("UMAP Cell type level 1 50 DEGs HC AUCell.png", width = 11, height = 7, units = 'in', res = 300)
FeaturePlot(HCsAUCell50, features = "AUC", label = F, repel = F,pt.size = 1, raster = T) +
  scale_colour_gradientn(colours = viridis(3), limits = c(0,0.3), breaks = seq(0, 0.3, 0.1), na.value = "grey", name = "AUC")+
  ggtitle("HC AUC")+
  guides(x = axis, y = axis)+ xlab('UMAP 1')+ylab('UMAP 2')+
  coord_equal(ylim =c(-16.37798,14.92378), xlim=c(-14.94919,14.49963))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        axis.line = element_line(arrow = arrow(type='closed', length = unit(3,'pt'))),
        axis.title = element_text(hjust = 0, size = 12, family = 'arial'),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks = element_blank())
dev.off()
