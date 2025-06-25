library(viridis)
library(colorspace)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice) 
library(dplyr)
library(viridis)
library(ggokabeito)
library(scCustomize)
library(colorBlindness)
library(groupdata2)

#############################################################
### OTULIN UMAPs: density/genotype and consensus for level 1 identities
#############################################################

setwd("/Volumes/Markle Lab/Data/Christine Mariskanish/10604-JM experiment/OTULIN/")

object.sct.filtered <- readRDS("object_postfiltering.rds")

df <- Embeddings(object.sct.filtered, reduction = "ref.umap")
df <- as.data.frame(df)
colnames(df) <- c("UMAP1", "UMAP2")
metadata <- object.sct.filtered@meta.data
# Include cell type predictions in the metadata
df <- cbind(df, metadata)

#downsample for UMAP visualization
df <- downsample(df, cat_col = 'newer.ident')
table(df$newer.ident)
#31885 cells per genotype (WT and R57C)

df$predicted.celltype <- NA
df$predicted.celltype[which(df$predicted.celltype.l1 == "B")] <- "B cells"
df$predicted.celltype[which(df$predicted.celltype.l1 == "NK")] <- "NK cells"
df$predicted.celltype[which(df$predicted.celltype.l1 == "CD4 T")] <- "CD4+ T cells"
df$predicted.celltype[which(df$predicted.celltype.l1 == "CD8 T")] <- "CD8+ T cells"
df$predicted.celltype[which(df$predicted.celltype.l1 == "Mono")] <- "Monocytes"
df$predicted.celltype[which(df$predicted.celltype.l1 == "DC")] <- "Dendritic cells"
df$predicted.celltype[which(df$predicted.celltype.l1 == "other")] <- "Others"
df$predicted.celltype[which(df$predicted.celltype.l1 == "other T")] <- "Other T cells"
table(df$predicted.celltype, df$predicted.celltype.l1)

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


cluster_colors <- c(
  "HC" = "blue",
  "R57C" = "red"
)


end_colors <- lighten(cluster_colors, amount = 0.8)


blend_colors <- function(cluster, density, colors, colors_end) {
  base_color <- colors[cluster]
  end_color <- colors_end[cluster]
  if (density > .999){blend_color = end_color} else {blend_color <- colorRampPalette(c(base_color, end_color))(100)[as.integer(density * 100) + 1]}
}

dfc1 <- df[df$newer.ident  == 'HC',]
dfp16 <- df[df$newer.ident  == 'R57C',]


### Controls normalized
dfc1$Density <- get_density(dfc1$UMAP1, dfc1$UMAP2, n =500)

dfc1 <- dfc1 %>%
  mutate(normalized_density = (Density - min(Density)) / (max(Density) - min(Density)))

# Apply blending function
dfc1$color <- mapply(blend_colors, dfc1$newer.ident, dfc1$normalized_density, MoreArgs = list(colors = end_colors, colors_end = cluster_colors))


### Patients normalized
dfp16$Density <- get_density(dfp16$UMAP1, dfp16$UMAP2, n =500)

dfp16 <- dfp16 %>%
  mutate(normalized_density = (Density - min(Density)) / (max(Density) - min(Density)))

# Apply blending function
dfp16$color <- mapply(blend_colors, dfp16$newer.ident, dfp16$normalized_density, MoreArgs = list(colors = end_colors, colors_end = cluster_colors))

axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(3, "cm")
)

safe_colorblind_palette <- c("#117733", "#4B0055", 
                             "#000000", "#888888", "orange",
                             "#CC79A7", "#28E2E5","#A65628")
colorspace::swatchplot(safe_colorblind_palette, cvd = TRUE)

umap_plot_comb <- ggplot(df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(colour = predicted.celltype), size = 3) + coord_equal(ylim =c(-16.61637,15.09251), xlim=c(-15.39274,14.75713))+
  scale_colour_manual(breaks = c("B cells","CD4+ T cells", "CD8+ T cells", "Monocytes", "Dendritic cells", "Other T cells", "NK cells",
                                 "Others"), values = safe_colorblind_palette[1:8])+
  theme_minimal(base_size = 14) +
  #annotate(geom = "label", label = "B cells",
          # x = -10, y = -10, size = 7, colour = "black", family = 'arial')+
  #annotate(geom = "label", label = "CD4+ T cells",
           #x = 3, y = 9, size = 7, colour = "black", family = 'arial')+
  #annotate(geom = "label", label = "CD8+ T cells",
           #x = 7, y = -1, size = 7, colour = "black", family = 'arial')+
  #annotate("segment", x = 1, xend = 1.5, y = -6, yend = -8, colour = "black") +
  #annotate(geom = "label", label = "NK cells",
           #x = -1, y = -6.5, size = 7, colour = "black", family = 'arial')+
  #annotate(geom = "label", label = "Monocytes",
           #x = -9, y = 7, size = 7, colour = "black", family = 'arial')+
  #annotate(geom = "label", label = "Other T cells",
           #x = 9, y = -10, size = 7, colour = "black", family = 'arial')+
  #annotate(geom = "label", label = "Others",
           #x = -2, y = 4, size = 7, colour = "black", family = 'arial')+
  #annotate("segment", x = -8, xend = -2, y = 12, yend = 13.5, colour = "black") +
  #annotate(geom = "label", label = "Dendritic cells",
           #x = -9, y = 13.5, size = 7, colour = "black", family = 'arial')+
  guides(x = axis, y = axis)+ xlab('UMAP 1')+ylab('UMAP 2')+
  #theme_set(theme_bw(18))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(arrow = arrow(type='closed', length = unit(3,'pt'))),
        axis.title = element_text(hjust = 0, size = 12),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")+ggtitle("Combined")+theme(panel.grid.major = element_blank(),
                                                            panel.grid.minor = element_blank())+theme(aspect.ratio = 1)
b <- ggplot_build(umap_plot_comb)
b$layout$panel_params[[1]]$x.range
b$layout$panel_params[[1]]$y.range

umap_plot_c1 <- ggplot(dfc1, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = color), alpha = 0.4, size = 3) + coord_equal(ylim =c(-16.61637,15.09251), xlim=c(-15.39274,14.75713))+
  scale_color_identity()+
  theme_minimal(base_size = 14) +
  guides(x = axis, y = axis)+ xlab('UMAP 1')+ylab('UMAP 2')+
  #theme_set(theme_bw(18))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(arrow = arrow(type='closed', length = unit(3,'pt'))),
        axis.title = element_text(hjust = 0, size = 12),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") + ggtitle('HC')+theme(panel.grid.major = element_blank(),
                                                         panel.grid.minor = element_blank())+theme(aspect.ratio = 1)


umap_plot_p16 <- ggplot(dfp16, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = color), alpha = 0.4, size = 3) + coord_equal(ylim =c(-16.61637,15.09251), xlim=c(-15.39274,14.75713))+
  scale_color_identity()+
  theme_minimal(base_size = 14) +
  guides(x = axis, y = axis)+ xlab('UMAP 1')+ylab('UMAP 2')+
  #theme_set(theme_bw(18))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(arrow = arrow(type='closed', length = unit(3,'pt'))),
        axis.title = element_text(hjust = 0, size = 12),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")+ggtitle("R57C") +theme(panel.grid.major = element_blank(),
                                                                panel.grid.minor = element_blank())+theme(aspect.ratio = 1)

umap_plot_comb$theme$text$family = 'arial'
umap_plot_c1$theme$text$family = 'arial'

umap_plot_p16$theme$text$family = 'arial'

umap_plot_comb$theme$text$size = 12
umap_plot_c1$theme$text$size = 12

umap_plot_p16$theme$text$size = 12

## Define function for legend
gradient <- function(cell_type, colors, colors_end){
  base_color = colors[cell_type]
  end_color = colors_end[cell_type]
  colorRampPalette(c(base_color, end_color))(10)
}
names <- c(
  "HCs" = "blue",
  "R57C" = "red"
)
## Create Legend dataframe
legend <- data.frame(
  Patient = rep(names(names), each = 10),
  density = rep(seq(0,1,length.out = 10), length(cluster_colors)),
  color = unlist(lapply(names(cluster_colors), function(ct) gradient(ct, end_colors, cluster_colors)))
)
legendHCs <- legend[legend$Patient  == 'HCs',]
legendHom <- legend[legend$Patient  == "R57C",]

legend_plot_HCs <- ggplot(legendHCs, aes(x = Patient, y = density, fill = color)) +
  geom_tile() +
  scale_fill_identity() +
  labs(title = "Density", x = "", y = "") +
  theme_minimal(base_size = 12, base_family = 'arial') +
  theme(axis.text.x = element_blank(),plot.title = element_text(family = 'arial', size = 12, hjust = .5),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none"
  )

legend_plot_Hom <- ggplot(legendHom, aes(x = Patient, y = density, fill = color)) +
  geom_tile() +
  scale_fill_identity() +
  labs(title = "Density", x = "", y = "") +
  theme_minimal(base_size = 12, base_family = 'arial') +
  theme(axis.text.x = element_blank(),plot.title = element_text(family = 'arial', size = 12,
                                                                hjust = .5),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none"
  )


legend_plot_HCs <- legend_plot_HCs + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())


legend_plot_Hom <- legend_plot_Hom + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
# Convert the custom legend plot to a grob
legend_grob_HCs <- ggplotGrob(legend_plot_HCs)

legend_grob_Hom <- ggplotGrob(legend_plot_Hom)

cowplot::plot_grid(umap_plot_comb, NULL, NULL,umap_plot_c1, legend_grob_HCs, NULL,umap_plot_p16,
                   legend_grob_Hom, labels=NULL, ncol = 3, rel_widths=c(6,0.5,.5),
                   rel_heights = c(6,6,6))
setwd('Figures')
png('Supervised_UMAP level1 no labels.png',height = 20,width = 20, units = 'in', res = 300)
cowplot::plot_grid(umap_plot_comb, NULL,NULL, NULL, umap_plot_c1, NULL,legend_grob_HCs, NULL, umap_plot_p16,NULL,
                   legend_grob_Hom, NULL, labels=NULL, ncol = 4, rel_widths=c(6,0,0.5,3),
                   rel_heights = c(6,6,6))
dev.off()
##########################
### stacked bar plots
##########################
df <- Embeddings(object.sct.filtered, reduction = "ref.umap")
df <- as.data.frame(df)
colnames(df) <- c("UMAP1", "UMAP2")
metadata <- object.sct.filtered@meta.data
# Include cell type predictions in the metadata
df <- cbind(df, metadata)
df$predicted.celltype <- NA
df$predicted.celltype[which(df$predicted.celltype.l1 == "B")] <- "B cells"
df$predicted.celltype[which(df$predicted.celltype.l1 == "NK")] <- "NK cells"
df$predicted.celltype[which(df$predicted.celltype.l1 == "CD4 T")] <- "CD4+ T cells"
df$predicted.celltype[which(df$predicted.celltype.l1 == "CD8 T")] <- "CD8+ T cells"
df$predicted.celltype[which(df$predicted.celltype.l1 == "Mono")] <- "Monocytes"
df$predicted.celltype[which(df$predicted.celltype.l1 == "DC")] <- "Dendritic cells"
df$predicted.celltype[which(df$predicted.celltype.l1 == "other")] <- "Others"
df$predicted.celltype[which(df$predicted.celltype.l1 == "other T")] <- "Other T cells"
table(df$predicted.celltype, df$predicted.celltype.l1)
cell_type_proportions <- df %>%
  group_by(new.ident, newer.ident, predicted.celltype) %>%
  summarise(celltype_count = n()) %>%
  group_by(new.ident) %>%
  mutate(total_cell_count = sum(celltype_count),
         proportion = celltype_count / total_cell_count)
# check if cluster_proportions is correct
x <- sum(cell_type_proportions$proportion[1:8])
x

bar_plot <- ggplot(cell_type_proportions, aes(x = new.ident, y = proportion, fill = predicted.celltype)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(#title = "Distribution of Patient Cells per Cluster",
    x = "",
    y = "Frequency",
    fill = "Cell type") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.location = "panel", axis.ticks.x = element_blank(),
        legend.title = element_text(size = 14, family = 'arial'),
        legend.text = element_text(size = 12, family = 'arial'),
        axis.title.y = element_text(size = 12, family = 'arial'),
        axis.text = element_text(size = 12, family = 'arial'))
bar_plot$theme$text$family = 'arial'

# Display the plot
safe_colorblind_palette <- c("#117733", "#4B0055", 
                             "#000000", "#888888", "orange",
                             "#28E2E5", "#CC79A7","#A65628")
bar <- bar_plot + scale_fill_manual(breaks = c("B cells","CD4+ T cells", "CD8+ T cells", "Monocytes", "Dendritic cells", "NK cells",
                                               "Other T cells","Others"), values = safe_colorblind_palette[1:8])
bar <- bar + theme(legend.key.size = unit(0.4, "cm"))
bar <- bar+theme(panel.grid.major.x = element_blank())
bar

##########################
### histograms by cell type
##########################

scaleFUN <- function(x) sprintf("%.2f", x)

png(filename = 'Frequency of L1 identity barplots.png', width = 11, height = 7, units = 'in', res = 300)
ggplot(data = cell_type_proportions, mapping = aes(x = newer.ident,
                                                   y = proportion,
                                                   fill= newer.ident)) +
  stat_summary(aes(fill = newer.ident),
               fun.data = mean_se, geom = "col", color ="black", position = position_dodge(1)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position = position_dodge(1)) +
  geom_jitter(
              position = position_jitterdodge(dodge.width = 1, jitter.width = .1),
              size = 1) +
  theme_minimal(base_size = 13) + ylab('Proportion')+ xlab('')+
  scale_fill_manual(values = c("blue", "red")) +
  theme(strip.background =element_rect(fill="#FFFFD8"),
              strip.text = element_text(size = 12, color = "black", family = 'arial'),
              axis.text.x  = element_blank(),
              axis.text.y = element_text(size = 10, color = "black", family = 'arial'),
              axis.title.y = element_text(size = 15, color = "black", family = 'arial'),
              axis.ticks.x = element_blank(),
              panel.grid.major.x = element_blank())+
  facet_wrap(~predicted.celltype)+
  scale_y_continuous(expand = c(0, 0), labels=scaleFUN)+
  NoLegend()+
  coord_cartesian(ylim = c(0, 0.75))
dev.off()


write.csv(cell_type_proportions, '/Volumes/Markle Lab/Data/Christine Mariskanish/10604-JM experiment/OTULIN/cell_type_proportions.csv')
