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
library(MASS)
#############################################################
### OTULIN feature density UMAPs
#############################################################

setwd("/Volumes/Markle Lab/Data/Christine Mariskanish/10604-JM experiment/OTULIN/")

object.sct.filtered <- readRDS("object_postfiltering.rds")

getwd()
setwd('Figures')

Patients<- subset(object.sct.filtered, subset= newer.ident =="R57C")
HCs<- subset(object.sct.filtered, subset= newer.ident =="HC")


axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(3, "cm")
)
setwd('Figures/Feature Density UMAPs/')
png("HC CXCL8 Feature Density.png", width = 11, height = 7, units = 'in', res = 300)
Nebulosa::plot_density(HCs, features = c("CXCL8"),size = 2) + 
  scale_colour_gradientn(colours = c("grey",viridis(1),"#FDE725FF"),limits = c(0.0,0.5), breaks = seq(0.0, 0.5, 0.1), na.value = "grey")+
  ggtitle("HC CXCL8 Feature Density")+guides(x = axis, y = axis)+ xlab('UMAP 1')+ylab('UMAP 2')+
  coord_equal(ylim =c(-16.61637,15.09251), xlim=c(-15.39274,14.75713))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        axis.line = element_line(arrow = arrow(type='closed', length = unit(3,'pt'))),
        axis.title = element_text(hjust = 0, size = 12, family = 'arial'),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks = element_blank())
dev.off()

png("R57C CXCL8 Feature Density.png", width = 11, height = 7, units = 'in', res = 300)
Nebulosa::plot_density(Patients, features = c("CXCL8"),size = 2) + 
  scale_colour_gradientn(colours = c("grey",viridis(1),"#FDE725FF"),limits = c(0.0,0.5), breaks = seq(0.0, 0.5, 0.1), na.value = "grey")+
  ggtitle("R57C CXCL8 Feature Density")+guides(x = axis, y = axis)+ xlab('UMAP 1')+ylab('UMAP 2')+
  coord_equal(ylim =c(-16.61637,15.09251), xlim=c(-15.39274,14.75713))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        axis.line = element_line(arrow = arrow(type='closed', length = unit(3,'pt'))),
        axis.title = element_text(hjust = 0, size = 12, family = 'arial'),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks = element_blank())
dev.off()
png("HC IL1B Feature Density.png", width = 11, height = 7, units = 'in', res = 300)
Nebulosa::plot_density(HCs, features = c("IL1B"), size = 2) + 
  scale_colour_gradientn(colours = c("grey",viridis(1),"#FDE725FF"),limits = c(0.0,0.5), breaks = seq(0.0, 0.5, 0.1), na.value = "grey")+
  ggtitle("HC IL1B Feature Density")+guides(x = axis, y = axis)+ xlab('UMAP 1')+ylab('UMAP 2')+
  coord_equal(ylim =c(-16.61637,15.09251), xlim=c(-15.39274,14.75713))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        axis.line = element_line(arrow = arrow(type='closed', length = unit(3,'pt'))),
        axis.title = element_text(hjust = 0, size = 12, family = 'arial'),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks = element_blank())
dev.off()
png("R57C IL1B Feature Density.png", width = 11, height = 7, units = 'in', res = 300)
Nebulosa::plot_density(Patients, features = c("IL1B"),size = 2) + 
  scale_colour_gradientn(colours = c("grey",viridis(1),"#FDE725FF"),limits = c(0.0,0.5), breaks = seq(0.0, 0.5, 0.1), na.value = "grey")+
  ggtitle("R57C IL1B Feature Density")+guides(x = axis, y = axis)+ xlab('UMAP 1')+ylab('UMAP 2')+
  coord_equal(ylim =c(-16.61637,15.09251), xlim=c(-15.39274,14.75713))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        axis.line = element_line(arrow = arrow(type='closed', length = unit(3,'pt'))),
        axis.title = element_text(hjust = 0, size = 12, family = 'arial'),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks = element_blank())
dev.off()






