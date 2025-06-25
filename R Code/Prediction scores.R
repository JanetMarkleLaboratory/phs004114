##########################
### load in the necessary libraries and set the working directory
##########################
install.packages("extrafont")
library(extrafont)
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

getwd()
setwd("/Volumes/pmi/Data/Research/Markle Lab/Data/Christine Mariskanish/10604-JM experiment/OTULIN/")


##########################
### Predicted cell type level 1 scores boxplots
##########################
object <- readRDS('prefiltered_object.rds')

object$newer.ident <- NA
object$newer.ident[which(object$orig.ident == "0014")] <- "HC"
object$newer.ident[which(object$orig.ident == "0015")] <- "HC"
object$newer.ident[which(object$orig.ident == "0016")] <- "HC"
object$newer.ident[which(object$orig.ident == "0012")] <- "R57C"
object$newer.ident[which(object$orig.ident == "0013")] <- "R57C"
object$newer.ident <- factor(x = object$newer.ident, levels = c("HC", "R57C"))
object$new.ident <- factor(x = object$new.ident, levels = c("HC 1", "HC 2", "HC 3", "R57C 1", "R57C 2"))

df <- data.frame(object$newer.ident)
df <- tibble::rownames_to_column(df, "newer.ident")
score1 <- data.frame(object$predicted.celltype.l1.score)
score1 <- tibble::rownames_to_column(score1, "newer.ident")
score2 <- data.frame(object$predicted.celltype.l2.score)
score2 <- tibble::rownames_to_column(score2, "newer.ident")
ident2 <- data.frame(object$predicted.celltype.l2)
ident2 <- tibble::rownames_to_column(ident2, "newer.ident")
ident1 <- data.frame(object$predicted.celltype.l1)
ident1 <- tibble::rownames_to_column(ident1, "newer.ident")
sample <- data.frame(object$new.ident)
sample <- tibble::rownames_to_column(sample, "newer.ident")
df <- merge(score1, df, by = "newer.ident")
df <- merge(score2, df, by = "newer.ident")
df <- merge(ident2, df, by= "newer.ident")
df <- merge(ident1, df, by= "newer.ident")
df<- merge(df, sample, by = "newer.ident")

names <- c(
  "HC" = "blue",
  "R57C" = "red"
)

boxplot.l1 <- df %>%
  ggplot( aes(x=factor(object.new.ident), y=object.predicted.celltype.l1.score, fill=factor(object.newer.ident))) +
  geom_violin()+ geom_point(size = 0)+
  geom_boxplot(width=0.25, color="black", alpha=0.01,outlier.colour="blue", outlier.size=0.1, outlier.shape = NA)+
  scale_fill_manual(values =names) +
  theme_bw()+coord_cartesian(ylim = c(0,1)) +
  theme(text=element_text(family="arial", size=12),
        legend.position="none"
  ) +
  xlab("") + ylab("Cell Type Level 1 Prediction Score")+ 
  theme(panel.grid.major.x = element_blank())+
  geom_hline(yintercept = .75, linetype = 2)
boxplot.l1
boxplot.l1.granular <- df %>%
  ggplot( aes(x=factor(object.new.ident), y=object.predicted.celltype.l1.score, fill=factor(object.newer.ident))) +
  geom_violin()+ geom_point(size = 0)+
  facet_wrap(~object.predicted.celltype.l1)+
  geom_boxplot(width=0.25, color="black", alpha=0.01,outlier.colour="blue", outlier.size=0.1, outlier.shape = NA)+
  scale_fill_manual(values =names) +
  theme_bw()+coord_cartesian(ylim = c(0,1)) +
  theme(text=element_text(family="arial", size=12),
        legend.position="none"
  ) +
  xlab("") + ylab("Cell Type Level 1 Prediction Score")+ 
  theme(panel.grid.major.x = element_blank())+
  geom_hline(yintercept = .75, linetype = 2)
boxplot.l1
boxplot.l1.granular
##########################
### Predicted cell type level 2 scores boxplots
##########################

boxplot.l2 <- df %>%
  ggplot( aes(x=object.new.ident, y=object.predicted.celltype.l2.score, fill=object.newer.ident)) +
  facet_wrap(~object.predicted.celltype.l2)+
  geom_violin() + geom_point(size = 0)+geom_boxplot(width=0.25, color="black", alpha=0.01, 
                                                    outlier.colour="blue", outlier.size=0.1, outlier.shape = NA)+
  scale_fill_manual(values=names) +
  theme_bw()+coord_cartesian(ylim = c(0,1)) +
  theme(text=element_text(family="arial", size=12),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.position="none"
  ) +
  xlab("") + ylab("Cell Type Level 2 Prediction Score")+ theme(panel.grid.major.x = element_blank())#+

boxplot.l2

