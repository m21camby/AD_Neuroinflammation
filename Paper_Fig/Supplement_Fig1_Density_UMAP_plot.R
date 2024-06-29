#! /bin/env RScript
# written by SJK at 1. Oct. 2020

# From the Samantha's paper in the figure legend:
# it is tumor-specific tSNE representation as indicated by 
# a color scale ranging from gray (no expression) to dark blue (high expression) with regions of high expression highlighted

# From the Samantha's paper in the figure text:
# This combined analysis enabled us to recapitulate all previously identified cell types without the need to resort to advanced 
# sample alignment methods21 as cells from both control and double-mutant mice were distributed evenly on the tSNE for many clusters 
# as shown in the local density plot

# From the Samantha's paper in the method text:
# In method part, The relative density of two sample groups on the t-SNE was plotted using the log2 ratio of two separate 
# 2D kernel density estimators interpolated on the t-SNE coordinates of each cell. 

library(Seurat)
library(ggplot2)
library(dplyr)

source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")
source("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200927_Benedikt_UMAP_Code.R")

load(file="/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200210_figures_for_meeting_report_cleaned_Seurat_object_res0.8.Robj")

# colors.use has to be changed by each figure in 20200927_Benedikt_UMAP_Code
# WT vs AD colors.use=c('red','blue','black')
# AD vs ADp40KO colors.use=c('black', 'red','blue')
# WT vs ADp40KO colors.use=c('red','black','blue')


##############
# WT vs AD
##############
data9set_cleaned_sub.SO <- subset(data9set_cleaned.SO, subset = sample %in% c("AD", "Ctrl"))

data9set_cleaned_sub.SO$sample <- data9set_cleaned_sub.SO$sample %>% as.factor

g1 <- smooth_UMAPPlot(data9set_cleaned_sub.SO, group.by = "sample", legend.size=12, ref = 1, alt = 2, label1 = "WT", label2 = "APPPS1")  


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig1_Density_UMAP_plot_WT_AD.png", 
       plot = g1, 
       scale = 1, width = 8, height = 8, units = "in", device = "png",
       dpi = 300)

print("5th")

###############
# AD vs ADp40KO
###############
data9set_cleaned_sub.SO <- subset(data9set_cleaned.SO, subset = sample %in% c("AD", "ADp40KO"))
data9set_cleaned_sub.SO$sample <- data9set_cleaned_sub.SO$sample %>% as.factor

g2 <- smooth_UMAPPlot(data9set_cleaned_sub.SO, group.by = "sample", legend.size=12, ref = 2, alt = 3, label1 = "APPPS1", label2 = "APPPS1.il12b-/-")  


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig1_Density_UMAP_plot_AD_ADp40KO.png", 
       plot = g2, 
       scale = 1, width = 8, height = 8, units = "in", device = "png",
       dpi = 300)

###############
# WT vs ADp40KO
###############
data9set_cleaned_sub.SO <- subset(data9set_cleaned.SO, subset = sample %in% c("ADp40KO", "Ctrl"))
data9set_cleaned_sub.SO$sample <- data9set_cleaned_sub.SO$sample %>% as.factor

g3 <- smooth_UMAPPlot(data9set_cleaned_sub.SO, group.by = "sample", legend.size=12, ref = 1, alt = 3, label1 = "WT", label2 = "APPPS1.il12b-/-")  

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig1_Density_UMAP_plot_WT_ADp40KO.png", 
       plot = g3, 
       scale = 1, width = 8, height = 8, units = "in", device = "png",
       dpi = 300)


##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig1_Density_UMAP_plot_session_info.txt")


