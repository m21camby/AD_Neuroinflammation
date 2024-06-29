#! /bin/env RScript
# written by SJK at 23. Mar. 2020

.libPaths(c("/home/skim/R/usr_lib", .libPaths()))

library(Seurat)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)
library(tidyr)
library(grid)

source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

load(file="/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200210_figures_for_meeting_report_cleaned_Seurat_object_res0.8.Robj")
data9set_cleaned_MG.SO <- subset(data9set_cleaned.SO, subset = seurat_clusters %in% c(3, 8))

data9set_cleaned_Ctrl_MG.SO <- subset(data9set_cleaned_MG.SO, subset = sample %in% c("Ctrl"))
data9set_cleaned_AD_MG.SO <- subset(data9set_cleaned_MG.SO, subset = sample %in% c("AD"))
data9set_cleaned_ADp40KO_MG.SO <- subset(data9set_cleaned_MG.SO, subset = sample %in% c("ADp40KO"))


###########################
# PAN Microglia markers
###########################

t1 <- textGrob("Hexb", gp=gpar(fontsize=20, col="black", font = "Helvetica"), hjust = 0.35)
t2 <- textGrob("Arhgap5", gp=gpar(fontsize=20, col="black", font = "Helvetica"), hjust = 0.35)


f1 <- FeaturePlot(data9set_cleaned_Ctrl_MG.SO, features = c("Hexb"), pt.size = 0.5) + 
  ggtitle(expression(~WT)) + 
  scale_x_continuous(limits = c(-5,5)) + scale_y_continuous(limits = c(18, 28)) + 
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", midpoint = 2) + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.6, size = 20), panel.border = element_blank(), axis.line = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) 

f2 <- FeaturePlot(data9set_cleaned_AD_MG.SO, features = c("Hexb"), pt.size = 0.5) + 
  ggtitle(expression(~APPPS1p40)) +
  scale_x_continuous(limits = c(-5,5)) + scale_y_continuous(limits = c(18, 28)) + 
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", midpoint = 2) + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.6, size = 20), panel.border = element_blank(), axis.line = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) 

f3 <- FeaturePlot(data9set_cleaned_ADp40KO_MG.SO, features = c("Hexb"), pt.size = 0.5) +  
  ggtitle(expression(~APPPS1.p40^{"-/-"})) + 
  scale_x_continuous(limits = c(-5,5)) + scale_y_continuous(limits = c(18, 28)) + 
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", midpoint = 2) + 
  theme(plot.title = element_text(hjust = 0.6, size = 20), panel.border = element_blank(), axis.line = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) 

f4 <- FeaturePlot(data9set_cleaned_Ctrl_MG.SO, features = c("Arhgap5"), pt.size = 0.5) + 
  scale_x_continuous(limits = c(-5,5)) + scale_y_continuous(limits = c(18, 28)) + 
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", midpoint = 2) + 
  theme(legend.position = "none", plot.title = element_blank(), panel.border = element_blank(), axis.line = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) 

f5 <- FeaturePlot(data9set_cleaned_AD_MG.SO, features = c("Arhgap5"), pt.size = 0.5) + 
  scale_x_continuous(limits = c(-5,5)) + scale_y_continuous(limits = c(18, 28)) + 
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", midpoint = 2) + 
  theme(legend.position = "none", plot.title = element_blank(), panel.border = element_blank(), axis.line = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) 

f6 <- FeaturePlot(data9set_cleaned_ADp40KO_MG.SO, features = c("Arhgap5"), pt.size = 0.5) +  
  scale_x_continuous(limits = c(-5,5)) + scale_y_continuous(limits = c(18, 28)) + 
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", midpoint = 2) + 
  theme(plot.title = element_blank(), panel.border = element_blank(), axis.line = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) 

g1 <- arrangeGrob(t1, f1,f2,f3,t2, f4,f5,f6, ncol = 4, widths = c(1, 1.6, 1.6, 2), heights = c(2.4,2))


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_DAM_Modified_PAN_MG_markers.png", 
       plot = g1, 
       scale = 1, width = 12, height = 6, units = "in", device = "png",
       dpi = 300)

#######################
# Homeostatic marker
#######################

t1 <- textGrob("C1qa", gp=gpar(fontsize=20, col="black", font = "Helvetica"), hjust = 0.35)
t2 <- textGrob("Cx3cr1", gp=gpar(fontsize=20, col="black", font = "Helvetica"), hjust = 0.35)

f1 <- FeaturePlot(data9set_cleaned_Ctrl_MG.SO, features = c("C1qa"), order = TRUE,pt.size = 0.5) + 
  ggtitle(expression(~WT)) + 
  scale_x_continuous(limits = c(-5,5)) + scale_y_continuous(limits = c(18, 28)) + 
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", midpoint = 2) + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.6, size = 20), panel.border = element_blank(), axis.line = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) 

f2 <- FeaturePlot(data9set_cleaned_AD_MG.SO, features = c("C1qa"), order = TRUE,pt.size = 0.5) + 
  ggtitle(expression(~APPPS1p40)) +
  scale_x_continuous(limits = c(-5,5)) + scale_y_continuous(limits = c(18, 28)) + 
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", midpoint = 2) + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.6, size = 20), panel.border = element_blank(), axis.line = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) 

f3 <- FeaturePlot(data9set_cleaned_ADp40KO_MG.SO, features = c("C1qa"), order = TRUE,pt.size = 0.5) +  
  ggtitle(expression(~APPPS1.p40^{"-/-"})) + 
  scale_x_continuous(limits = c(-5,5)) + scale_y_continuous(limits = c(18, 28)) + 
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", midpoint = 2) + 
  theme(plot.title = element_text(hjust = 0.6, face = "bold", size = 20), panel.border = element_blank(), axis.line = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) 

f4 <- FeaturePlot(data9set_cleaned_Ctrl_MG.SO, features = c("Cx3cr1"), pt.size = 0.5) + 
  scale_x_continuous(limits = c(-5,5)) + scale_y_continuous(limits = c(18, 28)) + 
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", midpoint = 2) + 
  theme(legend.position = "none", plot.title = element_blank(), panel.border = element_blank(), axis.line = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) 

f5 <- FeaturePlot(data9set_cleaned_AD_MG.SO, features = c("Cx3cr1"), pt.size = 0.5) + 
  scale_x_continuous(limits = c(-5,5)) + scale_y_continuous(limits = c(18, 28)) + 
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", midpoint = 2) + 
  theme(legend.position = "none", plot.title = element_blank(), panel.border = element_blank(), axis.line = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) 

f6 <- FeaturePlot(data9set_cleaned_ADp40KO_MG.SO, features = c("Cx3cr1"), pt.size = 0.5) +  
  scale_x_continuous(limits = c(-5,5)) + scale_y_continuous(limits = c(18, 28)) + 
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", midpoint = 2) + 
  theme(plot.title = element_blank(), panel.border = element_blank(), axis.line = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) 

g1 <- arrangeGrob(t1, f1,f2,f3,t2, f4,f5,f6, ncol = 4, widths = c(1, 1.6, 1.6, 2), heights = c(2.4,2))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_DAM_Modified_Homeo_MG_markers.png",
       plot = g1,
       scale = 1, width = 12, height = 6, units = "in", device = "png",
       dpi = 300)


#############################
# DAM markers
#############################

t1 <- textGrob("Cst7", gp=gpar(fontsize=20, col="black", font = "Helvetica"), hjust = 0.35)
t2 <- textGrob("Clec7a", gp=gpar(fontsize=20, col="black", font = "Helvetica"), hjust = 0.35)
t3 <- textGrob("Cd9", gp=gpar(fontsize=20, col="black", font = "Helvetica"), hjust = 0.35)

f1 <- FeaturePlot(data9set_cleaned_Ctrl_MG.SO, features = c("Cst7"), order = TRUE,pt.size = 0.5) + 
  ggtitle(expression(~WT)) + 
  scale_x_continuous(limits = c(-5,5)) + scale_y_continuous(limits = c(18, 28)) + 
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", midpoint = 2) + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.6, size = 20), panel.border = element_blank(), axis.line = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) 

f2 <- FeaturePlot(data9set_cleaned_AD_MG.SO, features = c("Cst7"), order = TRUE,pt.size = 0.5) + 
  ggtitle(expression(~APPPS1p40)) +
  scale_x_continuous(limits = c(-5,5)) + scale_y_continuous(limits = c(18, 28)) + 
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", midpoint = 2) + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.6, size = 20), panel.border = element_blank(), axis.line = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) 

f3 <- FeaturePlot(data9set_cleaned_ADp40KO_MG.SO, features = c("Cst7"), order = TRUE,pt.size = 0.5) +  
  ggtitle(expression(~APPPS1.p40^{"-/-"})) + 
  scale_x_continuous(limits = c(-5,5)) + scale_y_continuous(limits = c(18, 28)) + 
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", midpoint = 2) + 
  theme(plot.title = element_text(hjust = 0.6, face = "bold", size = 20), panel.border = element_blank(), axis.line = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) 

f4 <- FeaturePlot(data9set_cleaned_Ctrl_MG.SO, features = c("Clec7a"), order = TRUE,pt.size = 0.5) + 
  scale_x_continuous(limits = c(-5,5)) + scale_y_continuous(limits = c(18, 28)) + 
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", midpoint = 2) + 
  theme(legend.position = "none", plot.title = element_blank(), panel.border = element_blank(), axis.line = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) 

f5 <- FeaturePlot(data9set_cleaned_AD_MG.SO, features = c("Clec7a"), order = TRUE,pt.size = 0.5) + 
  scale_x_continuous(limits = c(-5,5)) + scale_y_continuous(limits = c(18, 28)) + 
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", midpoint = 2) + 
  theme(legend.position = "none", plot.title = element_blank(), panel.border = element_blank(), axis.line = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) 

f6 <- FeaturePlot(data9set_cleaned_ADp40KO_MG.SO, features = c("Clec7a"), pt.size = 0.5) +  
  scale_x_continuous(limits = c(-5,5)) + scale_y_continuous(limits = c(18, 28)) + 
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", midpoint = 2) + 
  theme(plot.title = element_blank(), panel.border = element_blank(), axis.line = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) 

f7 <- FeaturePlot(data9set_cleaned_Ctrl_MG.SO, features = c("Cd9"), pt.size = 0.5) +
  scale_x_continuous(limits = c(-5,5)) + scale_y_continuous(limits = c(18, 28)) +
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", midpoint = 2) +
  theme(legend.position = "none", plot.title = element_blank(), panel.border = element_blank(), axis.line = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

f8 <- FeaturePlot(data9set_cleaned_AD_MG.SO, features = c("Cd9"), pt.size = 0.5) +
  scale_x_continuous(limits = c(-5,5)) + scale_y_continuous(limits = c(18, 28)) +
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", midpoint = 2) +
  theme(legend.position = "none", plot.title = element_blank(), panel.border = element_blank(), axis.line = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

f9 <- FeaturePlot(data9set_cleaned_ADp40KO_MG.SO, features = c("Cd9"), order = TRUE, pt.size = 0.5) +
  scale_x_continuous(limits = c(-5,5)) + scale_y_continuous(limits = c(18, 28)) +
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", midpoint = 2) +
  theme(plot.title = element_blank(), panel.border = element_blank(), axis.line = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())





g1 <- arrangeGrob(t1, f1,f2,f3,t2, f4,f5,f6, t3, f7, f8, f9, ncol = 4, widths = c(1, 1.6, 1.6, 2), heights = c(2.4,2,2))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_DAM_Modified_DAM_MG_markers.png",
       plot = g1,
       scale = 1, width = 12, height = 9, units = "in", device = "png",
       dpi = 300)

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_DAM_markers_Modified_session_info.txt")



