#! /bin/env RScript
# written by SJK at 2. Mar. 2020
# modified at 04. Nov. 2020
# : libpaths removed and pdf and svg files created

#.libPaths(c("/home/skim/R/usr_lib", .libPaths()))

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

load(file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200221_figure1_for_meeting_report_Seurat_object.Robj")

d1 <- DimPlot(data9set_cleaned.SO, reduction = "umap", 
        label = TRUE, label.size = 6, 
        cols = c("#CC0000", "#FF6600", "#FF9900","#CC9900", "#FFCC66", "#99CC00","#003300", "#99CCCC", "#009999", "#3399FF","#000033")) + 
  theme(legend.position = "none") 

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_whole_UMAP.png", 
       plot = d1, 
  scale = 1, width = 10, height = 10, units = "in", device = "png",
  dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_whole_UMAP.pdf", 
       plot = d1, 
       scale = 1, width = 10, height = 10, units = "in", device = cairo_pdf,
       dpi = 300)
ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_whole_UMAP.svg", 
       plot = d1, 
       scale = 1, width = 10, height = 10, units = "in", device = svg,
       dpi = 300)

d2 <- DimPlot(data9set_cleaned.SO, reduction = "umap",
        cols = c("#CC0000", "#FF6600", "#FF9900","#CC9900", "#FFCC66", "#99CC00","#003300", "#99CCCC", "#009999", "#3399FF","#000033")) + theme(legend.position = "none")


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_whole_UMAP_No_label.png",
       plot = d2,
  scale = 1, width = 10, height = 10, units = "in", device = "png",
  dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_whole_UMAP_No_label.pdf",
       plot = d2,
       scale = 1, width = 10, height = 10, units = "in", device = cairo_pdf,
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_whole_UMAP_No_label.svg",
       plot = d2,
       scale = 1, width = 10, height = 10, units = "in", device = "svg",
       dpi = 300)


##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_whole_UMAP_session_info.txt")




