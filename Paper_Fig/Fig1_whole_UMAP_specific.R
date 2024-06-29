#! /bin/env RScript
# written by SJK at 14. Dec. 2020
# UMAP with more specific neuronal cell type

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

# In this plot, I remove zero counts and only retain expressing cell and show by violin plot

load(file="/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200210_figures_for_meeting_report_cleaned_Seurat_object_res0.8.Robj")

# New assigned cell type

meta.df <- data9set_cleaned.SO@meta.data %>% mutate(cell_type = case_when(seurat_clusters %in% c(0, 10) ~ "Dentate Gyrus",
                                                                          seurat_clusters %in% c(2, 13, 18, 40) ~ "CA1",
                                                                          seurat_clusters %in% c(9, 21, 27, 32) ~ "CA2/3",
                                                                          seurat_clusters %in% c(11, 14, 15, 17, 19, 20) ~ "subiculum",
                                                                          seurat_clusters %in% c(24, 29) ~ "Unidentified\nNeurons",
                                                                          seurat_clusters %in% c(7, 16, 22, 30) ~ "Inhibitory Neurons",
                                                                          seurat_clusters %in% c(6) ~ "MOL",
                                                                          seurat_clusters %in% c(1, 5) ~ "MFOL",
                                                                          seurat_clusters %in% c(38) ~ "NFOL",
                                                                          seurat_clusters %in% c(12) ~ "OPC",
                                                                          seurat_clusters %in% c(3, 8) ~ "Microglia",
                                                                          seurat_clusters %in% c(4) ~ "Astrocytes",
                                                                          seurat_clusters %in% c(36) ~ "Vascular",
                                                                          seurat_clusters %in% c(39) ~ "VLMC",
                                                                          seurat_clusters %in% c(26) ~ "Choroid",
                                                                          seurat_clusters %in% c(23) ~ "Fibroblast",
                                                                          seurat_clusters %in% c(28) ~ "Cajal",
                                                                          seurat_clusters %in% c(35) ~ "Pericyte",
                                                                          seurat_clusters %in% c(37) ~ "Macrophage"))


data9set_cleaned.SO$cell_type <- meta.df$cell_type

# give levels to large cell type
data9set_cleaned.SO$cell_type <- factor(data9set_cleaned.SO$cell_type, levels = c("Dentate Gyrus", "CA1", "CA2/3", "subiculum",  "Unidentified\nNeurons", 
                                                                                  "Inhibitory Neurons", "MOL", "MFOL", "OPC", "NFOL", "Microglia","Astrocytes",
                                                                                  "Vascular", "VLMC", "Choroid", "Fibroblast", "Cajal", "Pericyte", "Macrophage"))

data9set_cleaned.SO@active.ident <- data9set_cleaned.SO$cell_type

levels(data9set_cleaned.SO@active.ident)

DimPlot(data9set_cleaned.SO, reduction = "umap", 
        label = TRUE, label.size = 6)

d1 <- DimPlot(data9set_cleaned.SO, reduction = "umap", 
              label = TRUE, label.size = 6, 
              cols = c("#CC0000", "#FF0000", "#990000", "#660000", "#FF3366", 
                       "#FF6600",
                       "#006666", "#669999","#99CCCC","#9966FF", 
                       "#99CC00","#FFCC66", "#000033", "#000033","#CC9900", "#3399FF", "#FF9900", "#3399FF","#000033")) + 
  theme(legend.position = "none") 
d1

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_whole_UMAP_specific.pdf", 
       plot = d1, 
       scale = 1, width = 10, height = 10, units = "in", device = cairo_pdf,
       dpi = 300)

# New color
# EX CC0000
# IN FF6600
#  "MOL" (006666) "MFOL" (669999)   "OPC" (99CCCC) "NFOL" (9966FF)
# [11] "Microglia"(99CC00)   "Astrocytes" (FFCC66)   "Vascular" (000033)   "VLMC"  (000033)   "Choroid" (CC9900)           
# [16] "Fibroblast" (3399FF)  "Cajal" (FF9900)   "Pericyte" (3399FF)   "Macrophage" (003300)
# 
# Old color
# "Cajal Retzius"(FF9900)           "Choroid Plexus" (CC9900)        
# [5] "Astrocytes" (FFCC66)             "Microglia" (99CC00)              "Macrophages" (003300)            "Oligodendrocytes"(99CCCC)       
# [9] "OPC" (009999)                    "Fibroblast"  (3399FF)            "Vascular cells" (000033)
#  
# "Cajal Retzius"           "Astrocytes"             
# [9] "Microglia"               "Oligo"                   "OPC"                     "Fibroblast"             
# [13] "Vascular Endothelial"    "VLMC"                    "Pericytes"               "Choroid Plexus"
# 

d2 <- DimPlot(data9set_cleaned.SO, reduction = "umap", 
              label = FALSE, 
              cols = c("#CC0000", "#FF0000", "#990000", "#660000", "#FF3366", 
                       "#FF6600",
                       "#006666", "#669999","#99CCCC","#9966FF", 
                       "#99CC00","#FFCC66", "#000033", "#000033","#CC9900", "#3399FF", "#FF9900", "#3399FF","#000033")) + 
  theme(legend.position = "none") 

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_whole_UMAP_specific_no_label.pdf", 
       plot = d2, 
       scale = 1, width = 10, height = 10, units = "in", device = cairo_pdf,
       dpi = 300)

# Astrocytes
d3 <- DimPlot(data9set_cleaned.SO, reduction = "umap", 
              label = FALSE, 
              cols = c("#CCCCCC", "#CCCCCC", "#CCCCCC", "#CCCCCC", "#CCCCCC", 
                       "#CCCCCC",
                       "#CCCCCC", "#CCCCCC","#CCCCCC","#CCCCCC", 
                       "#CCCCCC","#FFCC66", "#CCCCCC", "#CCCCCC","#CCCCCC", "#CCCCCC", "#CCCCCC", "#CCCCCC","#CCCCCC")) + 
  theme(legend.position = "none", 
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) 
d3
ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_whole_UMAP_specific_no_label_Astrocytes.pdf", 
       plot = d3, 
       scale = 1, width = 10, height = 10, units = "in", device = cairo_pdf,
       dpi = 300)


