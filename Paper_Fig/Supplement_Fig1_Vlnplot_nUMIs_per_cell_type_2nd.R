#! /bin/env RScript
# written by SJK at 22. Oct. 2020

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

load(file="/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200210_figures_for_meeting_report_cleaned_Seurat_object_res0.8.Robj")

# New assigned cell type

meta.df <- data9set_cleaned.SO@meta.data %>% mutate(cell_type = case_when(seurat_clusters %in% c(0, 10) ~ "Dentate_Gyrus",
                                                                          seurat_clusters %in% c(2, 13, 18, 40) ~ "CA1",
                                                                          seurat_clusters %in% c(9, 21, 27, 32) ~ "CA2_3",
                                                                          seurat_clusters %in% c(11, 14, 15, 17, 19, 20) ~ "subiculum",
                                                                          seurat_clusters %in% c(24, 29) ~ "Unidentified_Neurons",
                                                                          seurat_clusters %in% c(7, 16, 22, 30) ~ "Inhibitory_Neurons",
                                                                          seurat_clusters %in% c(6) ~ "MOL",
                                                                          seurat_clusters %in% c(1, 5) ~ "MFOL",
                                                                          seurat_clusters %in% c(12) ~ "NFOL",
                                                                          seurat_clusters %in% c(38) ~ "OPC",
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
data9set_cleaned.SO$cell_type <- factor(data9set_cleaned.SO$cell_type, levels = c("Dentate_Gyrus", "CA1", "CA2_3", "subiculum",  "Unidentified_Neurons", 
                                                                                  "Inhibitory_Neurons", "MOL", "MFOL", "OPC", "NFOL", "Microglia","Astrocytes",
                                                                                  "Vascular", "VLMC", "Choroid", "Fibroblast", "Cajal", "Pericyte", "Macrophage"))

# Again with big cell types

data9set_cleaned.SO$cell_types <- ifelse(data9set_cleaned.SO$cell_type %in% c("CA1", "CA2_3", "Dentate_Gyrus","Unidentified_Neurons", "subiculum"), "Excitatory Neurons",
                                         ifelse(data9set_cleaned.SO$cell_type %in% "Inhibitory_Neurons", "Inhibitory Neurons",
                                                ifelse(data9set_cleaned.SO$cell_type %in% c("Microglia"), "Microglia",
                                                       ifelse(data9set_cleaned.SO$cell_type %in% c("Macrophage"), "Macrophages",
                                                          ifelse(data9set_cleaned.SO$cell_type %in% c("MOL", "MFOL", "NFOL"), "Oligodendrocytes",
                                                              ifelse(data9set_cleaned.SO$cell_type %in% c("OPC"), "OPC",
                                                                     ifelse(data9set_cleaned.SO$cell_type %in% c("Astrocytes"), "Asctrocytes",
                                                                            ifelse(data9set_cleaned.SO$cell_type %in% c("Cajal"), "Cajal Retzius",
                                                                                   ifelse(data9set_cleaned.SO$cell_type %in% c("Fibroblast"), "Fibroblast",
                                                                                          ifelse(data9set_cleaned.SO$cell_type %in% c("Vascular", "VLMC", "Pericyte"), "Vascular cells", "Choroid Plexus"))))))))))



data9set_cleaned.SO$cell_types <- factor(data9set_cleaned.SO$cell_types, levels = c("Excitatory Neurons", "Inhibitory Neurons", "Cajal Retzius", "Choroid Plexus", "Asctrocytes", "Microglia", "Macrophages", "Oligodendrocytes", "OPC", "Fibroblast", "Vascular cells"))



# meta data
data9set_cleaned.meta <- data9set_cleaned.SO@meta.data %>% as.data.frame

v1 <- VlnPlot(data9set_cleaned.SO, features = "nCount_RNA", group.by = "cell_types", pt.size = 0, cols = c("#CC0000", "#FF6600", "#FF9900","#CC9900", "#FFCC66", "#99CC00","#003300", "#99CCCC", "#009999", "#3399FF","#000033")) + 
  ggtitle("nUMIs per cell by cell type") + ylim(c(0, 31000)) +
  theme(legend.position = "none", axis.title.x = element_blank()) + 
  annotate("text", x = 1, y = 30500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Excitatory Neurons", ]$nCount_RNA)) +
  annotate("text", x = 2, y = 30500, label = round(median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Inhibitory Neurons", ]$nCount_RNA))) +
  annotate("text", x = 3, y = 30500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Cajal Retzius", ]$nCount_RNA)) +
  annotate("text", x = 4, y = 30500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Choroid Plexus", ]$nCount_RNA)) +
  annotate("text", x = 5, y = 30500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Asctrocytes", ]$nCount_RNA)) +
  annotate("text", x = 6, y = 30500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Microglia", ]$nCount_RNA)) +
  annotate("text", x = 7, y = 30500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Macrophages", ]$nCount_RNA)) +
  annotate("text", x = 8, y = 30500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Oligodendrocytes", ]$nCount_RNA)) +
  annotate("text", x = 9, y = 30500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "OPC", ]$nCount_RNA)) +
  annotate("text", x = 10, y = 30500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Fibroblast", ]$nCount_RNA)) +
  annotate("text", x = 11, y = 30500, label = round(median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Vascular cells", ]$nCount_RNA)))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig1_Vlnplot_nUMIs_per_cell_type_2nd.png", 
       plot = v1, 
       scale = 1, width = 9, height = 6, units = "in", device = "png",
       dpi = 300)

v1 <- VlnPlot(data9set_cleaned.SO, features = "nFeature_RNA", group.by = "cell_types", pt.size = 0, cols = c("#CC0000", "#FF6600", "#FF9900","#CC9900", "#FFCC66", "#99CC00","#003300", "#99CCCC", "#009999", "#3399FF","#000033")) + 
  ggtitle("nGenes per cell by cell type") + ylim(c(0, 8000)) +
  theme(legend.position = "none", axis.title.x = element_blank()) + 
  annotate("text", x = 1, y = 7500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Excitatory Neurons", ]$nFeature_RNA)) +
  annotate("text", x = 2, y = 7500, label = round(median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Inhibitory Neurons", ]$nFeature_RNA))) +
  annotate("text", x = 3, y = 7500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Cajal Retzius", ]$nFeature_RNA)) +
  annotate("text", x = 4, y = 7500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Choroid Plexus", ]$nFeature_RNA)) +
  annotate("text", x = 5, y = 7500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Asctrocytes", ]$nFeature_RNA)) +
  annotate("text", x = 6, y = 7500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Microglia", ]$nFeature_RNA)) +
  annotate("text", x = 7, y = 7500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Macrophages", ]$nFeature_RNA)) +
  annotate("text", x = 8, y = 7500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Oligodendrocytes", ]$nFeature_RNA)) +
  annotate("text", x = 9, y = 7500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "OPC", ]$nFeature_RNA)) +
  annotate("text", x = 10, y = 7500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Fibroblast", ]$nFeature_RNA)) +
  annotate("text", x = 11, y = 7500, label = round(median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Vascular cells", ]$nFeature_RNA)))


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig1_Vlnplot_nGenes_per_cell_type_2nd.png", 
       plot = v1, 
       scale = 1, width = 9, height = 6, units = "in", device = "png",
       dpi = 300)

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig1_Vlnplot_nUMIs_per_cell_type_2nd_session_info.txt")
