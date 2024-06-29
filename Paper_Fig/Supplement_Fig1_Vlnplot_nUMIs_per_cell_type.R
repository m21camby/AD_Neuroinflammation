#! /bin/env RScript
# written by SJK at 11. May. 2020


.libPaths(c("/home/skim/R/usr_lib", .libPaths()))


library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

load(file="/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200210_figures_for_meeting_report_cleaned_Seurat_object_res0.8.Robj")

data9set_cleaned.SO$cell_type <- data9set_cleaned.SO@active.ident

data9set_cleaned.SO$cell_types <- ifelse(data9set_cleaned.SO$cell_type %in% c("CA1 Neurons", "CA2/CA3 Neuron", "Dentate Gyrus","Neurons", "Subiculum"), "Excitatory Neurons",
                                      ifelse(data9set_cleaned.SO$cell_type %in% "Inhibitory Interneurons", "Inhibitory Interneurons",
                                             ifelse(data9set_cleaned.SO$cell_type %in% c("Microglia"), "Microglia",
                                                    ifelse(data9set_cleaned.SO$cell_type %in% c("Oligo"), "Oligodendrocytes",
                                                           ifelse(data9set_cleaned.SO$cell_type %in% c("OPC"), "OPC",
                                                                  ifelse(data9set_cleaned.SO$cell_type %in% c("Astrocytes"), "Asctrocytes",
                                                                         ifelse(data9set_cleaned.SO$cell_type %in% c("Cajal Retzius"), "Cajal Retzius",
                                                                                ifelse(data9set_cleaned.SO$cell_type %in% c("Fibroblast"), "Fibroblast",
                                                                                       ifelse(data9set_cleaned.SO$cell_type %in% c("Vascular Endothelial", "VLMC", "Pericytes"), "Vascular cells","Choroid Plexus")))))))))


data9set_cleaned.SO$cell_types <- ifelse(data9set_cleaned.SO$seurat_clusters %in% c(37), "Macrophages", data9set_cleaned.SO$cell_types)

data9set_cleaned.SO$cell_types <- factor(data9set_cleaned.SO$cell_types, levels = c("Excitatory Neurons", "Inhibitory Interneurons", "Cajal Retzius", "Choroid Plexus", "Asctrocytes", "Microglia", "Macrophages", "Oligodendrocytes", "OPC", "Fibroblast", "Vascular cells"))

# meta data
data9set_cleaned.meta <- data9set_cleaned.SO@meta.data %>% as.data.frame

median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Excitatory Neurons", ]$nCount_RNA)


v1 <- VlnPlot(data9set_cleaned.SO, features = "nCount_RNA", group.by = "cell_types", pt.size = 0, cols = c("#CC0000", "#FF6600", "#FF9900","#CC9900", "#FFCC66", "#99CC00","#003300", "#99CCCC", "#009999", "#3399FF","#000033")) + 
  ggtitle("nUMIs per cell by cell type") + ylim(c(0, 31000)) +
  theme(legend.position = "none", axis.title.x = element_blank()) + 
  annotate("text", x = 1, y = 30500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Excitatory Neurons", ]$nCount_RNA)) +
  annotate("text", x = 2, y = 30500, label = round(median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Inhibitory Interneurons", ]$nCount_RNA))) +
  annotate("text", x = 3, y = 30500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Cajal Retzius", ]$nCount_RNA)) +
  annotate("text", x = 4, y = 30500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Choroid Plexus", ]$nCount_RNA)) +
  annotate("text", x = 5, y = 30500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Asctrocytes", ]$nCount_RNA)) +
  annotate("text", x = 6, y = 30500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Microglia", ]$nCount_RNA)) +
  annotate("text", x = 7, y = 30500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Macrophages", ]$nCount_RNA)) +
  annotate("text", x = 8, y = 30500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Oligodendrocytes", ]$nCount_RNA)) +
  annotate("text", x = 9, y = 30500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "OPC", ]$nCount_RNA)) +
  annotate("text", x = 10, y = 30500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Fibroblast", ]$nCount_RNA)) +
  annotate("text", x = 11, y = 30500, label = round(median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Vascular cells", ]$nCount_RNA)))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig1_Vlnplot_nUMIs_per_cell_type.png", 
       plot = v1, 
  scale = 1, width = 9, height = 6, units = "in", device = "png",
  dpi = 300)

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig1_Vlnplot_nUMIs_per_cell_type_session_info.txt")



