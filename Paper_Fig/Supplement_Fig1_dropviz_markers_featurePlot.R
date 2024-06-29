#! /bin/env RScript
# written by SJK at 22. Oct. 2020

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(gridExtra)

source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

load(file="/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200210_figures_for_meeting_report_cleaned_Seurat_object_res0.8.Robj")

markers <- list.files("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DropViz_cluster_Markers/")

dropviz_featureplot <- function(markers = "Astrocytes_Gja11", titletext = "Astrocytes dropviz markers"){
  FeaturePlot(object = data9set_cleaned.SO, features = markers, cols = c("lightgray", "darkred"), order = TRUE) + 
    ggtitle(titletext) + theme_void() + theme(plot.title = element_text(size = 12, family = "helvetica", color = "black", hjust = 0.5),
                                                                 legend.text = element_text(size = 12, family = "helvetica", color = "black", hjust = 0.5))
}

genes <- read_csv(paste0("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DropViz_cluster_Markers/", markers[13]))


for(i in 1:17){
  genes <- read_csv(paste0("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DropViz_cluster_Markers/", markers[i]))
  genes_list <- genes$gene  
  sample <- str_remove(markers[i], "cluster-markers_")
  sample <- str_remove(sample, ".csv")
  data9set_cleaned.SO <- AddModuleScore(data9set_cleaned.SO, features =  list(genes_list), name = sample)
}


g1 <- dropviz_featureplot(markers = "Neurogenesis_Sox41", titletext = "Neurogenesis markers")
g2 <- dropviz_featureplot(markers = "Neuron_Dentate_C1ql21", titletext = "Dentate Gyrus markers")
g3 <- dropviz_featureplot(markers = "Neuron_CA2CA3_Pvrl31", titletext = "CA2/3 markers")
g4 <- dropviz_featureplot(markers = "Neuron_CA1_Subiculum_Postsubiculum_Entorhinal1", titletext = "CA1, subiculum, postsubiculum, Entorhinal markers")
g5 <- dropviz_featureplot(markers = "Subiculum_Slc17a61", titletext = "subiculum markers")
g6 <- dropviz_featureplot(markers = "Subiculum_Entorhinal_Nxph31", titletext = "subiculum, Entorhinal markers")
g7 <- dropviz_featureplot(markers = "Interneuron_Gad21", titletext = "Interneuron markers")
g8 <- dropviz_featureplot(markers = "Neuron_CajalRetzius_Lhx11", titletext = "Cajal Retzius  markers")
g9 <- dropviz_featureplot(markers = "Choroid_Plexus_Ttr1", titletext = "Choroid_Plexus markers")
g10 <- dropviz_featureplot(markers = "Oligodendrocyte_Tfr1", titletext = "Oligodendrocyte markers")
g11 <- dropviz_featureplot(markers = "Polydendrocyte_Tnr1", titletext = "Polydendrocyte markers")
g12 <- dropviz_featureplot(markers = "Microglia_Macrophage_C1qb1", titletext = "Microglia, Macrophage markers")
g13 <- dropviz_featureplot(markers = "Astrocytes_Gja11", titletext = "Astrocytes  markers")
g14 <- dropviz_featureplot(markers = "Endothelial_Flt11", titletext = "Endothelial markers")
g15 <- dropviz_featureplot(markers = "Ependyma1", titletext = "Ependyma markers")
g16 <- dropviz_featureplot(markers = "Fibroblast_Dcn1", titletext = "Fibroblast markers")
g17 <- dropviz_featureplot(markers = "Mural_Rgs5Acta21", titletext = "Mural markers")


v1 <- arrangeGrob(g1, g2, g3, g4, g5, g6, g7, g8, ncol = 2)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig1_dropviz_markers_featurePlot_1st.png", 
       plot = v1, 
       scale = 1, width = 9, height = 15, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig1_dropviz_markers_featurePlot_1st.pdf", 
       plot = v1, 
       scale = 1, width = 9, height = 15, units = "in", device = cairo_pdf,
       dpi = 300)

v1 <- arrangeGrob(g9, g10, g11, g12, g13, g14, g15, g16, ncol = 2)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig1_dropviz_markers_featurePlot_2nd.png", 
       plot = v1, 
       scale = 1, width = 9, height = 15, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig1_dropviz_markers_featurePlot_2nd.pdf", 
       plot = v1, 
       scale = 1, width = 9, height = 15, units = "in", device = cairo_pdf,
       dpi = 300)


##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig1_dropviz_markers_featurePlot_session_info.txt")

