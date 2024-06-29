

library(Seurat, lib.loc = "/data/rajewsky/home/skim/R/")
library(ggplot2)
library(dplyr)
library(tidyverse)
source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

# In this plot, I remove zero counts and only retain expressing cell and show by violin plot

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
data9set_cleaned.SO$cell_type <- factor(data9set_cleaned.SO$cell_type, levels = c("Dentate_Gyrus", "CA1", "CA2_3", "subiculum",  "Unidentified_Neurons", 
                                                                                  "Inhibitory_Neurons", "MOL", "MFOL", "OPC", "NFOL", "Microglia","Astrocytes",
                                                                                  "Vascular", "VLMC", "Choroid", "Fibroblast", "Cajal", "Pericyte", "Macrophage"))



meta.df <- meta.df %>% mutate(large_cell_type = case_when(cell_type %in% c("Dentate_Gyrus", "CA1", "CA2_3","subiculum","Unidentified_Neurons") ~ "Excitatory Neurons",
                                                          cell_type %in% c("Inhibitory_Neurons") ~ "Inhibitory Neurons",
                                                          cell_type %in% c("MOL", "MFOL","NFOL") ~ "Oligodendrocytes",
                                                          cell_type %in% c("Astrocytes") ~ "Astrocytes",
                                                          cell_type %in% c("Microglia") ~ "Microglia",
                                                          cell_type %in% c("OPC") ~ "OPC",
                                                          cell_type %in% c("Vascular", "VLMC","Choroid","Fibroblast","Cajal","Pericyte","Macrophage") ~ "Rest"))


meta.df <- transform(meta.df, pyscenic1 = paste0(cell_type, "_",sample))
meta.df <- transform(meta.df, pyscenic2 = paste0(cell_type, "_",sample, "_", gemgroup))
meta.df <- transform(meta.df, pyscenic3 = paste0(large_cell_type, "_",sample))
meta.df <- transform(meta.df, pyscenic4 = paste0(large_cell_type, "_",sample, "_", gemgroup))

write.csv(meta.df, file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/pyscenic/20220315_pyscenic_meta_Seurat_meta.csv")


