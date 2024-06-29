

library(Seurat, lib.loc = "/data/rajewsky/home/skim/R/usr_lib/")
library(ggplot2)
library(dplyr)
library(tidyverse)

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

meta.df$cell_ID <- names(data9set_cleaned.SO@active.ident)

#DimPlot(data9set_cleaned.SO, reduction = "umap", label = TRUE)

write.csv(meta.df, file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/pyscenic/20220310_pyscenic_pre_data_processing_meta_Seurat_pyscenic.csv")

write.table(meta.df$cell_ID, file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/pyscenic/20220310_pyscenic_pre_data_processing_meta_Seurat_pyscenic_CB.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
