library(Seurat, lib.loc = "/data/rajewsky/home/skim/R/")
library(ggplot2)
library(dplyr)
library(gridExtra)
source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

load(file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200221_figure1_for_meeting_report_Seurat_object.Robj")


data9set_cleaned.SO$cell_type <- data9set_cleaned.SO@active.ident

meta <- data9set_cleaned.SO@meta.data
meta$replicate <- paste0(meta$sample, "_", meta$gemgroup)

plot.data <- meta %>%
  dplyr::select(sample, replicate = replicate, cell_type = cell_type) %>%
  mutate(cell_type = cell_type) %>%
  group_by(replicate, cell_type) %>%
  summarise(count = n()) %>%
  mutate(sample_total = sum(count)) %>%
  mutate(cell_prop = count / sample_total) 

write.csv(plot.data, file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Ekaterina/cell_frequency_1st.csv")


DimPlot(object = data9set_cleaned.SO, reduction = 'umap', group.by = "seurat_clusters")

meta.df <- data9set_cleaned.SO@meta.data %>% mutate(cell_type = case_when(seurat_clusters %in% c(0, 10) ~ "Dentate_Gyrus",
                                                                          seurat_clusters %in% c(2, 13, 18, 40) ~ "CA1",
                                                                          seurat_clusters %in% c(9, 21, 27, 32) ~ "CA2_3",
                                                                          seurat_clusters %in% c(11, 14, 15, 17, 19, 20) ~ "subiculum",
                                                                          seurat_clusters %in% c(24, 29) ~ "Unidentified_Neurons",
                                                                          seurat_clusters %in% c(7, 16, 22, 30) ~ "Inhibitory_Neurons",
                                                                          seurat_clusters %in% c(6) ~ "MOL",
                                                                          seurat_clusters %in% c(1, 5) ~ "MFOL",
                                                                          seurat_clusters %in% c(12) ~ "OPC",
                                                                          seurat_clusters %in% c(38) ~ "NFOL",
                                                                          seurat_clusters %in% c(3, 8) ~ "Microglia",
                                                                          seurat_clusters %in% c(4) ~ "Astrocytes",
                                                                          seurat_clusters %in% c(36) ~ "Vascular",
                                                                          seurat_clusters %in% c(39) ~ "VLMC",
                                                                          seurat_clusters %in% c(26) ~ "Choroid",
                                                                          seurat_clusters %in% c(23) ~ "Fibroblast",
                                                                          seurat_clusters %in% c(28) ~ "Cajal",
                                                                          seurat_clusters %in% c(35) ~ "Pericyte",
                                                                          seurat_clusters %in% c(37) ~ "Macrophage"))

meta.df$replicate <- paste0(meta.df$sample, "_", meta.df$gemgroup)

plot.data2 <- meta.df %>%
  dplyr::select(sample, replicate = replicate, cell_type = cell_type) %>%
  mutate(cell_type = cell_type) %>%
  group_by(replicate, cell_type) %>%
  summarise(count = n()) %>%
  mutate(sample_total = sum(count)) %>%
  mutate(cell_prop = count / sample_total) 

write.csv(plot.data2, file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Ekaterina/cell_frequency_2nd.csv")
