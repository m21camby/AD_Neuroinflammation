#! /bin/env RScript
# written by SJK at 23. Mar. 2020

.libPaths(c("/home/skim/R/usr_lib", .libPaths()))

library(Seurat)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)
library(tidyr)

source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

load(file="/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200210_figures_for_meeting_report_cleaned_Seurat_object_res0.8.Robj")
data9set_cleaned_MG.SO <- subset(data9set_cleaned.SO, subset = seurat_clusters %in% c(3, 8))

# combine meta data and UMAP
data9set_cleaned.meta.data <- data9set_cleaned.SO@meta.data
data9set_cleaned_MG_umap <- as.data.frame(data9set_cleaned.SO@reductions$umap@cell.embeddings)
data9set_cleaned_MG_umap.df <- cbind(data9set_cleaned.meta.data, data9set_cleaned_MG_umap)

# Mark Microglia cluster
data9set_cleaned_MG_umap.df$MG <- ifelse(data9set_cleaned_MG_umap.df$seurat_clusters %in% c(3, 8), "MG", "others")

# subset Microglia cluster
data9set_cleaned_sub_MG_umap.df <- data9set_cleaned_MG_umap.df[data9set_cleaned_MG_umap.df$MG %in% "MG", ]

# Microglia whole UMAP
g1 <- ggplot(data9set_cleaned_MG_umap.df, aes(x = UMAP_1, y = UMAP_2)) + geom_point(color = "gray", size = 0.1) + 
  geom_point(data = data9set_cleaned_MG_umap.df[data9set_cleaned_MG_umap.df$MG %in% "MG", ], aes(x = UMAP_1, y = UMAP_2), color = "#99CC00", size = 0.1) + 
  theme_classic() + 
  theme(axis.title = element_text(color = "black", size =15),axis.text = element_text(color = "black", size =12))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_Microglia_whole.png",
       plot = g1,
       scale = 1, width = 10, height = 10, units = "in", device = "png",
       dpi = 300)


# Microglia subcluster
g2 <- ggplot(data9set_cleaned_sub_MG_umap.df, aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters)) + geom_point(size = 0.1) + 
  scale_color_manual(values=c("#E69F00", "#56B4E9"), label = paste0("cluster ", c(3, 8))) + 
  theme_classic() + 
  theme(axis.title = element_text(color = "black", size =15),axis.text = element_text(color = "black", size =12), legend.title = element_blank(), legend.position = c(0.8, 0.2)) + scale_x_continuous(limits = c(-5,5)) + scale_y_continuous(limits = c(18, 28))  + 
  guides(color = guide_legend(override.aes = list(size = 6), label.theme = element_text(size = 15)))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_Microglia_sub.png",
       plot = g2,
       scale = 1, width = 7, height = 7, units = "in", device = "png",
       dpi = 300)



##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_session_info.txt")

