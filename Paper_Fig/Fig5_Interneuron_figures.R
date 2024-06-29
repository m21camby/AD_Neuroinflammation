#! /bin/env RScript
# written by SJK at 2. April. 2020
# This file is for creating Inhibitory UMAP

.libPaths(c("/home/skim/R/usr_lib", .libPaths()))

library(Seurat)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)
library(tidyr)

source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

load(file="/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200210_figures_for_meeting_report_cleaned_Seurat_object_res0.8.Robj")
data9set_cleaned.SO@meta.data$sample <- ifelse((data9set_cleaned.SO$gemgroup == 1 | data9set_cleaned.SO$gemgroup == 4 | data9set_cleaned.SO$gemgroup == 7), 
                                               "Ctrl", ifelse((data9set_cleaned.SO$gemgroup == 2 | data9set_cleaned.SO$gemgroup == 5 | data9set_cleaned.SO$gemgroup == 8), "ADp40KO", "AD"))



data9set_IN.SO <- subset(data9set_cleaned.SO, subset = seurat_clusters %in% c(7, 16, 22, 30))

# combine meta data and UMAP
data9set.meta.data <- data9set_cleaned.SO@meta.data
data9set_IN_umap <- as.data.frame(data9set_cleaned.SO@reductions$umap@cell.embeddings)
data9set_IN_umap.df <- cbind(data9set.meta.data, data9set_IN_umap)

# Mark Inhibitory cluster
data9set_IN_umap.df$celltypes <- ifelse(data9set_IN_umap.df$seurat_clusters %in% c(7, 16, 22, 30), "Inhibitory Interneurons", "others")

# subset Inhibitory cluster
data9set_IN_umap_sub.df <- data9set_IN_umap.df[data9set_IN_umap.df$celltypes %in% "Inhibitory Interneurons", ]


# Inhibitory whole UMAP
g1 <- ggplot(data9set_IN_umap.df, aes(x = UMAP_1, y = UMAP_2)) + geom_point(color = "gray", size = 0.1) + 
  geom_point(data = data9set_IN_umap_sub.df[data9set_IN_umap_sub.df$celltypes %in% "Inhibitory Interneurons", ], aes(x = UMAP_1, y = UMAP_2), color = "#FF6600", size = 0.1) +
  theme_classic() + 
  theme(axis.title = element_text(color = "black", size =15),axis.text = element_text(color = "black", size =12))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Interneuron_figures_whole.png",
       plot = g1,
       scale = 1, width = 10, height = 10, units = "in", device = "png",
       dpi = 300)

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Interneuron_figures_session_info.txt")