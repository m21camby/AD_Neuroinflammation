#! /bin/env RScript
# written by SJK at 2. April. 2020
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

# Dentate gyrus: 0, 10
# CA1 Neurons: 2, 13, 18
# CA2/CA3 Neuron: 9, 21, 27
# Subiculum: 11, 14, 15, 17, 19, 20
# Neurons: 24, 29, 32, 40

data9set_EX.SO <- subset(data9set_cleaned.SO, subset = seurat_clusters %in% c(0, 10, 2, 13, 18, 9, 21, 27, 11, 14, 15, 17, 19, 20, 24, 29, 32, 40))

# combine meta data and UMAP
data9set.meta.data <- data9set_cleaned.SO@meta.data
data9set_EX_umap <- as.data.frame(data9set_cleaned.SO@reductions$umap@cell.embeddings)
data9set_EX_umap.df <- cbind(data9set.meta.data, data9set_EX_umap)

# Mark Excitatory cluster
data9set_EX_umap.df$celltypes <- ifelse(data9set_EX_umap.df$seurat_clusters %in% c(0, 10, 2, 13, 18, 9, 21, 27, 11, 14, 15, 17, 19, 20, 24, 29, 32, 40), "Excitatory Neurons", "others")

# subset Excitatory cluster
data9set_EX_umap_sub.df <- data9set_EX_umap.df[data9set_EX_umap.df$celltypes %in% "Excitatory Neurons", ]


# Excitatory whole UMAP
g1 <- ggplot(data9set_EX_umap.df, aes(x = UMAP_1, y = UMAP_2)) + geom_point(color = "gray", size = 0.1) + 
  geom_point(data = data9set_EX_umap_sub.df[data9set_EX_umap_sub.df$celltypes %in% "Excitatory Neurons", ], aes(x = UMAP_1, y = UMAP_2), color = "#CC0000", size = 0.1) +
  theme_classic() + 
  theme(axis.title = element_text(color = "black", size =15),axis.text = element_text(color = "black", size =12))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Excitatory_neuron_figures_whole.png",
       plot = g1,
       scale = 1, width = 10, height = 10, units = "in", device = "png",
       dpi = 300)

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Excitatory_neuron_figures_session_info.txt")