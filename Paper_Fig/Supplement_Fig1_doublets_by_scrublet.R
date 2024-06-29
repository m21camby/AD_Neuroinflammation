#! /bin/env RScript
# written by SJK at 4. Mar. 2020

.libPaths(c("/data/murphy/shared_libs/R", "/usr/local/lib/R/site-librar","/usr/local/lib/R/library","/usr/lib64/R/library","/usr/share/R/library", .libPaths() ))

library(Seurat)
library(ggplot2)
library(dplyr)
library(gridExtra)
source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

load("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200117_Clustering_DE_by_MAST_Seurat_object.Robj")

data9set.SO@meta.data$sample <- ifelse((data9set.SO$gemgroup == 1 | data9set.SO$gemgroup == 4 | data9set.SO$gemgroup == 7), "Ctrl", 
                                       ifelse((data9set.SO$gemgroup == 2 | data9set.SO$gemgroup == 5 | data9set.SO$gemgroup == 8), "ADp40KO", "AD"))



data9set.SO@meta.data$sample <- factor(data9set.SO@meta.data$sample, levels = c("Ctrl", "AD", "ADp40KO"))

scrublet_predicted.df <- read.table("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200205_doublets_by_scrublet_predicted_doublets_data.frame.txt", sep = "\t")

data9set.df <- as.data.frame(data9set.SO@meta.data)
data9set.df$cell_ID <- rownames(data9set.df)
data9set_scrublet.df <- merge(data9set.df, scrublet_predicted.df, by = "cell_ID")
data9set_scrublet.df$scrublet_prediction <- ifelse(data9set_scrublet.df$scrublet_predicted == 1, "doublets", "singlet")

data9set_umap <- as.data.frame(data9set.SO@reductions$umap@cell.embeddings)
data9set_umap$cell_ID <- rownames(data9set_umap)

data9set_umap <- merge(data9set_umap, data9set_scrublet.df, by = "cell_ID")



g1 <- ggplot(data9set_umap) + geom_point(aes(x = UMAP_1, y = UMAP_2, color = scrublet_prediction), size = 0.5) + 
  theme_classic() + ggtitle("doublets prediction by scrublet") +
  scale_color_manual(values=c("#FF9900", "#56B4E9")) + 
  guides(colour = guide_legend(override.aes = list(size=7))) + 
  geom_segment(aes(x = -16.5, y = 12, xend = -16.5, yend = 9), arrow = arrow(length=unit(0.3,'cm')),color='black',size=0.5) + 
  geom_segment(aes(x = -7, y = 17, xend = -5, yend = 15), arrow = arrow(length=unit(0.3,'cm')),color='black',size=0.5) + 
  geom_segment(aes(x = 8.6, y = 23.3, xend = 6.8, yend = 20.7), arrow = arrow(length=unit(0.3,'cm')),color='black',size=0.5) +
  geom_segment(aes(x = 3, y = -20, xend = 3, yend = -23), arrow = arrow(length=unit(0.3,'cm')),color='black',size=0.5) +
  geom_segment(aes(x = 10, y = -19, xend = 10, yend = -22), arrow = arrow(length=unit(0.3,'cm')),color='black',size=0.5) + 
  geom_segment(aes(x = 17, y = -15, xend = 19, yend = -17), arrow = arrow(length=unit(0.3,'cm')),color='black',size=0.5) +
  geom_segment(aes(x = 21, y = -12, xend = 21.3, yend = -15), arrow = arrow(length=unit(0.3,'cm')),color='black',size=0.5) +
  geom_segment(aes(x = 24.5, y = -13, xend = 22.7, yend = -15), arrow = arrow(length=unit(0.3,'cm')),color='black',size=0.5) +
  theme(legend.text = element_text(size = 15), legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 20, color = "black"), axis.title = element_blank(),
      axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank())
  
ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig1_doublets_by_scrublet.png", 
       plot = g1, 
       scale = 1, width = 8, height = 6, units = "in", device = "png",
       dpi = 300)

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig1_doublets_by_scrublet.png_session_info.txt")

