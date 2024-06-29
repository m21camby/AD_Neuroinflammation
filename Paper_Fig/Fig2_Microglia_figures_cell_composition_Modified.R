#! /bin/env RScript
# written by SJK at 4. Ma7. 2020

# This file is for not colored by genotype but same color as Fig2(A) 

.libPaths(c("/home/skim/R/usr_lib", .libPaths()))

library(Seurat)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)
library(tidyr)
library(grid)

source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

load(file="/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200210_figures_for_meeting_report_cleaned_Seurat_object_res0.8.Robj")
data9set_cleaned_MG.SO <- subset(data9set_cleaned.SO, subset = seurat_clusters %in% c(3, 8))

#########################
# distribution by sample
#########################
MG_meta <- as.data.frame(data9set_cleaned_MG.SO@meta.data)
MG_umap <- as.data.frame(data9set_cleaned_MG.SO@reductions$umap@cell.embeddings)
MG_merge <- cbind(MG_meta, MG_umap)

MG_merge_Ctrl <- MG_merge[MG_merge$sample %in% "Ctrl", ]

g1 <- ggplot(MG_merge_Ctrl, aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters)) + geom_point(size = 0.1) + 
  theme_classic() + scale_color_manual(values=c("#E69F00", "#56B4E9")) + 
  ggtitle("WT") + 
  scale_x_continuous(limits = c(-5,5)) + 
  scale_y_continuous(limits = c(18, 28)) + 
  theme(title = element_text(size = 18, family = "Helvetica"), legend.position = "none", plot.title = element_text(hjust = 0.6), panel.border = element_blank(), axis.line = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

MG_merge_AD <- MG_merge[MG_merge$sample %in% "AD", ]

g2 <- ggplot(MG_merge_AD, aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters)) + geom_point(size = 0.1) + 
  theme_classic() + scale_color_manual(values=c("#E69F00", "#56B4E9")) + 
  ggtitle("APPPS1") + 
  scale_x_continuous(limits = c(-5,5)) + 
  scale_y_continuous(limits = c(18, 28)) + 
  theme(title = element_text(size = 18, family = "Helvetica"), legend.position = "none", plot.title = element_text(hjust = 0.6), panel.border = element_blank(), axis.line = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

MG_merge_ADp40KO <- MG_merge[MG_merge$sample %in% "ADp40KO", ]

g3 <- ggplot(MG_merge_ADp40KO, aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters)) + geom_point(size = 0.1) + 
  theme_classic() + scale_color_manual(values=c("#E69F00", "#56B4E9")) + 
  ggtitle(expression(~APPPS1.p40^{"-/-"})) + 
  scale_x_continuous(limits = c(-5,5)) + 
  scale_y_continuous(limits = c(18, 28)) + 
  theme(title = element_text(size = 18, family = "Helvetica"), legend.position = "none", plot.title = element_text(hjust = 0.6), panel.border = element_blank(), axis.line = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

d1 <- arrangeGrob(g1, g2, g3,  ncol = 3)



ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_cell_composition_Modified_figure1.png",
       plot = d1,
       scale = 1, width = 12, height = 4, units = "in", device = "png",
       dpi = 300)



##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_cell_composition_Modified_session_info.txt")
