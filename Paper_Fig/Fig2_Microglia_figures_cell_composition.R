#! /bin/env RScript
# written by SJK at 23. Mar. 2020

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

d1 <- DimPlot(object = data9set_cleaned_MG.SO, group.by = "sample", pt.size = 0.05, split.by = "sample") 
d1 <- d1 + scale_colour_viridis_d() + 
  scale_x_continuous(limits = c(-5,5)) + 
  scale_y_continuous(limits = c(18, 28)) + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.6), panel.border = element_blank(), axis.line = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_cell_composition_figure1.png",
       plot = d1,
       scale = 1, width = 12, height = 4, units = "in", device = "png",
       dpi = 300)

##########################
# bar plot for composition
##########################

plot.data <- data9set_cleaned_MG.SO@meta.data %>%
  dplyr::select(sample, cluster = seurat_clusters) %>%
  mutate(cluster = cluster) %>%
  group_by(cluster, sample) %>%
  summarise(count = n()) %>%
  mutate(clust_total = sum(count)) %>%
  mutate(clust_prop = count / clust_total) %>%
  group_by(sample) %>%
  mutate(dataset_total = sum(count)) %>%
  ungroup() %>%
  mutate(dataset_prop = count / dataset_total)

plot.data$sample <- factor(plot.data$sample, levels = c("Ctrl", "AD", "ADp40KO"))

g1 <- ggplot(plot.data, aes(x = cluster, y = count, fill = sample)) +
  geom_col() + 
  theme_classic() + 
  scale_y_continuous(limits = c(0, 4000), expand = c(0,0)) + 
  scale_fill_viridis_d() + 
  geom_text(aes(label = count), color = "darkgray", size = 5,  position = position_stack(0.4)) + 
  theme(axis.text.x = element_text(size =14, color = "black"), axis.title.x = element_blank(), axis.text.y = element_text(size =14, color = "black"),
                                                                                                       axis.title.y = element_text(size =14), legend.title = element_blank(), legend.text = element_text(size = 12)) + guides(fill = guide_legend(override.aes = list(size = 8), label.theme = element_text(size = 15)))

g2 <- ggplot(plot.data, aes(x = cluster, y = clust_prop, fill = sample)) +
  geom_col() + theme_classic() + scale_y_continuous(expand = c(0,0)) + scale_fill_viridis_d() + 
  ylab("percentage") + 
  geom_text(aes(label = round(clust_prop, 3)), color = "darkgray", size = 5,  position = position_stack(0.4)) + 
  theme(axis.text.x = element_text(size =14, color = "black"), axis.title.x = element_blank(), axis.text.y = element_text(size =14, color = "black"),
                                                                                                                      axis.title.y = element_text(size =14), legend.position = "none")

g3 <- arrangeGrob(g1, g2, ncol = 2, widths = c(0.6, 0.4))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_cell_composition_figure2.png", 
       plot = g3, 
       scale = 1, width = 8, height = 5, units = "in", device = "png",
       dpi = 300)


##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_cell_composition_session_info.txt")

