#! /bin/env RScript
# written by SJK at 11. May. 2020

.libPaths(c("/home/skim/R/usr_lib", .libPaths()))

library(Seurat)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)
library(tidyr)

source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

load(file="/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200210_figures_for_meeting_report_cleaned_Seurat_object_res0.8.Robj")


########################
# only WT
########################

data9set_cleaned_Ctrl.SO <- subset(data9set_cleaned.SO, subset = sample %in% "Ctrl")


d1 <- DimPlot(data9set_cleaned_Ctrl.SO, reduction = "umap", group.by = "gemgroup", cols = c(rep(viridis(3))[1], "#CC9900", rep(viridis(3))[2]), pt.size = 0.00001) +
  scale_color_manual(labels = c("WT_1", "WT_2", "WT_3"), values = c("#FF0000", "#003366", "#FFCC33"))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig1_batch_effects_for_small_samples_WT_only.png",
       plot = d1,
       scale = 1, width = 7.2, height = 6, units = "in", device = "png",
       dpi = 300)


#############################
# by 1 sample in eah genotype
#############################
data9set_cleaned_123.SO <- subset(data9set_cleaned.SO, subset = gemgroup %in% c(1,2,3))

d2 <- DimPlot(data9set_cleaned_123.SO, reduction = "umap", group.by = "gemgroup", cols = c(rep(viridis(3))[1], "#CC9900", rep(viridis(3))[2]), pt.size = 0.00001) +
  scale_color_manual(labels = c("WT_1", "AD.p40KO_1", "AD_1"), values = c("#FF0000", "#003366", "#FFCC33"))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig1_batch_effects_for_small_samples_1sample_genotype.png",
       plot = d2,
       scale = 1, width = 7.5, height = 6, units = "in", device = "png",
       dpi = 300)

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig1_batch_effects_for_small_samples_session_info.txt")


