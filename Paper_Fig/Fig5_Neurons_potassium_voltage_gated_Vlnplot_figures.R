#! /bin/env RScript
# written by SJK at 2. April. 2020
# This file is for creating vlnplot for specific gene in neurons 

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


#########################
# cluster 13 CA1 Neurons
# Kcnh1
#########################

data9set_13.SO <- subset(data9set_cleaned.SO, subset = seurat_clusters %in% c(13))

data9set_13_Ctrl.SO <- subset(data9set_13.SO, subset = sample %in% "Ctrl")
data9set_13_AD.SO <- subset(data9set_13.SO, subset = sample %in% "AD")
data9set_13_ADp40KO.SO <- subset(data9set_13.SO, subset = sample %in% "ADp40KO")


v1 <- VlnPlot(data9set_13_Ctrl.SO, features = "Kcnh1", y.max =3, cols = viridis(3)[1], pt.size = 0) +
  xlab ("Cluster 13 (Ctrl)") + ggtitle(" ") + 
  theme(axis.title.x = element_text(size = 18, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"), 
        axis.text.x = element_blank(), 
        title = element_text(size = 10),
        legend.position = "none")

v2 <- VlnPlot(data9set_13_AD.SO, features = "Kcnh1", y.max =3, cols = viridis(3)[2], pt.size = 0) +
  xlab ("Cluster 13 (AD)") + ylab("  ") + 
  theme(axis.title.x = element_text(size = 18, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_blank(), 
        title = element_text(size = 20),
        legend.position = "none", 
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank())

v3 <- VlnPlot(data9set_13_ADp40KO.SO, features = "Kcnh1", y.max =3, cols = "#CC9966", pt.size = 0) +
  xlab ("Cluster 13 (ADp40KO)") + ylab("  ") + ggtitle(" ") + 
  theme(axis.title.x = element_text(size = 18, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_blank(), 
        title = element_text(size = 10),
        legend.position = "none", 
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank())

g1 <- arrangeGrob(v1, v2, v3, ncol = 3)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Neurons_potassium_voltage_gated_Vlnplot_Kcnh1_figures_whole.png",
       plot = g1,
       scale = 1, width = 9, height = 4.5, units = "in", device = "png",
       dpi = 300)


#########################
# cluster 20 subiculum
# Kcnc2
#########################

data9set_20.SO <- subset(data9set_cleaned.SO, subset = seurat_clusters %in% c(20))

data9set_20_Ctrl.SO <- subset(data9set_20.SO, subset = sample %in% "Ctrl")
data9set_20_AD.SO <- subset(data9set_20.SO, subset = sample %in% "AD")
data9set_20_ADp40KO.SO <- subset(data9set_20.SO, subset = sample %in% "ADp40KO")

v1 <- VlnPlot(data9set_20_Ctrl.SO, features = "Kcnc2", y.max =5, cols = viridis(3)[1], pt.size = 0) +
  xlab ("Cluster 20 (Ctrl)") + ggtitle(" ") + 
  theme(axis.title.x = element_text(size = 18, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"), 
        axis.text.x = element_blank(), 
        title = element_text(size = 10),
        legend.position = "none")

v2 <- VlnPlot(data9set_20_AD.SO, features = "Kcnc2", y.max =5, cols = viridis(3)[2], pt.size = 0) +
  xlab ("Cluster 20 (AD)") + ylab("  ") + 
  theme(axis.title.x = element_text(size = 18, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_blank(), 
        title = element_text(size = 20),
        legend.position = "none", 
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank())

v3 <- VlnPlot(data9set_20_ADp40KO.SO, features = "Kcnc2", y.max =5, cols = "#CC9966", pt.size = 0) +
  xlab ("Cluster 20 (ADp40KO)") + ylab("  ") + ggtitle(" ") + 
  theme(axis.title.x = element_text(size = 18, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_blank(), 
        title = element_text(size = 10),
        legend.position = "none", 
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank())

g1 <- arrangeGrob(v1, v2, v3, ncol = 3)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Neurons_potassium_voltage_gated_Vlnplot_Kcnc2_figures_whole.png",
       plot = g1,
       scale = 1, width = 9, height = 4.5, units = "in", device = "png",
       dpi = 300)

#########################
# cluster 20 subiculum
# Kcnd2
#########################

v1 <- VlnPlot(data9set_20_Ctrl.SO, features = "Kcnd2", y.max =6, cols = viridis(3)[1], pt.size = 0) +
  xlab ("Cluster 20 (Ctrl)") + ggtitle(" ") + 
  theme(axis.title.x = element_text(size = 18, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"), 
        axis.text.x = element_blank(), 
        title = element_text(size = 10),
        legend.position = "none")

v2 <- VlnPlot(data9set_20_AD.SO, features = "Kcnd2", y.max =6, cols = viridis(3)[2], pt.size = 0) +
  xlab ("Cluster 20 (AD)") + ylab("  ") + 
  theme(axis.title.x = element_text(size = 18, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_blank(), 
        title = element_text(size = 20),
        legend.position = "none", 
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank())

v3 <- VlnPlot(data9set_20_ADp40KO.SO, features = "Kcnd2", y.max =6, cols = "#CC9966", pt.size = 0) +
  xlab ("Cluster 20 (ADp40KO)") + ylab("  ") + ggtitle(" ") + 
  theme(axis.title.x = element_text(size = 18, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_blank(), 
        title = element_text(size = 10),
        legend.position = "none", 
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank())

g1 <- arrangeGrob(v1, v2, v3, ncol = 3)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Neurons_potassium_voltage_gated_Vlnplot_Kcnd2_figures_whole.png",
       plot = g1,
       scale = 1, width = 9, height = 4.5, units = "in", device = "png",
       dpi = 300)

#########################
# cluster 20 subiculum
# Kcnip4
#########################

v1 <- VlnPlot(data9set_20_Ctrl.SO, features = "Kcnip4", y.max =6, cols = viridis(3)[1], pt.size = 0) +
  xlab ("Cluster 20 (Ctrl)") + ggtitle(" ") + 
  theme(axis.title.x = element_text(size = 18, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"), 
        axis.text.x = element_blank(), 
        title = element_text(size = 10),
        legend.position = "none")

v2 <- VlnPlot(data9set_20_AD.SO, features = "Kcnip4", y.max =6, cols = viridis(3)[2], pt.size = 0) +
  xlab ("Cluster 20 (AD)") + ylab("  ") + 
  theme(axis.title.x = element_text(size = 18, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_blank(), 
        title = element_text(size = 20),
        legend.position = "none", 
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank())

v3 <- VlnPlot(data9set_20_ADp40KO.SO, features = "Kcnip4", y.max =6, cols = "#CC9966", pt.size = 0) +
  xlab ("Cluster 20 (ADp40KO)") + ylab("  ") + ggtitle(" ") + 
  theme(axis.title.x = element_text(size = 18, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.x = element_blank(), 
        title = element_text(size = 10),
        legend.position = "none", 
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank())

g1 <- arrangeGrob(v1, v2, v3, ncol = 3)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Neurons_potassium_voltage_gated_Vlnplot_Kcnip4_figures_whole.png",
       plot = g1,
       scale = 1, width = 9, height = 4.5, units = "in", device = "png",
       dpi = 300)

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Neurons_potassium_voltage_gated_Vlnplot_figures_session_info.txt")