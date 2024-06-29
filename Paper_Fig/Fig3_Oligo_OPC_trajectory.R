#! /bin/env RScript
# written by SJK at 20. Apr. 2020

# This file is for Figure 3 Oligodendrocytes & OPC trajectory feature plot

.libPaths(c("/home/skim/R/usr_lib", .libPaths()))

library(Seurat)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)
library(tidyr)
library(ggrepel)

source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

load(file="/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200210_figures_for_meeting_report_cleaned_Seurat_object_res0.8.Robj")

data9set_OL_OPC.SO <- subset(data9set_cleaned.SO, subset = seurat_clusters %in% c(1, 5, 6, 12, 38))


f1 <- FeaturePlot(data9set_OL_OPC.SO, features = c("Pdgfra")) + 
  scale_x_continuous(limits = c(10,25)) + 
  scale_y_continuous(limits = c(-20, 10)) + 
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", limits=c(0, 5), midpoint = 2.5, breaks=c(0,2.5,5)) + 
  ggtitle("OPC (Pdgfra)") + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

f2 <- FeaturePlot(data9set_OL_OPC.SO, features = c("Fyn")) + 
  scale_x_continuous(limits = c(10,25)) + 
  scale_y_continuous(limits = c(-20, 10)) + 
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", limits=c(0, 5), midpoint = 2.5, breaks=c(0,2.5,5)) + 
  ggtitle("COPs (Fyn)") + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

f3 <- FeaturePlot(data9set_OL_OPC.SO, features = c("Tcf7l2")) + 
  scale_x_continuous(limits = c(10,25)) + 
  scale_y_continuous(limits = c(-20, 10)) + 
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", limits=c(0, 5), midpoint = 2.5, breaks=c(0,2.5,5)) + 
  ggtitle("NFOL (Tcf7l2)") + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

f4 <- FeaturePlot(data9set_OL_OPC.SO, features = c("Opalin")) + 
  scale_x_continuous(limits = c(10,25)) + 
  scale_y_continuous(limits = c(-20, 10)) + 
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", limits=c(0, 5), midpoint = 2.5, breaks=c(0,2.5,5)) + 
  ggtitle("MFOL (Opalin)") + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

f5 <- FeaturePlot(data9set_OL_OPC.SO, features = c("Apod")) + 
  scale_x_continuous(limits = c(10,25)) + 
  scale_y_continuous(limits = c(-20, 10)) + 
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", limits=c(0, 5), midpoint = 2.5, breaks=c(0,2.5,5)) + 
  ggtitle("MOL (Apod)") + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank())

g1 <- arrangeGrob(f1, f2, f3, f4, f5, ncol = 5, widths= c(1,1,1,1,1.35))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory.png",
       plot = g1,
       scale = 1, width = 12, height = 4, units = "in", device = "png",
       dpi = 300)

################################
# Color change
################################

# Nikos recommended to gray90 and blue color

f1 <- FeaturePlot(data9set_OL_OPC.SO, features = c("Pdgfra"), order = FALSE, pt.size = 0.5) + 
  scale_x_continuous(limits = c(10,25)) + 
  scale_y_continuous(limits = c(-20, 10)) + 
  scale_color_gradient(low="gray90", high='blue', limits=c(0, 5)) + 
  ggtitle("OPC (Pdgfra)") + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

f2 <- FeaturePlot(data9set_OL_OPC.SO, features = c("Fyn"), order = FALSE, pt.size = 0.5) + 
  scale_x_continuous(limits = c(10,25)) + 
  scale_y_continuous(limits = c(-20, 10)) + 
  scale_color_gradient(low="gray90", high='blue', limits=c(0, 5)) + 
  ggtitle("COPs (Fyn)") + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

f3 <- FeaturePlot(data9set_OL_OPC.SO, features = c("Tcf7l2"), order = FALSE, pt.size = 0.5) + 
  scale_x_continuous(limits = c(10,25)) + 
  scale_y_continuous(limits = c(-20, 10)) + 
  scale_color_gradient(low="gray90", high='blue', limits=c(0, 5)) +
  ggtitle("NFOL (Tcf7l2)") + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

f4 <- FeaturePlot(data9set_OL_OPC.SO, features = c("Opalin"), order = FALSE, pt.size = 0.5) + 
  scale_x_continuous(limits = c(10,25)) + 
  scale_y_continuous(limits = c(-20, 10)) + 
  scale_color_gradient(low="gray90", high='blue', limits=c(0, 5)) +
  ggtitle("MFOL (Opalin)") + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

f5 <- FeaturePlot(data9set_OL_OPC.SO, features = c("Apod"), order = FALSE, pt.size = 0.5) + 
  scale_x_continuous(limits = c(10,25)) + 
  scale_y_continuous(limits = c(-20, 10)) + 
  scale_color_gradient(low="gray90", high='blue', limits=c(0, 5)) +
  ggtitle("MOL (Apod)") + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank())

g1 <- arrangeGrob(f1, f2, f3, f4, f5, ncol = 5, widths= c(1,1,1,1,1.35))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_color.pdf",
       plot = g1,
       scale = 1, width = 12, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)

g1 <- arrangeGrob(f1, f3, f4, f5, ncol =4, widths= c(1,1,1,1.3))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_color_2nd.pdf",
       plot = g1,
       scale = 1, width = 10, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_info.txt")
