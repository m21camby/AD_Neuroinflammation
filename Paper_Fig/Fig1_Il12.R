#! /bin/env RScript
# written by SJK at 2. Mar. 2020

.libPaths(c("/home/skim/R/usr_lib", .libPaths()))


library(Seurat)
library(ggplot2)
library(dplyr)
library(gridExtra)
source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

load(file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200221_figure1_for_meeting_report_Seurat_object.Robj")

data9set_cleaned_Ctrl.SO <- subset(data9set_cleaned.SO, subset = sample %in% "Ctrl")
data9set_cleaned_AD.SO <- subset(data9set_cleaned.SO, subset = sample %in% "AD")
data9set_cleaned_ADp40KO.SO <- subset(data9set_cleaned.SO, subset = sample %in% "ADp40KO")

f1 <- FeaturePlot(data9set_cleaned_Ctrl.SO, features = "Il12b", order = TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
f1 <- f1 + scale_colour_gradientn(colours = c("#003366","003366"), values = c(0.1, 0), na.value = "#003366") + ggtitle("Ctrl") + labs(y = "Il12b") + 
  theme(legend.position = "none", 
        axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title.x = element_blank())


f2 <- FeaturePlot(data9set_cleaned_AD.SO, features = "Il12b", sort.cell =  TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
f2 <- f2 + scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", midpoint = 1.5) + ggtitle("AD") + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank())

f3 <- FeaturePlot(data9set_cleaned_ADp40KO.SO, features = "Il12b", order = TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
f3 <- f3 + scale_color_gradient2(low="#003366", high='#003366', mid="#003366", midpoint = 1.5, na.value = "#003366") + ggtitle("AD.p40KO") + theme(legend.position = "none",axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank())


g1 <- arrangeGrob(f1, f3, f2, ncol = 3, widths= c(1,1,1.2))



ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Il12.png", 
       plot = g1, 
  scale = 1, width = 15, height = 5, units = "in", device = "png",
  dpi = 300)


#################
# Different color
#################

f1 <- FeaturePlot(data9set_cleaned_Ctrl.SO, features = "Il12b", order = TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
f1 <- f1 + scale_colour_gradientn(colours = c("#CCCCCC","#CCCCCC"), values = c(0.1, 0), na.value = "#CCCCCC") + ggtitle("WT") + labs(y = "Il12b") + 
  theme(legend.position = "none", 
        axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title.x = element_blank())


f2 <- FeaturePlot(data9set_cleaned_AD.SO, features = "Il12b", sort.cell =  TRUE, order = TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5, cols = c("#CCCCCC", "#0000FF"))
f2 <- f2 + ggtitle("AD") + scale_color_continuous(high = "blue", low = "gray", limit = c(0,3)) + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank())

f3 <- FeaturePlot(data9set_cleaned_ADp40KO.SO, features = "Il12b", order = TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
f3 <- f3 + scale_color_gradient2(low="#CCCCCC", high='#CCCCCC', mid="#CCCCCC", midpoint = 1.5, na.value = "#CCCCCC") + ggtitle("AD.p40KO") + 
  theme(legend.position = "none",axis.line=element_blank(),  axis.ticks=element_blank(), axis.text=element_blank(), axis.title = element_blank())


g1 <- arrangeGrob(f1, f3, f2, ncol = 3, widths= c(1,1,1.2))



ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Il12_blue_gray.png",
       plot = g1,
  scale = 1, width = 15, height = 5, units = "in", device = "png",
  dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Il12_blue_gray.pdf",
       plot = g1,
       scale = 1, width = 15, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Il12_session_info.txt")


