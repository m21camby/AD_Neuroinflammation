#! /bin/env RScript
# written by SJK at 8. May. 2020

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

#####################
# Il12rb1
#####################
f1 <- FeaturePlot(data9set_cleaned_Ctrl.SO, features = "Il12rb1", cols = c("#CCCCCC", "#0000FF"), order = TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
f1 <- f1 + ggtitle("WT") + labs(y = "Il12rb1") + scale_color_continuous(high = "blue", low = "gray", limit = c(0,3)) + 
  theme(legend.position = "none", 
        axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title.x = element_blank())


f2 <- FeaturePlot(data9set_cleaned_AD.SO, features = "Il12rb1", cols = c("#CCCCCC", "#0000FF"), sort.cell =  TRUE, order = TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
f2 <- f2 + ggtitle("AD") + scale_color_continuous(high = "blue", low = "gray", limit = c(0,3)) + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank())

f3 <- FeaturePlot(data9set_cleaned_ADp40KO.SO, features = "Il12rb1", cols = c("#CCCCCC", "#0000FF"), order = TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
f3 <- f3 + ggtitle("AD.p40KO") + scale_color_continuous(high = "blue", low = "gray", limit = c(0,3)) + 
  theme(legend.position = "none",axis.line=element_blank(),  axis.ticks=element_blank(), axis.text=element_blank(), axis.title = element_blank())


g1 <- arrangeGrob(f1, f3, f2, ncol = 3, widths= c(1.1,1.05,1.2))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Il12rb1_blue_gray.png", 
       plot = g1, 
  scale = 1, width = 15, height = 5, units = "in", device = "png",
  dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Il12rb1_blue_gray.pdf", 
       plot = g1, 
       scale = 1, width = 15, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

#####################
# Il12rb2
#####################
f1 <- FeaturePlot(data9set_cleaned_Ctrl.SO, features = "Il12rb2", cols = c("#CCCCCC", "#0000FF"), order = TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
f1 <- f1 + ggtitle("WT") + labs(y = "Il12rb2") + scale_color_continuous(high = "blue", low = "gray", limit = c(0,3)) + 
  theme(legend.position = "none", 
        axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title.x = element_blank())


f2 <- FeaturePlot(data9set_cleaned_AD.SO, features = "Il12rb2", cols = c("#CCCCCC", "#0000FF"), sort.cell =  TRUE, order = TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
f2 <- f2 + ggtitle("AD") + scale_color_continuous(high = "blue", low = "gray", limit = c(0,3)) + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank())

f3 <- FeaturePlot(data9set_cleaned_ADp40KO.SO, features = "Il12rb2", cols = c("#CCCCCC", "#0000FF"), order = TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
f3 <- f3 + ggtitle("AD.p40KO") + scale_color_continuous(high = "blue", low = "gray", limit = c(0,3)) + 
  theme(legend.position = "none",axis.line=element_blank(),  axis.ticks=element_blank(), axis.text=element_blank(), axis.title = element_blank())


g2 <- arrangeGrob(f1, f3, f2, ncol = 3, widths= c(1.1,1.05,1.2))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Il12rb2_blue_gray.png",
       plot = g2,
  scale = 1, width = 15, height = 5, units = "in", device = "png",
  dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Il12rb2_blue_gray.pdf",
       plot = g2,
       scale = 1, width = 15, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

########################
# Il23a
########################
f1 <- FeaturePlot(data9set_cleaned_Ctrl.SO, features = "Il23a", cols = c("#CCCCCC", "#0000FF"), order = TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
f1 <- f1 + ggtitle("WT") + labs(y = "Il23a") + 
  theme(legend.position = "none", 
        axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title.x = element_blank())


f2 <- FeaturePlot(data9set_cleaned_AD.SO, features = "Il23a", cols = c("#CCCCCC", "#0000FF"), sort.cell =  TRUE, order = TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
f2 <- f2 + ggtitle("AD") + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank())

f3 <- FeaturePlot(data9set_cleaned_ADp40KO.SO, features = "Il23a", cols = c("#CCCCCC", "#0000FF"), order = TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
f3 <- f3 + ggtitle("AD.p40KO") + 
  theme(legend.position = "none",axis.line=element_blank(),  axis.ticks=element_blank(), axis.text=element_blank(), axis.title = element_blank())


g3 <- arrangeGrob(f1, f3, f2, ncol = 3, widths= c(1.1,1.05,1.2))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Il23a_blue_gray.png",
       plot = g3,
  scale = 1, width = 15, height = 5, units = "in", device = "png",
  dpi = 300)


#######################
# Il12a
#######################

f1 <- FeaturePlot(data9set_cleaned_Ctrl.SO, features = "Il12a", cols = c("#CCCCCC", "#0000FF"), order = TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
f1 <- f1 + ggtitle("WT") + labs(y = "Il12a") + 
  theme(legend.position = "none", 
        axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title.x = element_blank())


f2 <- FeaturePlot(data9set_cleaned_AD.SO, features = "Il12a", cols = c("#CCCCCC", "#0000FF"), sort.cell =  TRUE, order = TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
f2 <- f2 + ggtitle("AD") + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank())

f3 <- FeaturePlot(data9set_cleaned_ADp40KO.SO, features = "Il12a", cols = c("#CCCCCC", "#0000FF"), order = TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
f3 <- f3 + ggtitle("AD.p40KO") + 
  theme(legend.position = "none",axis.line=element_blank(),  axis.ticks=element_blank(), axis.text=element_blank(), axis.title = element_blank())

g4 <- arrangeGrob(f1, f3, f2, ncol = 3, widths= c(1.1,1.05,1.2))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Il12a_blue_gray.png",
       plot = g4,
  scale = 1, width = 15, height = 5, units = "in", device = "png",
  dpi = 300)

#######################
# Il23r
#######################

f1 <- FeaturePlot(data9set_cleaned_Ctrl.SO, features = "Il23r", cols = c("#CCCCCC", "#0000FF"), order = TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
f1 <- f1 + ggtitle("WT") + labs(y = "Il23r") + scale_color_continuous(high = "blue", low = "gray", limit = c(0,3)) +
  theme(legend.position = "none", 
        axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title.x = element_blank())


f2 <- FeaturePlot(data9set_cleaned_AD.SO, features = "Il23r", cols = c("#CCCCCC", "#0000FF"), sort.cell =  TRUE, order = TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
f2 <- f2 + ggtitle("AD") + scale_color_continuous(high = "blue", low = "gray", limit = c(0,3)) + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank())

f3 <- FeaturePlot(data9set_cleaned_ADp40KO.SO, features = "Il23r", cols = c("#CCCCCC", "#0000FF"), order = TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
f3 <- f3 + ggtitle("AD.p40KO") + scale_color_continuous(high = "blue", low = "gray", limit = c(0,3)) + 
  theme(legend.position = "none",axis.line=element_blank(),  axis.ticks=element_blank(), axis.text=element_blank(), axis.title = element_blank())

g5 <- arrangeGrob(f1, f3, f2, ncol = 3, widths= c(1.1,1.05,1.2))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Il23r_blue_gray.png",
       plot = g5,
  scale = 1, width = 15, height = 5, units = "in", device = "png",
  dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Il23r_blue_gray_New_Scale.pdf",
       plot = g5,
       scale = 1, width = 15, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)


##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Il12Rb1_Il12Rb2_Il23a_Il12a_blue_gray_session_info.txt")

