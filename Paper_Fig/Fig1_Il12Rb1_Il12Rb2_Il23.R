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

f1 <- FeaturePlot(data9set_cleaned_Ctrl.SO, features = "Il23r", order = TRUE, pt.size = 0.5)
f1 <- f1 + scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", midpoint = 1.5, limits=c(0, 2.5), breaks=c(0,1.5,2.5))  + ggtitle("Ctrl") + labs(y = "Il23r") + 
  theme(legend.position = "none",
        axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15))


f2 <- FeaturePlot(data9set_cleaned_AD.SO, features = "Il23r", order = TRUE, pt.size = 0.5)
f2 <- f2 + scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", midpoint = 1.5, limits=c(0, 2.5), breaks=c(0,1.5,2.5)) + ggtitle("AD") + 
  theme(legend.position = "none",
        axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank())

f3 <- FeaturePlot(data9set_cleaned_ADp40KO.SO, features = "Il23r", order = TRUE, pt.size = 0.5)
f3 <- f3 + scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", limits=c(0, 2.5), midpoint = 1.5, breaks=c(0,1.5,2.5))  +   ggtitle("ADp40KO") + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank())

g1 <- arrangeGrob(f1, f2, f3, ncol = 3, widths= c(1.1,1.05,1.2))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Il23r.png", 
       plot = g1, 
  scale = 1, width = 15, height = 5, units = "in", device = "png",
  dpi = 300)

#####################
# Il12rb2
#####################
f1 <- FeaturePlot(data9set_cleaned_Ctrl.SO, features = "Il12rb2", order = TRUE, pt.size = 0.5)
f1 <- f1 + scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", limits=c(0, 3.5), midpoint = 2, breaks=c(0,2,3.5))  + 
  ggtitle("Ctrl") + labs(y = "Il12rb2") + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        legend.position = "none")


f2 <- FeaturePlot(data9set_cleaned_AD.SO, features = "Il12rb2", order = TRUE, pt.size = 0.5)
f2 <- f2 + scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", limits=c(0, 3.5), midpoint = 2, breaks=c(0,2,4)) + 
  ggtitle("AD") + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

f3 <- FeaturePlot(data9set_cleaned_ADp40KO.SO, features = "Il12rb2", order = TRUE, pt.size = 0.5)
f3 <- f3 + scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", limits=c(0, 3.5), midpoint = 2, breaks=c(0,2,3.5))  +   
  ggtitle("ADp40KO") + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank())

g2 <- arrangeGrob(f1, f2, f3, ncol = 3, widths= c(1.1,1.05,1.2))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Il12rb2.png",
       plot = g2,
  scale = 1, width = 15, height = 5, units = "in", device = "png",
  dpi = 300)


########################
# Il23a
########################
f1 <- FeaturePlot(data9set_cleaned_Ctrl.SO, features = "Il23a", order = TRUE, pt.size = 0.5)
f1 <- f1 + scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", limits=c(0, 3), midpoint = 1.5, breaks=c(0,1.5,3))  + 
  ggtitle("Ctrl") + labs(y = "Il23a") + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        legend.position = "none")


f2 <- FeaturePlot(data9set_cleaned_AD.SO, features = "Il23a", order = TRUE, pt.size = 0.5)
f2 <- f2 + scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", limits=c(0, 3), midpoint = 1.5, breaks=c(0,1.5,3)) + 
  ggtitle("AD") + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

f3 <- FeaturePlot(data9set_cleaned_ADp40KO.SO, features = "Il23a", order = TRUE, pt.size = 0.5)
f3 <- f3 + scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", limits=c(0, 3), midpoint = 1.5, breaks=c(0,1.5,3))  +   
  ggtitle("ADp40KO") + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank())

g3 <- arrangeGrob(f1, f2, f3, ncol = 3, widths= c(1.1,1.05,1.2))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Il23a.png",
       plot = g3,
  scale = 1, width = 15, height = 5, units = "in", device = "png",
  dpi = 300)


#######################
# Il12a
#######################
f1 <- FeaturePlot(data9set_cleaned_Ctrl.SO, features = "Il12a", order = TRUE, pt.size = 0.5)
f1 <- f1 + scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", limits=c(0, 3.5), midpoint = 2, breaks=c(0,2,3.5))  + 
  ggtitle("Ctrl") + labs(y = "Il12a") + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        legend.position = "none")


f2 <- FeaturePlot(data9set_cleaned_AD.SO, features = "Il12a", order = TRUE, pt.size = 0.5)
f2 <- f2 + scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", limits=c(0, 3.5), midpoint = 2, breaks=c(0,2,3.5)) + 
  ggtitle("AD") + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

f3 <- FeaturePlot(data9set_cleaned_ADp40KO.SO, features = "Il12a", order = TRUE, pt.size = 0.5)
f3 <- f3 + scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", limits=c(0, 3.5), midpoint = 2, breaks=c(0,2,3.5))  +   
  ggtitle("ADp40KO") + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank())


g4 <- arrangeGrob(f1, f2, f3, ncol = 3, widths= c(1.1,1.05,1.2))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Il12a.png",
       plot = g4,
  scale = 1, width = 15, height = 5, units = "in", device = "png",
  dpi = 300)

#######################
# Il12rb1
#######################
f1 <- FeaturePlot(data9set_cleaned_Ctrl.SO, features = "Il12rb1", order = TRUE, pt.size = 0.5)
f1 <- f1 + scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", limits=c(0, 3.5), midpoint = 2, breaks=c(0,2,3.5))  +
  ggtitle("Ctrl") + labs(y = "Il12rb1") +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        legend.position = "none")


f2 <- FeaturePlot(data9set_cleaned_AD.SO, features = "Il12rb1", order = TRUE, pt.size = 0.5)
f2 <- f2 + scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", limits=c(0, 3.5), midpoint = 2, breaks=c(0,2,3.5)) +
  ggtitle("AD") +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

f3 <- FeaturePlot(data9set_cleaned_ADp40KO.SO, features = "Il12rb1", order = TRUE, pt.size = 0.5)
f3 <- f3 + scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", limits=c(0, 3.5), midpoint = 2, breaks=c(0,2,3.5))  +
  ggtitle("AD.p40KO") +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title = element_blank())


g5 <- arrangeGrob(f1, f2, f3, ncol = 3, widths= c(1.1,1.05,1.2))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Il12rb1.png",
       plot = g5,
  scale = 1, width = 15, height = 5, units = "in", device = "png",
  dpi = 300)




##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Il12Rb1_Il12Rb2_Il23a_Il12a_session_info.txt")

