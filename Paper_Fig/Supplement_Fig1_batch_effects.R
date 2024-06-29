#! /bin/env RScript
# written by SJK at 3. Mar. 2020

.libPaths(c("/data/murphy/shared_libs/R", "/usr/local/lib/R/site-librar","/usr/local/lib/R/library","/usr/lib64/R/library","/usr/share/R/library", .libPaths() ))

library(Seurat)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)
source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

load(file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200117_Clustering_DE_by_MAST_Seurat_object.Robj")


data9set.SO@meta.data$sample <- ifelse((data9set.SO$gemgroup == 1 | data9set.SO$gemgroup == 4 | data9set.SO$gemgroup == 7), "Ctrl", ifelse((data9set.SO$gemgroup == 2 | data9set.SO$gemgroup == 5 | data9set.SO$gemgroup == 8), "ADp40KO", "AD"))


data9set.SO@meta.data$sample <- factor(data9set.SO@meta.data$sample, levels = c("Ctrl", "AD", "ADp40KO"))


d1 <- DimPlot(data9set.SO, cells =  rownames(data9set.SO@meta.data[data9set.SO@meta.data$gemgroup %in% "1", ]), reduction = "umap", cols = rep(viridis(3)[1], 45))

d2 <- DimPlot(data9set.SO, cells =  rownames(data9set.SO@meta.data[data9set.SO@meta.data$gemgroup %in% "2", ]), reduction = "umap", cols = rep(viridis(3)[2], 45))

d3 <- DimPlot(data9set.SO, cells =  rownames(data9set.SO@meta.data[data9set.SO@meta.data$gemgroup %in% "3", ]), reduction = "umap", cols = rep("#CC9966", 45))

d4 <- DimPlot(data9set.SO, cells =  rownames(data9set.SO@meta.data[data9set.SO@meta.data$gemgroup %in% "4", ]), reduction = "umap", cols = rep(viridis(3)[1], 45))

d5 <- DimPlot(data9set.SO, cells =  rownames(data9set.SO@meta.data[data9set.SO@meta.data$gemgroup %in% "5", ]), reduction = "umap", cols = rep(viridis(3)[2], 45))

d6 <- DimPlot(data9set.SO, cells =  rownames(data9set.SO@meta.data[data9set.SO@meta.data$gemgroup %in% "6", ]), reduction = "umap", cols = rep("#CC9966", 45))

d7 <- DimPlot(data9set.SO, cells =  rownames(data9set.SO@meta.data[data9set.SO@meta.data$gemgroup %in% "7", ]), reduction = "umap", cols = rep(viridis(3)[1], 45))

d8 <- DimPlot(data9set.SO, cells =  rownames(data9set.SO@meta.data[data9set.SO@meta.data$gemgroup %in% "8", ]), reduction = "umap", cols = rep(viridis(3)[2], 45))

d9 <- DimPlot(data9set.SO, cells =  rownames(data9set.SO@meta.data[data9set.SO@meta.data$gemgroup %in% "9", ]), reduction = "umap", cols = rep("#CC9966", 45))

d1 <- d1 + theme(legend.position = "none",
        axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank()) + 
  geom_text(geom = "text", aes(x = 18, y = 25), label = "Ctrl-1") +
  geom_text(geom = "text", aes(x = 18, y = 23), 
            label = paste0("cells: ",nrow(data9set.SO@meta.data[data9set.SO@meta.data$gemgroup %in% "1", ])))
d2 <- d2 + theme(legend.position = "none",
        axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank()) + 
  geom_text(geom = "text", aes(x = 18, y = 25), label = "ADp40KO-1") +
  geom_text(geom = "text", aes(x = 18, y = 23), 
            label = paste0("cells: ",nrow(data9set.SO@meta.data[data9set.SO@meta.data$gemgroup %in% "2", ])))
d3 <- d3 + theme(legend.position = "none",
        axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank()) + 
  geom_text(geom = "text", aes(x = 18, y = 25), label = "AD-1") +
  geom_text(geom = "text", aes(x = 18, y = 23), 
            label = paste0("cells: ",nrow(data9set.SO@meta.data[data9set.SO@meta.data$gemgroup %in% "3", ])))
d4 <- d4 + theme(legend.position = "none",
        axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank()) + 
  geom_text(geom = "text", aes(x = 18, y = 25), label = "Ctrl-2") +
  geom_text(geom = "text", aes(x = 18, y = 23), 
            label = paste0("cells: ",nrow(data9set.SO@meta.data[data9set.SO@meta.data$gemgroup %in% "4", ])))
d5 <- d5 + theme(legend.position = "none",
        axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank()) + 
  geom_text(geom = "text", aes(x = 18, y = 25), label = "ADp40KO-2") +
  geom_text(geom = "text", aes(x = 18, y = 23), 
            label = paste0("cells: ",nrow(data9set.SO@meta.data[data9set.SO@meta.data$gemgroup %in% "5", ])))
d6 <- d6 + theme(legend.position = "none",
        axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank()) + 
  geom_text(geom = "text", aes(x = 18, y = 25), label = "AD-2") +
  geom_text(geom = "text", aes(x = 18, y = 23), 
            label = paste0("cells: ",nrow(data9set.SO@meta.data[data9set.SO@meta.data$gemgroup %in% "6", ])))
d7 <- d7 + theme(legend.position = "none",
        axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank()) + 
  geom_text(geom = "text", aes(x = 18, y = 25), label = "Ctrl-3") +
  geom_text(geom = "text", aes(x = 18, y = 23), 
            label = paste0("cells: ",nrow(data9set.SO@meta.data[data9set.SO@meta.data$gemgroup %in% "7", ])))
d8 <- d8 + theme(legend.position = "none",
        axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank()) + 
  geom_text(geom = "text", aes(x = 18, y = 25), label = "ADp40KO-3") +
  geom_text(geom = "text", aes(x = 18, y = 23), 
            label = paste0("cells: ",nrow(data9set.SO@meta.data[data9set.SO@meta.data$gemgroup %in% "8", ])))
d9 <- d9 + theme(legend.position = "none",
        axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank()) + 
  geom_text(geom = "text", aes(x = 18, y = 25), label = "AD-3") +
  geom_text(geom = "text", aes(x = 18, y = 23), 
            label = paste0("cells: ",nrow(data9set.SO@meta.data[data9set.SO@meta.data$gemgroup %in% "9", ])))

g1 <- arrangeGrob(d1, d4, d7, d3, d6, d9, d2, d5, d8, ncol = 3)


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig1_batch_effects.png", 
       plot = g1, 
  scale = 1, width = 15, height = 15, units = "in", device = "png",
  dpi = 300)


##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig1_batch_effects_session_info.txt")




