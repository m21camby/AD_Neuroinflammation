#! /bin/env RScript
# written by SJK at 24. Mar. 2020

.libPaths(c("/home/skim/R/usr_lib", .libPaths()))

library(Seurat)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)
library(tidyr)
library(grid)
library(ggrepel)
source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

#################
# Cluster 3
#################
cluster3_AD_ADp40KO <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/DE_ADp40KO_AD/20200117_Clustering_DE_by_MAST_each_Cluster_3ADp40KO_AD_corrected.csv", row.names = 1)
cluster3_AD_ADp40KO$gene <- rownames(cluster3_AD_ADp40KO)
cluster3_AD_ADp40KO$avg_logFC <- -1*cluster3_AD_ADp40KO$avg_logFC

# Plus
cluster3_AD_ADp40KO_plus.df <- cluster3_AD_ADp40KO[which(cluster3_AD_ADp40KO$avg_logFC > 0.2  & cluster3_AD_ADp40KO$p_val_adj < .01),]
# Minus
cluster3_AD_ADp40KO_minus.df <- cluster3_AD_ADp40KO[which(cluster3_AD_ADp40KO$avg_logFC < -0.2 & cluster3_AD_ADp40KO$p_val_adj < .01),]
# Il12b gene
cluster3_AD_ADp40KO_Il12b.df <-  cluster3_AD_ADp40KO[cluster3_AD_ADp40KO$gene == "Il12b", ]

g1 <- ggplot(cluster3_AD_ADp40KO) + geom_point(aes(x = avg_logFC, y = -log10(p_val_adj))) +
  geom_point(data = cluster3_AD_ADp40KO_plus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), color = "red") +
  geom_point(data = cluster3_AD_ADp40KO_minus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), color = "blue") +
  ggtitle("AD versus ADp40KO in Microglia cluster 3") + 
  theme_classic() + 
  xlab("log2 fold change") + 
  ylab("-log10(adjusted p-value)") +
  scale_y_continuous(expand = c(0,0), limits = c(0,15)) + xlim(-0.5, 0.5) + 
  theme(legend.position = "none",
        plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12, color = "black")) +
  geom_text_repel(data = cluster3_AD_ADp40KO_plus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), label = cluster3_AD_ADp40KO_plus.df$gene, force = 10) + 
  geom_text_repel(data = cluster3_AD_ADp40KO_minus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), label = cluster3_AD_ADp40KO_minus.df$gene, force = 10) + 
  geom_text_repel(data = cluster3_AD_ADp40KO_Il12b.df, aes(x = avg_logFC, y = -log10(p_val_adj)), label = cluster3_AD_ADp40KO_Il12b.df$gene, nudge_y = 1)

#################
# Cluster 8
#################

cluster8_AD_ADp40KO <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/DE_ADp40KO_AD/20200117_Clustering_DE_by_MAST_each_Cluster_8ADp40KO_AD_corrected.csv", row.names = 1)
cluster8_AD_ADp40KO$gene <- rownames(cluster8_AD_ADp40KO)
cluster8_AD_ADp40KO$avg_logFC <- -1*cluster8_AD_ADp40KO$avg_logFC

# Plus
cluster8_AD_ADp40KO_plus.df <- cluster8_AD_ADp40KO[which(cluster8_AD_ADp40KO$avg_logFC > 0.2  & cluster8_AD_ADp40KO$p_val_adj < .01),]
# Minus
cluster8_AD_ADp40KO_minus.df <- cluster8_AD_ADp40KO[which(cluster8_AD_ADp40KO$avg_logFC < -0.2 & cluster8_AD_ADp40KO$p_val_adj < .01),]

g2 <- ggplot(cluster8_AD_ADp40KO) + geom_point(aes(x = avg_logFC, y = -log10(p_val_adj))) +
  geom_point(data = cluster8_AD_ADp40KO_plus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), color = "red") +
  geom_point(data = cluster8_AD_ADp40KO_minus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), color = "blue") +
  ggtitle("AD versus ADp40KO in DAM Microglia cluster 8") + 
  theme_classic() + 
  xlab("log2 fold change") + 
  ylab("-log10(adjusted p-value)") +
  scale_y_continuous(expand = c(0,0), limits = c(0,10)) + xlim(-0.75, 0.75) + 
  theme(legend.position = "none",
        plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12, color = "black")) +
  geom_text_repel(data = cluster8_AD_ADp40KO_plus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), label = cluster8_AD_ADp40KO_plus.df$gene, force = 10) + 
  geom_text_repel(data = cluster8_AD_ADp40KO_minus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), label = cluster8_AD_ADp40KO_minus.df$gene, force = 10)

g12 <- arrangeGrob(g1, g2, ncol = 2)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_Cluster_DE_volcano_plot.png",
       plot = g12,
       scale = 1, width = 11, height = 5, units = "in", device = "png",
       dpi = 300)



##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_Cluster_DE_volcano_plot_session_info.txt")

