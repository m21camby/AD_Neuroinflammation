#! /bin/env RScript
# written by SJK at 5. Mar. 2020

.libPaths(c("/data/murphy/shared_libs/R", "/usr/local/lib/R/site-librar","/usr/local/lib/R/library","/usr/lib64/R/library","/usr/share/R/library", .libPaths() ))

library(Seurat)
library(ggplot2)
library(dplyr)
library(gridExtra)
source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

load("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200117_Clustering_DE_by_MAST_Seurat_object.Robj")


data9set.SO@meta.data$sample <- ifelse((data9set.SO$gemgroup == 1 | data9set.SO$gemgroup == 4 | data9set.SO$gemgroup == 7), "Ctrl", 
                                       ifelse((data9set.SO$gemgroup == 2 | data9set.SO$gemgroup == 5 | data9set.SO$gemgroup == 8), "ADp40KO", "AD"))



data9set.SO@meta.data$sample <- factor(data9set.SO@meta.data$sample, levels = c("Ctrl", "AD", "ADp40KO"))

scrublet_predicted.df <- read.table("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200205_doublets_by_scrublet_predicted_doublets_data.frame.txt", sep = "\t")

data9set.df <- as.data.frame(data9set.SO@meta.data)
data9set.df$cell_ID <- rownames(data9set.df)
data9set_scrublet.df <- merge(data9set.df, scrublet_predicted.df, by = "cell_ID")
data9set_scrublet.df$scrublet_prediction <- ifelse(data9set_scrublet.df$scrublet_predicted == 1, "doublets", "singlet")

data9set_umap <- as.data.frame(data9set.SO@reductions$umap@cell.embeddings)
data9set_umap$cell_ID <- rownames(data9set_umap)

data9set_umap <- merge(data9set_umap, data9set_scrublet.df, by = "cell_ID")

g1 <- ggplot(data9set_umap, aes(x = nCount_RNA, y = nFeature_RNA, color = scrublet_prediction, alpha = scrublet_prediction)) + geom_point(size = 1) + 
  scale_color_manual(values = c("#FF9900","#3399FF")) + labs(x = "nUMIs", y = "nGenes") + theme_classic() + 
  scale_alpha_manual(values = c(singlet = 0.2, doublets = 1))  + ggtitle("nUMIs and nGenes in doublets") + 
  theme(plot.title = element_text(size = 18, hjust = 0.5), axis.text.x = element_text(size =14, color = "black"), axis.title.x = element_text(size =16, color = "black"), axis.text.y = element_text(size =14, color = "black"), axis.title.y = element_text(size =16), legend.position = c(0.8, 0.3), legend.title = element_blank(), legend.text = element_text(size = 16)) + 
  guides(colour = guide_legend(override.aes = list(size=5)))


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig1_doublets_by_scrublet_analysis_gene_plot.png", 
       plot = g1, 
       scale = 1, width = 7, height = 6, units = "in", device = "png",
       dpi = 300)

#############
# bar plot
#############

cluster_scrublet.df <- data.frame(cluster = data9set_umap$seurat_clusters, doublets = data9set_umap$scrublet_prediction)
cluster_scrublet.df <- cluster_scrublet.df %>% group_by(cluster) %>% count(doublets)
summary_cluster_scrublet.df <- data.frame(cluster = as.character(), scrublet_doublets_percent = as.numeric(), stringsAsFactors = FALSE)

j = 0
for(i in seq(1, by = 2, len = 41)){
  #print(i)
  #print(j)
  summary_cluster_scrublet.df[j+1,1] <- as.character(j)
  summary_cluster_scrublet.df[j+1,2] <- cluster_scrublet.df[i,3]/(cluster_scrublet.df[i,3]+ cluster_scrublet.df[i+1,3])
  j = j +1 
}

# cluster 41
summary_cluster_scrublet.df <- rbind(summary_cluster_scrublet.df, data.frame(cluster = "41", scrublet_doublets_percent = 1))

# 42 ~ 44
j = 42
for(i in seq(84, by = 2, len = 3)){
  #print(i)
  #print(j)
  summary_cluster_scrublet.df[j+1,1] <- as.character(j)
  summary_cluster_scrublet.df[j+1,2] <- cluster_scrublet.df[i,3]/(cluster_scrublet.df[i,3]+ cluster_scrublet.df[i+1,3])
  j = j +1 
}

summary_cluster_scrublet.df$cluster <- factor(summary_cluster_scrublet.df$cluster, levels = as.character(seq(0,44)))
summary_cluster_scrublet.df$scrublet_singlet_percent <- 1 - summary_cluster_scrublet.df$scrublet_doublets_percent

summary_cluster_scrublet.df$color <- ifelse(summary_cluster_scrublet.df$scrublet_doublets_percent > 0.5, "doublets cluster", "singlet cluster")


g2 <- ggplot(data=summary_cluster_scrublet.df, aes(x=cluster, y=scrublet_doublets_percent, fill = color)) +
  geom_bar(stat="identity", position=position_dodge())+ theme_classic()+ scale_y_continuous(expand = c(0,0)) + 
  ggtitle("percent of doublets by scrublet in each cluster") + 
  labs(y = "doublets percent") + 
  theme(plot.title = element_text(size = 18, hjust = 0.5), axis.text.x = element_text(size =12, color = "black"), axis.title.x = element_text(size =18, color = "black"), axis.text.y = element_text(size =14, color = "black"), axis.title.y = element_text(size =18), legend.title = element_blank(), legend.text = element_text(size = 16), legend.position = c(0.3,0.6)) + 
  scale_fill_manual(values = c("#FF9900","#3399FF"))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig1_doublets_by_scrublet_analysis_bar_plot.png", 
       plot = g2, 
       scale = 1, width = 10, height = 6, units = "in", device = "png",
       dpi = 300)


##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig1_doublets_by_scrublet_analysis.png_session_info.txt")




