#! /bin/env RScript
# written by SJK at 17. Apr. 2020

# This file is for Figure 3 Oligodendrocytes representative figures & volcano plots

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

######################### 
# UMAP 
#########################
# combine meta data and UMAP
data9set.meta.data <- data9set_cleaned.SO@meta.data
data9set_OL_OPC_umap <- as.data.frame(data9set_cleaned.SO@reductions$umap@cell.embeddings)
data9set_OL_OPC_umap.df <- cbind(data9set.meta.data, data9set_OL_OPC_umap)

# Mark Oligo cluster
data9set_OL_OPC_umap.df$celltypes <- ifelse(data9set_OL_OPC_umap.df$seurat_clusters %in% c(1, 5, 6), "Oligo", ifelse(data9set_OL_OPC_umap.df$seurat_clusters %in% c(12, 38, 41, 42, 44), "OPC", "others"))

# subset Oligo cluster
data9set_OL_OPC_sub_umap.df <- data9set_OL_OPC_umap.df[data9set_OL_OPC_umap.df$celltypes %in% "Oligo", ]

# subset OPC cluster
data9set_OL_OPC_sub2_umap.df <- data9set_OL_OPC_umap.df[data9set_OL_OPC_umap.df$celltypes %in% "OPC", ]


# Microglia whole UMAP
g1 <- ggplot(data9set_OL_OPC_umap.df, aes(x = UMAP_1, y = UMAP_2)) + geom_point(color = "gray", size = 0.1) + 
  geom_point(data = data9set_OL_OPC_sub_umap.df[data9set_OL_OPC_sub_umap.df$celltypes %in% "Oligo", ], aes(x = UMAP_1, y = UMAP_2), color = "#99CCCC", size = 0.1) +
  geom_point(data = data9set_OL_OPC_sub2_umap.df[data9set_OL_OPC_sub2_umap.df$celltypes %in% "OPC", ], aes(x = UMAP_1, y = UMAP_2), color = "#009999", size = 0.1) +
  theme_classic() + 
  theme(axis.title = element_text(color = "black", size =15),axis.text = element_text(color = "black", size =12))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_figures_whole.png",
       plot = g1,
       scale = 1, width = 10, height = 10, units = "in", device = "png",
       dpi = 300)


##############################
# volcano plot & scatter plot
##############################

#################
# 1. AD vs Ctrl
#################
OL_DE_AD_Ctrl.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/Cell_type/20200221_Cell_type_DE_by_MAST_OligodendrocytesAD_Ctrl.csv", row.names = 1)

OL_DE_AD_Ctrl.df$gene <- rownames(OL_DE_AD_Ctrl.df)


# log2FC > 1 or padj > 50
OL_DE_AD_Ctrl_plus.df <- OL_DE_AD_Ctrl.df[which(OL_DE_AD_Ctrl.df$avg_logFC > 0.5 | OL_DE_AD_Ctrl.df$avg_logFC > 0 & -log10(OL_DE_AD_Ctrl.df$p_val_adj) > 25),]
# log2FC < 0 and padj > 25
OL_DE_AD_Ctrl_minus.df <- OL_DE_AD_Ctrl.df[which(OL_DE_AD_Ctrl.df$avg_logFC < 0 & -log10(OL_DE_AD_Ctrl.df$p_val_adj) > 25),]
# Il12b gene
OL_DE_AD_Ctrl_Il12b.df <-  OL_DE_AD_Ctrl.df[OL_DE_AD_Ctrl.df$gene == "Il12b", ]

# volcano plot
g1 <- ggplot(OL_DE_AD_Ctrl.df) + geom_point(aes(x = avg_logFC, y = -log10(p_val_adj))) + 
  geom_point(data = OL_DE_AD_Ctrl_plus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), color = "red") +
  geom_point(data = OL_DE_AD_Ctrl_minus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), color = "blue") +
  ggtitle("AD versus Ctrl in Oligodendrocytes") + 
  theme_classic() + 
  xlab("log2 fold change") + 
  ylab("-log10(adjusted p-value)") +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) + 
  theme(legend.position = "none",
        plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12, color = "black")) +
  geom_text_repel(data = OL_DE_AD_Ctrl_plus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), label = OL_DE_AD_Ctrl_plus.df$gene, force = 10) + 
  geom_text_repel(data = OL_DE_AD_Ctrl_minus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), label = OL_DE_AD_Ctrl_minus.df$gene, force = 10) + 
  geom_text_repel(data = OL_DE_AD_Ctrl_Il12b.df, aes(x = avg_logFC, y = -log10(p_val_adj)), label = OL_DE_AD_Ctrl_Il12b.df$gene, nudge_y = 8)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_figures_DE_AD_Ctrl_volcano.png",
       plot = g1,
       scale = 1, width = 5.5, height = 5, units = "in", device = "png",
       dpi = 300)

######################
# 2. AD vs ADp40KO
######################
OL_DE_AD_ADp40KO.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/Cell_type/20200221_Cell_type_DE_by_MAST_OligodendrocytesAD_ADp40KO.csv", row.names = 1)

OL_DE_AD_ADp40KO.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/Cell_type/20200221_Cell_type_DE_by_MAST_OligodendrocytesAD_ADp40KO.csv", row.names = 1)

OL_DE_AD_ADp40KO.df$gene <- rownames(OL_DE_AD_ADp40KO.df)


# log2FC > 1 or padj > 50
OL_DE_AD_ADp40KO_plus.df <- OL_DE_AD_ADp40KO.df[which(OL_DE_AD_ADp40KO.df$avg_logFC > 0.5 | OL_DE_AD_ADp40KO.df$avg_logFC > 0 & -log10(OL_DE_AD_ADp40KO.df$p_val_adj) > 25),]
# log2FC < 0 and padj > 25
OL_DE_AD_ADp40KO_minus.df <- OL_DE_AD_ADp40KO.df[which(OL_DE_AD_ADp40KO.df$avg_logFC < 0 & -log10(OL_DE_AD_ADp40KO.df$p_val_adj) > 25),]
# Il12b gene
OL_DE_AD_ADp40KO_Il12b.df <-  OL_DE_AD_ADp40KO.df[OL_DE_AD_ADp40KO.df$gene == "Il12b", ]

# volcano plot
g2 <- ggplot(OL_DE_AD_ADp40KO.df) + geom_point(aes(x = avg_logFC, y = -log10(p_val_adj))) + 
  geom_point(data = OL_DE_AD_ADp40KO_plus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), color = "red") +
  geom_point(data = OL_DE_AD_ADp40KO_minus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), color = "blue") +
  ggtitle("AD versus ADp40KO in Oligodendrocytes") + 
  theme_classic() + 
  xlab("log2 fold change") + 
  ylab("-log10(adjusted p-value)") +
  scale_y_continuous(expand = c(0,0), limits = c(0,210)) + 
  theme(legend.position = "none",
        plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12, color = "black")) +
  geom_text_repel(data = OL_DE_AD_ADp40KO_plus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), label = OL_DE_AD_ADp40KO_plus.df$gene, force = 10) + 
  geom_text_repel(data = OL_DE_AD_ADp40KO_minus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), label = OL_DE_AD_ADp40KO_minus.df$gene, force = 10) + 
  geom_text_repel(data = OL_DE_AD_ADp40KO_Il12b.df, aes(x = avg_logFC, y = -log10(p_val_adj)), label = OL_DE_AD_ADp40KO_Il12b.df$gene, nudge_y = 8)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_figures_DE_AD_ADp40KO_volcano.png",
       plot = g2,
       scale = 1, width = 5.5, height = 5, units = "in", device = "png",
       dpi = 300)

######################
# 3. ADp40KO vs Ctrl
###################### 
OL_DE_ADp40KO_Ctrl.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/Cell_type/20200221_Cell_type_DE_by_MAST_OligodendrocytesADp40KO_Ctrl.csv", row.names = 1)

OL_DE_ADp40KO_Ctrl.df$gene <- rownames(OL_DE_ADp40KO_Ctrl.df)


# log2FC > 1 or padj > 50
OL_DE_ADp40KO_Ctrl_plus.df <- OL_DE_ADp40KO_Ctrl.df[which(OL_DE_ADp40KO_Ctrl.df$avg_logFC > 0.5 | OL_DE_ADp40KO_Ctrl.df$avg_logFC > 0 & -log10(OL_DE_ADp40KO_Ctrl.df$p_val_adj) > 25),]
# log2FC < 0 and padj > 25
OL_DE_ADp40KO_Ctrl_minus.df <- OL_DE_ADp40KO_Ctrl.df[which(OL_DE_ADp40KO_Ctrl.df$avg_logFC < -0.5 & -log10(OL_DE_ADp40KO_Ctrl.df$p_val_adj) > 25),]


g3 <- ggplot(OL_DE_ADp40KO_Ctrl.df) + geom_point(aes(x = avg_logFC, y = -log10(p_val_adj))) + 
  geom_point(data = OL_DE_ADp40KO_Ctrl_plus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), color = "red") +
  geom_point(data = OL_DE_ADp40KO_Ctrl_minus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), color = "blue") +
  ggtitle("ADp40KO versus Ctrl in Oligodendrocytes") + 
  theme_classic() + 
  xlab("log2 fold change") + 
  ylab("-log10(adjusted p-value)") +
  scale_y_continuous(expand = c(0,0), limits = c(0,150)) + 
  theme(legend.position = "none",
        plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12, color = "black")) +
  geom_text_repel(data = OL_DE_ADp40KO_Ctrl_plus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), label = OL_DE_ADp40KO_Ctrl_plus.df$gene, force = 10) + 
  geom_text_repel(data = OL_DE_ADp40KO_Ctrl_minus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), label = OL_DE_ADp40KO_Ctrl_minus.df$gene, force = 10)


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_figures_DE_ADp40KO_Ctrl_volcano.png",
       plot = g3,
       scale = 1, width = 5.5, height = 5, units = "in", device = "png",
       dpi = 300)

############################
# 4. AD vs Ctrl scatter plot
############################

data9set_OL.SO <- subset(data9set_cleaned.SO, subset = seurat_clusters %in% c(1, 5, 6))


# split by sample
data9set_cleaned_OL_Ctrl.SO <- subset(data9set_OL.SO, subset = sample %in% "Ctrl")
data9set_cleaned_OL_AD.SO <- subset(data9set_OL.SO, subset = sample %in% "AD")
data9set_cleaned_OL_ADp40KO.SO <- subset(data9set_OL.SO, subset = sample %in% "ADp40KO")

# extract UMI counts and calculate average per each gene
OL_Ctrl_count.df <- data.frame((GetAssayData(data9set_cleaned_OL_Ctrl.SO, slot = "counts")))
OL_AD_count.df <- data.frame((GetAssayData(data9set_cleaned_OL_AD.SO, slot = "counts")))
OL_ADp40KO_count.df <- data.frame((GetAssayData(data9set_cleaned_OL_ADp40KO.SO, slot = "counts")))

OL_Ctrl_count_cpm.df <- apply(OL_Ctrl_count.df, 2, function(x) (x/sum(x)) * 1000000)
OL_AD_count_cpm.df <- apply(OL_AD_count.df, 2, function(x) (x/sum(x)) * 1000000)
OL_ADp40KO_count_cpm.df <- apply(OL_ADp40KO_count.df, 2, function(x) (x/sum(x)) * 1000000)

OL_Ctrl_count_cpm.df <- data.frame(Matrix::rowMeans(OL_Ctrl_count_cpm.df))
OL_AD_count_cpm.df <- data.frame(Matrix::rowMeans(OL_AD_count_cpm.df))
OL_ADp40KO_count_cpm.df <- data.frame(Matrix::rowMeans(OL_ADp40KO_count_cpm.df))


OL_count_cpm.df <- cbind(OL_Ctrl_count_cpm.df, OL_AD_count_cpm.df, OL_ADp40KO_count_cpm.df)
OL_count_cpm.df <- as.data.frame(OL_count_cpm.df)

# remove genes where all zero in samples
OL_count_cpm.df <- OL_count_cpm.df[apply(OL_count_cpm.df == 0, 1, sum) != 3, ]
colnames(OL_count_cpm.df) <- c("Ctrl", "AD", "ADp40KO")
# add pseudocount 0.001 and logarithm 
OL_count_cpm.df <- log(OL_count_cpm.df + 1)
OL_count_cpm.df$gene <- rownames(OL_count_cpm.df)

# extract genes from volcano plots
OL_count_cpm_plus.df <- OL_count_cpm.df[OL_count_cpm.df$gene %in% OL_DE_AD_Ctrl_plus.df$gene, ]
OL_count_cpm_minus.df <- OL_count_cpm.df[OL_count_cpm.df$gene %in% OL_DE_AD_Ctrl_minus.df$gene, ]
OL_count_cpm_Il12b.df <- OL_count_cpm.df[OL_count_cpm.df$gene %in% OL_DE_AD_Ctrl_Il12b.df$gene, ]

# scatter plot for AD vs Ctrl
g4 <- ggplot(OL_count_cpm.df, aes(x = Ctrl, y = AD)) + geom_point(size = 0.5) +
  geom_point(data = OL_count_cpm_plus.df,  aes(x = Ctrl, y = AD), color = "red", size = 0.8) +
  geom_point(data = OL_count_cpm_minus.df,  aes(x = Ctrl, y = AD), color = "blue", size = 0.8) +
  ggtitle("AD versus Ctrl in Oligodendrocytes") +
  theme_classic() +
  xlab("log (average CPM + 1) Ctrl") +
  ylab("log (average CPM + 1) AD") +
  theme(plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12, color = "black")) +
  geom_abline(slope = 1, color="red", linetype="dashed", size=.5) +
  geom_text_repel(data = OL_count_cpm_plus.df, aes(x = Ctrl, y = AD), label = OL_count_cpm_plus.df$gene, force = 30, nudge_y = .8) +
  geom_text_repel(data = OL_count_cpm_minus.df, aes(x = Ctrl, y = AD), label = OL_count_cpm_minus.df$gene, force = 40, nudge_y = -.6) +
  geom_text_repel(data = OL_count_cpm_Il12b.df, aes(x = Ctrl, y = AD), label = OL_count_cpm_Il12b.df$gene, nudge_y = 1, nudge_x = -0.1)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_figures_DE_AD_Ctrl_scatter.png",
       plot = g4,
       scale = 1, width = 5.5, height = 5, units = "in", device = "png",
       dpi = 300)

############################
# 5. AD vs ADp40KO scatter plot
############################

OL_count_cpm_plus.df <- OL_count_cpm.df[OL_count_cpm.df$gene %in% OL_DE_AD_ADp40KO_plus.df$gene, ]
OL_count_cpm_minus.df <- OL_count_cpm.df[OL_count_cpm.df$gene %in% OL_DE_AD_ADp40KO_minus.df$gene, ]
OL_count_cpm_Il12b.df <- OL_count_cpm.df[OL_count_cpm.df$gene %in% OL_DE_AD_ADp40KO_Il12b.df$gene, ]

# scatter plot for AD vs Ctrl
g5 <- ggplot(OL_count_cpm.df, aes(x = ADp40KO, y = AD)) + geom_point(size = 0.5) +
  geom_point(data = OL_count_cpm_plus.df,  aes(x = ADp40KO, y = AD), color = "red", size = 0.8) +
  geom_point(data = OL_count_cpm_minus.df,  aes(x = ADp40KO, y = AD), color = "blue", size = 0.8) +
  ggtitle("AD versus ADp40KO in Oligodendrocytes") +
  theme_classic() +
  xlab("log (average CPM + 1) ADp40KO") +
  ylab("log (average CPM + 1) AD") +
  theme(plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12, color = "black")) +
  geom_abline(slope = 1, color="red", linetype="dashed", size=.5) +
  geom_text_repel(data = OL_count_cpm_plus.df, aes(x = ADp40KO, y = AD), label = OL_count_cpm_plus.df$gene, force = 10, nudge_y = .8) +
  geom_text_repel(data = OL_count_cpm_minus.df, aes(x = ADp40KO, y = AD), label = OL_count_cpm_minus.df$gene, force = 10, nudge_y = -.6) +
  geom_text_repel(data = OL_count_cpm_Il12b.df, aes(x = ADp40KO, y = AD), label = OL_count_cpm_Il12b.df$gene, nudge_y = 1, nudge_x = -0.1)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_figures_DE_AD_ADp40KO_scatter.png",
       plot = g5,
       scale = 1, width = 5.5, height = 5, units = "in", device = "png",
       dpi = 300)

############################
# 5. ADp40KO vs Ctrl scatter plot
############################

# extract genes from volcano plots
OL_count_cpm_plus.df <- OL_count_cpm.df[OL_count_cpm.df$gene %in% OL_DE_ADp40KO_Ctrl_plus.df$gene, ]
OL_count_cpm_minus.df <- OL_count_cpm.df[OL_count_cpm.df$gene %in% OL_DE_ADp40KO_Ctrl_minus.df$gene, ]


# scatter plot
g6 <- ggplot(OL_count_cpm.df, aes(x = Ctrl, y = ADp40KO)) + geom_point(size = 0.5) +
  geom_point(data = OL_count_cpm_plus.df,  aes(x = Ctrl, y = ADp40KO), color = "red", size = 0.8) +
  geom_point(data = OL_count_cpm_minus.df,  aes(x = Ctrl, y = ADp40KO), color = "blue", size = 0.8) +
  ggtitle("ADp40KO versus Ctrl in Oligodendrocytes") +
  theme_classic() +
  xlab("log (average CPM + 1) Ctrl") +
  ylab("log (average CPM + 1) ADp40KO") +
  theme(plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12, color = "black")) +
  geom_abline(slope = 1, color="red", linetype="dashed", size=.5) +
  geom_text_repel(data = OL_count_cpm_plus.df, aes(x = Ctrl, y = ADp40KO), label = OL_count_cpm_plus.df$gene, force = 10, nudge_y = .3) +
  geom_text_repel(data = OL_count_cpm_minus.df, aes(x = Ctrl, y = ADp40KO), label = OL_count_cpm_minus.df$gene, force = 10, nudge_y = -.3) 

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_figures_DE_ADp40KO_Ctrl_scatter.png",
       plot = g6,
       scale = 1, width = 5.5, height = 5, units = "in", device = "png",
       dpi = 300)

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_figures_session_info.txt")



