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

load(file="/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200210_figures_for_meeting_report_cleaned_Seurat_object_res0.8.Robj")
data9set_cleaned_MG.SO <- subset(data9set_cleaned.SO, subset = seurat_clusters %in% c(3, 8))

#############################
# loading AD vs Ctrl DE file
#############################
# In here, I subset absolute avg_logFC > 0.2 & p_val_adj < 0.01 for figures

MG_DE_AD_Ctrl.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/Cell_type/20200221_Cell_type_DE_by_MAST_MicrogliaAD_Ctrl.csv", row.names = 1)
MG_DE_AD_Ctrl.df <- MG_DE_AD_Ctrl.df[which(abs(MG_DE_AD_Ctrl.df$avg_logFC) > 0.2 & MG_DE_AD_Ctrl.df$p_val_adj < 0.01), ]
MG_DE_AD_Ctrl.df$gene <- rownames(MG_DE_AD_Ctrl.df)

# log2FC > 1 or padj > 50
MG_DE_AD_Ctrl_plus.df <- MG_DE_AD_Ctrl.df[which(MG_DE_AD_Ctrl.df$avg_logFC > 1 | MG_DE_AD_Ctrl.df$avg_logFC > 0 & -log10(MG_DE_AD_Ctrl.df$p_val_adj) > 50),]
# log2FC < 0 and padj > 25
MG_DE_AD_Ctrl_minus.df <- MG_DE_AD_Ctrl.df[which(MG_DE_AD_Ctrl.df$avg_logFC < 0 & -log10(MG_DE_AD_Ctrl.df$p_val_adj) > 25),]
# Il12b gene
MG_DE_AD_Ctrl_Il12b.df <-  MG_DE_AD_Ctrl.df[MG_DE_AD_Ctrl.df$gene == "Il12b", ]

# volcano plot
g1 <- ggplot(MG_DE_AD_Ctrl.df) + geom_point(aes(x = avg_logFC, y = -log10(p_val_adj))) + 
  geom_point(data = MG_DE_AD_Ctrl_plus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), color = "red") +
  geom_point(data = MG_DE_AD_Ctrl_minus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), color = "blue") +
  ggtitle("AD versus Ctrl in Microglia") + 
  theme_classic() + 
  xlab("log2 fold change") + 
  ylab("-log10(adjusted p-value)") +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) + 
  theme(legend.position = "none",
        plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12, color = "black")) +
  geom_text_repel(data = MG_DE_AD_Ctrl_plus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), label = MG_DE_AD_Ctrl_plus.df$gene, force = 10) + 
  geom_text_repel(data = MG_DE_AD_Ctrl_minus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), label = MG_DE_AD_Ctrl_minus.df$gene, force = 10) + 
  geom_text_repel(data = MG_DE_AD_Ctrl_Il12b.df, aes(x = avg_logFC, y = -log10(p_val_adj)), label = MG_DE_AD_Ctrl_Il12b.df$gene, nudge_y = 8)

#############################
# loading AD vs ADp40KO DE file
#############################
# In here, I subset absolute avg_logFC > 0.1 & p_val_adj < 0.01 for figures
MG_DE_AD_ADp40KO.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/Cell_type/20200221_Cell_type_DE_by_MAST_MicrogliaAD_ADp40KO.csv", row.names = 1)
MG_DE_AD_ADp40KO.df <- MG_DE_AD_ADp40KO.df[which(abs(MG_DE_AD_ADp40KO.df$avg_logFC) > 0.1 & MG_DE_AD_ADp40KO.df$p_val_adj < 0.01), ]
MG_DE_AD_ADp40KO.df$gene <- rownames(MG_DE_AD_ADp40KO.df)

# log2FC > 0.3 or padj > 50
MG_DE_AD_ADp40KO_plus.df <- MG_DE_AD_ADp40KO.df[which(MG_DE_AD_ADp40KO.df$avg_logFC > 0.3 | MG_DE_AD_ADp40KO.df$avg_logFC > 0 & -log10(MG_DE_AD_ADp40KO.df$p_val_adj) > 50),]
# padj > 10
MG_DE_AD_ADp40KO_minus.df <- MG_DE_AD_ADp40KO.df[which(MG_DE_AD_ADp40KO.df$avg_logFC < 0 & -log10(MG_DE_AD_ADp40KO.df$p_val_adj) > 10),]
# Il12b gene
MG_DE_AD_ADp40KO_Il12b.df <-  MG_DE_AD_ADp40KO.df[MG_DE_AD_ADp40KO.df$gene == "Il12b", ]

# volcano plot
# I used scale_y_log10 as Ttr is only high number of padj
g2 <- ggplot(MG_DE_AD_ADp40KO.df) + geom_point(aes(x = avg_logFC, y = -log10(p_val_adj))) + 
  geom_point(data = MG_DE_AD_ADp40KO_plus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), color = "red") +
  geom_point(data = MG_DE_AD_ADp40KO_minus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), color = "blue") +
  ggtitle("AD versus ADp40KO in Microglia") + 
  theme_classic() + 
  xlab("log2 fold change") + 
  ylab("-log10(adjusted p-value)") +
  scale_y_log10(expand = c(0,0), breaks = c(5,10,15,20, 50, 130), limits = c(1,130)) + 
  theme(legend.position = "none",
        plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12, color = "black")) +
  geom_text_repel(data = MG_DE_AD_ADp40KO_plus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), label = MG_DE_AD_ADp40KO_plus.df$gene, force = 20, nudge_y = -0.1) + 
  geom_text_repel(data = MG_DE_AD_ADp40KO_minus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), label = MG_DE_AD_ADp40KO_minus.df$gene, force = 10) + 
  geom_text_repel(data = MG_DE_AD_ADp40KO_Il12b.df, aes(x = avg_logFC, y = -log10(p_val_adj)), label = MG_DE_AD_ADp40KO_Il12b.df$gene, nudge_y = 0.1)

g12 <- arrangeGrob(g1, g2, ncol = 2)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_DE_figures_volcano.png",
       plot = g12,
       scale = 1, width = 11, height = 5, units = "in", device = "png",
       dpi = 300)

#############################
# loading ADp40KO vs Ctrl DE file
#############################
# In here, I subset absolute avg_logFC > 0.1 & p_val_adj < 0.01 for figures
MG_DE_ADp40KO_Ctrl.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/Cell_type/20200221_Cell_type_DE_by_MAST_MicrogliaADp40KO_Ctrl.csv", row.names = 1)
MG_DE_ADp40KO_Ctrl.df <- MG_DE_ADp40KO_Ctrl.df[which(abs(MG_DE_ADp40KO_Ctrl.df$avg_logFC) > 0.2 & MG_DE_ADp40KO_Ctrl.df$p_val_adj < 0.01), ]
MG_DE_ADp40KO_Ctrl.df$gene <- rownames(MG_DE_ADp40KO_Ctrl.df)

# log2FC > 0.3 or padj > 50
MG_DE_ADp40KO_Ctrl_plus.df <- MG_DE_ADp40KO_Ctrl.df[which(MG_DE_ADp40KO_Ctrl.df$avg_logFC > 1 | MG_DE_ADp40KO_Ctrl.df$avg_logFC > 0 & -log10(MG_DE_ADp40KO_Ctrl.df$p_val_adj) > 50),]
# padj > 10
MG_DE_ADp40KO_Ctrl_minus.df <- MG_DE_ADp40KO_Ctrl.df[which(MG_DE_ADp40KO_Ctrl.df$avg_logFC < 0 & -log10(MG_DE_ADp40KO_Ctrl.df$p_val_adj) > 25),]
# Il12b gene (no Il12b DE genes)
#MG_DE_ADp40KO_Ctrl_Il12b.df <-  MG_DE_ADp40KO_Ctrl.df[MG_DE_ADp40KO_Ctrl.df$gene == "Il12b", ]

g3 <- ggplot(MG_DE_ADp40KO_Ctrl.df) + geom_point(aes(x = avg_logFC, y = -log10(p_val_adj))) + 
  geom_point(data = MG_DE_ADp40KO_Ctrl_plus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), color = "red") +
  geom_point(data = MG_DE_ADp40KO_Ctrl_minus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), color = "blue") +
  ggtitle("ADp40KO versus Ctrl in Microglia") + 
  theme_classic() + 
  xlab("log2 fold change") + 
  ylab("-log10(adjusted p-value)") +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +   
  theme(legend.position = "none",
        plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12, color = "black")) +
  geom_text_repel(data = MG_DE_ADp40KO_Ctrl_plus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), label = MG_DE_ADp40KO_Ctrl_plus.df$gene, force = 20, nudge_y = -0.1) + 
  geom_text_repel(data = MG_DE_ADp40KO_Ctrl_minus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), label = MG_DE_ADp40KO_Ctrl_minus.df$gene, force = 10) 

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_DE_ADp40KO_Ctrl_figures_volcano.png",
       plot = g3,
       scale = 1, width = 5.5, height = 5, units = "in", device = "png",
       dpi = 300)

################################
# scatter plot of UMI count
################################
# Scatterplot showing the average molecules count (log scale) of sample compared with other sample (y axis).
# Pseudocount is 0.001

# split by sample
data9set_cleaned_MG_Ctrl.SO <- subset(data9set_cleaned_MG.SO, subset = sample %in% "Ctrl")
data9set_cleaned_MG_AD.SO <- subset(data9set_cleaned_MG.SO, subset = sample %in% "AD")
data9set_cleaned_MG_ADp40KO.SO <- subset(data9set_cleaned_MG.SO, subset = sample %in% "ADp40KO")

# extract UMI counts and calculate average per each gene
MG_Ctrl_count.df <- data.frame(Matrix::rowMeans(GetAssayData(data9set_cleaned_MG_Ctrl.SO, slot = "counts")))
MG_AD_count.df <- data.frame(Matrix::rowMeans(GetAssayData(data9set_cleaned_MG_AD.SO, slot = "counts")))
MG_ADp40KO_count.df <- data.frame(Matrix::rowMeans(GetAssayData(data9set_cleaned_MG_ADp40KO.SO, slot = "counts")))

MG_count.df <- cbind(MG_Ctrl_count.df, MG_AD_count.df, MG_ADp40KO_count.df)



# remove genes where all zero in samples
MG_count.df <- MG_count.df[apply(MG_count.df == 0, 1, sum) != 3, ]
colnames(MG_count.df) <- c("Ctrl", "AD", "ADp40KO")
# add pseudocount 0.001 and logarithm 
MG_count.df <- log(MG_count.df + 0.001)

MG_count.df$gene <- rownames(MG_count.df)
MG_count.df <- as.data.frame(MG_count.df)

# extract genes from volcano plots
MG_count_plus.df <- MG_count.df[MG_count.df$gene %in% MG_DE_AD_Ctrl_plus.df$gene, ]
MG_count_minus.df <- MG_count.df[MG_count.df$gene %in% MG_DE_AD_Ctrl_minus.df$gene, ]
MG_count_Il12b.df <- MG_count.df[MG_count.df$gene %in% MG_DE_AD_Ctrl_Il12b.df$gene, ]

# scatter plot for AD vs Ctrl
g4 <- ggplot(MG_count.df, aes(x = Ctrl, y = AD)) + geom_point(size = 0.5) + 
  geom_point(data = MG_count_plus.df,  aes(x = Ctrl, y = AD), color = "red", size = 0.8) +
  geom_point(data = MG_count_minus.df,  aes(x = Ctrl, y = AD), color = "blue", size = 0.8) +
  ggtitle("AD versus Ctrl in Microglia") + 
  theme_classic() + 
  xlab("log (average UMI counts) Ctrl") + 
  ylab("log (average UMI counts) AD") + 
  theme(plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12, color = "black")) +
  geom_abline(slope = 1, color="red", linetype="dashed", size=.5) + 
  geom_text_repel(data = MG_count_plus.df, aes(x = Ctrl, y = AD), label = MG_count_plus.df$gene, force = 30, nudge_y = .8) + 
  geom_text_repel(data = MG_count_minus.df, aes(x = Ctrl, y = AD), label = MG_count_minus.df$gene, force = 40, nudge_y = -.6) + 
  geom_text_repel(data = MG_count_Il12b.df, aes(x = Ctrl, y = AD), label = MG_count_Il12b.df$gene, nudge_y = 3, nudge_x = -0.1)

# extract genes from volcano plots
MG_count_plus.df2 <- MG_count.df[MG_count.df$gene %in% MG_DE_AD_ADp40KO_plus.df$gene, ]
MG_count_minus.df2 <- MG_count.df[MG_count.df$gene %in% MG_DE_AD_ADp40KO_minus.df$gene, ]
MG_count_Il12b.df2 <- MG_count.df[MG_count.df$gene %in% MG_DE_AD_ADp40KO_Il12b.df$gene, ]

# scatter plot for AD vs ADp40KO
g5 <- ggplot(MG_count.df, aes(x = ADp40KO, y = AD)) + geom_point(size = 0.5) + 
  geom_point(data = MG_count_plus.df2,  aes(x = ADp40KO, y = AD), color = "red", size = 0.8) +
  geom_point(data = MG_count_minus.df2,  aes(x = ADp40KO, y = AD), color = "blue", size = 0.8) +
  ggtitle("AD versus ADp40KO in Microglia") + 
  theme_classic() + 
  xlab("log (average UMI counts) ADp40KO") + 
  ylab("log (average UMI counts) AD") + 
  theme(plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12, color = "black")) +
  geom_abline(slope = 1, color="red", linetype="dashed", size=.5) + 
  geom_text_repel(data = MG_count_plus.df2, aes(x = ADp40KO, y = AD), label = MG_count_plus.df2$gene, force = 50, nudge_y = 1.5, nudge_x = -1.2, max.iter = 3000) + 
  geom_text_repel(data = MG_count_minus.df2, aes(x = ADp40KO, y = AD), label = MG_count_minus.df2$gene, force = 40, nudge_y = -.8) + 
  geom_text_repel(data = MG_count_Il12b.df2, aes(x = ADp40KO, y = AD), label = MG_count_Il12b.df2$gene, nudge_y = 1, nudge_x = -0.1)

g45 <- arrangeGrob(g4, g5, ncol = 2)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_DE_figures_scatter.png",
       plot = g45,
       scale = 1, width = 11, height = 5, units = "in", device = "png",
       dpi = 300)

# extract genes from volcano plots
MG_count_plus.df3 <- MG_count.df[MG_count.df$gene %in% MG_DE_ADp40KO_Ctrl_plus.df$gene, ]
MG_count_minus.df3 <- MG_count.df[MG_count.df$gene %in% MG_DE_ADp40KO_Ctrl_minus.df$gene, ]

g6 <- ggplot(MG_count.df, aes(x = Ctrl, y = ADp40KO)) + geom_point(size = 0.5) + 
  geom_point(data = MG_count_plus.df3,  aes(x = Ctrl, y = ADp40KO), color = "red", size = 0.8) +
  geom_point(data = MG_count_minus.df3,  aes(x = Ctrl, y = ADp40KO), color = "blue", size = 0.8) +
  ggtitle("ADp40KO versus Ctrl in Microglia") + 
  theme_classic() + 
  xlab("log (average UMI counts) Ctrl") + 
  ylab("log (average UMI counts) ADp40KO") + 
  theme(plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12, color = "black")) +
  geom_abline(slope = 1, color="red", linetype="dashed", size=.5) + 
  geom_text_repel(data = MG_count_plus.df3, aes(x = Ctrl, y = ADp40KO), label = MG_count_plus.df3$gene, force = 10, nudge_y = .8) + 
  geom_text_repel(data = MG_count_minus.df3, aes(x = Ctrl, y = ADp40KO), label = MG_count_minus.df3$gene, force = 10, nudge_y = -.6)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_DE_ADp40KO_Ctrl_figures_scatter.png",
       plot = g6,
       scale = 1, width = 5.5, height = 5, units = "in", device = "png",
       dpi = 300)

######################################
# scatter plot by logarithmic CPM + 1
######################################
# extract log(LogNormalize + 1) and calculate average per each gene
MG_Ctrl_cpm.df <- data.frame(Matrix::rowMeans(GetAssayData(data9set_cleaned_MG_Ctrl.SO, slot = "data")))
MG_AD_cpm.df <- data.frame(Matrix::rowMeans(GetAssayData(data9set_cleaned_MG_AD.SO, slot = "data")))
MG_ADp40KO_cpm.df <- data.frame(Matrix::rowMeans(GetAssayData(data9set_cleaned_MG_ADp40KO.SO, slot = "data")))

MG_cpm.df <- cbind(MG_Ctrl_cpm.df, MG_AD_cpm.df, MG_ADp40KO_cpm.df)

# remove genes where all zero in samples
MG_cpm.df <- MG_cpm.df[apply(MG_cpm.df == 0, 1, sum) != 3, ]
colnames(MG_cpm.df) <- c("Ctrl", "AD", "ADp40KO")
# No need to add pseudocount 0.001 and logarithm (Already done)

MG_cpm.df$gene <- rownames(MG_cpm.df)
MG_cpm.df <- as.data.frame(MG_cpm.df)

# AD vs Ctdrl
# extract genes from volcano plots
MG_cpm_plus.df <- MG_cpm.df[MG_cpm.df$gene %in% MG_DE_AD_Ctrl_plus.df$gene, ]
MG_cpm_minus.df <- MG_cpm.df[MG_cpm.df$gene %in% MG_DE_AD_Ctrl_minus.df$gene, ]
MG_cpm_Il12b.df <- MG_cpm.df[MG_cpm.df$gene %in% MG_DE_AD_Ctrl_Il12b.df$gene, ]

g7 <- ggplot(MG_cpm.df, aes(x = Ctrl, y = AD)) + geom_point(size = 0.5) + 
  geom_point(data = MG_cpm_plus.df,  aes(x = Ctrl, y = AD), color = "red", size = 0.8) +
  geom_point(data = MG_cpm_minus.df,  aes(x = Ctrl, y = AD), color = "blue", size = 0.8) +
  ggtitle("AD versus Ctrl in Microglia") + 
  theme_classic() + 
  xlab("log (average Normalize counts + 1) Ctrl") + 
  ylab("log (average Normalize counts + 1) AD") + 
  theme(plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12, color = "black")) +
  geom_abline(slope = 1, color="red", linetype="dashed", size=.5) + 
  geom_text_repel(data = MG_cpm_plus.df, aes(x = Ctrl, y = AD), label = MG_cpm_plus.df$gene, force = 5, nudge_y = .5) + 
  geom_text_repel(data = MG_cpm_minus.df, aes(x = Ctrl, y = AD), label = MG_cpm_minus.df$gene, force = 5, nudge_y = -.6) + 
  geom_text_repel(data = MG_cpm_Il12b.df, aes(x = Ctrl, y = AD), label = MG_cpm_Il12b.df$gene, nudge_y = 1, nudge_x = -0.1)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_DE_AD_Ctrl_figures_Normalize_scatter.png",
       plot = g7,
       scale = 1, width = 5.5, height = 5, units = "in", device = "png",
       dpi = 300)


# AD vs Ctrl
# extract genes from volcano plots
MG_cpm_plus.df2 <- MG_cpm.df[MG_cpm.df$gene %in% MG_DE_AD_ADp40KO_plus.df$gene, ]
MG_cpm_minus.df2 <- MG_cpm.df[MG_cpm.df$gene %in% MG_DE_AD_ADp40KO_minus.df$gene, ]
MG_cpm_Il12b.df2 <- MG_cpm.df[MG_cpm.df$gene %in% MG_DE_AD_ADp40KO_Il12b.df$gene, ]

g8 <- ggplot(MG_cpm.df, aes(x = ADp40KO, y = AD)) + geom_point(size = 0.5) + 
  geom_point(data = MG_cpm_plus.df2,  aes(x = ADp40KO, y = AD), color = "red", size = 0.8) +
  geom_point(data = MG_cpm_minus.df2,  aes(x = ADp40KO, y = AD), color = "blue", size = 0.8) +
  ggtitle("AD versus ADp40KO in Microglia") + 
  theme_classic() + 
  xlab("log (average Normalize counts + 1) ADp40KO") + 
  ylab("log (average Normalize counts + 1) AD") + 
  theme(plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12, color = "black")) +
  geom_abline(slope = 1, color="red", linetype="dashed", size=.5) + 
  geom_text_repel(data = MG_cpm_plus.df2, aes(x = ADp40KO, y = AD), label = MG_cpm_plus.df2$gene, force = 10, nudge_y = 0.5, nudge_x = -.2) + 
  geom_text_repel(data = MG_cpm_minus.df2, aes(x = ADp40KO, y = AD), label = MG_cpm_minus.df2$gene, force = 10, nudge_y = -.5) + 
  geom_text_repel(data = MG_cpm_Il12b.df2, aes(x = ADp40KO, y = AD), label = MG_cpm_Il12b.df2$gene, nudge_y = 1, nudge_x = -0.1)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_DE_AD_ADp40KO_figures_Normalize_scatter.png",
       plot = g8,
       scale = 1, width = 5.5, height = 5, units = "in", device = "png",
       dpi = 300)

# ADp40KO vs Ctrl
# extract genes from volcano plots
MG_cpm_plus.df3 <- MG_cpm.df[MG_cpm.df$gene %in% MG_DE_ADp40KO_Ctrl_plus.df$gene, ]
MG_cpm_minus.df3 <- MG_cpm.df[MG_cpm.df$gene %in% MG_DE_ADp40KO_Ctrl_minus.df$gene, ]

g9 <- ggplot(MG_cpm.df, aes(x = Ctrl, y = ADp40KO)) + geom_point(size = 0.5) + 
  geom_point(data = MG_cpm_plus.df3,  aes(x = Ctrl, y = ADp40KO), color = "red", size = 0.8) +
  geom_point(data = MG_cpm_minus.df3,  aes(x = Ctrl, y = ADp40KO), color = "blue", size = 0.8) +
  ggtitle("ADp40KO versus Ctrl in Microglia") + 
  theme_classic() + 
  xlab("log (average Normalize counts + 1) Ctrl") + 
  ylab("log (average Normalize counts + 1) ADp40KO") + 
  theme(plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12, color = "black")) +
  geom_abline(slope = 1, color="red", linetype="dashed", size=.5) + 
  geom_text_repel(data = MG_cpm_plus.df3, aes(x = Ctrl, y = ADp40KO), label = MG_cpm_plus.df3$gene, force = 10, nudge_y = .5) + 
  geom_text_repel(data = MG_cpm_minus.df3, aes(x = Ctrl, y = ADp40KO), label = MG_cpm_minus.df3$gene, force = 10, nudge_y = -.6)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_DE_ADp40KO_Ctrl_figures_Normalize_scatter.png",
       plot = g9,
       scale = 1, width = 5.5, height = 5, units = "in", device = "png",
       dpi = 300)

########################################
# scatter plot by logarithmic CPM + 1
########################################

# extract UMI counts and calculate average per each gene
MG_Ctrl_count.df <- data.frame((GetAssayData(data9set_cleaned_MG_Ctrl.SO, slot = "counts")))
MG_AD_count.df <- data.frame((GetAssayData(data9set_cleaned_MG_AD.SO, slot = "counts")))
MG_ADp40KO_count.df <- data.frame((GetAssayData(data9set_cleaned_MG_ADp40KO.SO, slot = "counts")))

MG_Ctrl_count_cpm.df <- apply(MG_Ctrl_count.df, 2, function(x) (x/sum(x)) * 1000000)
MG_AD_count_cpm.df <- apply(MG_AD_count.df, 2, function(x) (x/sum(x)) * 1000000)
MG_ADp40KO_count_cpm.df <- apply(MG_ADp40KO_count.df, 2, function(x) (x/sum(x)) * 1000000)

MG_Ctrl_count_cpm.df <- data.frame(Matrix::rowMeans(MG_Ctrl_count_cpm.df))
MG_AD_count_cpm.df <- data.frame(Matrix::rowMeans(MG_AD_count_cpm.df))
MG_ADp40KO_count_cpm.df <- data.frame(Matrix::rowMeans(MG_ADp40KO_count_cpm.df))


MG_count_cpm.df <- cbind(MG_Ctrl_count_cpm.df, MG_AD_count_cpm.df, MG_ADp40KO_count_cpm.df)
MG_count_cpm.df <- as.data.frame(MG_count_cpm.df)

# remove genes where all zero in samples
MG_count_cpm.df <- MG_count_cpm.df[apply(MG_count_cpm.df == 0, 1, sum) != 3, ]
colnames(MG_count_cpm.df) <- c("Ctrl", "AD", "ADp40KO")
# add pseudocount 0.001 and logarithm 
MG_count_cpm.df <- log(MG_count_cpm.df + 1)
MG_count_cpm.df$gene <- rownames(MG_count_cpm.df)


# extract genes from volcano plots
MG_count_cpm_plus.df <- MG_count_cpm.df[MG_count_cpm.df$gene %in% MG_DE_AD_Ctrl_plus.df$gene, ]
MG_count_cpm_minus.df <- MG_count_cpm.df[MG_count_cpm.df$gene %in% MG_DE_AD_Ctrl_minus.df$gene, ]
MG_count_cpm_Il12b.df <- MG_count_cpm.df[MG_count_cpm.df$gene %in% MG_DE_AD_Ctrl_Il12b.df$gene, ]

# scatter plot for AD vs Ctrl
g10 <- ggplot(MG_count_cpm.df, aes(x = Ctrl, y = AD)) + geom_point(size = 0.5) +
  geom_point(data = MG_count_cpm_plus.df,  aes(x = Ctrl, y = AD), color = "red", size = 0.8) +
  geom_point(data = MG_count_cpm_minus.df,  aes(x = Ctrl, y = AD), color = "blue", size = 0.8) +
  ggtitle("AD versus Ctrl in Microglia") +
  theme_classic() +
  xlab("log (average CPM + 1) Ctrl") +
  ylab("log (average CPM + 1) AD") +
  theme(plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12, color = "black")) +
  geom_abline(slope = 1, color="red", linetype="dashed", size=.5) +
  geom_text_repel(data = MG_count_cpm_plus.df, aes(x = Ctrl, y = AD), label = MG_count_cpm_plus.df$gene, force = 30, nudge_y = .8) +
  geom_text_repel(data = MG_count_cpm_minus.df, aes(x = Ctrl, y = AD), label = MG_count_cpm_minus.df$gene, force = 40, nudge_y = -.6) +
  geom_text_repel(data = MG_count_cpm_Il12b.df, aes(x = Ctrl, y = AD), label = MG_count_cpm_Il12b.df$gene, nudge_y = 3, nudge_x = -0.1)

# extract genes from volcano plots
MG_count_cpm_plus.df2 <- MG_count_cpm.df[MG_count_cpm.df$gene %in% MG_DE_AD_ADp40KO_plus.df$gene, ]
MG_count_cpm_minus.df2 <- MG_count_cpm.df[MG_count_cpm.df$gene %in% MG_DE_AD_ADp40KO_minus.df$gene, ]
MG_count_cpm_Il12b.df2 <- MG_count_cpm.df[MG_count_cpm.df$gene %in% MG_DE_AD_ADp40KO_Il12b.df$gene, ]

# scatter plot for AD vs ADp40KO
g11 <- ggplot(MG_count_cpm.df, aes(x = ADp40KO, y = AD)) + geom_point(size = 0.5) +
  geom_point(data = MG_count_cpm_plus.df2,  aes(x = ADp40KO, y = AD), color = "red", size = 0.8) +
  geom_point(data = MG_count_cpm_minus.df2,  aes(x = ADp40KO, y = AD), color = "blue", size = 0.8) +
  ggtitle("AD versus ADp40KO in Microglia") +
  theme_classic() +
  xlab("log (average CPM + 1) ADp40KO") +
  ylab("log (average CPM + 1) AD") +
  theme(plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12, color = "black")) +
  geom_abline(slope = 1, color="red", linetype="dashed", size=.5) +
  geom_text_repel(data = MG_count_cpm_plus.df2, aes(x = ADp40KO, y = AD), label = MG_count_cpm_plus.df2$gene, force = 50, nudge_y = 1.5, nudge_x = -1.2, max.iter = 3000) +
  geom_text_repel(data = MG_count_cpm_minus.df2, aes(x = ADp40KO, y = AD), label = MG_count_cpm_minus.df2$gene, force = 40, nudge_y = -.8) +
  geom_text_repel(data = MG_count_cpm_Il12b.df2, aes(x = ADp40KO, y = AD), label = MG_count_cpm_Il12b.df2$gene, nudge_y = 1, nudge_x = -0.1)

g1011 <- arrangeGrob(g10, g11, ncol = 2)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_DE_figures_CPM_scatter.png",
       plot = g1011,
       scale = 1, width = 11, height = 5, units = "in", device = "png",
       dpi = 300)

# extract genes from volcano plots
MG_count_cpm_plus.df3 <- MG_count_cpm.df[MG_count_cpm.df$gene %in% MG_DE_ADp40KO_Ctrl_plus.df$gene, ]
MG_count_cpm_minus.df3 <- MG_count_cpm.df[MG_count_cpm.df$gene %in% MG_DE_ADp40KO_Ctrl_minus.df$gene, ]

g12 <- ggplot(MG_count_cpm.df, aes(x = Ctrl, y = ADp40KO)) + geom_point(size = 0.5) +
  geom_point(data = MG_count_cpm_plus.df3,  aes(x = Ctrl, y = ADp40KO), color = "red", size = 0.8) +
  geom_point(data = MG_count_cpm_minus.df3,  aes(x = Ctrl, y = ADp40KO), color = "blue", size = 0.8) +
  ggtitle("ADp40KO versus Ctrl in Microglia") +
  theme_classic() +
  xlab("log (average CPM + 1) Ctrl") +
  ylab("log (average CPM + 1) ADp40KO") +
  theme(plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12, color = "black")) +
  geom_abline(slope = 1, color="red", linetype="dashed", size=.5) +
  geom_text_repel(data = MG_count_cpm_plus.df3, aes(x = Ctrl, y = ADp40KO), label = MG_count_cpm_plus.df3$gene, force = 10, nudge_y = .8) +
  geom_text_repel(data = MG_count_cpm_minus.df3, aes(x = Ctrl, y = ADp40KO), label = MG_count_cpm_minus.df3$gene, force = 10, nudge_y = -.6)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_DE_ADp40KO_Ctrl_figures_CPM_scatter.png",
       plot = g12,
       scale = 1, width = 5.5, height = 5, units = "in", device = "png",
       dpi = 300)


##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_DE_figures_session_info.txt")



