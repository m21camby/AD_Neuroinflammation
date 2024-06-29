#! /bin/env RScript
# written by SJK at 8. Ma7. 2020

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

# split by sample
data9set_cleaned_MG_Ctrl.SO <- subset(data9set_cleaned_MG.SO, subset = sample %in% "Ctrl")
data9set_cleaned_MG_AD.SO <- subset(data9set_cleaned_MG.SO, subset = sample %in% "AD")
data9set_cleaned_MG_ADp40KO.SO <- subset(data9set_cleaned_MG.SO, subset = sample %in% "ADp40KO")


# extract UMI counts and calculate average per each gene
MG_Ctrl_count.df <- data.frame(Matrix::rowMeans(GetAssayData(data9set_cleaned_MG_Ctrl.SO, slot = "counts")))
MG_AD_count.df <- data.frame(Matrix::rowMeans(GetAssayData(data9set_cleaned_MG_AD.SO, slot = "counts")))
MG_ADp40KO_count.df <- data.frame(Matrix::rowMeans(GetAssayData(data9set_cleaned_MG_ADp40KO.SO, slot = "counts")))

MG_Ctrl_count_cpm.df <- apply(MG_Ctrl_count.df, 2, function(x) (x/sum(x)) * 1000000)
MG_AD_count_cpm.df <- apply(MG_AD_count.df, 2, function(x) (x/sum(x)) * 1000000)
MG_ADp40KO_count_cpm.df <- apply(MG_ADp40KO_count.df, 2, function(x) (x/sum(x)) * 1000000)

MG_Ctrl_count_cpm.df <- data.frame(Matrix::rowMeans(MG_Ctrl_count_cpm.df))
MG_AD_count_cpm.df <- data.frame(Matrix::rowMeans(MG_AD_count_cpm.df))
MG_ADp40KO_count_cpm.df <- data.frame(Matrix::rowMeans(MG_ADp40KO_count_cpm.df))


MG_count_cpm.df <- cbind(MG_Ctrl_count_cpm.df, MG_AD_count_cpm.df, MG_ADp40KO_count_cpm.df)
MG_count_cpm.df <- as.data.frame(MG_count_cpm.df)
MG_count_cpm.df$gene <- rownames(MG_count_cpm.df)
# Remove Ttr
MG_count_cpm.df <- MG_count_cpm.df[MG_count_cpm.df$gene != "Ttr", ]
MG_count_cpm.df$gene <- NULL


# remove genes where all zero in samples
MG_count_cpm.df <- MG_count_cpm.df[apply(MG_count_cpm.df == 0, 1, sum) != 3, ]
colnames(MG_count_cpm.df) <- c("Ctrl", "AD", "ADp40KO")
# add pseudocount 0.001 and logarithm 
MG_count_cpm.df <- log(MG_count_cpm.df + 1)
MG_count_cpm.df$gene <- rownames(MG_count_cpm.df)


MG_DE_AD_ADp40KO.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/Cell_type/20200221_Cell_type_DE_by_MAST_MicrogliaAD_ADp40KO.csv", row.names = 1)
MG_DE_AD_ADp40KO.df <- MG_DE_AD_ADp40KO.df[which(abs(MG_DE_AD_ADp40KO.df$avg_logFC) > 0.1 & MG_DE_AD_ADp40KO.df$p_val_adj < 0.01), ]
MG_DE_AD_ADp40KO.df$gene <- rownames(MG_DE_AD_ADp40KO.df)

# log2FC > 0.25 or padj > 2
MG_DE_AD_ADp40KO_plus.df <- MG_DE_AD_ADp40KO.df[which(MG_DE_AD_ADp40KO.df$avg_logFC > 0.25 & -log10(MG_DE_AD_ADp40KO.df$p_val_adj) > 2),]
# padj > 2
MG_DE_AD_ADp40KO_minus.df <- MG_DE_AD_ADp40KO.df[which(MG_DE_AD_ADp40KO.df$avg_logFC < -0.25 & -log10(MG_DE_AD_ADp40KO.df$p_val_adj) > 2),]
# Il12b gene
MG_DE_AD_ADp40KO_Il12b.df <-  MG_DE_AD_ADp40KO.df[MG_DE_AD_ADp40KO.df$gene == "Il12b", ]



# extract genes from volcano plots
MG_count_cpm_plus.df2 <- MG_count_cpm.df[MG_count_cpm.df$gene %in% MG_DE_AD_ADp40KO_plus.df$gene, ]
# Remove Gm26520 for labeling
MG_count_cpm_plus.df2 <- MG_count_cpm_plus.df2[MG_count_cpm_plus.df2$gene != "Gm26520",]

MG_count_cpm_minus.df2 <- MG_count_cpm.df[MG_count_cpm.df$gene %in% MG_DE_AD_ADp40KO_minus.df$gene, ]
MG_count_cpm_Il12b.df2 <- MG_count_cpm.df[MG_count_cpm.df$gene %in% MG_DE_AD_ADp40KO_Il12b.df$gene, ]


# Gm26520
MG_count_cpm_plus.df3 <- MG_count_cpm.df[MG_count_cpm.df$gene %in% "Gm26520", ]


# scatter plot for AD vs ADp40KO
g1 <- ggplot(MG_count_cpm.df, aes(x = ADp40KO, y = AD)) + geom_point(size = 1, color = "#666666", alpha = 0.2) +
  geom_point(data = MG_count_cpm_plus.df2,  aes(x = ADp40KO, y = AD), color = "red", size = 2) +
  geom_point(data = MG_count_cpm_plus.df3,  aes(x = ADp40KO, y = AD), color = "red", size = 2) + 
  geom_point(data = MG_count_cpm_minus.df2,  aes(x = ADp40KO, y = AD), color = "blue", size = 2) +
  ggtitle("AD vs ADp40KO in Microglia") +
  theme_classic() +
  xlab("log (average CPM + 1) ADp40KO") +
  ylab("log (average CPM + 1) AD") +
  theme(plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12, color = "black")) +
  geom_abline(slope = 1, color="red", linetype="dashed", size=.5) +
  geom_text_repel(data = MG_count_cpm_plus.df2, aes(x = ADp40KO, y = AD), label = MG_count_cpm_plus.df2$gene, force = 15, nudge_y = 0.7, nudge_x = -1, max.iter = 1000) +
  geom_text_repel(data = MG_count_cpm_minus.df2, aes(x = ADp40KO, y = AD), label = MG_count_cpm_minus.df2$gene, force = 10, nudge_y = -.5) +
  geom_text_repel(data = MG_count_cpm_Il12b.df2, aes(x = ADp40KO, y = AD), label = MG_count_cpm_Il12b.df2$gene, nudge_y = 1, nudge_x = -0.1) +
  geom_text_repel(data = MG_count_cpm_plus.df3, aes(x = ADp40KO, y = AD), label = MG_count_cpm_plus.df3$gene, force = 15, nudge_y = -0.5, nudge_x = -1)


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_DE_scatter_AD_ADp40KO.png",
       plot = g1,
       scale = 1, width = 5.5, height = 5, units = "in", device = "png",
       dpi = 300)

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_DE_scatter_AD_ADp40KO_session_info.txt")







