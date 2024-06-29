#! /bin/env RScript
# written by SJK at 25. Mar. 2020

.libPaths(c("/home/skim/R/usr_lib", .libPaths()))

library(Seurat)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyr)
library(ggrepel)

source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

####################################################################
# This analysis is for comparison 
# with DAM cluster from Ido Amit's paper (figureS1 (D) & (E))
# Although some markers show similar expression, in here we double checked Genes and GO and compare with his data 
####################################################################


################################
# scatter plot
################################
# Scatterplot showing the average molecules count (log2 scale) of sample compared with other sample (y axis).
load(file="/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200210_figures_for_meeting_report_cleaned_Seurat_object_res0.8.Robj")

# split by sample
data9set_cleaned_MG_c3.SO <- subset(data9set_cleaned.SO, subset = seurat_clusters %in% 3)
data9set_cleaned_MG_c8.SO <- subset(data9set_cleaned.SO, subset = seurat_clusters %in% 8)

# extract UMI counts and calculate average per each gene
MG_c3_count.df <- data.frame(Matrix::rowMeans(GetAssayData(data9set_cleaned_MG_c3.SO, slot = "counts")))
MG_c8_count.df <- data.frame(Matrix::rowMeans(GetAssayData(data9set_cleaned_MG_c8.SO, slot = "counts")))

MG_c3_count_cpm.df <- apply(MG_c3_count.df, 2, function(x) (x/sum(x)) * 1000000)
MG_c8_count_cpm.df <- apply(MG_c8_count.df, 2, function(x) (x/sum(x)) * 1000000)

MG_c3_count_cpm.df <- data.frame(Matrix::rowMeans(MG_c3_count_cpm.df))
MG_c8_count_cpm.df <- data.frame(Matrix::rowMeans(MG_c8_count_cpm.df))


MG_count_cpm.df <- cbind(MG_c3_count_cpm.df, MG_c8_count_cpm.df)
MG_count_cpm.df <- as.data.frame(MG_count_cpm.df)

# remove genes where all zero in samples
MG_count_cpm.df <- MG_count_cpm.df[apply(MG_count_cpm.df == 0, 1, sum) != 1, ]
colnames(MG_count_cpm.df) <- c("cluster3", "cluster8")
# add pseudocount 0.001 and logarithm 
MG_count_cpm.df <- log(MG_count_cpm.df + 1)
MG_count_cpm.df$gene <- rownames(MG_count_cpm.df)

# DAM plus gene and minus genes list
DAM_plus <- c("Cst7", "Lpl", "Apoe", "Clec7a", "Ank",  "Spp1", "Igf1", "Csf1", "Gpnmb", "Gm11428")
DAM_plus2 <- "Axl"
DAM_plus3 <- "Fam20c"
DAM_plus4 <- "Itgax"
DAM_plus5 <- "Cybb"

DAM_minus <- c("Tmem119",  "Serinc3")
DAM_minus2 <- "Cx3cr1"
DAM_minus3 <- "P2ry12"

# extract genes from volcano plots
MG_count_cpm_plus.df <- MG_count_cpm.df[MG_count_cpm.df$gene %in% DAM_plus, ]
MG_count_cpm_minus.df <- MG_count_cpm.df[MG_count_cpm.df$gene %in% DAM_minus, ]

MG_count_cpm_plus2.df <- MG_count_cpm.df[MG_count_cpm.df$gene %in% DAM_plus2, ]
MG_count_cpm_plus3.df <- MG_count_cpm.df[MG_count_cpm.df$gene %in% DAM_plus3, ]
MG_count_cpm_minus2.df <- MG_count_cpm.df[MG_count_cpm.df$gene %in% DAM_minus2, ]
MG_count_cpm_plus4.df <- MG_count_cpm.df[MG_count_cpm.df$gene %in% DAM_plus4, ]
MG_count_cpm_minus3.df <- MG_count_cpm.df[MG_count_cpm.df$gene %in% DAM_minus3, ]
MG_count_cpm_plus5.df <- MG_count_cpm.df[MG_count_cpm.df$gene %in% DAM_plus5, ]

g1 <- ggplot(MG_count_cpm.df, aes(x = cluster3, y = cluster8)) + geom_point(size = 1, color = "#666666", alpha = 0.2) +
  geom_point(data = MG_count_cpm_plus.df,  aes(x = cluster3, y = cluster8), color = "red", size = 2) +
  geom_point(data = MG_count_cpm_plus2.df,  aes(x = cluster3, y = cluster8), color = "red", size = 2) +
  geom_point(data = MG_count_cpm_plus3.df,  aes(x = cluster3, y = cluster8), color = "red", size = 2) +
  geom_point(data = MG_count_cpm_plus4.df,  aes(x = cluster3, y = cluster8), color = "red", size = 2) +
  geom_point(data = MG_count_cpm_plus5.df,  aes(x = cluster3, y = cluster8), color = "red", size = 2) +
  geom_point(data = MG_count_cpm_minus.df,  aes(x = cluster3, y = cluster8), color = "blue", size = 2) +
  geom_point(data = MG_count_cpm_minus2.df,  aes(x = cluster3, y = cluster8), color = "blue", size = 2) +
  geom_point(data = MG_count_cpm_minus3.df,  aes(x = cluster3, y = cluster8), color = "blue", size = 2) +
  ggtitle("cluster3 vs cluster8 in Microglia") + 
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) + 
  theme_classic() +
  xlab("log (average CPM + 1) cluster3") +
  ylab("log (average CPM + 1) cluster8") +
  theme(plot.title = element_text(size = 16, hjust = 0.5, family = "helvetica"),
        axis.title = element_text(size = 16, family = "helvetica"),
        axis.text = element_text(size = 15, color = "black", family = "helvetica")) +
  geom_abline(slope = 1, color="red", linetype="dashed", size=.5) +
  geom_text_repel(data = MG_count_cpm_plus.df, aes(x = cluster3, y = cluster8), label = MG_count_cpm_plus.df$gene, force = 20, nudge_y = .5) +
  geom_text_repel(data = MG_count_cpm_minus.df, aes(x = cluster3, y = cluster8), label = MG_count_cpm_minus.df$gene, force = 5, nudge_y = -.6) +
  geom_text_repel(data = MG_count_cpm_plus2.df, aes(x = cluster3, y = cluster8), label = MG_count_cpm_plus2.df$gene, force = 5, nudge_y = -.2) +
  geom_text_repel(data = MG_count_cpm_plus3.df, aes(x = cluster3, y = cluster8), label = MG_count_cpm_plus3.df$gene, force = 5, nudge_x = .2, nudge_y = -.5) + 
  geom_text_repel(data = MG_count_cpm_minus2.df, aes(x = cluster3, y = cluster8), label = MG_count_cpm_minus2.df$gene, force = 5, nudge_y = 0.1, nudge_x = .2) +
  geom_text_repel(data = MG_count_cpm_plus4.df, aes(x = cluster3, y = cluster8), label = MG_count_cpm_plus4.df$gene, force = 5, nudge_x = 0.7) + 
  geom_text_repel(data = MG_count_cpm_minus3.df, aes(x = cluster3, y = cluster8), label = MG_count_cpm_minus3.df$gene, force = 5, nudge_x = 0.7, nudge_y = -.2) +
  geom_text_repel(data = MG_count_cpm_plus5.df, aes(x = cluster3, y = cluster8), label = MG_count_cpm_plus5.df$gene, force = 5, nudge_y = 1)
  

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_DAM_comparison_scatter_Modified.png",
       plot = g1,
       scale = 1, width = 5.5, height = 5, units = "in", device = "png",
       dpi = 300)

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_DAM_comparison_Modified_session_info.txt")


