library(gridExtra)
library(dplyr)
library(ggplot2)
library(tidyr)
library(Seurat)
library(ggrepel)
library(readxl)
library(tidyr)
library(scran)
library(edgeR)
library(topGO)
library(viridis)
library(grid)
library(biomaRt)
library(scales)

#source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

load(file = "/nfs/team292/sk27/tmp/20200210_figures_for_meeting_report_cleaned_Seurat_object_res0.8.Robj")


# --------------------------------------- #
# Fig 7.A
# --------------------------------------- #


# --------------------------------------- #
# Fig 7.B
# --------------------------------------- #

# --------------------------------------- #
# Data preparation 
data9set_cleaned_MG.SO <- subset(data9set_cleaned.SO, subset = seurat_clusters %in% c(3, 8))


# Cluster 8 AD AD vs Cluster 3 WT
data9set_cleaned_3.SO <- subset(data9set_cleaned.SO, subset = seurat_clusters %in% c(3))
data9set_cleaned_3_WT.SO <- subset(data9set_cleaned_3.SO, subset = sample %in% c("Ctrl"))

data9set_cleaned_8.SO <- subset(data9set_cleaned.SO, subset = seurat_clusters %in% c(8))
data9set_cleaned_8_AD.SO <- subset(data9set_cleaned_8.SO, subset = sample %in% c("AD"))

data9set_cleaned_3_WT_8_AD.SO <- merge(data9set_cleaned_3_WT.SO, data9set_cleaned_8_AD.SO)

# --------------------------------------- #

data9set_cleaned_3.SO <- subset(data9set_cleaned.SO, subset = seurat_clusters %in% c(3))
cluster3.df <- rowMeans(GetAssayData(data9set_cleaned_3.SO, slot = "data")) %>% as.data.frame
colnames(cluster3.df) <- "cluster3"

data9set_cleaned_8.SO <- subset(data9set_cleaned.SO, subset = seurat_clusters %in% c(8))
cluster8.df <- rowMeans(GetAssayData(data9set_cleaned_8.SO, slot = "data")) %>% as.data.frame
colnames(cluster8.df) <- "cluster8"

clusters.df <- cbind(cluster3.df, cluster8.df)
clusters.df$gene <- rownames(clusters.df)

# ---------------------------- # 
# DAM from Ido Amit
# ---------------------------- # 

DAM_markers <- read_xlsx("~/MDC/Final_revisions/mmc2_DAM_Ido_Amit.xlsx")
DAM_markers_UP <- DAM_markers[DAM_markers$`up/down` == 1, ]
DAM_markers_DOWN <- DAM_markers[DAM_markers$`up/down` == -1, ]

DAM_markers_UP$logFC <- log2(DAM_markers_UP$`Microglia3  (average UMI count)`/DAM_markers_UP$`Microglia1 (average UMI count)`)
colnames(DAM_markers_UP) <- c("gene", "-log10(p-value)", "MG1_avg_UMI", "MG2_avg_UMI", "MG3_avg_UMI", "up/down", "logFC_DAM")
DAM_markers_UP <- DAM_markers_UP[!is.infinite(DAM_markers_UP$logFC_DAM), ]


DAM_markers_DOWN$logFC <- log2(DAM_markers_DOWN$`Microglia3  (average UMI count)`/DAM_markers_DOWN$`Microglia1 (average UMI count)`)
colnames(DAM_markers_DOWN) <- c("gene", "-log10(p-value)", "MG1_avg_UMI", "MG2_avg_UMI", "MG3_avg_UMI", "up/down", "logFC_DAM")

DAM_All_markers <- rbind(DAM_markers_UP, DAM_markers_DOWN)

# DAM from Ido data
DAM_UP.df <- right_join(resNoFilt_final_sub, DAM_markers_UP, by = "gene")

DAM_DOWN.df <- right_join(resNoFilt_final_sub, DAM_markers_DOWN, by = "gene")


DAM_All.df <- rbind(DAM_UP.df, DAM_DOWN.df)
DAM_All.df <- left_join(DAM_All.df, clusters.df, by = "gene")
DAM_All.df$norm_exp <- ifelse(DAM_All.df$`up/down` == "1", DAM_All.df$cluster8, DAM_All.df$cluster3)
DAM_All.df_sig <- DAM_All.df[DAM_All.df$FDR < 0.05, ]

DAM_All.df_1 <- DAM_All.df[which(DAM_All.df$logFC_DAM > 6.5 & DAM_All.df$logFC < 0.2), ]
DAM_All.df_1 <- DAM_All.df_1[DAM_All.df_1$gene != "Atp6v0d2", ]
DAM_All.df_2 <- DAM_All.df[which(DAM_All.df$logFC_DAM > 6 & DAM_All.df$logFC > 1.2), ]
DAM_All.df_3 <- DAM_All.df[which(DAM_All.df$logFC_DAM > 2 & DAM_All.df$logFC > 2), ]
DAM_All.df_4 <- DAM_All.df[which(DAM_All.df$logFC_DAM > 6.5 & DAM_All.df$logFC > 0.2 & DAM_All.df$logFC < 1), ]
DAM_All.df_5 <- DAM_All.df[which(DAM_All.df$logFC_DAM < -1.5 & DAM_All.df$logFC < -0.7 ), ]
DAM_All.df_5 <- DAM_All.df_5[DAM_All.df_5$gene != "P2ry12", ]
DAM_All.df_6 <- DAM_All.df[DAM_All.df$gene == "P2ry12", ]


# ---------------------------- # 
# loading edgeR results
# ---------------------------- # 

resNoFilt_final <- read.csv("~/MDC/Final_revisions/Fig2_Microglia_figures_DAM_edgeR_DE_test_3_vs_8.csv", row.names = 1)

resNoFilt_final_up_0.5 <- resNoFilt_final[which(resNoFilt_final$logFC > 0.5 & resNoFilt_final$FDR < 0.01), ]
resNoFilt_final_up_0.5_sub <- resNoFilt_final_up_0.5[,c("logFC", "gene", "FDR")]

resNoFilt_final_down_0.5 <- resNoFilt_final[which(resNoFilt_final$logFC < -0.5 & resNoFilt_final$FDR < 0.01), ]
resNoFilt_final_down_0.5_sub <- resNoFilt_final_down_0.5[,c("logFC", "gene", "FDR")]

resNoFilt_final_0.5_sub <- rbind(resNoFilt_final_up_0.5_sub, resNoFilt_final_down_0.5_sub)

resNoFilt_final_DAM <- resNoFilt_final[resNoFilt_final$gene %in% DAM_markers3$`Gene name`, ]
resNoFilt_final_DAM_sub <- resNoFilt_final_DAM[,c("logFC", "gene", "FDR")]

resNoFilt_final_sub <- resNoFilt_final[,c("logFC", "gene", "FDR")]


# --------------------------------------- #
# extract UMI counts and calculate average per each gene
MG_c3_count.df <- data.frame(Matrix::rowMeans(GetAssayData(data9set_cleaned_3_WT.SO, slot = "counts")))
MG_c8_count.df <- data.frame(Matrix::rowMeans(GetAssayData(data9set_cleaned_8_AD.SO, slot = "counts")))

MG_c3_count_cpm.df <- apply(MG_c3_count.df, 2, function(x) (x/sum(x)) * 1000000)
MG_c8_count_cpm.df <- apply(MG_c8_count.df, 2, function(x) (x/sum(x)) * 1000000)

MG_c3_count_cpm.df <- data.frame(Matrix::rowMeans(MG_c3_count_cpm.df))
MG_c8_count_cpm.df <- data.frame(Matrix::rowMeans(MG_c8_count_cpm.df))

MG_count_cpm.df <- cbind(MG_c3_count_cpm.df, MG_c8_count_cpm.df)
MG_count_cpm.df <- as.data.frame(MG_count_cpm.df)

# remove genes where all zero in samples
MG_count_cpm.df <- MG_count_cpm.df[rowSums(MG_count_cpm.df) > 0, ]
colnames(MG_count_cpm.df) <- c("cluster3", "cluster8")

# add pseudocount 0.001 and logarithm 
MG_count_cpm.df <- log2(MG_count_cpm.df + 1)
MG_count_cpm.df$gene <- rownames(MG_count_cpm.df)

MG_count_cpm_plus.df <- MG_count_cpm.df[MG_count_cpm.df$gene %in% DAM_markers_UP$gene, ]
MG_count_cpm_minus.df <- MG_count_cpm.df[MG_count_cpm.df$gene %in% DAM_markers_DOWN$gene, ]

MG_count_plus_label1.df <- MG_count_cpm.df[MG_count_cpm.df$gene %in% c("Myo1e", "Mamdc2","Axl", "Cst7", "Lpl", "Apoe", "Clec7a", "Ank",  "Spp1", "Igf1", "Csf1", "Gm11428"), ]
MG_count_plus_label2.df <- MG_count_cpm.df[MG_count_cpm.df$gene %in% c("Il12b"), ]
MG_count_plus_label3.df <- MG_count_cpm.df[MG_count_cpm.df$gene %in% c("Fam20c", "Itgax", "Cybb"), ]

MG_count_mius_label1.df <- MG_count_cpm.df[MG_count_cpm.df$gene %in% c("Malat1","Tmem119", "P2ry12", "Serinc3", "Cx3cr1"), ]

g1 <- ggplot(MG_count_cpm.df, aes(x = cluster3, y = cluster8)) + geom_point(size = 1, color = "#666666", alpha = 0.2) + 
  geom_point(data = MG_count_cpm_plus.df,  aes(x = cluster3, y = cluster8), color = "red", size = 1.4, alpha = 0.5) +
  geom_point(data = MG_count_cpm_minus.df,  aes(x = cluster3, y = cluster8), color = "blue", size = 1.4, alpha = 0.5) +
  ggtitle("cluster3 WT vs cluster8 APPPS1 in Microglia") + 
  theme_classic() + xlab("log2 (average CPM + 1) cluster3 WT") +
  ylab("log2 (average CPM + 1) cluster8 APPPS1") +
  theme(plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
        axis.title.x = element_text(size = 15, family = "helvetica", color = "#E69F00"),
        axis.title.y = element_text(size = 15, family = "helvetica", color = "#56B4E9"),
        axis.text = element_text(size = 15, color = "black", family = "helvetica")) 

#  geom_abline(slope = 1, color="#FFFFCC", linetype="dashed", size=.8) + 
#  geom_text_repel(data = MG_count_plus_label1.df, aes(x = cluster3, y = cluster8), label = MG_count_plus_label1.df$gene, force = 20, nudge_y = 1) +
#  geom_text_repel(data = MG_count_plus_label2.df, aes(x = cluster3, y = cluster8), label = MG_count_plus_label2.df$gene, nudge_x = -0.5) +
#  geom_text_repel(data = MG_count_plus_label3.df, aes(x = cluster3, y = cluster8), label = MG_count_plus_label3.df$gene, force = 10, nudge_y = 0.5, nudge_x = 0.5) + 
#  geom_text_repel(data = MG_count_mius_label1.df, aes(x = cluster3, y = cluster8), label = MG_count_mius_label1.df$gene, force = 10, nudge_y = -0.7, nudge_x = 0.12) 

ggsave(filename = "~/MDC/Final_revisions/Fig7_B_no_label.png",
       plot = g1,
       scale = 1, width = 7.5, height = 6.5, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "~/MDC/Final_revisions/Fig7_B_no_label.pdf",
       plot = g1,
       scale = 1, width = 7.5, height = 6.5, units = "in", device = cairo_pdf,
       dpi = 300)


# --------------------------------------- #
# Fig 7.C
# --------------------------------------- #

# --------------------------------------- #
# Fig 7. D & E
# --------------------------------------- #

ensembl=useMart("ensembl")
#listDatasets(ensembl)
#mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://may2024.archive.ensembl.org/")
#human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://may2024.archive.ensembl.org/") 
#mouse = useEnsembl("ensembl","mmusculus_gene_ensembl", mirror = "useast")
#human = useEnsembl("ensembl","hsapiens_gene_ensembl", mirror = "useast")
#human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
#mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
human <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl', host = "www.ensembl.org")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "www.ensembl.org")
human2 <- useEnsembl("ensembl","hsapiens_gene_ensembl", mirror = "useast")
mouse2 <- useEnsembl("ensembl","mmusculus_gene_ensembl", mirror = "useast")




# load Microglia AD vs ADp40KO
AD_ADp40KO <- read.csv("~/MDC/Microglia_AD_ADp40KO.csv", row.names = 1)

# Depp genes
Depp <- read.csv("~/MDC/41586_2023_6120_MOESM3_ESM_tab.csv", row.names = 1)
Depp$gene <- rownames(Depp)
Depp_genes <- Depp[(Depp$avg_log2FC > 0.25 & Depp$p_val_adj < 0.01), ]$gene
Depp_genes2 <- Depp[(Depp$avg_log2FC < -0.25 & Depp$p_val_adj < 0.01), ]$gene

Depp_genes <- c(Depp_genes, Depp_genes2)

# Phagocytosis genes
Human_genes <- c("CLU", "APOJ", "TREM2", "SORL1", "CR1", "CD33", "PLCG2", "MS4A4A", "MS4A6A", "PILRA", "FERMT2",
                 "PTK2B", "CASS4", "ABI3", "BIN1", 'CD2AP', 'RIN3', 'RABEP1', 'APBB3', 
                 'ZKSCAN1', 'TP53INP1', 'ZYX', 'AP4E1', 'AP4M1', 'SPPL2A', 'GRN', 'PICALM', 'ABCA7', 'APOE')

AD_mouse_genes = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = Human_genes , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)

AD_mouse_genes <- AD_mouse_genes[!AD_mouse_genes$MGI.symbol %in% c("Ms4a4c", "Ms4a4b", "Ms4a4d", "Ms4a6c", "Ms4a6b", "Ms4a6d", "Siglecf"), ]
write.csv(AD_mouse_genes, "/nfs/team292/sk27/tmp/AD_Neuroinflammation/human_genes_to_mouse_conversion_Ekaterina.csv")

AD_mouse_genes <- read.csv( "/nfs/team292/sk27/tmp/AD_Neuroinflammation/human_genes_to_mouse_conversion_Ekaterina.csv")



# Myelin genes from Depp
# Plus
AD_ADp40KO_plus.df <- AD_ADp40KO[which(AD_ADp40KO$logFC > 0 & AD_ADp40KO$FDR < 0.05  & AD_ADp40KO$gene %in% Depp_genes),]
# Minus
AD_ADp40KO_minus.df <- AD_ADp40KO[which(AD_ADp40KO$logFC < 0 & AD_ADp40KO$FDR < 0.05 & AD_ADp40KO$gene %in% Depp_genes),]
# Il12b gene
AD_ADp40KO_Il12b.df <-  AD_ADp40KO[AD_ADp40KO$gene == "Il12b", ]

# plus
AD_ADp40KO_plus2.df <- AD_ADp40KO[which(AD_ADp40KO$logFC > 0 & AD_ADp40KO$FDR > 0.05  & AD_ADp40KO$gene %in% Depp_genes),]
# Minus
AD_ADp40KO_minus2.df <- AD_ADp40KO[which(AD_ADp40KO$logFC < 0 & AD_ADp40KO$FDR > 0.05 & AD_ADp40KO$gene %in% Depp_genes),]



# Myelin and amyloid genes from Depp
g1 <- ggplot(AD_ADp40KO) + geom_point(aes(x = logFC, y = -log10(FDR)), color = "gray", alpha = 0.5, shape = 1, size = 1) +
  geom_point(data = AD_ADp40KO_plus.df, aes(x = logFC, y = -log10(FDR)), color = "#d00000", shape = 19, size = 2.5) +
  geom_point(data = AD_ADp40KO_minus.df, aes(x = logFC, y = -log10(FDR)), color = "#023e8a", shape = 19, size = 2.5) +
  geom_point(data = AD_ADp40KO_plus2.df, aes(x = logFC, y = -log10(FDR)), color = "#d00000", alpha = 1, shape = 1, size = 2.5) +
  geom_point(data = AD_ADp40KO_minus2.df, aes(x = logFC, y = -log10(FDR)), color = "#023e8a", alpha = 1, shape = 1, size = 2.5) + 
  geom_point(data = AD_ADp40KO_Il12b.df, aes(x = logFC, y = -log10(FDR)), color = "#f72585", alpha = 1, shape = 19, size = 2.5) +
  
  ggtitle("AD vs ADp40KO in Microglia (Depp genes)") + 
  theme_classic() + 
  xlab("log2 fold change") + 
  ylab("-log10(adjusted p-value)") +
  scale_y_continuous(expand = c(0,0), limits = c(0,15), oob = squish) + xlim(-1, 1) + 
  theme(legend.position = "none",
        plot.title = element_text(size = 18, hjust = 0.5, family = "helvetica"),
        axis.title = element_text(size = 20, family = "helvetica"), 
        axis.text = element_text(size = 20, color = "black", family = "helvetica")) + 
  geom_text_repel(data = AD_ADp40KO_plus.df, aes(x = logFC, y = -log10(FDR)), label = AD_ADp40KO_plus.df$gene, force = 10, family= "helvetica", size = 5) + 
  #geom_text_repel(data = AD_ADp40KO_plus2.df, aes(x = logFC, y = -log10(FDR)), label = AD_ADp40KO_plus2.df$gene, force = 10) + 
  geom_text_repel(data = AD_ADp40KO_minus.df, aes(x = logFC, y = -log10(FDR)), label = AD_ADp40KO_minus.df$gene, nudge_y = 1, family= "helvetica", size = 5) +
  #geom_text_repel(data = AD_ADp40KO_minus2.df, aes(x = logFC, y = -log10(FDR)), label = AD_ADp40KO_minus2.df$gene, nudge_y = 1) + 
  geom_text_repel(data = AD_ADp40KO_Il12b.df, aes(x = logFC, y = -log10(FDR)), label = AD_ADp40KO_Il12b.df$gene, nudge_y = 1, family= "helvetica", size = 5)  

g1

g1 <- ggplot(AD_ADp40KO) + geom_point(aes(x = logFC, y = -log10(FDR)), color = "gray", alpha = 0.5, shape = 1, size = 1) +
  geom_point(data = AD_ADp40KO_plus.df, aes(x = logFC, y = -log10(FDR)), color = "#d00000", shape = 19, size = 2.5) +
  geom_point(data = AD_ADp40KO_minus.df, aes(x = logFC, y = -log10(FDR)), color = "#023e8a", shape = 19, size = 2.5) +
  geom_point(data = AD_ADp40KO_plus2.df, aes(x = logFC, y = -log10(FDR)), color = "#d00000", alpha = 1, shape = 1, size = 2.5) +
  geom_point(data = AD_ADp40KO_minus2.df, aes(x = logFC, y = -log10(FDR)), color = "#023e8a", alpha = 1, shape = 1, size = 2.5) + 
  geom_point(data = AD_ADp40KO_Il12b.df, aes(x = logFC, y = -log10(FDR)), color = "#f72585", alpha = 1, shape = 19, size = 2.5) +
  
  ggtitle("AD vs ADp40KO in Microglia (Depp genes)") + 
  theme_classic() + 
  xlab("log2 fold change") + 
  ylab("-log10(adjusted p-value)") +
  scale_y_continuous(expand = c(0,0), limits = c(0,15), oob = squish) + xlim(-1, 1) + 
  theme(legend.position = "none",
        plot.title = element_text(size = 18, hjust = 0.5, family = "helvetica"),
        axis.title = element_text(size = 20, family = "helvetica"), 
        axis.text = element_text(size = 20, color = "black", family = "helvetica")) 
g1

ggsave(filename = "~/MDC/Final_revisions/Fig7_E_no_label.png",
       plot = g1,
       scale = 1, width = 6.5, height = 6, units = "in", device = "png",
       dpi = 300)


ggsave(filename = "~/MDC/Final_revisions/Fig7_E_no_label.pdf",
       plot = g1,
       scale = 1, width = 6.5, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)


# -------------------- #
# Phagocytosis genes
# -------------------- #
Phagocytosis <- AD_mouse_genes$MGI.symbol

#AD_ADp40KO_plus.df <- AD_ADp40KO[which(AD_ADp40KO$gene %in% Phagocytosis),]

AD_ADp40KO_plus.df <- AD_ADp40KO[which(AD_ADp40KO$logFC > 0 & AD_ADp40KO$FDR < 0.05  & AD_ADp40KO$gene %in% Phagocytosis),]
# Minus
AD_ADp40KO_minus.df <- AD_ADp40KO[which(AD_ADp40KO$logFC < 0 & AD_ADp40KO$FDR < 0.05 & AD_ADp40KO$gene %in% Phagocytosis),]
# Il12b gene
AD_ADp40KO_Il12b.df <-  AD_ADp40KO[AD_ADp40KO$gene == "Il12b", ]
# plus
AD_ADp40KO_plus2.df <- AD_ADp40KO[which(AD_ADp40KO$logFC > 0 & AD_ADp40KO$FDR > 0.05  & AD_ADp40KO$gene %in% Phagocytosis),]
# Minus
AD_ADp40KO_minus2.df <- AD_ADp40KO[which(AD_ADp40KO$logFC < 0 & AD_ADp40KO$FDR > 0.05 & AD_ADp40KO$gene %in% Phagocytosis),]


g1 <- ggplot(AD_ADp40KO) + geom_point(aes(x = logFC, y = -log10(FDR)), color = "gray", alpha = 0.5, shape = 1, size = 1) +
  geom_point(data = AD_ADp40KO_plus.df, aes(x = logFC, y = -log10(FDR)), color = "#d00000", shape = 19, size = 2.5) +
  geom_point(data = AD_ADp40KO_minus.df, aes(x = logFC, y = -log10(FDR)), color = "#023e8a", shape = 19, size = 2.5) +
  geom_point(data = AD_ADp40KO_plus2.df, aes(x = logFC, y = -log10(FDR)), color = "#d00000", alpha = 1, shape = 1, size = 2.5) +
  geom_point(data = AD_ADp40KO_minus2.df, aes(x = logFC, y = -log10(FDR)), color = "#023e8a", alpha = 1, shape = 1, size = 2.5) + 
  geom_point(data = AD_ADp40KO_Il12b.df, aes(x = logFC, y = -log10(FDR)), color = "#f72585", alpha = 1, shape = 19, size = 2.5) +
  
  ggtitle("AD vs ADp40KO in Microglia (Phagocytosis genes)") + 
  theme_classic() + 
  xlab("log2 fold change") + 
  ylab("-log10(adjusted p-value)") +
  scale_y_continuous(expand = c(0,0), limits = c(0,15)) + xlim(-1, 1) + 
  theme(legend.position = "none",
        plot.title = element_text(size = 18, hjust = 0.5, family = "helvetica"),
        axis.title = element_text(size = 20, family = "helvetica"), 
        axis.text = element_text(size = 20, color = "black", family = "helvetica"))

g1


ggsave(filename = "~/MDC/Final_revisions/Fig7_D_no_label.png",
       plot = g1,
       scale = 1, width = 6.5, height = 6, units = "in", device = "png",
       dpi = 300)


ggsave(filename = "~/MDC/Final_revisions/Fig7_D_no_label.pdf",
       plot = g1,
       scale = 1, width = 6.5, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)