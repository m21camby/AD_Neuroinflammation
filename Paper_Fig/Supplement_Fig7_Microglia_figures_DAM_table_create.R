#! /bin/env RScript
# written by SJK at 20 Dec. 2020
# This file is for checking DAM vs DroNc MG DAM markers
# Total 853 genes
# Common 96 genes
# DAM specific: 461- 96 = 365
# DroNc MG specific: 488 - 96 = 392

library(gridExtra)
library(dplyr)
library(ggplot2)
library(tidyr)
library(Seurat)
library(ggrepel)
library(readxl)

source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")




DAM_markers <- read_xlsx("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DAM_markers/mmc2.xlsx")
DAM_markers_UP <- DAM_markers[DAM_markers$`up/down` == 1, ]
DAM_markers_DOWN <- DAM_markers[DAM_markers$`up/down` == -1, ]

DAM_markers_UP$logFC <- log2(DAM_markers_UP$`Microglia3  (average UMI count)`/DAM_markers_UP$`Microglia1 (average UMI count)`)
colnames(DAM_markers_UP) <- c("gene", "-log10(p-value)", "MG1_avg_UMI", "MG2_avg_UMI", "MG3_avg_UMI", "up/down", "logFC_DAM")

DAM_markers_UP <- DAM_markers_UP[!is.infinite(DAM_markers_UP$logFC_DAM), ]


DAM_markers_DOWN$logFC <- log2(DAM_markers_DOWN$`Microglia3  (average UMI count)`/DAM_markers_DOWN$`Microglia1 (average UMI count)`)
colnames(DAM_markers_DOWN) <- c("gene", "-log10(p-value)", "MG1_avg_UMI", "MG2_avg_UMI", "MG3_avg_UMI", "up/down", "logFC_DAM")

DAM_All_markers <- rbind(DAM_markers_UP, DAM_markers_DOWN)


############
resNoFilt_final <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cluster/Fig2_Microglia_figures_DAM_edgeR_DE_test_3_vs_8.csv", row.names = 1)

resNoFilt_final_up <- resNoFilt_final[which(resNoFilt_final$logFC > 0.25 & resNoFilt_final$FDR < 0.01), ]
resNoFilt_final_up <- resNoFilt_final_up[,c("logFC", "gene", "FDR")]

resNoFilt_final_down <- resNoFilt_final[which(resNoFilt_final$logFC < -0.25 & resNoFilt_final$FDR < 0.01), ]
resNoFilt_final_down <- resNoFilt_final_down[,c("logFC", "gene", "FDR")]

resNoFilt_final_sub <- rbind(resNoFilt_final_up, resNoFilt_final_down)


#############################
# Combine all
#############################

DAM_DroNc <- full_join(DAM_All_markers, resNoFilt_final_sub, by = "gene")


# DAM up Total
DAM_DroNc_UP <- full_join(DAM_markers_UP, resNoFilt_final_up, by = "gene")
DAM_DroNc_UP <- DAM_DroNc_UP[,c("gene", "-log10(p-value)", "logFC_DAM", "logFC","FDR" )]
colnames(DAM_DroNc_UP) <- c("gene", "-log10(p-value)", "logFC_DAM", "logFC_DroNc","FDR_DroNc")
write.csv(DAM_DroNc_UP, file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig7_Microglia_figures_DAM_table_create_DAM_DroNc_MG_UP_All.csv")


# DAM up common 
DAM_DroNc_UP_Common <- inner_join(DAM_markers_UP, resNoFilt_final_up, by = "gene")
DAM_DroNc_UP_Common <- DAM_DroNc_UP_Common[,c("gene", "-log10(p-value)", "logFC_DAM", "logFC","FDR" )]
colnames(DAM_DroNc_UP_Common) <- c("gene", "-log10(p-value)", "logFC_DAM", "logFC_DroNc","FDR_DroNc")
write.csv(DAM_DroNc_UP_Common, file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig7_Microglia_figures_DAM_table_create_DAM_DroNc_MG_UP_Common.csv")

# DAM up spec
DAM_DroNc_UP_DAM_spec <- left_join(DAM_markers_UP, resNoFilt_final_up, by = "gene")
DAM_DroNc_UP_DAM_spec <- DAM_DroNc_UP_DAM_spec[,c("gene", "-log10(p-value)", "logFC_DAM", "logFC","FDR" )]
colnames(DAM_DroNc_UP_DAM_spec) <- c("gene", "-log10(p-value)", "logFC_DAM", "logFC_DroNc","FDR_DroNc")
write.csv(DAM_DroNc_UP_DAM_spec, file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig7_Microglia_figures_DAM_table_create_DAM_UP_Specific.csv")

# DroNc up spec
DAM_DroNc_UP_DroNc_spec <- right_join(DAM_markers_UP, resNoFilt_final_up, by = "gene")
DAM_DroNc_UP_DroNc_spec <- DAM_DroNc_UP_DroNc_spec[,c("gene", "-log10(p-value)", "logFC_DAM", "logFC","FDR" )]
colnames(DAM_DroNc_UP_DroNc_spec) <- c("gene", "-log10(p-value)", "logFC_DAM", "logFC_DroNc","FDR_DroNc")
write.csv(DAM_DroNc_UP_DroNc_spec, file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig7_Microglia_figures_DAM_table_create_DroNc_MG_UP_Specific.csv")



