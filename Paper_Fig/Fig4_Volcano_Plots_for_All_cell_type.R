#! /bin/env RScript
# written by SJK at 22. Oct. 2020

library(gridExtra)
library(dplyr)
library(ggplot2)
library(tidyr)
library(Seurat)
library(ggrepel)

source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

#load(file="/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200210_figures_for_meeting_report_cleaned_Seurat_object_res0.8.Robj")

volcano_plot <- function(DE.df, DE.df_plus, DE.df_minus, DE.df_il12, title = "title", y_lim = c(0,20), x_lim = c(-1, 1), plus_threshold = 1E-5, minus_threshold = 1E-5, plus_force = 10, minus_force = 10, il12b = 3){
  ggplot(DE.df) + geom_point(aes(x = logFC, y = -log10(FDR)), alpha = 0.3) +
    geom_point(data = DE.df_plus, aes(x = logFC, y = -log10(FDR)), color = "red", alpha = 1) +
    geom_point(data = DE.df_minus, aes(x = logFC, y = -log10(FDR)), color = "blue", alpha = 1) +
    ggtitle(title) + 
    theme_classic() + 
    xlab("log fold change") + 
    ylab("-log10(adjusted p-value)") +
    scale_y_continuous(expand = c(0,0), limits = y_lim) + xlim(x_lim) + 
    theme(legend.position = "none",
          plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
          axis.title = element_text(size = 15), 
          axis.text = element_text(size = 12, color = "black")) +
    geom_text_repel(data = DE.df_plus[DE.df_plus$FDR < plus_threshold, ], aes(x = logFC, y = -log10(FDR)), label = DE.df_plus[DE.df_plus$FDR < plus_threshold, ]$gene, force = plus_force) + 
    geom_text_repel(data = DE.df_minus[DE.df_minus$FDR < minus_threshold, ], aes(x = logFC, y = -log10(FDR)), label = DE.df_minus[DE.df_minus$FDR < minus_threshold, ]$gene, force = minus_force) + 
    geom_text_repel(data = DE.df_il12, aes(x = logFC, y = -log10(FDR)), label = DE.df_il12$gene, nudge_y = il12b)
  
}

volcano_plot_no_labels <- function(DE.df, DE.df_plus, DE.df_minus, DE.df_il12, title = "title", y_lim = c(0,20), x_lim = c(-1, 1), plus_threshold = 1E-5, minus_threshold = 1E-5, plus_force = 10, minus_force = 10, il12b = 3){
  ggplot(DE.df) + geom_point(aes(x = logFC, y = -log10(FDR)), alpha = 0.3) +
    geom_point(data = DE.df_plus, aes(x = logFC, y = -log10(FDR)), color = "red", alpha = 1) +
    geom_point(data = DE.df_minus, aes(x = logFC, y = -log10(FDR)), color = "blue", alpha = 1) +
    ggtitle(title) + 
    theme_classic() + 
    xlab("log fold change") + 
    ylab("-log10(adjusted p-value)") +
    scale_y_continuous(expand = c(0,0), limits = y_lim) + xlim(x_lim) + 
    theme(legend.position = "none",
          plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
          axis.title = element_text(size = 15), 
          axis.text = element_text(size = 12, color = "black"))
  
}

volcano_plot_no_il12b <- function(DE.df, DE.df_plus, DE.df_minus, DE.df_il12, title = "title", y_lim = c(0,20), x_lim = c(-1, 1), plus_threshold = 1E-5, minus_threshold = 1E-5, plus_force = 10, minus_force = 10, il12b = 3){
  ggplot(DE.df) + geom_point(aes(x = logFC, y = -log10(FDR)), alpha = 0.3) +
    geom_point(data = DE.df_plus, aes(x = logFC, y = -log10(FDR)), color = "red", alpha = 1) +
    geom_point(data = DE.df_minus, aes(x = logFC, y = -log10(FDR)), color = "blue", alpha = 1) +
    ggtitle(title) + 
    theme_classic() + 
    xlab("log fold change") + 
    ylab("-log10(adjusted p-value)") +
    scale_y_continuous(expand = c(0,0), limits = y_lim) + xlim(x_lim) + 
    theme(legend.position = "none",
          plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
          axis.title = element_text(size = 15), 
          axis.text = element_text(size = 12, color = "black")) +
    geom_text_repel(data = DE.df_plus[DE.df_plus$FDR < plus_threshold, ], aes(x = logFC, y = -log10(FDR)), label = DE.df_plus[DE.df_plus$FDR < plus_threshold, ]$gene, force = plus_force) + 
    geom_text_repel(data = DE.df_minus[DE.df_minus$FDR < minus_threshold, ], aes(x = logFC, y = -log10(FDR)), label = DE.df_minus[DE.df_minus$FDR < minus_threshold, ]$gene, force = minus_force)
  
}

volcano_plot_no_il12b_no_labels <- function(DE.df, DE.df_plus, DE.df_minus, DE.df_il12, title = "title", y_lim = c(0,20), x_lim = c(-1, 1), plus_threshold = 1E-5, minus_threshold = 1E-5, plus_force = 10, minus_force = 10, il12b = 3){
  ggplot(DE.df) + geom_point(aes(x = logFC, y = -log10(FDR)), alpha = 0.3) +
    geom_point(data = DE.df_plus, aes(x = logFC, y = -log10(FDR)), color = "red", alpha = 1) +
    geom_point(data = DE.df_minus, aes(x = logFC, y = -log10(FDR)), color = "blue", alpha = 1) +
    ggtitle(title) + 
    theme_classic() + 
    xlab("log fold change") + 
    ylab("-log10(adjusted p-value)") +
    scale_y_continuous(expand = c(0,0), limits = y_lim) + xlim(x_lim) + 
    theme(legend.position = "none",
          plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
          axis.title = element_text(size = 15), 
          axis.text = element_text(size = 12, color = "black")) 
}

#########
# MOL
#########
AD_Ctrl.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/MOL_AD_Ctrl.csv", row.names = 1)

AD_Ctrl.df <- AD_Ctrl.df[AD_Ctrl.df$gene != "Ttr", ]

AD_Ctrl_plus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC > 0.5), ]
AD_Ctrl_minus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC < -0.5), ]
AD_Ctrl_il12b.df <- AD_Ctrl.df[which(AD_Ctrl.df$gene == "Il12b"), ]
v1 <- volcano_plot(AD_Ctrl.df, AD_Ctrl_plus.df, AD_Ctrl_minus.df, AD_Ctrl_il12b.df, title = "MOL APPPS1 vs WT", x_lim = c(-1, 1), y_lim = c(0,25), minus_force = 20, il12b = 1)

# Shirin request logFC +- 0.25 for all Oligodendrocytes cell type
AD_Ctrl_plus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC > 0.25), ]
AD_Ctrl_minus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC < -0.25), ]
AD_Ctrl_il12b.df <- AD_Ctrl.df[which(AD_Ctrl.df$gene == "Il12b"), ]
v1 <- volcano_plot(AD_Ctrl.df, AD_Ctrl_plus.df, AD_Ctrl_minus.df, AD_Ctrl_il12b.df, title = "MOL APPPS1 vs WT", x_lim = c(-1, 1), y_lim = c(0,25), minus_force = 20, il12b = 1)

# Shirin request no labels & logFC +- 0.25 for all Oligodendrocytes cell type
AD_Ctrl_plus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC > 0.25), ]
AD_Ctrl_minus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC < -0.25), ]
AD_Ctrl_il12b.df <- AD_Ctrl.df[which(AD_Ctrl.df$gene == "Il12b"), ]
v1 <- volcano_plot_no_labels(AD_Ctrl.df, AD_Ctrl_plus.df, AD_Ctrl_minus.df, AD_Ctrl_il12b.df, title = "MOL APPPS1 vs WT", x_lim = c(-1, 1), y_lim = c(0,25), minus_force = 20, il12b = 1)

# AD vs ADp40KO
AD_ADp40KO.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/MOL_AD_ADp40KO.csv", row.names = 1)
AD_ADp40KO.df <- AD_ADp40KO.df[AD_ADp40KO.df$gene != "Ttr", ]

AD_ADp40KO_plus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC > 0.5), ]
AD_ADp40KO_minus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC < -0.5), ]
AD_ADp40KO_il12b.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$gene == "Il12b"), ]
v2 <- volcano_plot(AD_ADp40KO.df, AD_ADp40KO_plus.df, AD_ADp40KO_minus.df, AD_ADp40KO_il12b.df, title = "MOL APPPS1 vs APPPS1.il12b-/-", x_lim = c(-1, 1), y_lim = c(0,25), minus_force = 40, il12b = 1)

# Shirin request logFC +- 0.25 for all Oligodendrocytes cell type
AD_ADp40KO_plus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC > 0.25), ]
AD_ADp40KO_minus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC < -0.25), ]
AD_ADp40KO_il12b.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$gene == "Il12b"), ]
v2 <- volcano_plot(AD_ADp40KO.df, AD_ADp40KO_plus.df, AD_ADp40KO_minus.df, AD_ADp40KO_il12b.df, title = "MOL APPPS1 vs APPPS1.il12b-/-", x_lim = c(-1, 1), y_lim = c(0,25), minus_force = 40, il12b = 1)

# Shirin request no labels & logFC +- 0.25 for all Oligodendrocytes cell type
AD_ADp40KO_plus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC > 0.25), ]
AD_ADp40KO_minus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC < -0.25), ]
AD_ADp40KO_il12b.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$gene == "Il12b"), ]
v2 <- volcano_plot_no_labels(AD_ADp40KO.df, AD_ADp40KO_plus.df, AD_ADp40KO_minus.df, AD_ADp40KO_il12b.df, title = "MOL APPPS1 vs APPPS1.il12b-/-", x_lim = c(-1, 1), y_lim = c(0,25), minus_force = 40, il12b = 1)


# ADp40KO vs Ctrl
ADp40KO_Ctrl.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/MOL_ADp40KO_Ctrl.csv", row.names = 1)
ADp40KO_Ctrl.df <- ADp40KO_Ctrl.df[ADp40KO_Ctrl.df$gene != "Ttr", ]

ADp40KO_Ctrl_plus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC > 0.5), ]
ADp40KO_Ctrl_minus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC < -0.5), ]
ADp40KO_Ctrl_il12b.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$gene == "Il12b"), ]
v3 <- volcano_plot(ADp40KO_Ctrl.df, ADp40KO_Ctrl_plus.df, ADp40KO_Ctrl_minus.df, ADp40KO_Ctrl_il12b.df, title = "MOL APPPS1.il12b-/- vs WT", x_lim = c(-1, 1), y_lim = c(0,25), plus_threshold = 1E-10, minus_threshold = 1E-10, minus_force = 40,
                   il12b = 2)

# Shirin request logFC +- 0.25 for all Oligodendrocytes cell type
ADp40KO_Ctrl_plus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC > 0.25), ]
ADp40KO_Ctrl_minus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC < -0.25), ]
ADp40KO_Ctrl_il12b.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$gene == "Il12b"), ]
v3 <- volcano_plot(ADp40KO_Ctrl.df, ADp40KO_Ctrl_plus.df, ADp40KO_Ctrl_minus.df, ADp40KO_Ctrl_il12b.df, title = "MOL APPPS1.il12b-/- vs WT", x_lim = c(-1, 1), y_lim = c(0,25), plus_threshold = 1E-10, minus_threshold = 1E-10, minus_force = 40,
                   il12b = 2)

# Shirin request no labels & logFC +- 0.25 for all Oligodendrocytes cell type
ADp40KO_Ctrl_plus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC > 0.25), ]
ADp40KO_Ctrl_minus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC < -0.25), ]
ADp40KO_Ctrl_il12b.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$gene == "Il12b"), ]
v3 <- volcano_plot_no_labels(ADp40KO_Ctrl.df, ADp40KO_Ctrl_plus.df, ADp40KO_Ctrl_minus.df, ADp40KO_Ctrl_il12b.df, title = "MOL APPPS1.il12b-/- vs WT", x_lim = c(-1, 1), y_lim = c(0,25), plus_threshold = 1E-10, minus_threshold = 1E-10, minus_force = 40,
                   il12b = 2)


g1 <- arrangeGrob(v1, v2, v3, ncol = 3)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig4_Volcano_Plots_for_All_MOL.png", 
       plot = g1, 
       scale = 1, width = 15, height = 4.5, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig4_Volcano_Plots_for_All_MOL_logFC_0.25.pdf", 
       plot = g1, 
       scale = 1, width = 15, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig4_Volcano_Plots_for_All_MOL_logFC_0.25_no_labels.pdf", 
       plot = g1, 
       scale = 1, width = 15, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

################
# MFOL
################

AD_Ctrl.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/MFOL_AD_Ctrl.csv", row.names = 1)
AD_Ctrl.df <- AD_Ctrl.df[AD_Ctrl.df$gene != "Ttr", ]

AD_Ctrl_plus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC > 0.5), ]
AD_Ctrl_minus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC < -0.5), ]
AD_Ctrl_il12b.df <- AD_Ctrl.df[which(AD_Ctrl.df$gene == "Il12b"), ]
v1 <- volcano_plot(AD_Ctrl.df, AD_Ctrl_plus.df, AD_Ctrl_minus.df, AD_Ctrl_il12b.df, title = "MFOL APPPS1 vs WT", x_lim = c(-1, 1), y_lim = c(0,90), minus_force = 20, il12b = 1)

# Shirin request logFC +- 0.25 for all Oligodendrocytes cell type
AD_Ctrl_plus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC > 0.25), ]
AD_Ctrl_minus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC < -0.25), ]
AD_Ctrl_il12b.df <- AD_Ctrl.df[which(AD_Ctrl.df$gene == "Il12b"), ]
v1 <- volcano_plot(AD_Ctrl.df, AD_Ctrl_plus.df, AD_Ctrl_minus.df, AD_Ctrl_il12b.df, title = "MFOL APPPS1 vs WT", x_lim = c(-1, 1), y_lim = c(0,90), plus_threshold = 1E-20, minus_threshold = 1E-20, minus_force = 20, il12b = 1)

# Shirin request no labels & logFC +- 0.25 for all Oligodendrocytes cell type
AD_Ctrl_plus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC > 0.25), ]
AD_Ctrl_minus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC < -0.25), ]
AD_Ctrl_il12b.df <- AD_Ctrl.df[which(AD_Ctrl.df$gene == "Il12b"), ]
v1 <- volcano_plot_no_labels(AD_Ctrl.df, AD_Ctrl_plus.df, AD_Ctrl_minus.df, AD_Ctrl_il12b.df, title = "MFOL APPPS1 vs WT", x_lim = c(-1, 1), y_lim = c(0,90), plus_threshold = 1E-20, minus_threshold = 1E-20, minus_force = 20, il12b = 1)



# AD vs ADp40KO
AD_ADp40KO.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/MFOL_AD_ADp40KO.csv", row.names = 1)
AD_ADp40KO.df <- AD_ADp40KO.df[AD_ADp40KO.df$gene != "Ttr", ]

AD_ADp40KO_plus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC > 0.5), ]
AD_ADp40KO_minus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC < -0.5), ]
AD_ADp40KO_il12b.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$gene == "Il12b"), ]
v2 <- volcano_plot(AD_ADp40KO.df, AD_ADp40KO_plus.df, AD_ADp40KO_minus.df, AD_ADp40KO_il12b.df, title = "MFOL APPPS1 vs APPPS1.il12b-/-", x_lim = c(-1, 1), y_lim = c(0,90), plus_threshold = 1E-20, minus_threshold = 1E-20, minus_force = 40, il12b = 1)

# Shirin request logFC +- 0.25 for all Oligodendrocytes cell type
AD_ADp40KO_plus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC > 0.25), ]
AD_ADp40KO_minus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC < -0.25), ]
AD_ADp40KO_il12b.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$gene == "Il12b"), ]
v2 <- volcano_plot(AD_ADp40KO.df, AD_ADp40KO_plus.df, AD_ADp40KO_minus.df, AD_ADp40KO_il12b.df, title = "MFOL APPPS1 vs APPPS1.il12b-/-", x_lim = c(-1, 1), y_lim = c(0,90), plus_threshold = 1E-20, minus_threshold = 1E-20, minus_force = 40, il12b = 1)

# Shirin request no labels logFC +- 0.25 for all Oligodendrocytes cell type
AD_ADp40KO_plus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC > 0.25), ]
AD_ADp40KO_minus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC < -0.25), ]
AD_ADp40KO_il12b.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$gene == "Il12b"), ]
v2 <- volcano_plot_no_labels(AD_ADp40KO.df, AD_ADp40KO_plus.df, AD_ADp40KO_minus.df, AD_ADp40KO_il12b.df, title = "MFOL APPPS1 vs APPPS1.il12b-/-", x_lim = c(-1, 1), y_lim = c(0,90), plus_threshold = 1E-20, minus_threshold = 1E-20, minus_force = 40, il12b = 1)


# ADp40KO vs Ctrl
ADp40KO_Ctrl.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/MFOL_ADp40KO_Ctrl.csv", row.names = 1)
ADp40KO_Ctrl.df <- ADp40KO_Ctrl.df[ADp40KO_Ctrl.df$gene != "Ttr", ]

ADp40KO_Ctrl_plus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC > 0.5), ]
ADp40KO_Ctrl_minus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC < -0.5), ]
ADp40KO_Ctrl_il12b.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$gene == "Il12b"), ]
v3 <- volcano_plot(ADp40KO_Ctrl.df, ADp40KO_Ctrl_plus.df, ADp40KO_Ctrl_minus.df, ADp40KO_Ctrl_il12b.df, title = "MFOL APPPS1.il12b-/- vs WT", x_lim = c(-1, 1), y_lim = c(0,90), plus_threshold = 1E-10, minus_threshold = 1E-10, minus_force = 40,
                   il12b = 10)

# Shirin request logFC +- 0.25 for all Oligodendrocytes cell type
ADp40KO_Ctrl_plus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC > 0.25), ]
ADp40KO_Ctrl_minus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC < -0.25), ]
ADp40KO_Ctrl_il12b.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$gene == "Il12b"), ]
v3 <- volcano_plot(ADp40KO_Ctrl.df, ADp40KO_Ctrl_plus.df, ADp40KO_Ctrl_minus.df, ADp40KO_Ctrl_il12b.df, title = "MFOL APPPS1.il12b-/- vs WT", x_lim = c(-1, 1), y_lim = c(0,90), plus_threshold = 1E-25, minus_threshold = 1E-30, minus_force = 40,
                   il12b = 10)

# Shirin request no labels & logFC +- 0.25 for all Oligodendrocytes cell type
ADp40KO_Ctrl_plus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC > 0.25), ]
ADp40KO_Ctrl_minus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC < -0.25), ]
ADp40KO_Ctrl_il12b.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$gene == "Il12b"), ]
v3 <- volcano_plot_no_labels(ADp40KO_Ctrl.df, ADp40KO_Ctrl_plus.df, ADp40KO_Ctrl_minus.df, ADp40KO_Ctrl_il12b.df, title = "MFOL APPPS1.il12b-/- vs WT", x_lim = c(-1, 1), y_lim = c(0,90), plus_threshold = 1E-25, minus_threshold = 1E-30, minus_force = 40,
                   il12b = 10)

g1 <- arrangeGrob(v1, v2, v3, ncol = 3)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig4_Volcano_Plots_for_All_MFOL.png", 
       plot = g1, 
       scale = 1, width = 15, height = 4.5, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig4_Volcano_Plots_for_All_MFOL_logFC_0.25.pdf", 
       plot = g1, 
       scale = 1, width = 15, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig4_Volcano_Plots_for_All_MFOL_logFC_0.25_no_labels.pdf", 
       plot = g1, 
       scale = 1, width = 15, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

###########################################
# NFOL -> OPC (change file names)
# Mistaken by NFOL and OPC in DE analysis
###########################################

AD_Ctrl.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/OPC_AD_Ctrl.csv", row.names = 1)
AD_Ctrl.df <- AD_Ctrl.df[AD_Ctrl.df$gene != "Ttr", ]

AD_Ctrl_plus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC > 0.5), ]
AD_Ctrl_minus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC < -0.5), ]
AD_Ctrl_il12b.df <- AD_Ctrl.df[which(AD_Ctrl.df$gene == "Il12b"), ]
v1 <- volcano_plot(AD_Ctrl.df, AD_Ctrl_plus.df, AD_Ctrl_minus.df, AD_Ctrl_il12b.df, title = "OPC APPPS1 vs WT", x_lim = c(-1.5, 1.5), y_lim = c(0,30), minus_force = 20, il12b = 1)

# Shirin request logFC +- 0.25 for all Oligodendrocytes cell type
AD_Ctrl_plus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC > 0.25), ]
AD_Ctrl_minus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC < -0.25), ]
AD_Ctrl_il12b.df <- AD_Ctrl.df[which(AD_Ctrl.df$gene == "Il12b"), ]
v1 <- volcano_plot(AD_Ctrl.df, AD_Ctrl_plus.df, AD_Ctrl_minus.df, AD_Ctrl_il12b.df, title = "OPC APPPS1 vs WT", x_lim = c(-1.5, 1.5), y_lim = c(0,30), minus_force = 20, il12b = 1)

# Shirin request no labels & logFC +- 0.25 for all Oligodendrocytes cell type
AD_Ctrl_plus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC > 0.25), ]
AD_Ctrl_minus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC < -0.25), ]
AD_Ctrl_il12b.df <- AD_Ctrl.df[which(AD_Ctrl.df$gene == "Il12b"), ]
v1 <- volcano_plot_no_labels(AD_Ctrl.df, AD_Ctrl_plus.df, AD_Ctrl_minus.df, AD_Ctrl_il12b.df, title = "OPC APPPS1 vs WT", x_lim = c(-1.5, 1.5), y_lim = c(0,30), minus_force = 20, il12b = 1)


# AD vs ADp40KO
AD_ADp40KO.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/OPC_AD_ADp40KO.csv", row.names = 1)
AD_ADp40KO.df <- AD_ADp40KO.df[AD_ADp40KO.df$gene != "Ttr", ]

AD_ADp40KO_plus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC > 0.5), ]
AD_ADp40KO_minus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC < -0.5), ]
AD_ADp40KO_il12b.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$gene == "Il12b"), ]
v2 <- volcano_plot(AD_ADp40KO.df, AD_ADp40KO_plus.df, AD_ADp40KO_minus.df, AD_ADp40KO_il12b.df, title = "OPC APPPS1 vs APPPS1.il12b-/-", x_lim = c(-1.5, 1.5), y_lim = c(0,30), minus_force = 10, il12b = 2)

# Shirin request logFC +- 0.25 for all Oligodendrocytes cell type
AD_ADp40KO_plus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC > 0.25), ]
AD_ADp40KO_minus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC < -0.25), ]
AD_ADp40KO_il12b.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$gene == "Il12b"), ]
v2 <- volcano_plot_no_il12b(AD_ADp40KO.df, AD_ADp40KO_plus.df, AD_ADp40KO_minus.df, AD_ADp40KO_il12b.df, title = "OPC APPPS1 vs APPPS1.il12b-/-", x_lim = c(-1.5, 1.5), y_lim = c(0,30), minus_force = 10, il12b = 2)

# Shirin request no labels & logFC +- 0.25 for all Oligodendrocytes cell type
AD_ADp40KO_plus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC > 0.25), ]
AD_ADp40KO_minus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC < -0.25), ]
AD_ADp40KO_il12b.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$gene == "Il12b"), ]
v2 <- volcano_plot_no_il12b_no_labels(AD_ADp40KO.df, AD_ADp40KO_plus.df, AD_ADp40KO_minus.df, AD_ADp40KO_il12b.df, title = "OPC APPPS1 vs APPPS1.il12b-/-", x_lim = c(-1.5, 1.5), y_lim = c(0,30), minus_force = 10, il12b = 2)


# ADp40KO vs Ctrl
ADp40KO_Ctrl.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/OPC_ADp40KO_Ctrl.csv", row.names = 1)
ADp40KO_Ctrl.df <- ADp40KO_Ctrl.df[ADp40KO_Ctrl.df$gene != "Ttr", ]

ADp40KO_Ctrl_plus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC > 0.5), ]
ADp40KO_Ctrl_minus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC < -0.5), ]
ADp40KO_Ctrl_il12b.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$gene == "Il12b"), ]
v3 <- volcano_plot(ADp40KO_Ctrl.df, ADp40KO_Ctrl_plus.df, ADp40KO_Ctrl_minus.df, ADp40KO_Ctrl_il12b.df, title = "OPC APPPS1.il12b-/- vs WT", x_lim = c(-1.5, 1.5), y_lim = c(0,30), plus_threshold = 1E-5, minus_threshold = 1E-5, minus_force = 20,
                   il12b = 2)

# Shirin request logFC +- 0.25 for all Oligodendrocytes cell type
ADp40KO_Ctrl_plus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC > 0.25), ]
ADp40KO_Ctrl_minus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC < -0.25), ]
ADp40KO_Ctrl_il12b.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$gene == "Il12b"), ]
v3 <- volcano_plot(ADp40KO_Ctrl.df, ADp40KO_Ctrl_plus.df, ADp40KO_Ctrl_minus.df, ADp40KO_Ctrl_il12b.df, title = "OPC APPPS1.il12b-/- vs WT", x_lim = c(-1.5, 1.5), y_lim = c(0,30), plus_threshold = 1E-5, minus_threshold = 1E-5, minus_force = 20,
                   il12b = 2)

# Shirin request logFC +- 0.25 for all Oligodendrocytes cell type
ADp40KO_Ctrl_plus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC > 0.25), ]
ADp40KO_Ctrl_minus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC < -0.25), ]
ADp40KO_Ctrl_il12b.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$gene == "Il12b"), ]
v3 <- volcano_plot_no_labels(ADp40KO_Ctrl.df, ADp40KO_Ctrl_plus.df, ADp40KO_Ctrl_minus.df, ADp40KO_Ctrl_il12b.df, title = "OPC APPPS1.il12b-/- vs WT", x_lim = c(-1.5, 1.5), y_lim = c(0,30), plus_threshold = 1E-5, minus_threshold = 1E-5, minus_force = 20,
                   il12b = 2)

g1 <- arrangeGrob(v1, v2, v3, ncol = 3)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig4_Volcano_Plots_for_All_OPC_New.png", 
       plot = g1, 
       scale = 1, width = 15, height = 4.5, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig4_Volcano_Plots_for_All_OPC_New_logFC_0.25.pdf", 
       plot = g1, 
       scale = 1, width = 15, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig4_Volcano_Plots_for_All_OPC_New_logFC_0.25_no_labels.pdf", 
       plot = g1, 
       scale = 1, width = 15, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

##################################
# NFOL
# OPC -> NFOL (change file names)
##################################

AD_Ctrl.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/NFOL_AD_Ctrl.csv", row.names = 1)
AD_Ctrl.df <- AD_Ctrl.df[AD_Ctrl.df$gene != "Ttr", ]

AD_Ctrl_plus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC > 0.5), ]
AD_Ctrl_minus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC < -0.5), ]
AD_Ctrl_il12b.df <- AD_Ctrl.df[which(AD_Ctrl.df$gene == "Il12b"), ]
v1 <- volcano_plot(AD_Ctrl.df, AD_Ctrl_plus.df, AD_Ctrl_minus.df, AD_Ctrl_il12b.df, title = "NFOL APPPS1 vs WT", x_lim = c(-2, 2), y_lim = c(0,7), minus_force = 20, minus_threshold = 1E-3, il12b = 0.5)

# Shirin request logFC +- 0.25 for all Oligodendrocytes cell type
AD_Ctrl_plus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC > 0.25), ]
AD_Ctrl_minus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC < -0.25), ]
AD_Ctrl_il12b.df <- AD_Ctrl.df[which(AD_Ctrl.df$gene == "Il12b"), ]
v1 <- volcano_plot(AD_Ctrl.df, AD_Ctrl_plus.df, AD_Ctrl_minus.df, AD_Ctrl_il12b.df, title = "NFOL APPPS1 vs WT", x_lim = c(-2, 2), y_lim = c(0,7), minus_force = 20, minus_threshold = 1E-3, il12b = 0.5)

# Shirin request no labels & logFC +- 0.25 for all Oligodendrocytes cell type
AD_Ctrl_plus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC > 0.25), ]
AD_Ctrl_minus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC < -0.25), ]
AD_Ctrl_il12b.df <- AD_Ctrl.df[which(AD_Ctrl.df$gene == "Il12b"), ]
v1 <- volcano_plot_no_labels(AD_Ctrl.df, AD_Ctrl_plus.df, AD_Ctrl_minus.df, AD_Ctrl_il12b.df, title = "NFOL APPPS1 vs WT", x_lim = c(-2, 2), y_lim = c(0,7), minus_force = 20, minus_threshold = 1E-3, il12b = 0.5)


# AD vs ADp40KO
AD_ADp40KO.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/NFOL_AD_ADp40KO.csv", row.names = 1)
AD_ADp40KO.df <- AD_ADp40KO.df[AD_ADp40KO.df$gene != "Ttr", ]

AD_ADp40KO_plus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC > 0.5), ]
AD_ADp40KO_minus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC < -0.5), ]
AD_ADp40KO_il12b.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$gene == "Il12b"), ]
v2 <- volcano_plot(AD_ADp40KO.df, AD_ADp40KO_plus.df, AD_ADp40KO_minus.df, AD_ADp40KO_il12b.df, title = "NFOL APPPS1 vs APPPS1.il12b-/-", x_lim = c(-2, 2), y_lim = c(0,7), minus_force = 10, il12b = 0.5)

# Shirin request logFC +- 0.25 for all Oligodendrocytes cell type
AD_ADp40KO_plus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC > 0.25), ]
AD_ADp40KO_minus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC < -0.25), ]
AD_ADp40KO_il12b.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$gene == "Il12b"), ]
v2 <- volcano_plot(AD_ADp40KO.df, AD_ADp40KO_plus.df, AD_ADp40KO_minus.df, AD_ADp40KO_il12b.df, title = "NFOL APPPS1 vs APPPS1.il12b-/-", x_lim = c(-2, 2), y_lim = c(0,7), minus_force = 10, il12b = 0.5)

# Shirin request no labels & logFC +- 0.25 for all Oligodendrocytes cell type
AD_ADp40KO_plus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC > 0.25), ]
AD_ADp40KO_minus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC < -0.25), ]
AD_ADp40KO_il12b.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$gene == "Il12b"), ]
v2 <- volcano_plot_no_labels(AD_ADp40KO.df, AD_ADp40KO_plus.df, AD_ADp40KO_minus.df, AD_ADp40KO_il12b.df, title = "NFOL APPPS1 vs APPPS1.il12b-/-", x_lim = c(-2, 2), y_lim = c(0,7), minus_force = 10, il12b = 0.5)


# ADp40KO vs Ctrl
ADp40KO_Ctrl.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/NFOL_ADp40KO_Ctrl.csv", row.names = 1)
ADp40KO_Ctrl.df <- ADp40KO_Ctrl.df[ADp40KO_Ctrl.df$gene != "Ttr", ]

ADp40KO_Ctrl_plus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC > 0.5), ]
ADp40KO_Ctrl_minus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC < -0.5), ]
ADp40KO_Ctrl_il12b.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$gene == "Il12b"), ]
v3 <- volcano_plot(ADp40KO_Ctrl.df, ADp40KO_Ctrl_plus.df, ADp40KO_Ctrl_minus.df, ADp40KO_Ctrl_il12b.df, title = "NFOL APPPS1.il12b-/- vs WT", x_lim = c(-2, 2), y_lim = c(0,7), plus_threshold = 1E-5, minus_threshold = 1E-5, minus_force = 20,
                   il12b = 0.5)

# Shirin request logFC +- 0.25 for all Oligodendrocytes cell type
ADp40KO_Ctrl_plus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC > 0.25), ]
ADp40KO_Ctrl_minus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC < -0.25), ]
ADp40KO_Ctrl_il12b.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$gene == "Il12b"), ]
v3 <- volcano_plot(ADp40KO_Ctrl.df, ADp40KO_Ctrl_plus.df, ADp40KO_Ctrl_minus.df, ADp40KO_Ctrl_il12b.df, title = "NFOL APPPS1.il12b-/- vs WT", x_lim = c(-2, 2), y_lim = c(0,7), plus_threshold = 1E-5, minus_threshold = 1E-5, minus_force = 20,
                   il12b = 0.5)

# Shirin request logFC +- 0.25 for all Oligodendrocytes cell type
ADp40KO_Ctrl_plus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC > 0.25), ]
ADp40KO_Ctrl_minus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC < -0.25), ]
ADp40KO_Ctrl_il12b.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$gene == "Il12b"), ]
v3 <- volcano_plot_no_labels(ADp40KO_Ctrl.df, ADp40KO_Ctrl_plus.df, ADp40KO_Ctrl_minus.df, ADp40KO_Ctrl_il12b.df, title = "NFOL APPPS1.il12b-/- vs WT", x_lim = c(-2, 2), y_lim = c(0,7), plus_threshold = 1E-5, minus_threshold = 1E-5, minus_force = 20,
                   il12b = 0.5)


g1 <- arrangeGrob(v1, v2, v3, ncol = 3)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig4_Volcano_Plots_for_All_NFOL_New.png", 
       plot = g1, 
       scale = 1, width = 15, height = 4.5, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig4_Volcano_Plots_for_All_NFOL_New_logFC_0.25.pdf", 
       plot = g1, 
       scale = 1, width = 15, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig4_Volcano_Plots_for_All_NFOL_New_logFC_0.25_no_labels.pdf", 
       plot = g1, 
       scale = 1, width = 15, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

################
# Astrocytes
################

AD_Ctrl.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/Astrocytes_AD_Ctrl.csv", row.names = 1)

AD_Ctrl.df <- AD_Ctrl.df[AD_Ctrl.df$gene != "Ttr", ]

AD_Ctrl_plus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC > 0.5), ]
AD_Ctrl_minus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC < -0.5), ]
AD_Ctrl_il12b.df <- AD_Ctrl.df[which(AD_Ctrl.df$gene == "Il12b"), ]

v1 <- volcano_plot(AD_Ctrl.df, AD_Ctrl_plus.df, AD_Ctrl_minus.df, AD_Ctrl_il12b.df, title = "Astrocytes APPPS1 vs WT", x_lim = c(-1.5, 1.5), y_lim = c(0,50), minus_force = 20, plus_threshold = 1E-12, minus_threshold = 1E-10, il12b = 3)

# AD vs ADp40KO
AD_ADp40KO.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/Astrocytes_AD_ADp40KO.csv", row.names = 1)
AD_ADp40KO.df <- AD_ADp40KO.df[AD_ADp40KO.df$gene != "Ttr", ]

AD_ADp40KO_plus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC > 0.5), ]
AD_ADp40KO_minus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC < -0.5), ]
AD_ADp40KO_il12b.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$gene == "Il12b"), ]

v2 <- volcano_plot(AD_ADp40KO.df, AD_ADp40KO_plus.df, AD_ADp40KO_minus.df, AD_ADp40KO_il12b.df, title = "Astrocytes APPPS1 vs APPPS1.il12b-/-", x_lim = c(-1.5, 1.5), y_lim = c(0,50), plus_threshold = 1E-8, minus_threshold = 1E-8, plus_force = 20, minus_force = 10, il12b = 3)

# ADp40KO vs Ctrl
ADp40KO_Ctrl.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/Astrocytes_ADp40KO_Ctrl.csv", row.names = 1)
ADp40KO_Ctrl.df <- ADp40KO_Ctrl.df[ADp40KO_Ctrl.df$gene != "Ttr", ]

ADp40KO_Ctrl_plus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC > 0.5), ]
ADp40KO_Ctrl_minus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC < -0.5), ]
ADp40KO_Ctrl_il12b.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$gene == "Il12b"), ]

v3 <- volcano_plot(ADp40KO_Ctrl.df, ADp40KO_Ctrl_plus.df, ADp40KO_Ctrl_minus.df, ADp40KO_Ctrl_il12b.df, title = "Astrocytes APPPS1.il12b-/- vs WT", x_lim = c(-1.5, 1.5), y_lim = c(0,50), plus_threshold = 1E-11, minus_threshold = 1E-10, minus_force = 20,
                   il12b = 3)

g1 <- arrangeGrob(v1, v2, v3, ncol = 3)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig4_Volcano_Plots_for_All_Astrocytes.png", 
       plot = g1, 
       scale = 1, width = 15, height = 4.5, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig4_Volcano_Plots_for_All_Astrocytes.pdf", 
       plot = g1, 
       scale = 1, width = 15, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

################
# Microglia
################

AD_Ctrl.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/Microglia_AD_Ctrl.csv", row.names = 1)
AD_Ctrl.df <- AD_Ctrl.df[AD_Ctrl.df$gene != "Ttr", ]

AD_Ctrl_plus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC > 0.5), ]
AD_Ctrl_minus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC < -0.5), ]
AD_Ctrl_il12b.df <- AD_Ctrl.df[which(AD_Ctrl.df$gene == "Il12b"), ]

v1 <- volcano_plot(AD_Ctrl.df, AD_Ctrl_plus.df, AD_Ctrl_minus.df, AD_Ctrl_il12b.df, title = "Microglia APPPS1 vs WT", x_lim = c(-2.5, 2.5), y_lim = c(0,150), minus_force = 20, plus_threshold = 1E-70, minus_threshold = 1E-40, il12b = 10)

# AD vs ADp40KO
AD_ADp40KO.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/Microglia_AD_ADp40KO.csv", row.names = 1)
AD_ADp40KO.df <- AD_ADp40KO.df[AD_ADp40KO.df$gene != "Ttr", ]

AD_ADp40KO_plus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC > 0.5), ]
AD_ADp40KO_minus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC < -0.5), ]
AD_ADp40KO_il12b.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$gene == "Il12b"), ]

v2 <- volcano_plot(AD_ADp40KO.df, AD_ADp40KO_plus.df, AD_ADp40KO_minus.df, AD_ADp40KO_il12b.df, title = "Microglia APPPS1 vs APPPS1.il12b-/-", x_lim = c(-2.5, 2.5), y_lim = c(0,150), plus_threshold = 1E-8, minus_threshold = 1E-8, plus_force = 20, minus_force = 10, il12b = 3)

# ADp40KO vs Ctrl
ADp40KO_Ctrl.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/Microglia_ADp40KO_Ctrl.csv", row.names = 1)
ADp40KO_Ctrl.df <- ADp40KO_Ctrl.df[ADp40KO_Ctrl.df$gene != "Ttr", ]

ADp40KO_Ctrl_plus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC > 0.5), ]
ADp40KO_Ctrl_minus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC < -0.5), ]
ADp40KO_Ctrl_il12b.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$gene == "Il12b"), ]

v3 <- volcano_plot(ADp40KO_Ctrl.df, ADp40KO_Ctrl_plus.df, ADp40KO_Ctrl_minus.df, ADp40KO_Ctrl_il12b.df, title = "Microglia APPPS1.il12b-/- vs WT", x_lim = c(-2.5, 2.5), y_lim = c(0,150), plus_threshold = 1E-65, minus_threshold = 1E-40, minus_force = 20,
                   il12b = 10)

g1 <- arrangeGrob(v1, v2, v3, ncol = 3)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig4_Volcano_Plots_for_All_Microglia.png", 
       plot = g1, 
       scale = 1, width = 15, height = 4.5, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig4_Volcano_Plots_for_All_Microglia.pdf", 
       plot = g1, 
       scale = 1, width = 15, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)


################
# Dentate_Gyrus
################

AD_Ctrl.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/Dentate_Gyrus_AD_Ctrl.csv", row.names = 1)
AD_Ctrl.df <- AD_Ctrl.df[AD_Ctrl.df$gene != "Ttr", ]

AD_Ctrl_plus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC > 0.5), ]
AD_Ctrl_minus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC < -0.5), ]
AD_Ctrl_il12b.df <- AD_Ctrl.df[which(AD_Ctrl.df$gene == "Il12b"), ]

v1 <- volcano_plot(AD_Ctrl.df, AD_Ctrl_plus.df, AD_Ctrl_minus.df, AD_Ctrl_il12b.df, title = "Dentate Gyrus APPPS1 vs WT", x_lim = c(-1, 1), y_lim = c(0,160), minus_force = 20, plus_threshold = 1E-30, minus_threshold = 1E-30, il12b = 2)

# AD vs ADp40KO
AD_ADp40KO.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/Dentate_Gyrus_AD_ADp40KO.csv", row.names = 1)
AD_ADp40KO.df <- AD_ADp40KO.df[AD_ADp40KO.df$gene != "Ttr", ]

AD_ADp40KO_plus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC > 0.5), ]
AD_ADp40KO_minus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC < -0.5), ]
AD_ADp40KO_il12b.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$gene == "Il12b"), ]

v2 <- volcano_plot(AD_ADp40KO.df, AD_ADp40KO_plus.df, AD_ADp40KO_minus.df, AD_ADp40KO_il12b.df, title = "Dentate Gyrus APPPS1 vs APPPS1.il12b-/-", x_lim = c(-1, 1), y_lim = c(0,160), plus_threshold = 1E-30, minus_threshold = 1E-30, plus_force = 20, minus_force = 10, il12b = 2)

# ADp40KO vs Ctrl
ADp40KO_Ctrl.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/Dentate_Gyrus_ADp40KO_Ctrl.csv", row.names = 1)
ADp40KO_Ctrl.df <- ADp40KO_Ctrl.df[ADp40KO_Ctrl.df$gene != "Ttr", ]

ADp40KO_Ctrl_plus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC > 0.5), ]
ADp40KO_Ctrl_minus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC < -0.5), ]
ADp40KO_Ctrl_il12b.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$gene == "Il12b"), ]

v3 <- volcano_plot(ADp40KO_Ctrl.df, ADp40KO_Ctrl_plus.df, ADp40KO_Ctrl_minus.df, ADp40KO_Ctrl_il12b.df, title = "Dentate Gyrus APPPS1.il12b-/- vs WT", x_lim = c(-1, 1), y_lim = c(0,160), plus_threshold = 1E-30, minus_threshold = 1E-30, minus_force = 20,
                   il12b = 15)

g1 <- arrangeGrob(v1, v2, v3, ncol = 3)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig4_Volcano_Plots_for_All_Dentate_Gyrus.png", 
       plot = g1, 
       scale = 1, width = 15, height = 4.5, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig4_Volcano_Plots_for_All_Dentate_Gyrus.pdf", 
       plot = g1, 
       scale = 1, width = 15, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

################
# CA1
################

AD_Ctrl.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/CA1_AD_Ctrl.csv", row.names = 1)
AD_Ctrl.df <- AD_Ctrl.df[AD_Ctrl.df$gene != "Ttr", ]

AD_Ctrl_plus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC > 0.5), ]
AD_Ctrl_minus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC < -0.5), ]
AD_Ctrl_il12b.df <- AD_Ctrl.df[which(AD_Ctrl.df$gene == "Il12b"), ]

v1 <- volcano_plot_no_il12b(AD_Ctrl.df, AD_Ctrl_plus.df, AD_Ctrl_minus.df, AD_Ctrl_il12b.df, title = "CA1 APPPS1 vs WT", x_lim = c(-1.5, 1.5), y_lim = c(0,160), minus_force = 20, plus_threshold = 1E-30, minus_threshold = 1E-30, il12b = 2)

# AD vs ADp40KO
AD_ADp40KO.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/CA1_AD_ADp40KO.csv", row.names = 1)
AD_ADp40KO.df <- AD_ADp40KO.df[AD_ADp40KO.df$gene != "Ttr", ]

AD_ADp40KO_plus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC > 0.5), ]
AD_ADp40KO_minus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC < -0.5), ]
AD_ADp40KO_il12b.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$gene == "Il12b"), ]

v2 <- volcano_plot_no_il12b(AD_ADp40KO.df, AD_ADp40KO_plus.df, AD_ADp40KO_minus.df, AD_ADp40KO_il12b.df, title = "CA1 APPPS1 vs APPPS1.il12b-/-", x_lim = c(-1.5, 1.5), y_lim = c(0,160), plus_threshold = 1E-30, minus_threshold = 1E-30, plus_force = 20, minus_force = 10, il12b = 2)

# ADp40KO vs Ctrl
ADp40KO_Ctrl.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/CA1_ADp40KO_Ctrl.csv", row.names = 1)
ADp40KO_Ctrl.df <- ADp40KO_Ctrl.df[ADp40KO_Ctrl.df$gene != "Ttr", ]

ADp40KO_Ctrl_plus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC > 0.5), ]
ADp40KO_Ctrl_minus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC < -0.5), ]
ADp40KO_Ctrl_il12b.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$gene == "Il12b"), ]

v3 <- volcano_plot(ADp40KO_Ctrl.df, ADp40KO_Ctrl_plus.df, ADp40KO_Ctrl_minus.df, ADp40KO_Ctrl_il12b.df, title = "CA1 APPPS1.il12b-/- vs WT", x_lim = c(-1.5, 1.5), y_lim = c(0,160), plus_threshold = 1E-35, minus_threshold = 1E-37, minus_force = 10,
                   il12b = 15)

g1 <- arrangeGrob(v1, v2, v3, ncol = 3)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig4_Volcano_Plots_for_All_CA1.png", 
       plot = g1, 
       scale = 1, width = 15, height = 4.5, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig4_Volcano_Plots_for_All_CA1.pdf", 
       plot = g1, 
       scale = 1, width = 15, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

################
# CA2/3
################

AD_Ctrl.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/CA2_3_AD_Ctrl.csv", row.names = 1)
AD_Ctrl.df <- AD_Ctrl.df[AD_Ctrl.df$gene != "Ttr", ]

AD_Ctrl_plus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC > 0.5), ]
AD_Ctrl_minus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC < -0.5), ]
AD_Ctrl_il12b.df <- AD_Ctrl.df[which(AD_Ctrl.df$gene == "Il12b"), ]

v1 <- volcano_plot_no_il12b(AD_Ctrl.df, AD_Ctrl_plus.df, AD_Ctrl_minus.df, AD_Ctrl_il12b.df, title = "CA2/3 APPPS1 vs WT", x_lim = c(-1.5, 1.5), y_lim = c(0,60), minus_force = 20, plus_threshold = 1E-15, minus_threshold = 1E-15, il12b = 2)

# AD vs ADp40KO
AD_ADp40KO.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/CA2_3_AD_ADp40KO.csv", row.names = 1)
AD_ADp40KO.df <- AD_ADp40KO.df[AD_ADp40KO.df$gene != "Ttr", ]

AD_ADp40KO_plus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC > 0.5), ]
AD_ADp40KO_minus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC < -0.5), ]
AD_ADp40KO_il12b.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$gene == "Il12b"), ]

v2 <- volcano_plot_no_il12b(AD_ADp40KO.df, AD_ADp40KO_plus.df, AD_ADp40KO_minus.df, AD_ADp40KO_il12b.df, title = "CA2/3 APPPS1 vs APPPS1.il12b-/-", x_lim = c(-1.5, 1.5), y_lim = c(0,60), plus_threshold = 1E-15, minus_threshold = 1E-15, plus_force = 20, minus_force = 10, il12b = 2)

# ADp40KO vs Ctrl
ADp40KO_Ctrl.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/CA2_3_ADp40KO_Ctrl.csv", row.names = 1)
ADp40KO_Ctrl.df <- ADp40KO_Ctrl.df[ADp40KO_Ctrl.df$gene != "Ttr", ]

ADp40KO_Ctrl_plus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC > 0.5), ]
ADp40KO_Ctrl_minus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC < -0.5), ]
ADp40KO_Ctrl_il12b.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$gene == "Il12b"), ]

v3 <- volcano_plot(ADp40KO_Ctrl.df, ADp40KO_Ctrl_plus.df, ADp40KO_Ctrl_minus.df, ADp40KO_Ctrl_il12b.df, title = "CA2/3 APPPS1.il12b-/- vs WT", x_lim = c(-1.5, 1.5), y_lim = c(0,60), plus_threshold = 1E-15, minus_threshold = 1E-15, minus_force = 10,
                   il12b = 7)

g1 <- arrangeGrob(v1, v2, v3, ncol = 3)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig4_Volcano_Plots_for_All_CA2_3.png", 
       plot = g1, 
       scale = 1, width = 15, height = 4.5, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig4_Volcano_Plots_for_All_CA2_3.pdf", 
       plot = g1, 
       scale = 1, width = 15, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

################
# subiculum
################

AD_Ctrl.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/subiculum_AD_Ctrl.csv", row.names = 1)
AD_Ctrl.df <- AD_Ctrl.df[AD_Ctrl.df$gene != "Ttr", ]

AD_Ctrl_plus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC > 0.5), ]
AD_Ctrl_minus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC < -0.5), ]
AD_Ctrl_il12b.df <- AD_Ctrl.df[which(AD_Ctrl.df$gene == "Il12b"), ]

v1 <- volcano_plot(AD_Ctrl.df, AD_Ctrl_plus.df, AD_Ctrl_minus.df, AD_Ctrl_il12b.df, title = "subiculum APPPS1 vs WT", x_lim = c(-1.5, 1.5), y_lim = c(0,160), minus_force = 20, plus_threshold = 1E-40, minus_threshold = 1E-35, plus_force = 20, il12b = 7)

# AD vs ADp40KO
AD_ADp40KO.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/subiculum_AD_ADp40KO.csv", row.names = 1)
AD_ADp40KO.df <- AD_ADp40KO.df[AD_ADp40KO.df$gene != "Ttr", ]

AD_ADp40KO_plus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC > 0.5), ]
AD_ADp40KO_minus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC < -0.5), ]
AD_ADp40KO_il12b.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$gene == "Il12b"), ]

v2 <- volcano_plot(AD_ADp40KO.df, AD_ADp40KO_plus.df, AD_ADp40KO_minus.df, AD_ADp40KO_il12b.df, title = "subiculum APPPS1 vs APPPS1.il12b-/-", x_lim = c(-1.5, 1.5), y_lim = c(0,160), plus_threshold = 1E-35, minus_threshold = 1E-32, plus_force = 20, minus_force = 10, il12b = 12)

# ADp40KO vs Ctrl
ADp40KO_Ctrl.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/subiculum_ADp40KO_Ctrl.csv", row.names = 1)
ADp40KO_Ctrl.df <- ADp40KO_Ctrl.df[ADp40KO_Ctrl.df$gene != "Ttr", ]

ADp40KO_Ctrl_plus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC > 0.5), ]
ADp40KO_Ctrl_minus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC < -0.5), ]
ADp40KO_Ctrl_il12b.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$gene == "Il12b"), ]

v3 <- volcano_plot(ADp40KO_Ctrl.df, ADp40KO_Ctrl_plus.df, ADp40KO_Ctrl_minus.df, ADp40KO_Ctrl_il12b.df, title = "subiculum APPPS1.il12b-/- vs WT", x_lim = c(-1.5, 1.5), y_lim = c(0,160), plus_threshold = 1E-30, minus_threshold = 1E-25, plus_force = 20, minus_force = 10,
                   il12b = 15)

g1 <- arrangeGrob(v1, v2, v3, ncol = 3)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig4_Volcano_Plots_for_All_subiculum.png", 
       plot = g1, 
       scale = 1, width = 15, height = 4.5, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig4_Volcano_Plots_for_All_subiculum.pdf", 
       plot = g1, 
       scale = 1, width = 15, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

################
# Inhibitory_Neurons
################

AD_Ctrl.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/Inhibitory_Neurons_AD_Ctrl.csv", row.names = 1)
AD_Ctrl.df <- AD_Ctrl.df[AD_Ctrl.df$gene != "Ttr", ]

AD_Ctrl_plus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC > 0.5), ]
AD_Ctrl_minus.df <- AD_Ctrl.df[which(AD_Ctrl.df$FDR < 0.01 & AD_Ctrl.df$logFC < -0.5), ]
AD_Ctrl_il12b.df <- AD_Ctrl.df[which(AD_Ctrl.df$gene == "Il12b"), ]

v1 <- volcano_plot_no_il12b(AD_Ctrl.df, AD_Ctrl_plus.df, AD_Ctrl_minus.df, AD_Ctrl_il12b.df, title = "Inhibitory Neurons APPPS1 vs WT", x_lim = c(-1.5, 1.5), y_lim = c(0,75), minus_force = 20, plus_threshold = 1E-40, minus_threshold = 1E-35, plus_force = 20, il12b = 7)

# AD vs ADp40KO
AD_ADp40KO.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/Inhibitory_Neurons_AD_ADp40KO.csv", row.names = 1)
AD_ADp40KO.df <- AD_ADp40KO.df[AD_ADp40KO.df$gene != "Ttr", ]

AD_ADp40KO_plus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC > 0.5), ]
AD_ADp40KO_minus.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$FDR < 0.01 & AD_ADp40KO.df$logFC < -0.5), ]
AD_ADp40KO_il12b.df <- AD_ADp40KO.df[which(AD_ADp40KO.df$gene == "Il12b"), ]

v2 <- volcano_plot_no_il12b(AD_ADp40KO.df, AD_ADp40KO_plus.df, AD_ADp40KO_minus.df, AD_ADp40KO_il12b.df, title = "Inhibitory Neurons APPPS1 vs APPPS1.il12b-/-", x_lim = c(-1.5, 1.5), y_lim = c(0,75), plus_threshold = 1E-15, minus_threshold = 1E-15, plus_force = 20, minus_force = 10, il12b = 12)

# ADp40KO vs Ctrl
ADp40KO_Ctrl.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/Inhibitory_Neurons_ADp40KO_Ctrl.csv", row.names = 1)
ADp40KO_Ctrl.df <- ADp40KO_Ctrl.df[ADp40KO_Ctrl.df$gene != "Ttr", ]

ADp40KO_Ctrl_plus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC > 0.5), ]
ADp40KO_Ctrl_minus.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$FDR < 0.01 & ADp40KO_Ctrl.df$logFC < -0.5), ]
ADp40KO_Ctrl_il12b.df <- ADp40KO_Ctrl.df[which(ADp40KO_Ctrl.df$gene == "Il12b"), ]

v3 <- volcano_plot(ADp40KO_Ctrl.df, ADp40KO_Ctrl_plus.df, ADp40KO_Ctrl_minus.df, ADp40KO_Ctrl_il12b.df, title = "Inhibitory Neurons APPPS1.il12b-/- vs WT", x_lim = c(-1.5, 1.5), y_lim = c(0,75), plus_threshold = 1E-18, minus_threshold = 1E-18, plus_force = 20, minus_force = 10,
                   il12b = 5)

g1 <- arrangeGrob(v1, v2, v3, ncol = 3)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig4_Volcano_Plots_for_All_Inhibitory_Neurons.png", 
       plot = g1, 
       scale = 1, width = 15, height = 4.5, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig4_Volcano_Plots_for_All_Inhibitory_Neurons.pdf", 
       plot = g1, 
       scale = 1, width = 15, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig4_Volcano_Plots_for_All_cell_type_session_info.txt")


