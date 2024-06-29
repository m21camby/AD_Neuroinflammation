#! /bin/env RScript
# written by SJK at 30. Sep. 2020
# This file is for Fig 6. lncRNAs analysis part figures
# What Shirin asked: Pvt1 and Neat 1 in microglia & Rmst1 in astrocytes, subiculum and dentate gyrus

library(gridExtra)
library(dplyr)
library(ggplot2)
library(tidyr)
library(Seurat)
library(ggrepel)
source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")



load(file="/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200210_figures_for_meeting_report_cleaned_Seurat_object_res0.8.Robj")

# New assigned cell type

meta.df <- data9set_cleaned.SO@meta.data %>% mutate(cell_type = case_when(seurat_clusters %in% c(0, 10) ~ "Dentate_Gyrus",
                                                                          seurat_clusters %in% c(2, 13, 18, 40) ~ "CA1",
                                                                          seurat_clusters %in% c(9, 21, 27, 32) ~ "CA2_3",
                                                                          seurat_clusters %in% c(11, 14, 15, 17, 19, 20) ~ "subiculum",
                                                                          seurat_clusters %in% c(24, 29) ~ "Unidentified_Neurons",
                                                                          seurat_clusters %in% c(7, 16, 22, 30) ~ "Inhibitory_Neurons",
                                                                          seurat_clusters %in% c(6) ~ "MOL",
                                                                          seurat_clusters %in% c(1, 5) ~ "MFOL",
                                                                          seurat_clusters %in% c(12) ~ "NFOL",
                                                                          seurat_clusters %in% c(38) ~ "OPC",
                                                                          seurat_clusters %in% c(3, 8) ~ "Microglia",
                                                                          seurat_clusters %in% c(4) ~ "Astrocytes",
                                                                          seurat_clusters %in% c(36) ~ "Vascular",
                                                                          seurat_clusters %in% c(39) ~ "VLMC",
                                                                          seurat_clusters %in% c(26) ~ "Choroid",
                                                                          seurat_clusters %in% c(23) ~ "Fibroblast",
                                                                          seurat_clusters %in% c(28) ~ "Cajal",
                                                                          seurat_clusters %in% c(35) ~ "Pericyte",
                                                                          seurat_clusters %in% c(37) ~ "Macrophage"))


data9set_cleaned.SO$cell_type <- meta.df$cell_type

# give levels to large cell type
data9set_cleaned.SO$cell_type <- factor(data9set_cleaned.SO$cell_type, levels = c("Dentate_Gyrus", "CA1", "CA2_3", "subiculum",  "Unidentified_Neurons", 
                                                                                  "Inhibitory_Neurons", "MOL", "MFOL", "OPC", "NFOL", "Microglia","Astrocytes",
                                                                                  "Vascular", "VLMC", "Choroid", "Fibroblast", "Cajal", "Pericyte", "Macrophage"))



CPM_calculation <- function(data9set_sub.SO){
  data9set_sub_ctrl.SO <- subset(data9set_sub.SO, subset = sample %in% "Ctrl")
  data9set_sub_AD.SO <- subset(data9set_sub.SO, subset = sample %in% "AD")
  data9set_sub_ADp40KO.SO <- subset(data9set_sub.SO, subset = sample %in% "ADp40KO")
  
  
  # 1. for ctrl
  data9set_cleaned_sub_ctrl.df <- as.data.frame(data9set_sub_ctrl.SO@assays$RNA@counts)
  cpm <- apply(data9set_cleaned_sub_ctrl.df, 2, function(x) (x/sum(x))*1000000)
  cpm.df <- rowMeans(cpm) %>% as.data.frame
  
  # add pseudocount 1
  cpm.df <- cpm.df + 1
  
  # log2
  cpm_log2.df <- as.data.frame(lapply(cpm.df, log2), row.names = rownames(cpm.df))
  
  # 2. for AD
  data9set_cleaned_sub_AD.df <- as.data.frame(data9set_sub_AD.SO@assays$RNA@counts)
  cpm2 <- apply(data9set_cleaned_sub_AD.df, 2, function(x) (x/sum(x))*1000000)
  cpm2.df <- rowMeans(cpm2) %>% as.data.frame
  
  # add pseudocount 1
  cpm2.df <- cpm2.df + 1
  
  # log2
  cpm2_log2.df <- as.data.frame(lapply(cpm2.df, log2), row.names = rownames(cpm2.df))
  
  # 3. for ADp40KO
  data9set_cleaned_sub_ADp40KO.df <- as.data.frame(data9set_sub_ADp40KO.SO@assays$RNA@counts)
  cpm3 <- apply(data9set_cleaned_sub_ADp40KO.df, 2, function(x) (x/sum(x))*1000000)
  cpm3.df <- rowMeans(cpm3) %>% as.data.frame
  
  # add pseudocount 1
  cpm3.df <- cpm3.df + 1
  
  # log2
  cpm3_log2.df <- as.data.frame(lapply(cpm3.df, log2), row.names = rownames(cpm3.df))
  
  cpm_final <- cbind(cpm_log2.df, cpm2_log2.df, cpm3_log2.df)
  
  colnames(cpm_final) <- c("WT", "AD", "ADp40KO")
  
  return(cpm_final)
}

scatter_plot <- function(cpm_final, x_axis = WT, y_axis = AD, xlabel = "x label", ylabel = "y label", cpm_final_sub){
  ggplot(cpm_final, aes_string(x = x_axis, y = y_axis)) + geom_point(size = 1, color = "#666666", alpha = 0.2 ) +
    geom_point(data = cpm_final_sub, aes_string(x = x_axis, y = y_axis), color = "red", size = 1.2) + 
    theme_classic() + scale_x_continuous(expand = c(0,0), limits = c(0, 15)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,15)) + 
    xlab(xlabel) +
    ylab(ylabel) +
    theme(plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
          axis.title = element_text(size = 12, family = "helvetica"),
          axis.text = element_text(size = 12, color = "black", family = "helvetica")) +
    geom_abline(slope = 1, color="red", linetype="dashed", size=.5) +
    geom_text_repel(data = cpm_final_sub, aes_string(x = x_axis, y = y_axis), label = rownames(cpm_final_sub), size = 5)
}

#########################
# Microglia (Pvt1, Neat1)
#########################

Seurat_sub.SO <- subset(data9set_cleaned.SO, subset = cell_type %in% c("Microglia"))

cpm_final <- CPM_calculation(Seurat_sub.SO)

cpm_final_sub <- cpm_final[rownames(cpm_final) %in% c("Pvt1", "Neat1", "Actb"), ]

s1 <- scatter_plot(cpm_final, x_axis = "WT", y_axis = "AD", xlabel = "log2 avg CPM + 1 Microglia WT", ylabel = "log2 avg CPM + 1 Microglia APPPS1 ", cpm_final_sub = cpm_final_sub)
s2 <- scatter_plot(cpm_final, x_axis = "ADp40KO", y_axis = "AD", xlabel = "log2 avg CPM + 1 Microglia APPPS1.il12b-/-", ylabel = "log2 avg CPM + 1 Microglia APPPS1 ", cpm_final_sub = cpm_final_sub)
s3 <- scatter_plot(cpm_final, x_axis = "WT", y_axis = "ADp40KO", xlabel = "log2 avg CPM + 1 Microglia WT", ylabel = "log2 avg CPM + 1 Microglia APPPS1.il12b-/- ", cpm_final_sub = cpm_final_sub)

g1 <- arrangeGrob(s1, s2, s3, ncol = 3)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig6_lncRNAs_Pvt1_Neat1_Microglia_AD_WT_scatter.png",
       plot = s1,
       scale = 1, width = 5, height = 4.5, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig6_lncRNAs_Pvt1_Neat1_Microglia_AD_WT_scatter.pdf",
       plot = s1,
       scale = 1, width = 5, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig6_lncRNAs_Pvt1_Neat1_Microglia_AD_ADp40KO_scatter.png",
       plot = s2,
       scale = 1, width = 5, height = 4.5, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig6_lncRNAs_Pvt1_Neat1_Microglia_AD_ADp40KO_scatter.pdf",
       plot = s2,
       scale = 1, width = 5, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig6_lncRNAs_Pvt1_Neat1_Microglia_ADp40KO_WT_scatter.png",
       plot = s3,
       scale = 1, width = 5, height = 4.5, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig6_lncRNAs_Pvt1_Neat1_Microglia_ADp40KO_WT_scatter.pdf",
       plot = s3,
       scale = 1, width = 5, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

#########################
# Astrocytes (Rmst) 
#########################

Seurat_sub.SO <- subset(data9set_cleaned.SO, subset = cell_type %in% c("Astrocytes"))

cpm_final <- CPM_calculation(Seurat_sub.SO)

cpm_final_sub <- cpm_final[rownames(cpm_final) %in% c("Rmst", "Actb"), ]

s1 <- scatter_plot(cpm_final, x_axis = "WT", y_axis = "AD", xlabel = "log2 avg CPM + 1 Astrocytes WT", ylabel = "log2 avg CPM + 1 Astrocytes APPPS1 ", cpm_final_sub = cpm_final_sub)
s2 <- scatter_plot(cpm_final, x_axis = "ADp40KO", y_axis = "AD", xlabel = "log2 avg CPM + 1 Astrocytes APPPS1.il12b-/-", ylabel = "log2 avg CPM + 1 Astrocytes APPPS1 ", cpm_final_sub = cpm_final_sub)
s3 <- scatter_plot(cpm_final, x_axis = "WT", y_axis = "ADp40KO", xlabel = "log2 avg CPM + 1 Astrocytes WT", ylabel = "log2 avg CPM + 1 Astrocytes APPPS1.il12b-/- ", cpm_final_sub = cpm_final_sub)

g1 <- arrangeGrob(s1, s2, s3, ncol = 3)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig6_lncRNAs_Rmst_Astrocytes_AD_WT_scatter.png",
       plot = s1,
       scale = 1, width = 5, height = 4.5, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig6_lncRNAs_Rmst_Astrocytes_AD_ADp40KO_scatter.png",
       plot = s2,
       scale = 1, width = 5, height = 4.5, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig6_lncRNAs_Rmst_Astrocytes_ADp40KO_WT_scatter.png",
       plot = s3,
       scale = 1, width = 5, height = 4.5, units = "in", device = "png",
       dpi = 300)

#########################
# subiculum (Rmst) 
#########################

Seurat_sub.SO <- subset(data9set_cleaned.SO, subset = cell_type %in% c("subiculum"))

cpm_final <- CPM_calculation(Seurat_sub.SO)

cpm_final_sub <- cpm_final[rownames(cpm_final) %in% c("Rmst", "Actb"), ]

s1 <- scatter_plot(cpm_final, x_axis = "WT", y_axis = "AD", xlabel = "log2 avg CPM + 1 subiculum WT", ylabel = "log2 avg CPM + 1 subiculum APPPS1 ", cpm_final_sub = cpm_final_sub)
s2 <- scatter_plot(cpm_final, x_axis = "ADp40KO", y_axis = "AD", xlabel = "log2 avg CPM + 1 subiculum APPPS1.il12b-/-", ylabel = "log2 avg CPM + 1 subiculum APPPS1 ", cpm_final_sub = cpm_final_sub)
s3 <- scatter_plot(cpm_final, x_axis = "WT", y_axis = "ADp40KO", xlabel = "log2 avg CPM + 1 subiculum WT", ylabel = "log2 avg CPM + 1 subiculum APPPS1.il12b-/- ", cpm_final_sub = cpm_final_sub)

g1 <- arrangeGrob(s1, s2, s3, ncol = 3)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig6_lncRNAs_Rmst_subiculum_AD_WT_scatter.png",
       plot = s1,
       scale = 1, width = 5, height = 4.5, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig6_lncRNAs_Rmst_subiculum_AD_ADp40KO_scatter.png",
       plot = s2,
       scale = 1, width = 5, height = 4.5, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig6_lncRNAs_Rmst_subiculum_ADp40KO_WT_scatter.png",
       plot = s3,
       scale = 1, width = 5, height = 4.5, units = "in", device = "png",
       dpi = 300)


#########################
# Dentate Gyrus (Rmst) 
#########################

Seurat_sub.SO <- subset(data9set_cleaned.SO, subset = cell_type %in% c("Dentate_Gyrus"))

cpm_final <- CPM_calculation(Seurat_sub.SO)

cpm_final_sub <- cpm_final[rownames(cpm_final) %in% c("Rmst", "Actb"), ]

s1 <- scatter_plot(cpm_final, x_axis = "WT", y_axis = "AD", xlabel = "log2 avg CPM + 1 Dentate Gyrus WT", ylabel = "log2 avg CPM + 1 Dentate Gyrus APPPS1 ", cpm_final_sub = cpm_final_sub)
s2 <- scatter_plot(cpm_final, x_axis = "ADp40KO", y_axis = "AD", xlabel = "log2 avg CPM + 1 Dentate Gyrus APPPS1.il12b-/-", ylabel = "log2 avg CPM + 1 Dentate Gyrus APPPS1 ", cpm_final_sub = cpm_final_sub)
s3 <- scatter_plot(cpm_final, x_axis = "WT", y_axis = "ADp40KO", xlabel = "log2 avg CPM + 1 Dentate Gyrus WT", ylabel = "log2 avg CPM + 1 Dentate Gyrus APPPS1.il12b-/- ", cpm_final_sub = cpm_final_sub)

g1 <- arrangeGrob(s1, s2, s3, ncol = 3)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig6_lncRNAs_Rmst_Dentate_Gyrus_AD_WT_scatter.png",
       plot = s1,
       scale = 1, width = 5, height = 4.5, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig6_lncRNAs_Rmst_Dentate_Gyrus_AD_ADp40KO_scatter.png",
       plot = s2,
       scale = 1, width = 5, height = 4.5, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig6_lncRNAs_Rmst_Dentate_Gyrus_ADp40KO_WT_scatter.png",
       plot = s3,
       scale = 1, width = 5, height = 4.5, units = "in", device = "png",
       dpi = 300)

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig6_lncRNAs_Pvt1_Neat1_Rmst_Scatter_AND_volcano_plot_session_info.txt")
