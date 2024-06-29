#! /bin/env RScript
# written by SJK at 31. Oct. 2020

library(gridExtra)
library(dplyr)
library(ggplot2)
library(tidyr)
library(Seurat)
library(ggrepel)
library(readxl)
library(tidyr)
library(Cairo)
library(destiny)
library(SCORPIUS)
library(viridis)

source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")


#####################
# Seurat Oligo object
#####################
load(file="/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200210_figures_for_meeting_report_cleaned_Seurat_object_res0.8.Robj")

data9set_cleaned.SO$cell_type <- data9set_cleaned.SO@active.ident

data9set_sub.SO <- subset(data9set_cleaned.SO, subset = cell_type %in% c("Oligo", "OPC"))

data9set_sub.meta <- data9set_sub.SO@meta.data %>% as.data.frame

data9set_sub.meta$cell_barcode <- rownames(data9set_sub.meta)

data9set_sub.meta <- data9set_sub.meta %>% mutate(cell_type_precise = case_when(seurat_clusters %in% 6 ~ "MOL",
                                                                                seurat_clusters %in% c(1,5) ~ "MFOL",
                                                                                seurat_clusters %in% 12 ~ "OPC",
                                                                                seurat_clusters %in% 38 ~ "NFOL"))

rownames(data9set_sub.meta) <- data9set_sub.meta$cell_barcode

group_name <- data9set_sub.SO@meta.data$seurat_clusters %>% as.character

group_name <- factor(group_name, levels = c(6, 1, 5, 38, 12))

space <- readRDS("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20201027_IL12b_downstream_analysis_scorpius_space.rda")
traj <- readRDS("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20201027_IL12b_downstream_analysis_scorpius_trajectory_for.rda")

traj.df <- traj$time %>% as.data.frame
colnames(traj.df) <- "pseudotime"

traj.df <- cbind(traj.df, data9set_sub.meta)

traj.df <- traj.df[, c("pseudotime", "orig.ident", "sample", "cell_type", "cell_barcode", "cell_type_precise")]

traj.df$cell_type_precise <- factor(traj.df$cell_type_precise, levels = c("OPC", "NFOL", "MFOL", "MOL"))


expression.df <- as.data.frame(data9set_sub.SO@assays$RNA@data)

#################
# Graph function
#################

pseudograph <- function(gene = "Pdgfra"){
  gene.df <- expression.df[rownames(expression.df) %in% gene, ] %>% t 
  traj_gene.df <- cbind(traj.df, gene.df)
  colnames(traj_gene.df) <- c("pseudotime", "orig.ident", "sample", "cell_type", "cell_barcode", "cell_type_precise", "gene")
  return(traj_gene.df)  
}

pseudograph2 <- function(gene = c("Pdgfra", "Mbp")){
  gene.df <- expression.df[rownames(expression.df) %in% gene, ] %>% t %>% as.data.frame
  gene.df$exp <- Matrix::rowMeans(gene.df)
  gene.df <- gene.df[,c("exp"), drop = FALSE] %>% as.data.frame
  
  traj_gene.df <- cbind(traj.df, gene.df)
  colnames(traj_gene.df) <- c("pseudotime", "orig.ident", "sample", "cell_type", "cell_barcode", "cell_type_precise", "gene")
  return(traj_gene.df)  
}

gene_plot <- function(data = pseudo_graph, plot_title = "Bmp4 gene", coord = c(0,1)){
  ggplot(data, aes(x = pseudotime, y = gene, color = sample)) + 
    geom_point(size = 0.1, alpha = 0.1) + ylab("gene expression") +
    stat_smooth(method = "gam", formula = y ~ s(x), size = 1, se = FALSE) + 
    ggtitle(plot_title) +  scale_color_manual(labels = c("WT", "APPPS1", "APPPS1.il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) + 
    theme_classic() + 
    theme(axis.text = element_text(size = 12, color = "black", family = "helvetica"),
          axis.title = element_text(size = 12, color = "black", family = "helvetica"),
          plot.title = element_text(hjust = 0.5, size = 12, color = "black", family = "helvetica"),
          legend.position = "none",
          legend.title = element_blank(),
          legend.text = element_text(size = 11, color = "black", family = "helvetica")) +
    coord_cartesian(ylim = coord)
}

#########################################################
# negative regulation of oligodendrocyte differentiation
#########################################################

pseudo_graph <- pseudograph2(gene = c("Bmp4","Ctnnb1","Daam2","Dusp10","Nf1","Notch1","Sirt2","Tmem98"))

g1 <- ggplot(pseudo_graph, aes(x = pseudotime, y = gene, color = sample)) + 
  geom_point(size = 0.1, alpha = 0.2) + 
  stat_smooth(method = "gam", formula = y ~ s(x), size = 1, se = FALSE) + 
  ggtitle("negative regulation of oligodendrocyte differentiation") + 
  ylab("gene expression") + 
  coord_cartesian(ylim = c(0,1)) + 
  scale_color_manual(labels = c("WT", "APPPS1", "APPPS1.il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 12, color = "black", family = "helvetica"),
        axis.title = element_text(size = 12, color = "black", family = "helvetica"),
        plot.title = element_text(hjust = 0.5, size = 12, color = "black", family = "helvetica"),
        legend.title = element_blank(),
        legend.text = element_text(size = 11, color = "black", family = "helvetica"))

g2 <- ggplot(data = traj.df) + geom_segment(
  mapping=aes(x=pseudotime, y=-0.8, xend=pseudotime, yend=-.5, color = cell_type_precise),
  size=0.2, show.legend = FALSE) + theme_void() + 
  scale_color_manual(values = c("#99CCCC","#9966FF",  "#669999","#006666"))

g2_grob = ggplotGrob(g2)
g3 <- g1 + annotation_custom(grob = g2_grob, xmin = -0.05, xmax = 1.05, 
                             ymin = -0.25, ymax = 0) + coord_cartesian(ylim = c(-0.25,1))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_supple_trajectory_Negative_Diff.png",
       plot = g3,
       scale = 1, width = 6, height = 3.8, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_supple_trajectory_Negative_Diff.svg",
       plot = g3,
       scale = 1, width = 6, height = 3.8, units = "in", device = "svg",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_supple_trajectory_Negative_Diff.pdf",
       plot = g3,
       scale = 1, width = 6, height = 3.8, units = "in", device = cairo_pdf,
       dpi = 300)

###############
# separate gene
###############
pseudo_graph <- pseudograph(gene = "Bmp4")
g1 <- gene_plot(data = pseudo_graph, plot_title = "Bmp4 gene", coord = c(0,2.5))
pseudo_graph <- pseudograph(gene = "Ctnnb1")
g2 <- gene_plot(data = pseudo_graph, plot_title = "Ctnnb1 gene", coord = c(0,2.5))
pseudo_graph <- pseudograph(gene = "Daam2")
g3 <- gene_plot(data = pseudo_graph, plot_title = "Daam2 gene", coord = c(0,2.5))
pseudo_graph <- pseudograph(gene = "Dusp10")
g4 <- gene_plot(data = pseudo_graph, plot_title = "Dusp10 gene", coord = c(0,2.5))
pseudo_graph <- pseudograph(gene = "Nf1")
g5 <- gene_plot(data = pseudo_graph, plot_title = "Nf1 gene", coord = c(0,2.5))
pseudo_graph <- pseudograph(gene = "Notch1")
g6 <- gene_plot(data = pseudo_graph, plot_title = "Notch1 gene", coord = c(0,2.5))
pseudo_graph <- pseudograph(gene = "Sirt2")
g7 <- gene_plot(data = pseudo_graph, plot_title = "Sirt2 gene", coord = c(0,2.5))
pseudo_graph <- pseudograph(gene = "Tmem98")
g8 <- gene_plot(data = pseudo_graph, plot_title = "Tmem98 gene", coord = c(0,2.5))

g_all <- arrangeGrob(g1, g2, g3, g4, g5, g6, g7, g8, ncol = 3)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_supple_trajectory_Negative_Diff_genes.png",
       plot = g_all,
       scale = 1, width = 7, height = 6, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_supple_trajectory_Negative_Diff_genes.pdf",
       plot = g_all,
       scale = 1, width = 7, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_supple_trajectory_Negative_Diff_genes.svg",
       plot = g_all,
       scale = 1, width = 7, height = 6, units = "in", device = "svg",
       dpi = 300)


#########################################################
# postive regulation of oligodendrocyte differentiation
#########################################################

pseudo_graph <- pseudograph2(gene = c("Aspa","Enpp2",
                                      "Ptn", "Nkx2-2", 
                                      "Ptpra","Ptprz1","Qk","Tenm4", "Zfp365"))


g1 <- ggplot(pseudo_graph, aes(x = pseudotime, y = gene, color = sample)) + 
  geom_point(size = 0.1, alpha = 0.2) + 
  stat_smooth(method = "gam", formula = y ~ s(x), size = 1, se = FALSE) + 
  ggtitle("postive regulation of oligodendrocyte differentiation") + 
  ylab("gene expression") + 
  coord_cartesian(ylim = c(0,1.5)) + 
  scale_color_manual(labels = c("WT", "APPPS1", "APPPS1.il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 12, color = "black", family = "helvetica"),
        axis.title = element_text(size = 12, color = "black", family = "helvetica"),
        plot.title = element_text(hjust = 0.5, size = 12, color = "black", family = "helvetica"),
        legend.title = element_blank(),
        legend.text = element_text(size = 11, color = "black", family = "helvetica"))

g2 <- ggplot(data = traj.df) + geom_segment(
  mapping=aes(x=pseudotime, y=-0.8, xend=pseudotime, yend=-.5, color = cell_type_precise),
  size=0.2, show.legend = FALSE) + theme_void() + 
  scale_color_manual(values = c("#99CCCC","#9966FF",  "#669999","#006666"))

g2_grob = ggplotGrob(g2)
g3 <- g1 + annotation_custom(grob = g2_grob, xmin = -0.05, xmax = 1.05, 
                             ymin = -0.25, ymax = 0) + coord_cartesian(ylim = c(-0.25,1.5))


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_supple_trajectory_Positive_Diff.png",
       plot = g3,
       scale = 1, width = 6, height = 3.8, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_supple_trajectory_Positive_Diff.svg",
       plot = g3,
       scale = 1, width = 6, height = 3.8, units = "in", device = "svg",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_supple_trajectory_Positive_Diff.pdf",
       plot = g3,
       scale = 1, width = 6, height = 3.8, units = "in", device = cairo_pdf,
       dpi = 300)

###############
# separate gene
###############

pseudo_graph <- pseudograph(gene = "Aspa")
g1 <- gene_plot(data = pseudo_graph, plot_title = "Aspa gene", coord = c(0,4))
pseudo_graph <- pseudograph(gene = "Enpp2")
g2 <- gene_plot(data = pseudo_graph, plot_title = "Enpp2 gene", coord = c(0,4))
pseudo_graph <- pseudograph(gene = "Nkx2-2")
g3 <- gene_plot(data = pseudo_graph, plot_title = "Nkx2-2 gene", coord = c(0,4))
pseudo_graph <- pseudograph(gene = "Ptn")
g4 <- gene_plot(data = pseudo_graph, plot_title = "Ptn gene", coord = c(0,4))
pseudo_graph <- pseudograph(gene = "Ptpra")
g5 <- gene_plot(data = pseudo_graph, plot_title = "Ptpra gene", coord = c(0,4))
pseudo_graph <- pseudograph(gene = "Ptprz1")
g6 <- gene_plot(data = pseudo_graph, plot_title = "Ptprz1 gene", coord = c(0,4))
pseudo_graph <- pseudograph(gene = "Qk")
g7 <- gene_plot(data = pseudo_graph, plot_title = "Qk gene", coord = c(0,4))
pseudo_graph <- pseudograph(gene = "Tenm4")
g8 <- gene_plot(data = pseudo_graph, plot_title = "Tenm4 gene", coord = c(0,4))
pseudo_graph <- pseudograph(gene = "Zfp365")
g9 <- gene_plot(data = pseudo_graph, plot_title = "Zfp365 gene", coord = c(0,4))


g_all <- arrangeGrob(g1, g2, g3, g4, g5, g6, g7, g8, g9, ncol = 3)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_supple_trajectory_Positive_Diff_genes.png",
       plot = g_all,
       scale = 1, width = 7, height = 6, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_supple_trajectory_Positive_Diff_genes.pdf",
       plot = g_all,
       scale = 1, width = 7, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_supple_trajectory_Positive_Diff_genes.svg",
       plot = g_all,
       scale = 1, width = 7, height = 6, units = "in", device = "svg",
       dpi = 300)

#########################################################
# negative regulation of myelination
#########################################################

pseudo_graph <- pseudograph2(gene = c("Jam2","Mtmr2","Pten","Tmem98", "Epha4"))


g1 <- ggplot(pseudo_graph, aes(x = pseudotime, y = gene, color = sample)) + 
  geom_point(size = 0.1, alpha = 0.2) + 
  stat_smooth(method = "gam", formula = y ~ s(x), size = 1, se = FALSE) + 
  ggtitle("negative regulation of myelination") + 
  ylab("gene expression") + 
  coord_cartesian(ylim = c(0,1)) + 
  scale_color_manual(labels = c("WT", "APPPS1", "APPPS1.il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 12, color = "black", family = "helvetica"),
        axis.title = element_text(size = 12, color = "black", family = "helvetica"),
        plot.title = element_text(hjust = 0.5, size = 12, color = "black", family = "helvetica"),
        legend.title = element_blank(),
        legend.text = element_text(size = 11, color = "black", family = "helvetica"))

g2 <- ggplot(data = traj.df) + geom_segment(
  mapping=aes(x=pseudotime, y=-0.8, xend=pseudotime, yend=-.5, color = cell_type_precise),
  size=0.2, show.legend = FALSE) + theme_void() + 
  scale_color_manual(values = c("#99CCCC","#9966FF",  "#669999","#006666"))

g2_grob = ggplotGrob(g2)
g3 <- g1 + annotation_custom(grob = g2_grob, xmin = -0.05, xmax = 1.05, 
                             ymin = -0.25, ymax = 0) + coord_cartesian(ylim = c(-0.25,1))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_supple_trajectory_Negative_myel.png",
       plot = g3,
       scale = 1, width = 6, height = 3.8, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_supple_trajectory_Negative_myel.svg",
       plot = g3,
       scale = 1, width = 6, height = 3.8, units = "in", device = "svg",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_supple_trajectory_Negative_myel.pdf",
       plot = g3,
       scale = 1, width = 6, height = 3.8, units = "in", device = cairo_pdf,
       dpi = 300)

###############
# separate gene
###############
pseudo_graph <- pseudograph(gene = "Jam2")
g1 <- gene_plot(data = pseudo_graph, plot_title = "Jam2 gene", coord = c(0,1))
pseudo_graph <- pseudograph(gene = "Mtmr2")
g2 <- gene_plot(data = pseudo_graph, plot_title = "Mtmr2 gene", coord = c(0,1))
pseudo_graph <- pseudograph(gene = "Pten")
g3 <- gene_plot(data = pseudo_graph, plot_title = "Pten gene", coord = c(0,1))
pseudo_graph <- pseudograph(gene = "Tmem98")
g4 <- gene_plot(data = pseudo_graph, plot_title = "Tmem98 gene", coord = c(0,1))
pseudo_graph <- pseudograph(gene = "Epha4")
g5 <- gene_plot(data = pseudo_graph, plot_title = "Epha4 gene", coord = c(0,1))

g_all <- arrangeGrob(g1, g2, g3, g4, g5, ncol = 3)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_supple_trajectory_Negative_myel_genes.png",
       plot = g_all,
       scale = 1, width = 7, height = 4, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_supple_trajectory_Negative_myel_genes.pdf",
       plot = g_all,
       scale = 1, width = 7, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_supple_trajectory_Negative_myel_genes.svg",
       plot = g_all,
       scale = 1, width = 7, height = 4, units = "in", device = "svg",
       dpi = 300)


#########################################################
# positive regulation of myelination
#########################################################

pseudo_graph <- pseudograph2(gene = c("Dicer1", "Mag","Myrf","Nrg1","Pard3",
                                      "Sox10","Tppp","Trf"))


g1 <- ggplot(pseudo_graph, aes(x = pseudotime, y = gene, color = sample)) + 
  geom_point(size = 0.1, alpha = 0.2) + 
  stat_smooth(method = "gam", formula = y ~ s(x), size = 1, se = FALSE) + 
  ggtitle("positive regulation of myelination") + 
  ylab("gene expression") + 
  coord_cartesian(ylim = c(0,1)) + 
  scale_color_manual(labels = c("WT", "APPPS1", "APPPS1.il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 12, color = "black", family = "helvetica"),
        axis.title = element_text(size = 12, color = "black", family = "helvetica"),
        plot.title = element_text(hjust = 0.5, size = 12, color = "black", family = "helvetica"),
        legend.title = element_blank(),
        legend.text = element_text(size = 11, color = "black", family = "helvetica"))

g2 <- ggplot(data = traj.df) + geom_segment(
  mapping=aes(x=pseudotime, y=-0.8, xend=pseudotime, yend=-.5, color = cell_type_precise),
  size=0.2, show.legend = FALSE) + theme_void() + 
  scale_color_manual(values = c("#99CCCC","#9966FF",  "#669999","#006666"))

g2_grob = ggplotGrob(g2)
g3 <- g1 + annotation_custom(grob = g2_grob, xmin = -0.05, xmax = 1.05, 
                             ymin = -0.25, ymax = 0) + coord_cartesian(ylim = c(-0.25,1))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_supple_trajectory_Positive_myel.png",
       plot = g3,
       scale = 1, width = 6, height = 3.8, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_supple_trajectory_Positive_myel.svg",
       plot = g3,
       scale = 1, width = 6, height = 3.8, units = "in", device = "svg",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_supple_trajectory_Positive_myel.pdf",
       plot = g3,
       scale = 1, width = 6, height = 3.8, units = "in", device = cairo_pdf,
       dpi = 300)


###############
# separate gene
###############

pseudo_graph <- pseudograph(gene = "Dicer1")
g1 <- gene_plot(data = pseudo_graph, plot_title = "Dicer1 gene", coord = c(0,2.5))
pseudo_graph <- pseudograph(gene = "Mag")
g2 <- gene_plot(data = pseudo_graph, plot_title = "Mag gene", coord = c(0,2.5))
pseudo_graph <- pseudograph(gene = "Myrf")
g3 <- gene_plot(data = pseudo_graph, plot_title = "Myrf gene", coord = c(0,2.5))
pseudo_graph <- pseudograph(gene = "Nrg1")
g4 <- gene_plot(data = pseudo_graph, plot_title = "Nrg1 gene", coord = c(0,2.5))
pseudo_graph <- pseudograph(gene = "Pard3")
g5 <- gene_plot(data = pseudo_graph, plot_title = "Pard3 gene", coord = c(0,2.5))
pseudo_graph <- pseudograph(gene = "Sox10")
g6 <- gene_plot(data = pseudo_graph, plot_title = "Sox10 gene", coord = c(0,2.5))
pseudo_graph <- pseudograph(gene = "Tppp")
g7 <- gene_plot(data = pseudo_graph, plot_title = "Tppp gene", coord = c(0,2.5))
pseudo_graph <- pseudograph(gene = "Trf")
g8 <- gene_plot(data = pseudo_graph, plot_title = "Trf gene", coord = c(0,2.5))

g_all <- arrangeGrob(g1, g2, g3, g4, g5, g6, g7, g8, ncol = 3)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_supple_trajectory_Positive_myel_genes.png",
       plot = g_all,
       scale = 1, width = 7, height = 6, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_supple_trajectory_Positive_myel_genes.pdf",
       plot = g_all,
       scale = 1, width = 7, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_supple_trajectory_Positive_myel_genes.svg",
       plot = g_all,
       scale = 1, width = 7, height = 6, units = "in", device = "svg",
       dpi = 300)



