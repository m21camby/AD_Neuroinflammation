#! /bin/env RScript
# written by SJK at 30. Oct. 2020

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
traj2 <- readRDS("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20201027_IL12b_downstream_analysis_scorpius_trajectory_rev.rda")


palette <- 
  setNames(
    c("#999999", "#009999", "#006666", "#99CCCC", "red"),
    c(1,5,6,12,38)
  )


d1 <- draw_trajectory_plot(
  space,
  progression_group = group_name,
  path = traj2$path,
  point_size = 2,
  point_alpha = 0.8,
  path_size = 1, 
  contour = FALSE,
  progression_group_palette = palette) + 
  ggtitle("Oligo differentiation trajectory") + 
  scale_color_manual(name = "Cell type", labels = c("MOL(cluster 6)", "MFOL(cluster 1)", "MFOL(cluster 5)",  "NFOL(cluster 38)", "OPC(cluster 12)"),
                       values = c("#006666", "#669999", "#009999", "#9966FF",  "#99CCCC")) + 
  theme(axis.text = element_text(size = 12, color = "black", family = "helvetica"),
        axis.title = element_text(size = 12, color = "black", family = "helvetica"),
        plot.title = element_text(hjust = 0.5, size = 15, color = "black", family = "helvetica"),
        legend.title = element_text(size = 12, color = "black", family = "helvetica"),
        legend.text = element_text(size = 11, color = "black", family = "helvetica"))


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_trajectory.png",
       plot = d1,
       scale = 1, width = 7, height = 5.5, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_trajectory.svg",
       plot = d1,
       scale = 1, width = 7, height = 5.5, units = "in", device = "svg",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_trajectory.pdf",
       plot = d1,
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

######################
# genes trajectory
######################

traj <- readRDS("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20201027_IL12b_downstream_analysis_scorpius_trajectory_for.rda")


traj.df <- traj$time %>% as.data.frame
colnames(traj.df) <- "pseudotime"

traj.df <- cbind(traj.df, data9set_sub.meta)

traj.df <- traj.df[, c("pseudotime", "orig.ident", "sample", "cell_type", "cell_barcode", "cell_type_precise")]

traj.df$cell_type_precise <- factor(traj.df$cell_type_precise, levels = c("OPC", "NFOL", "MFOL", "MOL"))

expression.df <- as.data.frame(data9set_sub.SO@assays$RNA@data)

pseudograph2 <- function(gene = c("Pdgfra", "Mbp")){
  gene.df <- expression.df[rownames(expression.df) %in% gene, ] %>% t %>% as.data.frame
  gene.df$exp <- Matrix::rowMeans(gene.df)
  gene.df <- gene.df[,c("exp"), drop = FALSE] %>% as.data.frame
  
  traj_gene.df <- cbind(traj.df, gene.df)
  colnames(traj_gene.df) <- c("pseudotime", "orig.ident", "sample", "cell_type", "cell_barcode", "cell_type_precise", "gene")
  return(traj_gene.df)  
}

pseudo_graph <- pseudograph2(gene = c("Cdo1", "Pdgfra", "Cspg4"))

g1 <- ggplot(pseudo_graph, aes(x = pseudotime, y = gene, color = sample)) + 
  geom_point(size = 0.1, alpha = 0.2) + 
  stat_smooth(method = "gam", formula = y ~ s(x), size = 1, se = FALSE) + 
  ggtitle("OPC marker genes (Cdo1,Pdgfra,Cspg4)") + 
  ylab("gene expression") + 
  coord_cartesian(ylim = c(0,1)) + 
  scale_color_manual(labels = c("WT", "APPPS1", "APPPS1.il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 12, color = "black", family = "helvetica"),
        axis.title = element_text(size = 12, color = "black", family = "helvetica"),
        plot.title = element_text(hjust = 0.5, size = 15, color = "black", family = "helvetica"),
        legend.title = element_blank(),
        legend.text = element_text(size = 11, color = "black", family = "helvetica"))

g2 <- ggplot(data = traj.df) + geom_segment(
  mapping=aes(x=pseudotime, y=-0.8, xend=pseudotime, yend=-.5, color = cell_type_precise),
  size=0.2, show.legend = FALSE) + theme_void() + 
  scale_color_manual(values = c("#99CCCC","#9966FF",  "#669999","#006666"))

g2_grob = ggplotGrob(g2)
g3 <- g1 + annotation_custom(grob = g2_grob, xmin = -0.05, xmax = 1.05, 
                       ymin = -0.25, ymax = 0) + coord_cartesian(ylim = c(-0.25,1.5))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_trajectory_OPC.png",
       plot = g3,
       scale = 1, width = 6, height = 3.8, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_trajectory_OPC.svg",
       plot = g3,
       scale = 1, width = 6, height = 3.8, units = "in", device = "svg",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_trajectory_OPC.pdf",
       plot = g3,
       scale = 1, width = 6, height = 3.8, units = "in", device = cairo_pdf,
       dpi = 300)

###############
# NFOL markers
###############

pseudo_graph <- pseudograph2(gene = c("Fyn", "Enpp6", "Tmem163"))

g1 <- ggplot(pseudo_graph, aes(x = pseudotime, y = gene, color = sample)) + 
  geom_point(size = 0.1, alpha = 0.2) + 
  stat_smooth(method = "gam", formula = y ~ s(x), size = 1, se = FALSE) + 
  ggtitle("NFOL marker genes (Fyn,Enpp6,Tmem163)") + 
  ylab("gene expression") + 
  coord_cartesian(ylim = c(0,1)) + 
  scale_color_manual(labels = c("WT", "APPPS1", "APPPS1.il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 12, color = "black", family = "helvetica"),
        axis.title = element_text(size = 12, color = "black", family = "helvetica"),
        plot.title = element_text(hjust = 0.5, size = 15, color = "black", family = "helvetica"),
        legend.title = element_blank(),
        legend.text = element_text(size = 11, color = "black", family = "helvetica"))

g2 <- ggplot(data = traj.df) + geom_segment(
  mapping=aes(x=pseudotime, y=-0.8, xend=pseudotime, yend=-.5, color = cell_type_precise),
  size=0.2, show.legend = FALSE) + theme_void() + 
  scale_color_manual(values = c("#99CCCC","#9966FF",  "#669999","#006666"))

g2_grob = ggplotGrob(g2)

g4 <- g1 + annotation_custom(grob = g2_grob, xmin = -0.05, xmax = 1.05, 
                             ymin = -0.25, ymax = 0) + coord_cartesian(ylim = c(-0.25,2.5))


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_trajectory_NFOL.png",
       plot = g4,
       scale = 1, width = 6, height = 3.8, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_trajectory_NFOL.svg",
       plot = g4,
       scale = 1, width = 6, height = 3.8, units = "in", device = "svg",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_trajectory_NFOL.pdf",
       plot = g4,
       scale = 1, width = 6, height = 3.8, units = "in", device = cairo_pdf,
       dpi = 300)

# MOL markers
pseudo_graph <- pseudograph2(gene = c("Mbp", "Mal", "Mdrg1", "Apod"))

g1 <- ggplot(pseudo_graph, aes(x = pseudotime, y = gene, color = sample)) + 
  geom_point(size = 0.1, alpha = 0.2) + 
  stat_smooth(method = "gam", formula = y ~ s(x), size = 1, se = FALSE) + 
  ggtitle("MOL marker genes (Mbp,Mal,Mdrg1,Apod)") + 
  ylab("gene expression") + 
  coord_cartesian(ylim = c(0,1)) + 
  scale_color_manual(labels = c("WT", "APPPS1", "APPPS1.il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 12, color = "black", family = "helvetica"),
        axis.title = element_text(size = 12, color = "black", family = "helvetica"),
        plot.title = element_text(hjust = 0.5, size = 15, color = "black", family = "helvetica"),
        legend.title = element_blank(),
        legend.text = element_text(size = 11, color = "black", family = "helvetica"))

g2 <- ggplot(data = traj.df) + geom_segment(
  mapping=aes(x=pseudotime, y=-0.8, xend=pseudotime, yend=-.5, color = cell_type_precise),
  size=0.2, show.legend = FALSE) + theme_void() + 
  scale_color_manual(values = c("#99CCCC","#9966FF",  "#669999","#006666"))

g2_grob = ggplotGrob(g2)


g5 <- g1 + annotation_custom(grob = g2_grob, xmin = -0.05, xmax = 1.05, 
                             ymin = -0.25, ymax = 0) + coord_cartesian(ylim = c(-0.25,2.5))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_trajectory_MOL.png",
       plot = g5,
       scale = 1, width = 6, height = 3.8, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_trajectory_MOL.svg",
       plot = g5,
       scale = 1, width = 6, height = 3.8, units = "in", device = "svg",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_trajectory_MOL.pdf",
       plot = g5,
       scale = 1, width = 6, height = 3.8, units = "in", device = cairo_pdf,
       dpi = 300)

###################
# OPC & NFOL & MOL
###################


g3_2nd <- g3 + theme(legend.position = "none") + coord_cartesian(ylim = c(-0.25,2.5)) + 
  ggtitle("OPC marker genes")
g4_2nd <- g4 + theme(legend.position = "none") + 
  ggtitle("NFOL marker genes")
g5_2nd <- g5 + ggtitle("MOL marker genes")

g6 <- arrangeGrob(g3_2nd, g4_2nd, g5_2nd, ncol = 3, widths = c(1,1,1.5))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_trajectory_All_cell_type.png",
       plot = g6,
       scale = 1, width = 12, height = 4, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_trajectory_All_cell_type.svg",
       plot = g6,
       scale = 1, width = 12, height = 4, units = "in", device = "svg",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_trajectory_All_cell_type.pdf",
       plot = g6,
       scale = 1, width = 12, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)

########################
# Density of pseudotime
########################
traj.df$merged_time <- round(traj.df$pseudotime,2)

density.data <- traj.df %>%
  dplyr::select(pseudotime, orig.ident, sample, merged_time) %>%
  group_by(sample, merged_time) %>%
  summarise(count = n()) %>%
  mutate(total = sum(count)) %>%
  mutate(cell_prop = count / total)


g7 <- ggplot(density.data, aes(x = merged_time, y = cell_prop, color = sample)) + 
  geom_point(size = 0.2, alpha = 0.2) + 
  stat_smooth(method = "loess", formula = y ~ x, size = 1, se = FALSE) +
#  stat_smooth(method = "gam", formula = y ~ s(x), size = 1, se = FALSE)
  ylab("cell proportion density") + xlab("pseudotime") + 
  coord_cartesian(ylim = c(0, 0.03)) + 
  scale_color_manual(labels = c("WT", "APPPS1", "APPPS1.il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) + 
  theme_cowplot() + 
  theme(axis.text = element_text(size = 12, color = "black", family = "helvetica"),
        axis.title = element_text(size = 12, color = "black", family = "helvetica"),
        plot.title = element_text(hjust = 0.5, size = 15, color = "black", family = "helvetica"),
        legend.title = element_blank(),
        legend.text = element_text(size = 11, color = "black", family = "helvetica"))

g8 <- ggplot(data = traj.df) + geom_segment(
  mapping=aes(x=pseudotime, y=-0.8, xend=pseudotime, yend=-.5, color = cell_type_precise),
  size=0.2, show.legend = FALSE) + theme_void() + 
  scale_color_manual(values = c("#99CCCC","#9966FF",  "#669999","#006666"))


g8_grob = ggplotGrob(g8)
g9 <- g7 + annotation_custom(grob = g8_grob, xmin = -0.05, xmax = 1.05, 
                             ymin = -0.01, ymax = -0.002) + coord_cartesian(ylim = c(-0.01,0.03))

g9

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_trajectory_cell_density.pdf",
       plot = g9,
       scale = 1, width = 8, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_trajectory_New_trajectory_session_info.txt")


