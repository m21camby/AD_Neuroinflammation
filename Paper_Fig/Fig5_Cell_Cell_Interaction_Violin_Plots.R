
#library(Seurat)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)
library(tidyr)
library(grid)
library(cowplot)
library(ggrepel)
library(viridis)
library(xlsx)
library(RColorBrewer)
#.libPaths(c("/data/murphy/shared_libs/R", "/usr/local/lib/R/site-librar","/usr/local/lib/R/library","/usr/lib64/R/library","/usr/share/R/library", .libPaths() ))
library(Seurat, lib.loc = "/data/rajewsky/home/skim/R")

source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

load(file="/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200210_figures_for_meeting_report_cleaned_Seurat_object_res0.8.Robj")

meta.df <- data9set_cleaned.SO@meta.data %>% mutate(cell_type = case_when(seurat_clusters %in% c(0, 10) ~ "Dentate_Gyrus",
                                                                          seurat_clusters %in% c(2, 13, 18, 40) ~ "CA1",
                                                                          seurat_clusters %in% c(9, 21, 27, 32) ~ "CA2_3",
                                                                          seurat_clusters %in% c(11, 14, 15, 17, 19, 20) ~ "subiculum",
                                                                          seurat_clusters %in% c(24, 29) ~ "Unidentified_Neurons",
                                                                          seurat_clusters %in% c(7, 16, 22, 30) ~ "Inhibitory_Neurons",
                                                                          seurat_clusters %in% c(6) ~ "MOL",
                                                                          seurat_clusters %in% c(1, 5) ~ "MFOL",
                                                                          seurat_clusters %in% c(12) ~ "OPC",
                                                                          seurat_clusters %in% c(38) ~ "NFOL",
                                                                          seurat_clusters %in% c(3, 8) ~ "Microglia",
                                                                          seurat_clusters %in% c(4) ~ "Astrocytes",
                                                                          seurat_clusters %in% c(36) ~ "Vascular",
                                                                          seurat_clusters %in% c(39) ~ "VLMC",
                                                                          seurat_clusters %in% c(26) ~ "Choroid",
                                                                          seurat_clusters %in% c(23) ~ "Fibroblast",
                                                                          seurat_clusters %in% c(28) ~ "Cajal",
                                                                          seurat_clusters %in% c(35) ~ "Pericyte",
                                                                          seurat_clusters %in% c(37) ~ "Macrophage"))


#data9set_cleaned_updated.SO = UpdateSeuratObject(object = data9set_cleaned.SO)

data9set_cleaned.SO$cell_type <- meta.df$cell_type

# give levels to large cell type
data9set_cleaned.SO$cell_type <- factor(data9set_cleaned.SO$cell_type, levels = c("Dentate_Gyrus", "CA1", "CA2_3", "subiculum",  "Unidentified_Neurons", 
                                                                                  "Inhibitory_Neurons", "MOL", "MFOL", "OPC", "NFOL", "Microglia","Astrocytes",
                                                                                  "Vascular", "VLMC", "Choroid", "Fibroblast", "Cajal", "Pericyte", "Macrophage"))

table(meta.df$cell_type)

Cell_type <- c("Dentate_Gyrus", "CA1", "CA2_3", "subiculum", "Inhibitory_Neurons", "MOL", "MFOL", "OPC", "NFOL", "Microglia", "Astrocytes")

########################
# Extract data function
########################
Vlnplot_gene_subset <- function(Cell_type = c("Inhibitory_Neurons"), gene = "Vegfa"){
  data9set_sub.SO <- subset(data9set_cleaned.SO, subset = cell_type %in% Cell_type)
  
  data9set_sub_ctrl.SO <- subset(data9set_sub.SO, subset = sample %in% c("Ctrl"))
  data9set_sub_AD.SO <- subset(data9set_sub.SO, subset = sample %in% c("AD"))
  data9set_sub_ADp40KO.SO <- subset(data9set_sub.SO, subset = sample %in% c("ADp40KO"))
  
  IN_Ctrl_Vegfa <- data9set_sub_ctrl.SO@assays$RNA@data[rownames(data9set_sub_ctrl.SO@assays$RNA@data) %in% gene, ] %>% as.data.frame
  colnames(IN_Ctrl_Vegfa) <- "gene"
  IN_Ctrl_Vegfa$celltype <- "WT"
  
  IN_AD_Vegfa <- data9set_sub_AD.SO@assays$RNA@data[rownames(data9set_sub_AD.SO@assays$RNA@data) %in% gene, ] %>% as.data.frame
  colnames(IN_AD_Vegfa) <- "gene"
  IN_AD_Vegfa$celltype <- "APPPS1"
  
  IN_ADp40KO_Vegfa <- data9set_sub_ADp40KO.SO@assays$RNA@data[rownames(data9set_sub_ADp40KO.SO@assays$RNA@data) %in% gene, ] %>% as.data.frame
  colnames(IN_ADp40KO_Vegfa) <- "gene"
  IN_ADp40KO_Vegfa$celltype <- "APPPS1.Il12b.-/-"
  
  IN_Vegfa <- rbind(IN_Ctrl_Vegfa, IN_AD_Vegfa, IN_ADp40KO_Vegfa)
  
  IN_Vegfa$celltype <- factor(IN_Vegfa$celltype, levels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"))
 
  
  return(IN_Vegfa)
   
}


######################
# 1. Vegfa IN
######################
Vegfa_IN <- Vlnplot_gene_subset(Cell_type = c("Inhibitory_Neurons"), gene = "Vegfa")

g1 <- ggplot(Vegfa_IN, aes_string(x="celltype", y="gene", color = "celltype", fill="celltype")) + 
  geom_violin(trim=TRUE, adjust = 3) + theme_cowplot() + ylab("Vegfa gene expression (Inhibitory)") + 
  theme(legend.position = "None", axis.title.x = element_blank()) + ylim(c(0,3)) + 
  scale_fill_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) + 
  scale_color_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) +
  geom_text(x=1, y=3, label=round(mean(Vegfa_IN[Vegfa_IN$celltype %in% "WT", ]$gene),3), size = 5, color = "black") +
  geom_text(x=2, y=3, label=round(mean(Vegfa_IN[Vegfa_IN$celltype %in% "APPPS1", ]$gene),3), size = 5, color = "black") +
  geom_text(x=3, y=3, label=round(mean(Vegfa_IN[Vegfa_IN$celltype %in% "APPPS1.Il12b.-/-", ]$gene),3), size = 5, color = "black")

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_Violin_Plots_Vegfa_IN.pdf",
       plot = g1,
       scale = 1, width = 6, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)

######################
# 2. Nrp1 CA2/3
######################

Nrp1_CA23 <- Vlnplot_gene_subset(Cell_type = c("CA2_3"), gene = "Nrp1")

g2 <- ggplot(Nrp1_CA23, aes_string(x="celltype", y="gene", color = "celltype", fill="celltype")) + 
  geom_violin(trim=TRUE, adjust = 3) + theme_cowplot() + ylab("Nrp1 gene expression (CA2/3)") + 
  theme(legend.position = "None", axis.title.x = element_blank()) + ylim(c(0,4)) + 
  scale_fill_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) + 
  scale_color_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) +
  geom_text(x=1, y=4, label=round(mean(Nrp1_CA23[Nrp1_CA23$celltype %in% "WT", ]$gene),3), size = 5, color = "black") +
  geom_text(x=2, y=4, label=round(mean(Nrp1_CA23[Nrp1_CA23$celltype %in% "APPPS1", ]$gene),3), size = 5, color = "black") +
  geom_text(x=3, y=4, label=round(mean(Nrp1_CA23[Nrp1_CA23$celltype %in% "APPPS1.Il12b.-/-", ]$gene),3), size = 5, color = "black")

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_Violin_Plots_Nrp1_CA23.pdf",
       plot = g2,
       scale = 1, width = 6, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)

######################
# 3. Grin2b CA1
######################

Grin2b_CA1 <- Vlnplot_gene_subset(Cell_type = c("CA1"), gene = "Grin2b")

g3 <- ggplot(Grin2b_CA1, aes_string(x="celltype", y="gene", color = "celltype", fill="celltype")) + 
  geom_violin(trim=TRUE, adjust = 3) + theme_cowplot() + ylab("Grin2b gene expression (CA1)") + 
  theme(legend.position = "None", axis.title.x = element_blank()) + ylim(c(0,5)) + 
  scale_fill_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) + 
  scale_color_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) +
  geom_text(x=1, y=5, label=round(mean(Grin2b_CA1[Grin2b_CA1$celltype %in% "WT", ]$gene),3), size = 5, color = "black") +
  geom_text(x=2, y=5, label=round(mean(Grin2b_CA1[Grin2b_CA1$celltype %in% "APPPS1", ]$gene),3), size = 5, color = "black") +
  geom_text(x=3, y=5, label=round(mean(Grin2b_CA1[Grin2b_CA1$celltype %in% "APPPS1.Il12b.-/-", ]$gene),3), size = 5, color = "black")

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_Violin_Plots_Grin2b_CA1.pdf",
       plot = g3,
       scale = 1, width = 6, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)

######################
# 4. Nrp2 DG
######################

Nrp2_DG <- Vlnplot_gene_subset(Cell_type = c("Dentate_Gyrus"), gene = "Nrp2")
  
g4 <- ggplot(Nrp2_DG, aes_string(x="celltype", y="gene", color = "celltype", fill="celltype")) + 
  geom_violin(trim=TRUE, adjust = 3) + theme_cowplot() + ylab("Nrp2 gene expression (Dentate Gyrus)") + 
  theme(legend.position = "None", axis.title.x = element_blank()) + ylim(c(0,4)) + 
  scale_fill_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) + 
  scale_color_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) +
  geom_text(x=1, y=4, label=round(mean(Nrp2_DG[Nrp2_DG$celltype %in% "WT", ]$gene),3), size = 5, color = "black") +
  geom_text(x=2, y=4, label=round(mean(Nrp2_DG[Nrp2_DG$celltype %in% "APPPS1", ]$gene),3), size = 5, color = "black") +
  geom_text(x=3, y=4, label=round(mean(Nrp2_DG[Nrp2_DG$celltype %in% "APPPS1.Il12b.-/-", ]$gene),3), size = 5, color = "black")

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_Violin_Plots_Nrp2_DG.pdf",
       plot = g4,
       scale = 1, width = 6, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)
######################
# 5. Nrg1 SB
######################

Nrg1_SB <- Vlnplot_gene_subset(Cell_type = c("subiculum"), gene = "Nrg1")

g5 <- ggplot(Nrg1_SB, aes_string(x="celltype", y="gene", color = "celltype", fill="celltype")) + 
  geom_violin(trim=TRUE, adjust = 3) + theme_cowplot() + ylab("Nrg1 gene expression (subiculum)") + 
  theme(legend.position = "None", axis.title.x = element_blank()) + ylim(c(0,6)) + 
  scale_fill_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) + 
  scale_color_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) +
  geom_text(x=1, y=6, label=round(mean(Nrg1_SB[Nrg1_SB$celltype %in% "WT", ]$gene),3), size = 5, color = "black") +
  geom_text(x=2, y=6, label=round(mean(Nrg1_SB[Nrg1_SB$celltype %in% "APPPS1", ]$gene),3), size = 5, color = "black") +
  geom_text(x=3, y=6, label=round(mean(Nrg1_SB[Nrg1_SB$celltype %in% "APPPS1.Il12b.-/-", ]$gene),3), size = 5, color = "black")

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_Violin_Plots_Nrg1_SB.pdf",
       plot = g5,
       scale = 1, width = 6, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)

######################
# 6. Erbb4 MOL
######################

Erbb4_MOL <- Vlnplot_gene_subset(Cell_type = c("MOL"), gene = "Erbb4")

g6 <- ggplot(Erbb4_MOL, aes_string(x="celltype", y="gene", color = "celltype", fill="celltype")) + 
  geom_violin(trim=TRUE, adjust = 3) + theme_cowplot() + ylab("Erbb4 gene expression (MOL)") + 
  theme(legend.position = "None", axis.title.x = element_blank()) + ylim(c(0,6)) + 
  scale_fill_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) + 
  scale_color_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) +
  geom_text(x=1, y=6, label=round(mean(Erbb4_MOL[Erbb4_MOL$celltype %in% "WT", ]$gene),3), size = 5, color = "black") +
  geom_text(x=2, y=6, label=round(mean(Erbb4_MOL[Erbb4_MOL$celltype %in% "APPPS1", ]$gene),3), size = 5, color = "black") +
  geom_text(x=3, y=6, label=round(mean(Erbb4_MOL[Erbb4_MOL$celltype %in% "APPPS1.Il12b.-/-", ]$gene),3), size = 5, color = "black")

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_Violin_Plots_Erbb4_MOL.pdf",
       plot = g6,
       scale = 1, width = 6, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)

######################
# Erbb3 OPC
######################

Erbb3_OPC <- Vlnplot_gene_subset(Cell_type = c("OPC"), gene = "Erbb3")


g7 <- ggplot(Erbb3_OPC, aes_string(x="celltype", y="gene", color = "celltype", fill="celltype")) + 
  geom_violin(trim=TRUE, adjust = 3) + theme_cowplot() + ylab("Erbb3 gene expression (OPC)") + 
  theme(legend.position = "None", axis.title.x = element_blank()) + ylim(c(0,4)) + 
  scale_fill_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) + 
  scale_color_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) +
  geom_text(x=1, y=4, label=round(mean(Erbb3_OPC[Erbb3_OPC$celltype %in% "WT", ]$gene),3), size = 5, color = "black") +
  geom_text(x=2, y=4, label=round(mean(Erbb3_OPC[Erbb3_OPC$celltype %in% "APPPS1", ]$gene),3), size = 5, color = "black") +
  geom_text(x=3, y=4, label=round(mean(Erbb3_OPC[Erbb3_OPC$celltype %in% "APPPS1.Il12b.-/-", ]$gene),3), size = 5, color = "black")


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_Violin_Plots_Erbb3_OPC.pdf",
       plot = g7,
       scale = 1, width = 6, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)



Erbb3_NFOL <- Vlnplot_gene_subset(Cell_type = c("NFOL"), gene = "Erbb3")

g8 <- ggplot(Erbb3_NFOL, aes_string(x="celltype", y="gene", color = "celltype", fill="celltype")) + 
  geom_violin(trim=TRUE, adjust = 3) + theme_cowplot() + ylab("Erbb3 gene expression (NFOL)") + 
  theme(legend.position = "None", axis.title.x = element_blank()) + ylim(c(0,3)) + 
  scale_fill_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) + 
  scale_color_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) +
  geom_text(x=1, y=3, label=round(mean(Erbb3_NFOL[Erbb3_NFOL$celltype %in% "WT", ]$gene),3), size = 5, color = "black") +
  geom_text(x=2, y=3, label=round(mean(Erbb3_NFOL[Erbb3_NFOL$celltype %in% "APPPS1", ]$gene),3), size = 5, color = "black") +
  geom_text(x=3, y=3, label=round(mean(Erbb3_NFOL[Erbb3_NFOL$celltype %in% "APPPS1.Il12b.-/-", ]$gene),3), size = 5, color = "black")

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_Violin_Plots_Erbb3_NFOL.pdf",
       plot = g8,
       scale = 1, width = 6, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)



Erbb4_NFOL <- Vlnplot_gene_subset(Cell_type = c("NFOL"), gene = "Erbb4")

g9 <- ggplot(Erbb4_NFOL, aes_string(x="celltype", y="gene", color = "celltype", fill="celltype")) + 
  geom_violin(trim=TRUE, adjust = 3) + theme_cowplot() + ylab("Erbb4 gene expression (NFOL)") + 
  theme(legend.position = "None", axis.title.x = element_blank()) + ylim(c(0,5)) + 
  scale_fill_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) + 
  scale_color_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) +
  geom_text(x=1, y=5, label=round(mean(Erbb4_NFOL[Erbb4_NFOL$celltype %in% "WT", ]$gene),3), size = 5, color = "black") +
  geom_text(x=2, y=5, label=round(mean(Erbb4_NFOL[Erbb4_NFOL$celltype %in% "APPPS1", ]$gene),3), size = 5, color = "black") +
  geom_text(x=3, y=5, label=round(mean(Erbb4_NFOL[Erbb4_NFOL$celltype %in% "APPPS1.Il12b.-/-", ]$gene),3), size = 5, color = "black")


Fgfr1_OPC <- Vlnplot_gene_subset(Cell_type = c("OPC"), gene = "Fgfr1")

ggplot(Fgfr1_OPC, aes_string(x="celltype", y="gene", color = "celltype", fill="celltype")) + 
  geom_violin(trim=TRUE, adjust = 3) + theme_cowplot()



Epha4_OPC <- Vlnplot_gene_subset(Cell_type = c("OPC"), gene = "Epha4")

g11 <- ggplot(Epha4_OPC, aes_string(x="celltype", y="gene", color = "celltype", fill="celltype")) + 
  geom_violin(trim=TRUE, adjust = 3) + theme_cowplot() + ylab("Epha4 gene expression (OPC)") + 
  theme(legend.position = "None", axis.title.x = element_blank()) + ylim(c(0,4)) + 
  scale_fill_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) + 
  scale_color_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) +
  geom_text(x=1, y=4, label=round(mean(Epha4_OPC[Epha4_OPC$celltype %in% "WT", ]$gene),3), size = 5, color = "black") +
  geom_text(x=2, y=4, label=round(mean(Epha4_OPC[Epha4_OPC$celltype %in% "APPPS1", ]$gene),3), size = 5, color = "black") +
  geom_text(x=3, y=4, label=round(mean(Epha4_OPC[Epha4_OPC$celltype %in% "APPPS1.Il12b.-/-", ]$gene),3), size = 5, color = "black")

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_Violin_Plots_Epha4_OPC.pdf",
       plot = g11,
       scale = 1, width = 6, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)

######################
# Erbb4 MFOL
######################
Erbb4_MFOL <- Vlnplot_gene_subset(Cell_type = c("MFOL"), gene = "Erbb4")

g12 <- ggplot(Erbb4_MFOL, aes_string(x="celltype", y="gene", color = "celltype", fill="celltype")) + 
  geom_violin(trim=TRUE, adjust = 3) + theme_cowplot() + ylab("Erbb4 gene expression (MFOL)") + 
  theme(legend.position = "None", axis.title.x = element_blank()) + ylim(c(0,6)) + 
  scale_fill_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) + 
  scale_color_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) +
  geom_text(x=1, y=6, label=round(mean(Erbb4_MFOL[Erbb4_MFOL$celltype %in% "WT", ]$gene),3), size = 5, color = "black") +
  geom_text(x=2, y=6, label=round(mean(Erbb4_MFOL[Erbb4_MFOL$celltype %in% "APPPS1", ]$gene),3), size = 5, color = "black") +
  geom_text(x=3, y=6, label=round(mean(Erbb4_MFOL[Erbb4_MFOL$celltype %in% "APPPS1.Il12b.-/-", ]$gene),3), size = 5, color = "black")

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_Violin_Plots_Erbb4_MFOL.pdf",
       plot = g12,
       scale = 1, width = 6, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)

######################
# Erbb3 MFOL
######################
Erbb3_MFOL <- Vlnplot_gene_subset(Cell_type = c("MFOL"), gene = "Erbb3")

g13 <- ggplot(Erbb3_MFOL, aes_string(x="celltype", y="gene", color = "celltype", fill="celltype")) + 
  geom_violin(trim=TRUE, adjust = 3) + theme_cowplot() + ylab("Erbb3 gene expression (MFOL)") + 
  theme(legend.position = "None", axis.title.x = element_blank()) + ylim(c(0,5)) + 
  scale_fill_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) + 
  scale_color_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) +
  geom_text(x=1, y=5, label=round(mean(Erbb3_MFOL[Erbb3_MFOL$celltype %in% "WT", ]$gene),3), size = 5, color = "black") +
  geom_text(x=2, y=5, label=round(mean(Erbb3_MFOL[Erbb3_MFOL$celltype %in% "APPPS1", ]$gene),3), size = 5, color = "black") +
  geom_text(x=3, y=5, label=round(mean(Erbb3_MFOL[Erbb3_MFOL$celltype %in% "APPPS1.Il12b.-/-", ]$gene),3), size = 5, color = "black")

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_Violin_Plots_Erbb3_MFOL.pdf",
       plot = g13,
       scale = 1, width = 6, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)

######################
# Erbb4 NFOL
######################
Erbb4_NFOL <- Vlnplot_gene_subset(Cell_type = c("NFOL"), gene = "Erbb4")

g14 <- ggplot(Erbb4_NFOL, aes_string(x="celltype", y="gene", color = "celltype", fill="celltype")) + 
  geom_violin(trim=TRUE, adjust = 3) + theme_cowplot() + ylab("Erbb4 gene expression (NFOL)") + 
  theme(legend.position = "None", axis.title.x = element_blank()) + ylim(c(0,6)) + 
  scale_fill_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) + 
  scale_color_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) +
  geom_text(x=1, y=6, label=round(mean(Erbb4_NFOL[Erbb4_NFOL$celltype %in% "WT", ]$gene),3), size = 5, color = "black") +
  geom_text(x=2, y=6, label=round(mean(Erbb4_NFOL[Erbb4_NFOL$celltype %in% "APPPS1", ]$gene),3), size = 5, color = "black") +
  geom_text(x=3, y=6, label=round(mean(Erbb4_NFOL[Erbb4_NFOL$celltype %in% "APPPS1.Il12b.-/-", ]$gene),3), size = 5, color = "black")

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_Violin_Plots_Erbb4_NFOL.pdf",
       plot = g14,
       scale = 1, width = 6, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)


Erbb3_MOL <- Vlnplot_gene_subset(Cell_type = c("MOL"), gene = "Erbb3")

g15 <- ggplot(Erbb3_MOL, aes_string(x="celltype", y="gene", color = "celltype", fill="celltype")) + 
  geom_violin(trim=TRUE, adjust = 3) + theme_cowplot() + ylab("Erbb3 gene expression (MOL)") + 
  theme(legend.position = "None", axis.title.x = element_blank()) + ylim(c(0,6)) + 
  scale_fill_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) + 
  scale_color_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) +
  geom_text(x=1, y=6, label=round(mean(Erbb3_MOL[Erbb3_MOL$celltype %in% "WT", ]$gene),3), size = 5, color = "black") +
  geom_text(x=2, y=6, label=round(mean(Erbb3_MOL[Erbb3_MOL$celltype %in% "APPPS1", ]$gene),3), size = 5, color = "black") +
  geom_text(x=3, y=6, label=round(mean(Erbb3_MOL[Erbb3_MOL$celltype %in% "APPPS1.Il12b.-/-", ]$gene),3), size = 5, color = "black")

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_Violin_Plots_Erbb3_MOL.pdf",
       plot = g15,
       scale = 1, width = 6, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)


Nrg1_DG <- Vlnplot_gene_subset(Cell_type = c("Dentate_Gyrus"), gene = "Nrg1")

g16 <- ggplot(Nrg1_DG, aes_string(x="celltype", y="gene", color = "celltype", fill="celltype")) + 
  geom_violin(trim=TRUE, adjust = 3) + theme_cowplot() + ylab("Nrg1 gene expression (Dentate Gyrus)") + 
  theme(legend.position = "None", axis.title.x = element_blank()) + ylim(c(0,6)) + 
  scale_fill_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) + 
  scale_color_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) +
  geom_text(x=1, y=6, label=round(mean(Nrg1_DG[Nrg1_DG$celltype %in% "WT", ]$gene),3), size = 5, color = "black") +
  geom_text(x=2, y=6, label=round(mean(Nrg1_DG[Nrg1_DG$celltype %in% "APPPS1", ]$gene),3), size = 5, color = "black") +
  geom_text(x=3, y=6, label=round(mean(Nrg1_DG[Nrg1_DG$celltype %in% "APPPS1.Il12b.-/-", ]$gene),3), size = 5, color = "black")

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_Violin_Plots_Nrg1_DG.pdf",
       plot = g16,
       scale = 1, width = 6, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)

Nrg1_CA1 <- Vlnplot_gene_subset(Cell_type = c("CA1"), gene = "Nrg1")

g17 <- ggplot(Nrg1_CA1, aes_string(x="celltype", y="gene", color = "celltype", fill="celltype")) + 
  geom_violin(trim=TRUE, adjust = 3) + theme_cowplot() + ylab("Nrg1 gene expression (CA1)") + 
  theme(legend.position = "None", axis.title.x = element_blank()) + ylim(c(0,6)) + 
  scale_fill_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) + 
  scale_color_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) +
  geom_text(x=1, y=6, label=round(mean(Nrg1_CA1[Nrg1_CA1$celltype %in% "WT", ]$gene),3), size = 5, color = "black") +
  geom_text(x=2, y=6, label=round(mean(Nrg1_CA1[Nrg1_CA1$celltype %in% "APPPS1", ]$gene),3), size = 5, color = "black") +
  geom_text(x=3, y=6, label=round(mean(Nrg1_CA1[Nrg1_CA1$celltype %in% "APPPS1.Il12b.-/-", ]$gene),3), size = 5, color = "black")

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_Violin_Plots_Nrg1_CA1.pdf",
       plot = g17,
       scale = 1, width = 6, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)

Nrg2_DG <- Vlnplot_gene_subset(Cell_type = c("Dentate_Gyrus"), gene = "Nrg2")

g18 <- ggplot(Nrg2_DG, aes_string(x="celltype", y="gene", color = "celltype", fill="celltype")) + 
  geom_violin(trim=TRUE, adjust = 3) + theme_cowplot() + ylab("Nrg2 gene expression (Dentate Gyrus)") + 
  theme(legend.position = "None", axis.title.x = element_blank()) + ylim(c(0,6)) + 
  scale_fill_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) + 
  scale_color_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) +
  geom_text(x=1, y=6, label=round(mean(Nrg2_DG[Nrg2_DG$celltype %in% "WT", ]$gene),3), size = 5, color = "black") +
  geom_text(x=2, y=6, label=round(mean(Nrg2_DG[Nrg2_DG$celltype %in% "APPPS1", ]$gene),3), size = 5, color = "black") +
  geom_text(x=3, y=6, label=round(mean(Nrg2_DG[Nrg2_DG$celltype %in% "APPPS1.Il12b.-/-", ]$gene),3), size = 5, color = "black")

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_Violin_Plots_Nrg2_DG.pdf",
       plot = g18,
       scale = 1, width = 6, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)

Nrg4_CA23 <- Vlnplot_gene_subset(Cell_type = c("CA2_3"), gene = "Nrg4")

g19 <- ggplot(Nrg4_CA23, aes_string(x="celltype", y="gene", color = "celltype", fill="celltype")) + 
  geom_violin(trim=TRUE, adjust = 3) + theme_cowplot() + ylab("Nrg4 gene expression (CA2/3)") + 
  theme(legend.position = "None", axis.title.x = element_blank()) + ylim(c(0,6)) + 
  scale_fill_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) + 
  scale_color_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) +
  geom_text(x=1, y=6, label=round(mean(Nrg4_CA23[Nrg4_CA23$celltype %in% "WT", ]$gene),3), size = 5, color = "black") +
  geom_text(x=2, y=6, label=round(mean(Nrg4_CA23[Nrg4_CA23$celltype %in% "APPPS1", ]$gene),3), size = 5, color = "black") +
  geom_text(x=3, y=6, label=round(mean(Nrg4_CA23[Nrg4_CA23$celltype %in% "APPPS1.Il12b.-/-", ]$gene),3), size = 5, color = "black")

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_Violin_Plots_Nrg4_CA23.pdf",
       plot = g19,
       scale = 1, width = 6, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)

Nrg2_SB <- Vlnplot_gene_subset(Cell_type = c("subiculum"), gene = "Nrg2")

g20 <- ggplot(Nrg2_SB, aes_string(x="celltype", y="gene", color = "celltype", fill="celltype")) + 
  geom_violin(trim=TRUE, adjust = 3) + theme_cowplot() + ylab("Nrg2 gene expression (subiculum)") + 
  theme(legend.position = "None", axis.title.x = element_blank()) + ylim(c(0,6)) + 
  scale_fill_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) + 
  scale_color_manual(labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) +
  geom_text(x=1, y=6, label=round(mean(Nrg2_SB[Nrg2_SB$celltype %in% "WT", ]$gene),3), size = 5, color = "black") +
  geom_text(x=2, y=6, label=round(mean(Nrg2_SB[Nrg2_SB$celltype %in% "APPPS1", ]$gene),3), size = 5, color = "black") +
  geom_text(x=3, y=6, label=round(mean(Nrg2_SB[Nrg2_SB$celltype %in% "APPPS1.Il12b.-/-", ]$gene),3), size = 5, color = "black")

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_Violin_Plots_Nrg2_SB.pdf",
       plot = g20,
       scale = 1, width = 6, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)

################################
# VlnPlot
################################

data9set_cleaned.SO@active.ident <- data9set_cleaned.SO$cell_type

VlnPlot(data9set_cleaned.SO, features = "Nrg1", pt.size = 0) + theme(legend.position = "None")

VlnPlot(data9set_cleaned.SO, features = "Nrg2", pt.size = 0) + theme(legend.position = "None")

VlnPlot(data9set_cleaned.SO, features = "Erbb4", pt.size = 0) + theme(legend.position = "None")

VlnPlot(data9set_cleaned.SO, features = "Erbb3", pt.size = 0) + theme(legend.position = "None")

VlnPlot(data9set_cleaned.SO, features = "Zbtb16", pt.size = 0, group.by = "seurat_clusters") + theme(legend.position = "None")


DE_1_5 <- FindMarkers(data9set_cleaned.SO, ident.1 = "1", ident.2 = "5", group.by = "seurat_clusters", test.use = "MAST", logfc.threshold = 0.1, min.pct = 0.1)
DE_1_5$gene <- rownames(DE_1_5)

FeaturePlot(data9set_cleaned.SO, features = "Robo4")

DE_MOL <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/MOL_AD_ADp40KO.csv")

DE_DG <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/Dentate_Gyrus_AD_ADp40KO.csv")
DE_DG2 <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/Dentate_Gyrus_AD_Ctrl.csv")


DE_SB <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/subiculum_AD_ADp40KO.csv")
DE_SB2 <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/subiculum_AD_Ctrl.csv")

DE_CA1 <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/CA1_AD_ADp40KO.csv")


DE_NFOL <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/NFOL_AD_ADp40KO.csv")

DE_MFOL <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/MFOL_AD_ADp40KO.csv")


data9set_sub.SO <- subset(data9set_cleaned.SO, subset = cell_type %in% "MOL")

VlnPlot(data9set_sub.SO, features = "Erbb4", pt.size = 0, group.by = "gemgroup") + ylim(c(0,6)) + theme(legend.position = "None") +
  geom_text(x=1, y=5, label=round(mean(data9set_sub.df2[data9set_sub.df2$gemgroup %in% "1", ]$Erbb4),3), size = 4, color = "black") +
  geom_text(x=2, y=5, label=round(mean(data9set_sub.df2[data9set_sub.df2$gemgroup %in% "2", ]$Erbb4),3), size = 4, color = "black") +
  geom_text(x=3, y=5, label=round(mean(data9set_sub.df2[data9set_sub.df2$gemgroup %in% "3", ]$Erbb4),3), size = 4, color = "black") +
  geom_text(x=4, y=5, label=round(mean(data9set_sub.df2[data9set_sub.df2$gemgroup %in% "4", ]$Erbb4),3), size = 4, color = "black") +
  geom_text(x=5, y=5, label=round(mean(data9set_sub.df2[data9set_sub.df2$gemgroup %in% "5", ]$Erbb4),3), size = 4, color = "black") +
  geom_text(x=6, y=5, label=round(mean(data9set_sub.df2[data9set_sub.df2$gemgroup %in% "6", ]$Erbb4),3), size = 4, color = "black") +
  geom_text(x=7, y=5, label=round(mean(data9set_sub.df2[data9set_sub.df2$gemgroup %in% "7", ]$Erbb4),3), size = 4, color = "black") +
  geom_text(x=8, y=5, label=round(mean(data9set_sub.df2[data9set_sub.df2$gemgroup %in% "8", ]$Erbb4),3), size = 4, color = "black") +
  geom_text(x=9, y=5, label=round(mean(data9set_sub.df2[data9set_sub.df2$gemgroup %in% "9", ]$Erbb4),3), size = 4, color = "black")


summary(data9set_sub.df2[data9set_sub.df2$gemgroup %in% "6", ]$Erbb4)

table(data9set_sub.df2$gemgroup)

data9set_sub.meta <- data9set_sub.SO@meta.data  

data9set_sub.df <- GetAssayData(data9set_sub.SO, slot = "data") %>% as.data.frame

data9set_sub.df2 <- data9set_sub.df[rownames(data9set_sub.df) %in% c("Erbb4", "Mbp", "Nab2", "Reg2", "Zbtb16", "Eif2b4", "Gpc1", "Clu", "Dicer1"), ] %>% t %>% as.data.frame
data9set_sub.df2 <- cbind(data9set_sub.df2, gemgroup = data9set_sub.SO$gemgroup)


ggplot(data9set_sub.df2, aes(x=Erbb4, y=Mbp)) + geom_point()
cor(data9set_sub.df2$Erbb4, data9set_sub.df2$Mbp)


data9set_sub2.SO <- subset(data9set_cleaned.SO, subset = cell_type %in% "MFOL")
VlnPlot(data9set_sub2.SO, features = "Erbb4", pt.size = 0, group.by = "gemgroup") + theme(legend.position = "None")




