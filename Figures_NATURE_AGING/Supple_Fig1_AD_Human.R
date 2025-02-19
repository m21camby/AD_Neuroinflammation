
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(readr)
library(tidyverse)

Human <- readRDS("/nfs/team292/sk27/tmp/20211020_Seurat_create_from_matrix_seurat_3rd.rds")

Human_meta <- Human@meta.data %>% as.data.frame()

Human$subclass_label <- Human$subclass_label %>% as.factor()
Human@active.ident <- Human$subclass_label

meta.df <- Human@meta.data %>% mutate(broad_cell_type = case_when(subclass_label %in% c("Astro") ~ "Astrocytes",
                                                                     subclass_label %in% c("L2/3 IT", "L5 ET", "L5 IT", "L5/6 NP", "L6 CT", "L6 IT", "L6 IT Car3", "L6b") ~ "Excitatory Neurons",
                                                                     subclass_label %in% c("Lamp5", "Pvalb", "Sncg", "Sst", "Sst Chodl", "Vip") ~ "Inhibitory Neurons",
                                                                     subclass_label %in% c("OPC") ~ "OPC",
                                                                     subclass_label %in% c("Micro-PVM") ~ "Microglia",
                                                                     subclass_label %in% c("Oligo") ~ "Oligodendrocytes",
                                                                     subclass_label %in% c("Endo", "VLMC") ~ "Rest"))

Human$broad_cell_type <- meta.df$broad_cell_type

table(Human$broad_cell_type)

nonzero_counts <- function(gene = "Il12B"){
  sub.SO <- subset(Human, subset = broad_cell_type %in% c("Excitatory Neurons"))
  EX_Il12b <- sub.SO@assays$RNA@data[rownames(sub.SO@assays$RNA@data) %in% gene, ] %>% as.data.frame
  colnames(EX_Il12b) <- gene
  EX_Il12b$celltype <- "Excitatory Neurons"
  
  sub.SO <- subset(Human, subset = broad_cell_type %in% c("Inhibitory Neurons"))
  IN_Il12b <- sub.SO@assays$RNA@data[rownames(sub.SO@assays$RNA@data) %in% gene, ] %>% as.data.frame
  colnames(IN_Il12b) <- gene
  IN_Il12b$celltype <- "Inhibitory Neurons"
  
  sub.SO <- subset(Human, subset = broad_cell_type %in% c("Oligodendrocytes"))
  OL_Il12b <- sub.SO@assays$RNA@data[rownames(sub.SO@assays$RNA@data) %in% gene, ] %>% as.data.frame
  colnames(OL_Il12b) <- gene
  OL_Il12b$celltype <- "Oligodendrocytes"
  
  sub.SO <- subset(Human, subset = broad_cell_type %in% c("OPC"))
  OPC_Il12b <- sub.SO@assays$RNA@data[rownames(sub.SO@assays$RNA@data) %in% gene, ] %>% as.data.frame
  colnames(OPC_Il12b) <- gene
  OPC_Il12b$celltype <- "OPC"
  
  sub.SO <- subset(Human, subset = broad_cell_type %in% c("Microglia"))
  MG_Il12b <- sub.SO@assays$RNA@data[rownames(sub.SO@assays$RNA@data) %in% gene, ] %>% as.data.frame
  colnames(MG_Il12b) <- gene
  MG_Il12b$celltype <- "Microglia"
  
  sub.SO <- subset(Human, subset = broad_cell_type %in% c("Astrocytes"))
  AS_Il12b <- sub.SO@assays$RNA@data[rownames(sub.SO@assays$RNA@data) %in% gene, ] %>% as.data.frame
  colnames(AS_Il12b) <- gene
  AS_Il12b$celltype <- "Astrocytes"
  
  All <- rbind(EX_Il12b, IN_Il12b, OL_Il12b, OPC_Il12b, MG_Il12b, AS_Il12b)
  
  return(All)
  
}

# -----------
# IL12RB1
# -----------
IL12RB1 <- nonzero_counts(gene = "IL12RB1")
temp <- data.frame(IL12RB1 = c(-0.2), celltype = c("OPC"))
IL12RB1 <- rbind(IL12RB1, temp)
IL12RB1$celltype <- factor(IL12RB1$celltype, levels = c("Excitatory Neurons", "Inhibitory Neurons", "Oligodendrocytes", "OPC", "Microglia", "Astrocytes"))

v1 <- ggplot(IL12RB1, aes(x=celltype, y=IL12RB1, color=celltype))  +
  geom_violin(trim=FALSE, width = 1) + scale_color_manual(values=c("#CC0000", "#FF6600", "#99CCCC", "#009999", "#99CC00","#FFCC66")) +
  theme_classic() +
  geom_jitter(position=position_jitter(0.2), size = 3, alpha = 2) + scale_y_continuous(expand = c(0,0), limits = c(-0.1,4)) +  
  ggtitle("IL12RB1 expression") +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15, family = "helvetica", color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, size = 15, family = "helvetica", color = "black"),
        axis.text.y = element_text(size = 15, family = "helvetica", color = "black"),
        plot.title = element_text(family = "helvetica", color = "black", hjust = 0.5)) + expand_limits(y=0)

ggsave(filename = "~/MDC/Fig_Violin_Human_IL12RB1.png",
       plot = v1,
       scale = 1, width = 9, height = 6, units = "in", device = "png",
       dpi = 600)


ggsave(filename = "~/MDC/Final_revisions/Supple_Fig1_A.pdf",
       plot = v1,
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 600)

# -----------
# IL12RB2
# -----------
IL12RB2 <- nonzero_counts(gene = "IL12RB2")
#temp <- data.frame(IL12RB1 = c(-0.2), celltype = c("OPC"))
#IL12RB2 <- rbind(IL12RB2, temp)
IL12RB2$celltype <- factor(IL12RB2$celltype, levels = c("Excitatory Neurons", "Inhibitory Neurons", "Oligodendrocytes", "OPC", "Microglia", "Astrocytes"))

v1 <- ggplot(IL12RB2, aes(x=celltype, y=IL12RB2, color=celltype))  +
  geom_violin(trim=FALSE, width = 1) + scale_color_manual(values=c("#CC0000", "#FF6600", "#99CCCC", "#009999", "#99CC00","#FFCC66")) +
  theme_classic() +
  geom_jitter(position=position_jitter(0.2), size = 3, alpha = 2) + scale_y_continuous(expand = c(0,0), limits = c(-0.1,4)) +  
  ggtitle("IL12RB2 expression") +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15, family = "helvetica", color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, size = 15, family = "helvetica", color = "black"),
        axis.text.y = element_text(size = 15, family = "helvetica", color = "black"),
        plot.title = element_text(family = "helvetica", color = "black", hjust = 0.5)) + expand_limits(y=0)
v1
ggsave(filename = "~/MDC/Fig_Violin_Human_IL12RB2.png",
       plot = v1,
       scale = 1, width = 9, height = 6, units = "in", device = "png",
       dpi = 600)

ggsave(filename = "~/MDC/Final_revisions/Supple_Fig1_B.pdf",
       plot = v1,
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 600)

# -----------
# IL23R
# -----------
IL23R <- nonzero_counts(gene = "IL23R")
temp <- data.frame(IL23R = c(-0.2,-0.2,-0.2), celltype = c("Microglia", "Astrocytes", "OPC"))
IL23R <- rbind(IL23R, temp)
IL23R$celltype <- factor(IL23R$celltype, levels = c("Excitatory Neurons", "Inhibitory Neurons", "Oligodendrocytes", "OPC", "Microglia", "Astrocytes"))

v1 <- ggplot(IL23R, aes(x=celltype, y=IL23R, color=celltype))  +
  geom_violin(trim=FALSE, width = 1) + scale_color_manual(values=c("#CC0000", "#FF6600", "#99CCCC", "#009999", "#99CC00","#FFCC66")) +
  theme_classic() +
  geom_jitter(position=position_jitter(0.2), size = 3, alpha = 2) + scale_y_continuous(expand = c(0,0), limits = c(-0.1,4)) +  
  ggtitle("IL23R expression") +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15, family = "helvetica", color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, size = 15, family = "helvetica", color = "black"),
        axis.text.y = element_text(size = 15, family = "helvetica", color = "black"),
        plot.title = element_text(family = "helvetica", color = "black", hjust = 0.5)) + expand_limits(y=0)
v1
ggsave(filename = "~/MDC/Fig_Violin_Human_IL23R.png",
       plot = v1,
       scale = 1, width = 9, height = 6, units = "in", device = "png",
       dpi = 600)


ggsave(filename = "~/MDC/Final_revisions/Supple_Fig1_C.pdf",
       plot = v1,
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 600)
