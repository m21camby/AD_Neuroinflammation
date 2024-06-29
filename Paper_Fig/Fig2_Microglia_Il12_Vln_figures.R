#! /bin/env RScript
# written by SJK at 23. Oct. 2020

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

# In this plot, I remove zero counts and only retain expressing cell and show by violin plot

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
                                                                          seurat_clusters %in% c(38) ~ "NFOL",
                                                                          seurat_clusters %in% c(12) ~ "OPC",
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



nonzero_counts <- function(gene = "Il12b"){
  data9set_sub.SO <- subset(data9set_cleaned.SO, subset = cell_type %in% c("Dentate_Gyrus","subiculum","CA1", "CA2_3","Unidentified_Neurons"))
  EX_Il12b <- data9set_sub.SO@assays$RNA@data[rownames(data9set_sub.SO@assays$RNA@data) %in% gene, ] %>% as.data.frame
  colnames(EX_Il12b) <- gene
  EX_Il12b$celltype <- "Excitatory Neurons"
  
  data9set_sub.SO <- subset(data9set_cleaned.SO, subset = cell_type %in% c("Inhibitory_Neurons"))
  IN_Il12b <- data9set_sub.SO@assays$RNA@data[rownames(data9set_sub.SO@assays$RNA@data) %in% gene, ] %>% as.data.frame
  colnames(IN_Il12b) <- gene
  IN_Il12b$celltype <- "Inhibitory Neurons"
  
  data9set_sub.SO <- subset(data9set_cleaned.SO, subset = cell_type %in% c("MOL", "MFOL", "NFOL"))
  OL_Il12b <- data9set_sub.SO@assays$RNA@data[rownames(data9set_sub.SO@assays$RNA@data) %in% gene, ] %>% as.data.frame
  colnames(OL_Il12b) <- gene
  OL_Il12b$celltype <- "Oligodendrocytes"
  
  data9set_sub.SO <- subset(data9set_cleaned.SO, subset = cell_type %in% c("OPC"))
  OPC_Il12b <- data9set_sub.SO@assays$RNA@data[rownames(data9set_sub.SO@assays$RNA@data) %in% gene, ] %>% as.data.frame
  colnames(OPC_Il12b) <- gene
  OPC_Il12b$celltype <- "OPC"
  
  data9set_sub.SO <- subset(data9set_cleaned.SO, subset = cell_type %in% c("Microglia"))
  MG_Il12b <- data9set_sub.SO@assays$RNA@data[rownames(data9set_sub.SO@assays$RNA@data) %in% gene, ] %>% as.data.frame
  colnames(MG_Il12b) <- gene
  MG_Il12b$celltype <- "Microglia"
  
  data9set_sub.SO <- subset(data9set_cleaned.SO, subset = cell_type %in% c("Astrocytes"))
  AS_Il12b <- data9set_sub.SO@assays$RNA@data[rownames(data9set_sub.SO@assays$RNA@data) %in% gene, ] %>% as.data.frame
  colnames(AS_Il12b) <- gene
  AS_Il12b$celltype <- "Astrocytes"
  
  data9set_sub.SO <- subset(data9set_cleaned.SO, subset = cell_type %in% c("Vascular", "VLMC", "Choroid", "Fibroblast", "Cajal", "Pericyte", "Macrophage"))
  R_Il12b <- data9set_sub.SO@assays$RNA@data[rownames(data9set_sub.SO@assays$RNA@data) %in% gene, ] %>% as.data.frame
  colnames(R_Il12b) <- gene
  R_Il12b$celltype <- "Rest"
  
  
  All <- rbind(EX_Il12b, IN_Il12b, OL_Il12b, OPC_Il12b, MG_Il12b, AS_Il12b, R_Il12b)
  
  return(All)
  
}

#######################
# Il12b
#######################


Il12b <- nonzero_counts(gene = "Il12b")

Il12b2 <- Il12b[Il12b$Il12b != 0, ]
Il12b2$celltype <- factor(Il12b2$celltype, levels = c("Excitatory Neurons", "Inhibitory Neurons", "Oligodendrocytes", "OPC", "Microglia", "Astrocytes", "Rest"))

v1 <- ggplot(Il12b2, aes(x=celltype, y=Il12b, color=celltype))  + 
  geom_violin(trim=FALSE, width = 1) + scale_color_manual(values=c("#CC0000", "#FF6600", "#99CCCC", "#009999", "#99CC00","#FFCC66", "purple")) + 
  theme_classic() + 
  geom_jitter(position=position_jitter(0.2), size = 1, alpha = 0.3) + scale_y_continuous(expand = c(0,0)) +  
  ggtitle("Il12b expression excluding non-expressing cells") + 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12, family = "helvetica", color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, size = 12, family = "helvetica", color = "black"),
        axis.text.y = element_text(size = 12, family = "helvetica", color = "black"),
        plot.title = element_text(family = "helvetica", color = "black", hjust = 0.5))
v1
ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_Il12_Vln_figures_Il12b.png", 
       plot = v1, 
       scale = 1, width = 9, height = 6, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_Il12_Vln_figures_Il12b.pdf", 
       plot = v1, 
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)

#######################
# Il12a
#######################
Il12a <- nonzero_counts(gene = "Il12a")

Il12a2 <- Il12a[Il12a$Il12a != 0, ]
temp <- data.frame(Il12a = c(-0.2, -0.2, -0.2), celltype = c("Astrocytes", "Microglia", "OPC"))
Il12a2 <- rbind(Il12a2, temp)

Il12a2$celltype <- factor(Il12a2$celltype, levels = c("Excitatory Neurons", "Inhibitory Neurons", "Oligodendrocytes", "OPC", "Microglia", "Astrocytes","Rest"))

v1 <- ggplot(Il12a2, aes(x=celltype, y=Il12a, color=celltype)) + 
  geom_violin(trim=FALSE, width = 1) + 
  geom_jitter(position=position_jitter(0.2), size = 1, alpha = 0.3) + scale_y_continuous(expand = c(0,0), limits = c(0,4)) +  
  scale_color_manual(values=c("#FFCC66", "#CC0000", "#FF6600", "#99CC00", "#99CCCC", "#009999", "purple")) + theme_classic() + 
  ggtitle("Il12a expression excluding non-expressing cells") + 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12, family = "helvetica", color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, size = 12, family = "helvetica", color = "black"),
        axis.text.y = element_text(size = 12, family = "helvetica", color = "black"),
        plot.title = element_text(family = "helvetica", color = "black", hjust = 0.5))
v1
ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_Il12_Vln_figures_Il12a_New_scale.pdf", 
       plot = v1, 
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)

#######################
# Il23a
#######################

Il23a <- nonzero_counts(gene = "Il23a")

Il23a2 <- Il23a[Il23a$Il23a != 0, ]
temp <- data.frame(Il23a = c(-0.2, -0.2, -0.2), celltype = c("Astrocytes", "Oligodendrocytes", "OPC"))
Il23a2 <- rbind(Il23a2, temp)

Il23a2$celltype <- factor(Il23a2$celltype, levels = c("Excitatory Neurons", "Inhibitory Neurons", "Oligodendrocytes", "OPC", "Microglia", "Astrocytes","Rest"))

v1 <- ggplot(Il23a2, aes(x=celltype, y=Il23a, color=celltype)) + 
  geom_violin(trim=FALSE, width = 1) + 
  geom_jitter(position=position_jitter(0.2), size = 1, alpha = 0.3) + 
  scale_color_manual(values=c("#FFCC66", "#CC0000", "#FF6600", "#99CC00", "#99CCCC", "#009999", "purple")) + theme_classic() + 
  scale_y_continuous(expand = c(0,0), limits = c(0,4)) +  
  ggtitle("Il23a expression excluding non-expressing cells") + 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12, family = "helvetica", color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, size = 12, family = "helvetica", color = "black"),
        axis.text.y = element_text(size = 12, family = "helvetica", color = "black"),
        plot.title = element_text(family = "helvetica", color = "black", hjust = 0.5))
v1
ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_Il12_Vln_figures_Il23a.png", 
       plot = v1, 
       scale = 1, width = 9, height = 6, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_Il12_Vln_figures_Il23a_New_Scale.pdf", 
       plot = v1, 
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)

#######################
# Snap25
#######################

Snap25 <- nonzero_counts(gene = "Snap25")
Snap25 <- Snap25[Snap25$Snap25 != 0, ]
Snap25$celltype <- factor(Snap25$celltype, levels = c("Excitatory Neurons", "Inhibitory Neurons", "Oligodendrocytes", "OPC", "Microglia", "Astrocytes", "Rest"))


v1 <- ggplot(Snap25, aes(x=celltype, y=Snap25, color=celltype)) + 
  geom_violin(trim=FALSE, width = 1) + 
  geom_jitter(position=position_jitter(0.2), size = 1, alpha = 0.3) +  
  scale_color_manual(values=c("#CC0000", "#FF6600", "#99CCCC", "#009999", "#99CC00","#FFCC66", "purple")) + 
  theme_classic() + scale_y_continuous(expand = c(0,0)) +  
  ggtitle("Snap25 expression excluding non-expressing cells") + 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12, family = "helvetica", color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, size = 12, family = "helvetica", color = "black"),
        axis.text.y = element_text(size = 12, family = "helvetica", color = "black"),
        plot.title = element_text(family = "helvetica", color = "black", hjust = 0.5))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_Il12_Vln_figures_Snap25.pdf", 
       plot = v1, 
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)

#######################
# Mbp
#######################

Mbp <- nonzero_counts(gene = "Mbp")
Mbp <- Mbp[Mbp$Mbp != 0, ]
Mbp$celltype <- factor(Mbp$celltype, levels = c("Excitatory Neurons", "Inhibitory Neurons", "Oligodendrocytes", "OPC", "Microglia", "Astrocytes", "Rest"))


v1 <- ggplot(Mbp, aes(x=celltype, y=Mbp, color=celltype)) + 
  geom_violin(trim=FALSE, width = 1) + 
  geom_jitter(position=position_jitter(0.2), size = 1, alpha = 0.3) +  
  scale_color_manual(values=c("#CC0000", "#FF6600", "#99CCCC", "#009999", "#99CC00","#FFCC66", "purple")) + 
  theme_classic() + scale_y_continuous(expand = c(0,0)) +  
  ggtitle("Mbp expression excluding non-expressing cells") + 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12, family = "helvetica", color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, size = 12, family = "helvetica", color = "black"),
        axis.text.y = element_text(size = 12, family = "helvetica", color = "black"),
        plot.title = element_text(family = "helvetica", color = "black", hjust = 0.5))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_Il12_Vln_figures_Mbp.pdf", 
       plot = v1, 
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)



#######################
# Il12rb1
#######################

Il12rb1 <- nonzero_counts(gene = "Il12rb1")
Il12rb1 <- Il12rb1[Il12rb1$Il12rb1 != 0, ]
table(Il12rb1$celltype)
temp <- data.frame(Il12rb1 = c(-0.2), celltype = c("OPC"))
Il12rb1 <- rbind(Il12rb1, temp)


Il12rb1$celltype <- factor(Il12rb1$celltype, levels = c("Excitatory Neurons", "Inhibitory Neurons", "Oligodendrocytes", "OPC", "Microglia", "Astrocytes", "Rest"))

v1 <- ggplot(Il12rb1, aes(x=celltype, y=Il12rb1, color=celltype)) + 
  geom_violin(trim=FALSE, width = 1) + 
  geom_jitter(position=position_jitter(0.2), size = 1, alpha = 0.3) +  
  scale_color_manual(breaks = c("Excitatory Neurons", "Inhibitory Neurons", "Oligodendrocytes", "OPC", "Microglia", "Astrocytes", "Rest"),
                     values=c("#CC0000", "#FF6600", "#99CCCC", "#009999", "#99CC00","#FFCC66", "purple")) +
  theme_classic() + scale_y_continuous(expand = c(0,0), limits = c(0,4)) +  
  ggtitle("Il12rb1 expression excluding non-expressing cells") + 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12, family = "helvetica", color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, size = 12, family = "helvetica", color = "black"),
        axis.text.y = element_text(size = 12, family = "helvetica", color = "black"),
        plot.title = element_text(family = "helvetica", color = "black", hjust = 0.5))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_Il12_Vln_figures_Il12rb1.pdf", 
       plot = v1, 
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)

#######################
# Il12rb2
#######################

Il12rb2 <- nonzero_counts(gene = "Il12rb2")
Il12rb2 <- Il12rb2[Il12rb2$Il12rb2 != 0, ]
Il12rb2$celltype <- factor(Il12rb2$celltype, levels = c("Excitatory Neurons", "Inhibitory Neurons", "Oligodendrocytes", "OPC", "Microglia", "Astrocytes", "Rest"))


v1 <- ggplot(Il12rb2, aes(x=celltype, y=Il12rb2, color=celltype)) + 
  geom_violin(trim=FALSE, width = 1) + 
  geom_jitter(position=position_jitter(0.2), size = 1, alpha = 0.3) +  
  scale_color_manual(values=c("#CC0000", "#FF6600", "#99CCCC", "#009999", "#99CC00","#FFCC66", "purple")) + 
  theme_classic() + scale_y_continuous(expand = c(0,0), limits = c(0,4)) +  
  ggtitle("Il12rb2 expression excluding non-expressing cells") + 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12, family = "helvetica", color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, size = 12, family = "helvetica", color = "black"),
        axis.text.y = element_text(size = 12, family = "helvetica", color = "black"),
        plot.title = element_text(family = "helvetica", color = "black", hjust = 0.5))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_Il12_Vln_figures_Il12rb2.pdf", 
       plot = v1, 
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)


#######################
# Il23r
#######################

Il23r <- nonzero_counts(gene = "Il23r")
Il23r <- Il23r[Il23r$Il23r != 0, ]
temp <- data.frame(Il23r = c(-0.2, -0.2), celltype = c("Microglia", "OPC"))
Il23r <- rbind(Il23r, temp)
Il23r$celltype <- factor(Il23r$celltype, levels = c("Excitatory Neurons", "Inhibitory Neurons", "Oligodendrocytes", "OPC", "Microglia", "Astrocytes", "Rest"))


v1 <- ggplot(Il23r, aes(x=celltype, y=Il23r, color=celltype)) + 
  geom_violin(trim=FALSE, width = 1) + 
  geom_jitter(position=position_jitter(0.2), size = 1, alpha = 0.3) + 
  scale_color_manual(breaks = c("Excitatory Neurons", "Inhibitory Neurons", "Oligodendrocytes", "OPC", "Microglia", "Astrocytes", "Rest"),
                     values=c("#CC0000", "#FF6600", "#99CCCC", "#009999", "#99CC00","#FFCC66", "purple")) +
  theme_classic() +  scale_y_continuous(expand = c(0,0), limits = c(0,4)) +  
  ggtitle("Il23r expression excluding non-expressing cells") + 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12, family = "helvetica", color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, size = 12, family = "helvetica", color = "black"),
        axis.text.y = element_text(size = 12, family = "helvetica", color = "black"),
        plot.title = element_text(family = "helvetica", color = "black", hjust = 0.5))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_Il12_Vln_figures_Il23r.pdf", 
       plot = v1, 
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_Il12_Vln_figures_session_info.txt")