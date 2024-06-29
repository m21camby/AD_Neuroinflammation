#! /bin/env RScript
# written by SJK at 9. Dec. 2020

library(cowplot)
library(viridis)
library(gridExtra)
library(dplyr)
library(ggplot2)
library(tidyr)
library(Seurat)
library(ggrepel)
library(readxl)
library(tidyr)
library(Cairo)


source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")


##########################
# pie chart of biotype
##########################
# legend: Biotypes from 337 pri-miRNAs expressed in snRNA-seq

SP064_pri_sub.sum.df <- readRDS("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/miRNA_lncRNA/R_Scripts/20200505_pri_miRNA_dgCMatrix_Create_SP064_pri_sub.sum.rda")


g1 <- ggplot(SP064_pri_sub.sum.df, aes(x = 2, y= percent, fill = factor(percent, levels = c(0.567, 0.246, 0.187)))) +
  geom_bar(width = 1, stat = "identity") + 
  coord_polar(theta = "y", start=0) +
  theme_void() + theme(legend.text = element_text(size = 12), legend.title = element_blank()) +
  scale_fill_manual(values=c("#CC0000", "#FFCC66", "#006600"),
                    labels = c("Protein coding", "intergenic", "lncRNA")) +
  geom_text(aes(y = 0.7, label = "56.7%"), color = "white", size = 6) + 
  geom_text(aes(y = 0.32, label = "24.6%"), color = "white", size = 6) + 
  geom_text(aes(y = 0.1, label = "18.7%"), color = "white", size = 6)

#####################################
# barplot of read mapped to biotype
#####################################

SP064_pri_sub.sum.df2 <- readRDS("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/miRNA_lncRNA/R_Scripts/20200505_pri_miRNA_dgCMatrix_Create_SP064_pri_sub.sum_2nd.rda")

g2 <- ggplot(SP064_pri_sub.sum.df2, aes(fill=biotype, y=percent, x="gene_biotype")) + 
  geom_bar(position="fill", stat="identity") + 
  theme_classic() + 
  geom_text(aes(y=0.45, label="lncRNA: 69.2%"), vjust=1.6, color="white", size=4) + 
  geom_text(aes(y=0.90, label="protein coding:"), vjust=1.6, color="white", size=4) + 
  geom_text(aes(y=0.81, label="26.6%"), vjust=1.6, color="white", size=4) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,1.2), breaks = c(0,1)) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.title = element_blank(), 
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_blank()) + 
  annotate(geom="text", x=1, y=1.08, label="intergenic: 4.2%", color="black", size = 4) + 
  scale_fill_manual(values=c("#FFCC66", "#CC0000", "#006600")) +
  geom_segment(aes(x = 1, y = 0.999, xend = 1, yend = 1.05),color='black',size=0.5) 

blankPlot <- ggplot()+geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

g3 <- arrangeGrob(g1, blankPlot, g2, blankPlot, ncol = 2, widths = c(3, 1), heights = c(0.5, 2.5, 0.7), layout_matrix = cbind(c(1,1,1), c(2,3,4)))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig10_pri_miRNA_plots.pdf",
       plot = g3,
       scale = 1, width = 8, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)


################################
# mir146a
################################
# Seurat loading

load(file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/miRNA_lncRNA/R_Scripts/20200505_pri_miRNA_dgCMatrix_Create_Seurat_with_intergenic_pri_miRNA.obj")

f1 <- FeaturePlot(data9set.SO, features = "mir-146a", pt.size = 1, order = TRUE) + ggtitle("pri-mir-146a")

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig10_pri_miRNA_plots_mir-146a_feature_plot.pdf",
       plot = f1,
       scale = 1, width = 7, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)


data9set_MG.SO <- subset(data9set.SO, subset = seurat_clusters %in% c(3,8))

data9set_MG.meta <- data9set_MG.SO@meta.data

MG_mir_146a <- data9set_MG.SO@assays$RNA@data[rownames(data9set_MG.SO@assays$RNA@data) %in% "mir-146a", ] %>% as.data.frame
colnames(MG_mir_146a) <- "mir_146a"

MG_mir_146a <- cbind(MG_mir_146a, sample = data9set_MG.meta$sample)

mean(MG_mir_146a[MG_mir_146a$sample %in% "Ctrl", ]$mir_146a)
mean(MG_mir_146a[MG_mir_146a$sample %in% "AD", ]$mir_146a)
mean(MG_mir_146a[MG_mir_146a$sample %in% "ADp40KO", ]$mir_146a)

# remove zero counts but not used in here
MG_mir_146a2 <- MG_mir_146a[MG_mir_146a$mir_146a != 0, ]

v1 <- ggplot(MG_mir_146a, aes(x=sample, y=mir_146a, color=sample)) + coord_cartesian(ylim =  c(1,5)) + 
  theme_cowplot() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15, family = "helvetica"),
        axis.text = element_text(size = 15, family = "helvetica"),
        legend.position = "none") + 
  ylab("normalized expression") + 
  geom_violin(trim=FALSE, width = 1) + geom_jitter(position=position_jitter(0.2), size = 1, alpha = 0.3)  + 
  geom_text(x = 1, y = 5, label = paste0("avg: ", round(mean(MG_mir_146a[MG_mir_146a$sample %in% "Ctrl", ]$mir_146a), 2)), color = "black", family = "helvetica", size = 5) + 
  geom_text(x = 2, y = 5, label = paste0("avg: ", round(mean(MG_mir_146a[MG_mir_146a$sample %in% "AD", ]$mir_146a), 2)), color = "black", family = "helvetica", size = 5) +
  geom_text(x = 3, y = 5, label = paste0("avg: ", round(mean(MG_mir_146a[MG_mir_146a$sample %in% "ADp40KO", ]$mir_146a), 2)), color = "black", family = "helvetica", size = 5) +
  scale_color_manual(values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) + 
  scale_x_discrete(labels=c("Ctrl" = "WT", "AD" = "APPPS1", "ADp40KO" = "APPPS1.il12b"))


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig10_pri_miRNA_plots_mir-146a_violin_plot.pdf", 
       plot = v1, 
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)


# cell type specific
All_mir_146a <- data9set.SO@assays$RNA@data[rownames(data9set.SO@assays$RNA@data) %in% "mir-146a", ] %>% as.data.frame
colnames(All_mir_146a) <- "mir_146a"

data9set.meta <- data9set.SO@meta.data

data9set_cleaned.meta <- data9set_cleaned.SO@meta.data


All_mir_146a_merged <- cbind(All_mir_146a, data9set_cleaned.meta)

ggplot(All_mir_146a_merged, aes(x=cell_type, y=mir_146a, color=cell_type)) + 
  geom_violin(trim=FALSE, width = 1) + geom_jitter(position=position_jitter(0.2), size = 1, alpha = 0.3) + 
  coord_cartesian(ylim =  c(1,5)) + 
  theme_cowplot()

################################
# mir7-2
################################

f2 <- FeaturePlot(data9set.SO, features = "mir-7a-2", pt.size = 1, order = TRUE) + ggtitle("pri-mir-7a-2")

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig10_pri_miRNA_plots_mir-7a-2_feature_plot.pdf",
       plot = f2,
       scale = 1, width = 7, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)


# cell type specific
All_mir_7a_2 <- data9set.SO@assays$RNA@data[rownames(data9set.SO@assays$RNA@data) %in% "mir-7a-2", ] %>% as.data.frame
colnames(All_mir_7a_2) <- "mir_7a_2"

All_mir_7a_2_merged <- cbind(All_mir_7a_2, data9set_cleaned.meta)


ggplot(All_mir_7a_2_merged, aes(x=cell_type, y=mir_7a_2, color=cell_type)) + 
  geom_violin(trim=FALSE, width = 1) + geom_jitter(position=position_jitter(0.2), size = 1, alpha = 0.3) + 
  coord_cartesian(ylim =  c(0.5,3.5)) + 
  theme_cowplot()

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig10_pri_miRNA_plots_session_info.txt")
