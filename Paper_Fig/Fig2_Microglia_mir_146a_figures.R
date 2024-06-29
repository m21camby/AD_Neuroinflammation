#! /bin/env RScript
# written by SJK at 29. 08. 2020
# Creating z-score of mir-146a in Microglia
# File name: Fig2_Microglia_mir_146a_figures.R


##################
# library loading
##################


library(Seurat)
library(tidyverse)
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(scales)

source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

load(file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/miRNA_lncRNA/R_Scripts/20200505_pri_miRNA_dgCMatrix_Create_Seurat_with_intergenic_pri_miRNA.obj")
data9set.SO$cell_type <- data9set.SO@active.ident

# checking interesting pri-miRNA RNAs
MG_DE_AD_Ctrl <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/miRNA_lncRNA/R_Scripts/DE_csv_files/20200513_pri_miRNA_Cell_type_DE_by_MAST_MicrogliaAD_Ctrl.csv")
MG_DE_ADp40KO_Ctrl <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/miRNA_lncRNA/R_Scripts/DE_csv_files/20200513_pri_miRNA_Cell_type_DE_by_MAST_MicrogliaADp40KO_Ctrl.csv")
MG_DE_AD_ADp40KO <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/miRNA_lncRNA/R_Scripts/DE_csv_files/20200513_pri_miRNA_Cell_type_DE_by_MAST_MicrogliaAD_ADp40KO.csv")

data9set_MG.SO <- subset(data9set.SO, subset = cell_type %in% "Microglia")

# z-score calculation
z_score_DE_genes <- function(Seurat_object, list_genes){
  for(i in c("Ctrl", "AD", "ADp40KO")){
    if(i == "Ctrl"){
      cell_avg.df <- as.data.frame(Matrix::rowMeans(GetAssayData(object = subset(Seurat_object, subset = sample %in% i), slot = 'data')))
    }
    else{
      tempcell_avg.df <- as.data.frame(Matrix::rowMeans(GetAssayData(object = subset(Seurat_object, subset = sample %in% i), slot = 'data')))
      cell_avg.df <- cbind(cell_avg.df, tempcell_avg.df)
    }
  }
  #print("1st")
  cell_avg.df$gene <- rownames(cell_avg.df)
  #print("2nd")
  cell_avg.df <- cell_avg.df[cell_avg.df$gene %in% list_genes, ]
  #print("3rd")
  cell_avg.df$gene <- NULL
  
  # calculate z-score (high)
  cell_avg_z_score.df <- as.data.frame(t(apply(cell_avg.df, 1, function(x) (x - mean(x)) / sd(x))))
  return(cell_avg_z_score.df)
} 

list_genes <- c("mir-146a", "mir-466i", "Mir6236")
MG_z_score <- z_score_DE_genes(data9set_MG.SO, list_genes)
colnames(MG_z_score) <- c("WT", "APPPS1", "APPPS1.Il12-/-")

MG_z_score$gene <- rownames(MG_z_score)

MG_z_score2 <- gather(MG_z_score, sample, z_score, "WT":"APPPS1.il12-/-")

MG_z_score2$sample <- factor(MG_z_score2$sample, levels = rev(c("WT", "APPPS1", "APPPS1.il12-/-")))

g1 <- ggplot(data = MG_z_score2, mapping = aes(x = sample, y = gene, fill = z_score)) +
  ggtitle("pri-miRNA in Microglia") + 
  geom_tile()+  theme_classic() + theme(plot.title = element_text(size = 15, hjust = 0.5), 
                                        axis.title.y = element_blank(),
                                        axis.title.x = element_blank(),
                                        axis.text.x = element_text(size = 12, color = "black"),
                                        axis.text.y = element_text(size = 12, color = "black"),
                                        axis.line = element_line(color = "white"),
                                        axis.ticks = element_line(color = "white"),
                                        legend.title = element_blank(),
                                        legend.text = element_text(size = 11, color = "black")) +
  scale_fill_viridis_c(limits=c(-1, 1), oob=squish) + coord_flip() 



ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_mir_146a_figures_z_score.png", 
       plot = g1, 
       scale = 1, width = 6, height = 2.5, units = "in", device = "png",
       dpi = 300)

# mir-146a only

list_genes <- c("mir-146a")
MG_z_score <- z_score_DE_genes(data9set_MG.SO, list_genes)
colnames(MG_z_score) <- c("WT", "APPPS1", "APPPS1.Il12-/-")

MG_z_score$gene <- rownames(MG_z_score)

MG_z_score2 <- gather(MG_z_score, sample, z_score, "WT":"APPPS1.Il12-/-")

MG_z_score2$sample <- factor(MG_z_score2$sample, levels = c("WT", "APPPS1", "APPPS1.Il12-/-"))

g2 <- ggplot(data = MG_z_score2, mapping = aes(x = sample, y = gene, fill = z_score)) +
  ggtitle("pri-miR-146a in microglia") + 
  geom_tile()+  theme_classic() + theme(plot.title = element_text(size = 15, hjust = 0.5), 
                                        axis.title.y = element_blank(),
                                        axis.title.x = element_blank(),
                                        axis.text.x = element_text(angle = 60, vjust = 0.5, size = 12, color = "black"),
                                        axis.text.y = element_text(size = 12, color = "black"),
                                        axis.line = element_line(color = "white"),
                                        axis.ticks = element_line(color = "white"),
                                        legend.title = element_blank(),
                                        legend.text = element_text(size = 11, color = "black")) +
  scale_fill_viridis_c(limits=c(-1.5, 1.5), oob=squish)
g2
ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_mir_146a_figures_z_score_mir_146a.png", 
       plot = g2, 
       scale = 1, width = 6, height = 2, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_mir_146a_figures_z_score_mir_146a.pdf", 
       plot = g2, 
       scale = 1, width = 6, height = 2.2, units = "in", device = cairo_pdf,
       dpi = 300)

f1 <- FeaturePlot(data9set.SO, features = "mir-146a", order = TRUE, pt.size = 1) + 
  ggtitle("pri-mir-146a") + 
  theme(axis.line = element_blank(), 
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_mir_146a_figures_feature_plot.png", 
       plot = f1, 
       scale = 1, width = 7, height = 6, units = "in", device = "png",
       dpi = 300)

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_mir_146a_figures_session_info.txt")

