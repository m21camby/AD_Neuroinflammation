#! /bin/env RScript
# written by SJK at 15. May. 2020

.libPaths(c("/home/skim/R/usr_lib", .libPaths()))

library(Seurat)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)
library(tidyr)
library(grid)
library(ggrepel)
source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

avg_lncRNA_top10_z_score.df2 <- readRDS(file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200430_lncRNAs_2nd_analysis_avg_lncRNA_top20_z_score.rda")



g1 <- ggplot(data = avg_lncRNA_top10_z_score.df2, mapping = aes(x = Cell_type, y = gene, fill = z_score)) +
  geom_tile()+  ggtitle("Cell-type specific lncRNAs") + 
  theme_classic() + theme(axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5),
                                        axis.title.x = element_blank(), 
                                        axis.text.x = element_text(angle = 90, vjust = 0, size = 12, color = "black"), 
                                        axis.text.y = element_text(size = 8, color = "black"),
                                        axis.line = element_line(color = "white"),
                                        axis.ticks = element_line(color = "white"),
                                        legend.title = element_blank(),
                                        legend.text = element_text(size = 11, color = "black")) + 
  scale_fill_viridis_c(limits=c(-2, 2), na.value = "#FDE725FF", labels = c("-2", " ", " "," ", "2"))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig7_lncRNA_cell_type_specific_plot.png",
       plot = g1,
       scale = 1, width = 5, height = 7.5, units = "in", device = "png",
       dpi = 300)


##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig7_lncRNA_cell_type_specific_session_info.txt")

