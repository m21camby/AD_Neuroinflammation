#! /bin/env RScript
# written by SJK at 8. Mar. 2020

.libPaths(c("/home/skim/R/usr_lib", .libPaths()))


library(Seurat)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)
library(tidyr)
library(ggcorrplot)
library(Hmisc)
library(grid)

source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

merged.df <- read.table(file = "/data/rajewsky/home/skim/Microglia_Heppner/Bulk_transcriptome/R_Scripts/20200312_Bulk_Correlation_bulk_tpm_scRNA_data_cpm.txt")


scatter_corr_plot <- function(merged.df, x_input = x_axis, y_input = y_axis){
  g1 <- ggplot(merged.df, aes_string(x = x_input, y = y_input)) + geom_point()
  g1 <- g1 + theme_classic()+ scale_y_continuous(expand = c(0,0), limits = c(0,14))+ scale_x_continuous(expand = c(0,0), limits = c(0,14)) + theme(axis.text.x = element_text(size =14, color = "black", family = "Helvetica"), axis.title.x = element_blank(), axis.text.y = element_text(size =14, color = "black", family = "Helvetica"), axis.title.y = element_blank()) 

  return(g1)
}

merged.df2 <- data.frame(Ctrl_1_bulk = rowMeans(merged.df[,c(1:3)]), 
                         AD_1_bulk = rowMeans(merged.df[,c(4:6)]),
                         ADp40KO_1_bulk = rowMeans(merged.df[,c(7:9)]),
                         Ctrl_1_sc = rowMeans(merged.df[,c(10:12)]),
                         AD_1_sc = rowMeans(merged.df[,c(13:15)]),
                         ADp40KO_1_sc = rowMeans(merged.df[,c(16:18)]))



# text
t1 <- textGrob("bulk Ctrl", gp=gpar(fontsize=15, col="black", font = "Helvetica"), hjust = 0.2)
t2 <- textGrob(round(cor(merged.df2$Ctrl_1_bulk, merged.df2$AD_1_bulk, method = c("pearson")),2), 
              gp=gpar(fontsize=20, col="black", font = "Helvetica"), hjust = 0)
t3 <- textGrob(round(cor(merged.df2$Ctrl_1_bulk, merged.df2$ADp40KO_1_bulk, method = c("pearson")),2), 
              gp=gpar(fontsize=20, col="black", font = "Helvetica"), hjust = 0)
t4 <- textGrob(round(cor(merged.df2$Ctrl_1_bulk, merged.df2$Ctrl_1_sc, method = c("pearson")),2), 
              gp=gpar(fontsize=20, col="black", font = "Helvetica"), hjust = 0)
t5 <- textGrob(round(cor(merged.df2$Ctrl_1_bulk, merged.df2$AD_1_sc, method = c("pearson")),2), 
              gp=gpar(fontsize=20, col="black", font = "Helvetica"), hjust = 0)
t6 <- textGrob(round(cor(merged.df2$Ctrl_1_bulk, merged.df2$ADp40KO_1_sc, method = c("pearson")),2), 
              gp=gpar(fontsize=20, col="black", font = "Helvetica"), hjust = 0)

t7 <- textGrob("bulk AD", gp=gpar(fontsize=15, col="black", font = "Helvetica"), hjust = 0.2)
t8 <- textGrob(round(cor(merged.df2$AD_1_bulk, merged.df2$ADp40KO_1_bulk, method = c("pearson")),2), 
              gp=gpar(fontsize=20, col="black", font = "Helvetica"), hjust = 0)
t9 <- textGrob(round(cor(merged.df2$AD_1_bulk, merged.df2$Ctrl_1_sc, method = c("pearson")),2), 
              gp=gpar(fontsize=20, col="black", font = "Helvetica"), hjust = 0)
t10 <- textGrob(round(cor(merged.df2$AD_1_bulk, merged.df2$AD_1_sc, method = c("pearson")),2), 
              gp=gpar(fontsize=20, col="black", font = "Helvetica"), hjust = 0)
t11 <- textGrob(round(cor(merged.df2$AD_1_bulk, merged.df2$ADp40KO_1_sc, method = c("pearson")),2), 
              gp=gpar(fontsize=20, col="black", font = "Helvetica"), hjust = 0)

t12 <- textGrob("bulk ADp40KO", gp=gpar(fontsize=15, col="black", font = "Helvetica"), hjust = 0.4)
t13 <- textGrob(round(cor(merged.df2$ADp40KO_1_bulk, merged.df2$Ctrl_1_sc, method = c("pearson")),2), 
                gp=gpar(fontsize=20, col="black", font = "Helvetica"), hjust = 0)
t14 <- textGrob(round(cor(merged.df2$ADp40KO_1_bulk, merged.df2$AD_1_sc, method = c("pearson")),2), 
                gp=gpar(fontsize=20, col="black", font = "Helvetica"), hjust = 0)
t15 <- textGrob(round(cor(merged.df2$ADp40KO_1_bulk, merged.df2$ADp40KO_1_sc, method = c("pearson")),2), 
                gp=gpar(fontsize=20, col="black", font = "Helvetica"), hjust = 0)

t16 <- textGrob("sc Ctrl", gp=gpar(fontsize=15, col="black", font = "Helvetica"), hjust = 0.1)
t17 <- textGrob(round(cor(merged.df2$Ctrl_1_sc, merged.df2$AD_1_sc, method = c("pearson")),2), 
                gp=gpar(fontsize=20, col="black", font = "Helvetica"), hjust = 0)
t18 <- textGrob(round(cor(merged.df2$Ctrl_1_sc, merged.df2$ADp40KO_1_sc, method = c("pearson")),2), 
                gp=gpar(fontsize=20, col="black", font = "Helvetica"), hjust = 0)

t19 <- textGrob("sc AD", gp=gpar(fontsize=15, col="black", font = "Helvetica"), hjust = 0.1)
t20 <- textGrob(round(cor(merged.df2$AD_1_sc, merged.df2$ADp40KO_1_sc, method = c("pearson")),2), 
                gp=gpar(fontsize=20, col="black", font = "Helvetica"), hjust = 0)

t21 <- textGrob("sc ADp40KO", gp=gpar(fontsize=15, col="black", font = "Helvetica"), hjust = 0.5)

# graphs
g1 <- scatter_corr_plot(merged.df2, x = "AD_1_bulk", y = "Ctrl_1_bulk")

g2 <- scatter_corr_plot(merged.df2, x = "ADp40KO_1_bulk", y = "Ctrl_1_bulk")
g3 <- scatter_corr_plot(merged.df2, x = "ADp40KO_1_bulk", y = "AD_1_bulk")

g4 <- scatter_corr_plot(merged.df2, x = "Ctrl_1_sc", y = "Ctrl_1_bulk")
g5 <- scatter_corr_plot(merged.df2, x = "Ctrl_1_sc", y = "AD_1_bulk")
g6 <- scatter_corr_plot(merged.df2, x = "Ctrl_1_sc", y = "ADp40KO_1_bulk")

g7 <- scatter_corr_plot(merged.df2, x = "AD_1_sc", y = "Ctrl_1_bulk")
g8 <- scatter_corr_plot(merged.df2, x = "AD_1_sc", y = "AD_1_bulk")
g9 <- scatter_corr_plot(merged.df2, x = "AD_1_sc", y = "ADp40KO_1_bulk")
g10 <- scatter_corr_plot(merged.df2, x = "AD_1_sc", y = "Ctrl_1_sc")

g11 <- scatter_corr_plot(merged.df2, x = "ADp40KO_1_sc", y = "Ctrl_1_bulk")
g12 <- scatter_corr_plot(merged.df2, x = "ADp40KO_1_sc", y = "AD_1_bulk")
g13 <- scatter_corr_plot(merged.df2, x = "ADp40KO_1_sc", y = "ADp40KO_1_bulk")
g14 <- scatter_corr_plot(merged.df2, x = "ADp40KO_1_sc", y = "Ctrl_1_sc")
g15 <- scatter_corr_plot(merged.df2, x = "ADp40KO_1_sc", y = "AD_1_sc")

g1 <- arrangeGrob(t1, t2, t3, t4, t5, t6,
             g1, t7, t8, t9, t10, t11,
             g2, g3, t12, t13, t14, t15, 
             g4, g5, g6, t16, t17, t18,
             g7, g8, g9, g10, t19, t20,
             g11, g12, g13, g14, g15, t21, ncol = 6)





ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig2_bulk_sc_correlation_plot1.png", 
       plot = g1, 
       scale = 1, width = 10, height = 10, units = "in", device = "png",
       dpi = 300)



res5 <- rcorr(as.matrix(merged.df))

res5.df <- as.data.frame(res5$r)

g2 <- ggcorrplot(res5.df, method = "circle", colors = c("#6D9EC1", "white", "#E46726"), type = "lower")
g2 <- g2 + theme(axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.x = element_text(angle = 90, color="black", family = "Helvetica"),
  axis.text.y = element_text(color="black", family = "Helvetica"),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  legend.title = element_text(size = 15, color="black", family = "Helvetica"), 
  legend.text = element_text(size = 13, color="black", family = "Helvetica"),
  legend.justification = c(1, 0),
  legend.position = c(0.4, 0.8),
  legend.direction = "horizontal") +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                title.position = "top", title.hjust = 0.5)) 


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig2_bulk_sc_correlation_plot2.png",
       plot = g2,
       scale = 1, width = 8, height = 8, units = "in", device = "png",
       dpi = 300)


##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig2_bulk_sc_correlation_session_info.txt")

















