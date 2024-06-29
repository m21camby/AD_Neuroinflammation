#! /bin/env RScript
# written by SJK at 20. 08. 2020
# Creating Bulk deconvolution results by barplot

##################
# library loading
##################

.libPaths(c("/data/rajewsky/home/skim/R/usr_lib_Seurat/", "/data/rajewsky/home/skim/R/usr_lib_v4/"))
#.libPaths()
unloadNamespace("mgcv")
unloadNamespace("Matrix")

library(gridExtra)
library(dplyr)
library(ggplot2)
library(tidyr)
library(viridis)
library(readr)
library(tibble)
library(ggrepel)
library(forcats)
library(ggpubr)
library(viridis)
library(multcomp)

source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")
###############################
# decovolution exp comparison
###############################

deconv <- readRDS(file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200818_deconvolution_MuSic_reports/deconv.rda")

deconv$Cell_type <- factor(deconv$Cell_type, levels = c("Neurons", "Astrocytes", "Microglia","Oligo","OPC","Vascular", "Rest"))

g1 <- ggplot(deconv, aes(x=Cell_type, y=proportion, fill = exp)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) + 
  theme_classic() + 
  ylab("proportion (%)") + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_viridis_d() + 
  geom_errorbar(aes(ymin=proportion-sd, ymax=proportion+sd), width=.2,
                position=position_dodge(.9)) + ggtitle("deconvolution of bulk transcriptomes")  + 
  theme(axis.text.x = element_text(size =14, color = "black"), 
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size =14, color = "black"),
        axis.title.y = element_text(size = 16), 
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = c(0.85, 0.85), 
        legend.title = element_blank(),
        legend.text = element_text(size =14, color = "black"))


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig2_bulk_deconvolution_bar_plot_bulk_deconvolution_figure.png", 
       plot = g1, 
       scale = 1, width = 8.8, height = 5.5, units = "in", device = "png",
       dpi = 300)

###########
# box plot
###########
deconv_raw <- readRDS(file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200818_deconvolution_MuSic_reports/deconv_raw.rda")

deconv_raw$Cell_type <- factor(deconv_raw$Cell_type, levels = c("Neurons", "Astrocytes", "Microglia", "Oligo", "OPC", "Vascular", "Rest"))

deconv_raw$genotype <- factor(deconv_raw$genotype, levels = c("WT", "APPPS1", "APPPS1.il12b-/-"))

g3 <- ggplot(deconv_raw, aes(x=Cell_type, y=proportion, fill = genotype, color=genotype)) +
  geom_boxplot(fill="white") +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize= 0.5) +
  scale_color_manual(values = c("#440154FF", viridis(3)[2], "#CC9900"), labels = c("WT", "APPPS1", bquote(APPPS1.il12b^"-/-"))) +
  scale_fill_manual(values = c("#440154FF", viridis(3)[2], "#CC9900"), labels = c("WT", "APPPS1", bquote(APPPS1.il12b^"-/-"))) +
  theme_classic() +
  ggtitle("deconvolution of bulk transcriptomes") + 
  ylab("percent") +
  theme(axis.text.x = element_text(size =15, angle = 60, color = "black", vjust = 0.5, family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size =15, color = "black", family = "Helvetica"),
        axis.title.y = element_text(size =15, family = "Helvetica"),
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 15, family = "Helvetica"),
        legend.key.size = unit(1, "cm")) +
  ylim(0, 0.6) +
  stat_compare_means(aes(group = genotype), label = "p.signif", method = "anova", show.legend = FALSE, size = 6, family = "Helvetica") +
  guides(color = FALSE)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig2_bulk_deconvolution_bar_plot_bulk_deconvolution_figure_box_plot.png", 
       plot = g3, 
       scale = 1, width = 9, height = 6, units = "in", device = "png",
       dpi = 300)

g3 <- ggplot(deconv_raw, aes(x=Cell_type, y=proportion, fill = genotype, color=genotype)) +
  geom_boxplot(fill="white") +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize= 0.5) +
  scale_color_manual(values = c("#440154FF", viridis(3)[2], "#CC9900"), labels = c("WT", "APPPS1", bquote(APPPS1.il12b^"-/-"))) +
  scale_fill_manual(values = c("#440154FF", viridis(3)[2], "#CC9900"), labels = c("WT", "APPPS1", bquote(APPPS1.il12b^"-/-"))) +
  theme_classic() +
  ggtitle("deconvolution of bulk transcriptomes") + 
  ylab("percent") +
  theme(axis.text.x = element_text(size =15, angle = 60, color = "black", vjust = 0.5, family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size =15, color = "black", family = "Helvetica"),
        axis.title.y = element_text(size =15, family = "Helvetica"),
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 15, family = "Helvetica"),
        legend.key.size = unit(1, "cm")) +
  ylim(0, 0.6) +
  stat_compare_means(aes(group = genotype), label = "p.format", method = "anova", show.legend = FALSE, size = 4, family = "Helvetica") +
  guides(color = FALSE)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig2_bulk_deconvolution_bar_plot_bulk_deconvolution_figure_box_plot_p_format.png", 
       plot = g3, 
       scale = 1, width = 9, height = 6, units = "in", device = "png",
       dpi = 300)


######################
# F-value calculation
######################

df <- data.frame(celltype = character(),
                 Df = numeric(),
                 sum_sq = numeric(),
                 mean_sq = numeric(),
                 f_value = numeric(),
                 p_value = numeric())

celltypes <- deconv_raw$Cell_type %>% unique()

for(i in celltypes){
  print(i)
  res.aov <- aov(proportion ~ genotype, data = deconv_raw[deconv_raw$Cell_type %in% i, ])
  # Summary of the analysis
  res.df <- summary(res.aov)
  
  tmp <- data.frame(celltype = i,
                    Df = res.df[[1]]$Df[1],
                    sum_sq = res.df[[1]]$`Sum Sq`[1],
                    mean_sq = res.df[[1]]$`Mean Sq`[1],
                    f_value = res.df[[1]]$`F value`[1],
                    p_value = res.df[[1]]$`Pr(>F)`[1])
  
  res.df[[1]]$`F value`[1]
  res.df[[1]]$Df[1]
  res.df[[1]]$`Pr(>F)`[1]
  res.df[[1]]$`Sum Sq`[1]
  res.df[[1]]$`Mean Sq`[1]
  
  df <- rbind(df, tmp)
  
}

write.csv(df, file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig2_bulk_deconvolution_bar_plot_bulk_deconvolution_figure_box_plot_raw_data.csv")


###########################
# multiple comparison test
###########################

glht(res.aov, linfct = mcp(genotype = "Tukey"))

test <- deconv_raw[deconv_raw$Cell_type %in% celltypes[1], ]

model1 <- lm(proportion ~ genotype, data = test[, c(3,4)])
glht(model1, linfct = mcp(genotype = "Tukey"))

summary(glht(model1, linfct = mcp(genotype = "Tukey")))

#################################
# Post-hoc Tukey test
#################################
# ref: https://rpubs.com/aaronsc32/post-hoc-analysis-tukey

for(i in unique(deconv_raw$Cell_type)){
  print(i)
  plot.data.lm <- lm(proportion ~ genotype, data = deconv_raw[deconv_raw$Cell_type %in% i,])
  plot.data.av <- aov(plot.data.lm)
  summary(plot.data.av)
  
  tukey.test <- TukeyHSD(plot.data.av)
  write.csv(tukey.test$genotype, file = paste0("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig2_bulk_deconvolution_bar_plot_Tukey_", i,".csv"))
  
}

###############################################
# decovolution, ref, real snRNA-seq comparison
###############################################

final_cell_type <- readRDS(file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200818_deconvolution_MuSic_reports/final_cell_type.rda")

final_cell_type$Cell_type <- factor(final_cell_type$Cell_type, levels = c("Neurons", "Astrocytes", "Microglia","Oligo","OPC","Vascular", "Rest"))

g2 <- ggplot(final_cell_type, aes(x=Cell_type, y=proportion, fill  = dataset)) + 
  theme_classic() + scale_y_continuous(expand = c(0,0)) + 
  geom_bar(aes(x = Cell_type, fill = Cell_type, alpha = fct_rev(dataset)),
             stat = "identity", position=position_dodge()) + 
  geom_errorbar(aes(ymin=proportion-sd, ymax=proportion+sd, fill  = fct_rev(dataset)), width=.2,
                                                                    position=position_dodge(0.9)) +
  scale_alpha_discrete("dataset", range = c(.2, 1),
                         guide = guide_legend(reverse = TRUE)) + 
  guides(fill = FALSE) +  
  ylab("proportion (%)") + 
  scale_x_discrete("cell type", limits = rev(c("Neurons", "Astrocytes", "Microglia","Oligo","OPC","Vascular", "Rest"))) + 
  coord_flip() + 
  ggtitle("Comparison cell type composition") + 
  theme(axis.text.x = element_text(size =14, color = "black"), 
        axis.title.y = element_blank(), 
        axis.text.y = element_text(size =14, color = "black"),
        axis.title.x = element_text(size = 16), 
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = c(0.85, 0.5), 
        legend.title = element_blank(),
        legend.text = element_text(size =14, color = "black"))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig2_bulk_deconvolution_bar_plot_bulk_deconvolution_all_figure2.png", 
       plot = g2, 
       scale = 1, width = 9, height = 6, units = "in", device = "png",
       dpi = 300)

###############
# box plot
###############
final_cell_type_2nd <- readRDS(file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200818_deconvolution_MuSic_reports/final_cell_type_2nd.rda")
 
final_cell_type_2nd$dataset <- factor(final_cell_type_2nd$dataset, levels = c("reference", "snRNA-seq", "deconvolution"))

g4 <- ggplot(final_cell_type_2nd, aes(x=Cell_type, y=proportion, fill = dataset, color=dataset)) +
  geom_boxplot(fill="white") + 
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize= 0.3) +
    scale_color_manual(values = c("#CC3300", "#33CC33", "#003366"), labels = c("reference", "snRNA-seq", "deconvolution")) +
    scale_fill_manual(values = c("#CC3300", "#33CC33", "#003366"), labels = c("reference", "snRNA-seq", "deconvolution")) +
    theme_classic() +
    ggtitle("Comparison cell type composition") + 
    ylab("percent") +
    theme(axis.text.x = element_text(size =15, angle = 60, color = "black", vjust = 0.5, family = "Helvetica"),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size =15, color = "black", family = "Helvetica"),
          axis.title.y = element_text(size =15, family = "Helvetica"),
          plot.title = element_text(size = 20, hjust = 0.5),
          legend.title = element_blank(),
          legend.text = element_text(size = 15, family = "Helvetica"),
          legend.key.size = unit(1, "cm")) +
    ylim(0, 70) +
    guides(color = FALSE)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig2_bulk_deconvolution_bar_plot_bulk_deconvolution_all_figure2_boxplot.png", 
       plot = g4, 
       scale = 1, width = 9, height = 6, units = "in", device = "png",
       dpi = 300)

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig2_bulk_deconvolution_bar_plot_session_info.txt")



VlnPlot(data9set_cleaned.SO, features = "Hexb", pt.size = 0.1) + ggtitle("change") + theme(legend.position = "none") + geom_boxplot()


