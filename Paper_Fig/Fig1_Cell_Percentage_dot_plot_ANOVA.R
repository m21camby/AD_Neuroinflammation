#! /bin/env RScript
# written by SJK at 3. Mar. 2020

# ref: http://www.sthda.com/english/wiki/ggplot2-box-plot-quick-start-guide-r-software-and-data-visualization
# ref: http://www.sthda.com/english/wiki/ggplot2-dot-plot-quick-start-guide-r-software-and-data-visualization
# ref: http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/
# ref: http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/

.libPaths(c("/data/rajewsky/home/skim/R/usr_lib_Seurat/", "/data/rajewsky/home/skim/R/usr_lib_v4/"))
#.libPaths()
unloadNamespace("mgcv")
unloadNamespace("Matrix")

#.libPaths(c("/home/skim/R/usr_lib", .libPaths()))
#.libPaths(c("/data/murphy/shared_libs/R", "/usr/local/lib/R/site-librar","/usr/local/lib/R/library","/usr/lib64/R/library","/usr/share/R/library", .libPaths() ))

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(viridis)
library(car)
library(fBasics)
source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

load(file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200221_figure1_for_meeting_report_Seurat_object.Robj")

data9set_cleaned.SO$cell_type <- data9set_cleaned.SO@active.ident

plot.data <- data9set_cleaned.SO@meta.data %>%
    dplyr::select(sample, gemgroup, cell_type = cell_type) %>%
  dplyr::mutate(cell_type = cell_type) %>%
  dplyr::group_by(gemgroup, cell_type, sample) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::group_by(gemgroup) %>%
  dplyr::mutate(gemgroup_total = sum(count)) %>%
  dplyr::mutate(cell_prop = count / gemgroup_total)

plot.data <- plot.data[which(plot.data$cell_type %in% c("Excitatory Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", "Inhibitory Interneurons", "OPC")),]


g1 <- ggplot(plot.data, aes(x=cell_type, y=cell_prop, fill = sample, color=sample)) +
  geom_point(position=position_dodge(0.8), size = 2) + scale_color_viridis_d() + scale_fill_viridis_d() + theme_classic() +  
  ylab("percent") +
  theme(axis.text.x = element_text(size =15, angle = 90, color = "black"), 
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size =15, color = "black"), 
        axis.title.y = element_text(size =15), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 15)) + 
  ylim(0, 0.6) + 
  stat_compare_means(aes(group = sample), label = "p.signif", method = "anova", show.legend = FALSE)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Cell_Percentage_dot_plot_ANOVA.png", 
       plot = g1, 
       scale = 1, width = 9, height = 6, units = "in", device = "png",
       dpi = 300)

g1 <- ggplot(plot.data, aes(x=cell_type, y=cell_prop, fill = sample, color=sample)) +
  geom_point(position=position_dodge(0.8), size = 2) + scale_color_viridis_d() + scale_fill_viridis_d() + theme_classic() +  
  ylab("percent") +
  theme(axis.text.x = element_text(size =15, angle = 90, color = "black"), 
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size =15, color = "black"), 
        axis.title.y = element_text(size =15), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 15)) + 
  ylim(0, 0.6) + 
  stat_compare_means(aes(group = sample), label = "p.format", method = "anova", show.legend = FALSE)


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Cell_Percentage_dot_plot_ANOVA_p_format.png", 
       plot = g1, 
  scale = 1, width = 9, height = 6, units = "in", device = "png",
  dpi = 300)


#####################
# summary static save
#####################
df <- data.frame(celltype = character(),
                 Df = numeric(),
                 sum_sq = numeric(),
                 mean_sq = numeric(),
                 f_value = numeric(),
                 p_value = numeric())

celltypes <- plot.data$cell_type %>% unique()

for(i in celltypes){
  print(i)
  res.aov <- aov(cell_prop ~ sample, data = plot.data[plot.data$cell_type %in% i, ])
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

write.csv(df, file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Cell_Percentage_dot_plot_ANOVA_raw_data.csv")




g2 <- ggplot(plot.data, aes(x=cell_type, y=cell_prop, fill = sample, color=sample)) + 
  geom_boxplot(fill="white") +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize= 0.5) + 
  scale_color_viridis_d() + 
  scale_fill_viridis_d() + 
  theme_classic() +  
  ylab("percent") +
  theme(axis.text.x = element_text(size =15, angle = 90, color = "black"), 
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size =15, color = "black"), 
        axis.title.y = element_text(size =15), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 15)) + 
  ylim(0, 0.6) + 
  stat_compare_means(aes(group = sample), label = "p.signif", method = "anova", show.legend = FALSE, size = 7)
#  stat_compare_means(method = "anova", label.y = 0.5)
#  geom_text(geom = "text", aes(x = 4.5, y = 0.5), label = "* < 0.05", color = "black", size = 7)


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Cell_Percentage_box_dot_plot_ANOVA.png", 
       plot = g2, 
       scale = 1, width = 9, height = 6, units = "in", device = "png",
       dpi = 300)


#################################
# Box & dot plot for publication 
#################################

g3 <- ggplot(plot.data, aes(x=cell_type, y=cell_prop, fill = sample, color=sample)) +
  geom_boxplot(fill="white") +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize= 0.5) +
  scale_color_manual(values = c("#440154FF", viridis(3)[2], "#CC9900"), labels = c("WT", "APPPS1", bquote(APPPS1.p40^"-/-"))) +
  scale_fill_manual(values = c("#440154FF", viridis(3)[2], "#CC9900"), labels = c("WT", "APPPS1", bquote(APPPS1.p40^"-/-"))) +
  theme_classic() +
  ylab("percent") +
  theme(axis.text.x = element_text(size =15, angle = 60, color = "black", vjust = 0.5, family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size =15, color = "black", family = "Helvetica"),
        axis.title.y = element_text(size =15, family = "Helvetica"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15, family = "Helvetica"),
        legend.key.size = unit(1, "cm")) +
  ylim(0, 0.6) +
  stat_compare_means(aes(group = sample), label = "p.signif", method = "anova", show.legend = FALSE, size = 6, family = "Helvetica") +
  guides(color = FALSE)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Cell_Percentage_box_dot_plot_ANOVA_Modified.png",
       plot = g3,
       scale = 1, width = 9, height = 6, units = "in", device = "png",
       dpi = 300)

#################################
# Post-hoc Tukey test
#################################
# ref: https://rpubs.com/aaronsc32/post-hoc-analysis-tukey

for(i in unique(plot.data$cell_type)){
  print(i)
  plot.data.lm <- lm(cell_prop ~ sample, data = plot.data[plot.data$cell_type %in% i,])
  plot.data.av <- aov(plot.data.lm)
  summary(plot.data.av)
  
  tukey.test <- TukeyHSD(plot.data.av)
  write.csv(tukey.test$sample, file = paste0("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Cell_Percentage_dot_plot_ANOVA_Tukey_", i,".csv"))
  
}

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Cell_Percentage_dot_plot_ANOVA_session_info.txt")





