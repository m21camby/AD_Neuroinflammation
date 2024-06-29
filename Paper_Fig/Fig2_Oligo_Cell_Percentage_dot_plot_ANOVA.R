#! /bin/env RScript
# written by SJK at 10. Dec. 2020
# Email from Shirin Do you think we could get something like Figure 1e) 
# but only for OPC/NFOL/MFOL/MOL so we can also add some statistics?
# I initially implemented Kruskal-Wallis test but I do agin with anova (changed 2021.Dec.20)
 

.libPaths(c("/data/rajewsky/home/skim/R/usr_lib_Seurat/", "/data/rajewsky/home/skim/R/usr_lib_v4/"))
#.libPaths()
unloadNamespace("mgcv")
unloadNamespace("Matrix")

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
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


plot.data <- data9set_sub.meta %>%
  dplyr::select(sample, gemgroup, cell_type = cell_type_precise) %>%
  mutate(cell_type = cell_type) %>%
  group_by(gemgroup, cell_type, sample) %>%
  summarise(count = n()) %>%
  group_by(gemgroup) %>%
  mutate(gemgroup_total = sum(count)) %>%
  mutate(cell_prop = count / gemgroup_total)

plot.data$cell_type <- factor(plot.data$cell_type, levels = c("MOL", "MFOL", "NFOL", "OPC"))

g1 <- ggplot(plot.data, aes(x=cell_type, y=cell_prop, fill = sample, color=sample)) +
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
  ylim(0, 0.8) +
  stat_compare_means(aes(group = sample), label = "p.signif", method = "kruskal.test", show.legend = FALSE, size = 6, family = "Helvetica") +
  guides(color = FALSE)
g1

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Oligo_Cell_Percentage_dot_plot_ANOVA_box_dot_plot_ANOVA.pdf",
       plot = g1,
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)


#####################
# by all cells
#####################
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

plot.data <- meta.df %>%
  dplyr::select(sample, gemgroup, cell_type = cell_type) %>%
  mutate(cell_type = cell_type) %>%
  group_by(gemgroup, cell_type, sample) %>%
  summarise(count = n()) %>%
  group_by(gemgroup) %>%
  mutate(gemgroup_total = sum(count)) %>%
  mutate(cell_prop = count / gemgroup_total)


plot.data_OL <- plot.data[plot.data$cell_type %in% c("MOL", "MFOL", "OPC", "NFOL"), ]

plot.data_OL$cell_type <- factor(plot.data_OL$cell_type, levels = c("OPC", "NFOL", "MFOL", "MOL"))


# I do again from initial kruskal.test to anova (changed 2021.Dec.20)
g1 <- ggplot(plot.data_OL, aes(x=cell_type, y=cell_prop, fill = sample, color=sample)) +
  geom_boxplot(fill="white") +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize= 0.5) +
  scale_color_manual(values = c("#440154FF", viridis(3)[2], "#CC9900"), labels = c("WT", "APPPS1", bquote(APPPS1.Il12b^"-/-"))) +
  scale_fill_manual(values = c("#440154FF", viridis(3)[2], "#CC9900"), labels = c("WT", "APPPS1", bquote(APPPS1.Il12b^"-/-"))) +
  theme_classic() +
  ylab("percent") +
  theme(axis.text.x = element_text(size =15, angle = 60, color = "black", vjust = 0.5, family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size =15, color = "black", family = "Helvetica"),
        axis.title.y = element_text(size =15, family = "Helvetica"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15, family = "Helvetica"),
        legend.key.size = unit(1, "cm")) +
  ylim(0, 0.3) +
  stat_compare_means(aes(group = sample), label = "p.signif", method = "anova", show.legend = FALSE, size = 6, family = "Helvetica") +
  guides(color = FALSE)
g1

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Oligo_Cell_Percentage_dot_plot_ANOVA_box_dot_plot_ANOVA_Final_anova.pdf",
       plot = g1,
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)

# I do again from initial kruskal.test to anova (changed 2021.Dec.20)
g1 <- ggplot(plot.data_OL, aes(x=cell_type, y=cell_prop, fill = sample, color=sample)) +
  geom_boxplot(fill="white") +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize= 0.5) +
  scale_color_manual(values = c("#440154FF", viridis(3)[2], "#CC9900"), labels = c("WT", "APPPS1", bquote(APPPS1.Il12b^"-/-"))) +
  scale_fill_manual(values = c("#440154FF", viridis(3)[2], "#CC9900"), labels = c("WT", "APPPS1", bquote(APPPS1.Il12b^"-/-"))) +
  theme_classic() +
  ylab("percent") +
  theme(axis.text.x = element_text(size =15, angle = 60, color = "black", vjust = 0.5, family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size =15, color = "black", family = "Helvetica"),
        axis.title.y = element_text(size =15, family = "Helvetica"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15, family = "Helvetica"),
        legend.key.size = unit(1, "cm")) +
  ylim(0, 0.3) +
  stat_compare_means(aes(group = sample), label = "p.format", method = "anova", show.legend = FALSE, size = 6, family = "Helvetica") +
  guides(color = FALSE)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Oligo_Cell_Percentage_dot_plot_ANOVA_box_dot_plot_ANOVA_Final_p_format_anova.pdf",
       plot = g1,
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
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

celltypes <- plot.data_OL$cell_type %>% unique()

for(i in celltypes){
  print(i)
  res.aov <- aov(cell_prop ~ sample, data = plot.data_OL[plot.data_OL$cell_type %in% i, ])
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

write.csv(df, file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Oligo_Cell_Percentage_dot_plot_ANOVA_box_dot_plot_ANOVA_Final_raw_data.csv")




##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Oligo_Cell_Percentage_dot_plot_ANOVA_session_info.txt")


