

library(Seurat)
library(ggplot2)
library(gridExtra)
library(grid)
library(ggrepel)
library(DT)
library(biomaRt)
library(dplyr)
library(viridis)
library(tidyr)
library(scran)
library(edgeR)
library(topGO)
library(scales)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(viridis)

load(file="/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200210_figures_for_meeting_report_cleaned_Seurat_object_res0.8.Robj")


meta.df <- data9set_cleaned.SO@meta.data %>% mutate(cell_type = case_when(seurat_clusters %in% c(0, 10) ~ "Dentate_Gyrus",
                                                                          seurat_clusters %in% c(2, 13, 18, 40) ~ "CA1",
                                                                          seurat_clusters %in% c(9, 21, 27, 32) ~ "CA2_3",
                                                                          seurat_clusters %in% c(11, 14, 15, 17, 19, 20) ~ "subiculum",
                                                                          seurat_clusters %in% c(24, 29) ~ "Unidentified_Neurons",
                                                                          seurat_clusters %in% c(7, 16, 22, 30) ~ "Inhibitory_Neurons",
                                                                          seurat_clusters %in% c(6) ~ "MOL",
                                                                          seurat_clusters %in% c(1, 5) ~ "MFOL",
                                                                          seurat_clusters %in% c(12) ~ "NFOL",
                                                                          seurat_clusters %in% c(38) ~ "OPC",
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

plot.data <- data9set_cleaned.SO@meta.data %>%
  dplyr::select(sample, gemgroup, cell_type = cell_type) %>%
  mutate(cell_type = cell_type) %>%
  group_by(gemgroup, cell_type, sample) %>%
  summarise(count = n()) %>%
  group_by(gemgroup) %>%
  mutate(gemgroup_total = sum(count)) %>%
  mutate(cell_prop = count / gemgroup_total)


plot.data2 <- plot.data[plot.data$cell_type %in% c("CA1", "CA2_3", "Dentate_Gyrus", "subiculum", "Inhibitory_Neurons"), ]

plot.data2$cell_type <- factor(plot.data2$cell_type, levels = rev(c("CA1", "CA2_3", "Dentate_Gyrus", "subiculum", "Inhibitory_Neurons")))

g1 <- ggplot(plot.data2, aes(x=cell_type, y=cell_prop,  color=sample)) +
  geom_boxplot(outlier.shape=NA, position = position_dodge(-.9)) +
  geom_point(position = position_dodge(-.9), pch=20) + 
  scale_x_discrete(breaks=c("CA1", "CA2_3", "Dentate_Gyrus", "subiculum", "Inhibitory_Neurons"), labels=c("CA1", "CA2/3", "Dentate Gyrus", "subiculum", "Inhibitory Neurons")) + 
  scale_color_viridis_d(labels = c("WT", "APPPS1", "APPPS1.il12b-/-"))  + theme_classic() +  
  ylab("percent") +
  theme(axis.text.x = element_text(size =15, color = "black"), 
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size =15, color = "black"), 
        axis.title.y = element_blank(), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 15),
        legend.position = c(0.68, 0.2)) + 
  ylim(0, 0.4) + coord_flip() + 
  stat_compare_means(aes(group = sample), label = "p.signif", method = "anova", show.legend = FALSE)


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Neurons_cell_type_proportion_cell_type_boxplot_2nd.png",
       plot = g1,
       scale = 1, width = 7, height = 4.5, units = "in", device = "png",
       dpi = 300)


