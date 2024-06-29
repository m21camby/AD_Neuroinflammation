#! /bin/env RScript
# written by SJK at 27. 08. 2020
# Creating cell type proportion of Neuronal cell type by barplot
# File name: Fig5_Neurons_cell_type_proportion.R

##################
# library loading
##################

library(Seurat)
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

source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

# data loading
load(file="/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200210_figures_for_meeting_report_cleaned_Seurat_object_res0.8.Robj")

# Meta data for new Neuronal cell type assignment
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

# give levels to large cell type
data9set_cleaned.SO$cell_type <- factor(data9set_cleaned.SO$cell_type, levels = c("Dentate_Gyrus", "CA1", "CA2_3", "subiculum",  "Unidentified_Neurons", 
                                                                                  "Inhibitory_Neurons", "MOL", "MFOL", "OPC", "NFOL", "Microglia","Astrocytes",
                                                                                  "Vascular", "VLMC", "Choroid", "Fibroblast", "Cajal", "Pericyte", "Macrophage"))

plot.data2 <- data9set_cleaned.SO@meta.data %>%
  dplyr::select(sample, cell_type, gemgroup, cluster = seurat_clusters) %>%
  mutate(cluster = cluster) %>%
  group_by(cell_type, gemgroup, sample) %>% 
  dplyr::summarise(count = n()) %>% mutate(total = sum(count)) 

plot.data2 <- plot.data2 %>% dplyr::group_by(gemgroup) %>% dplyr::mutate(sample_total = sum(count)) 

# add total cell number
plot.data2$total <- 82298

# as data frame
plot.data2 <- as.data.frame(plot.data2)

# calculate by cell type
plot.data2$proportion <- plot.data2$count / plot.data2$sample_total

# levels to sample
plot.data2$sample <- factor(plot.data2$sample, levels = c("Ctrl", "AD", "ADp40KO"))

# percentage
plot.data2$proportion2 <- round(100*plot.data2$proportion,2)

plot.data2 <- plot.data2[c(1:54),]

g1 <- ggplot(plot.data2, aes(x=cell_type, y=proportion2, fill = sample, color=sample)) +
  geom_boxplot(fill="white") + 
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize= 0.3) +
  scale_color_manual(values = c("#440154FF", viridis(3)[2], "#CC9900"), labels = c("WT", "APPPS1", bquote(APPPS1.il12b^"-/-"))) +
  scale_fill_manual(values = c("#440154FF", viridis(3)[2], "#CC9900"), labels = c("WT", "APPPS1", bquote(APPPS1.il12b^"-/-"))) +
  theme_classic() +
  ggtitle("Neuronal cell types proportions") + 
  ylab("percent (%)") +
  theme(axis.text.x = element_text(size =15, angle = 45, color = "black", vjust = 0.5, family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size =15, color = "black", family = "Helvetica"),
        axis.title.y = element_text(size =15, family = "Helvetica"),
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 15, family = "Helvetica"),
        legend.key.size = unit(1, "cm")) +
  ylim(0, 35) +
  stat_compare_means(aes(group = sample), label = "p.signif", method = "anova", show.legend = FALSE, size = 6, family = "Helvetica") +
  guides(color = FALSE) + scale_x_discrete(breaks=c("Dentate_Gyrus", "CA1", "CA2_3", "subiculum", "Unidentified_Neurons", "Inhibitory_Neurons"),
                                           labels=c("Dentate Gyrus", "CA1", "CA2/3", "subiculum", "unidentified Neurons", "Inhibitory Neurons"))


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Neurons_cell_type_proportion_cell_type_boxplot.png", 
       plot = g1, 
       scale = 1, width = 9, height = 6, units = "in", device = "png",
       dpi = 300)

g1 <- ggplot(plot.data2, aes(x=cell_type, y=proportion2, fill = sample, color=sample)) +
  geom_boxplot(fill="white") + 
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize= 0.3) +
  scale_color_manual(values = c("#440154FF", viridis(3)[2], "#CC9900"), labels = c("WT", "APPPS1", bquote(APPPS1.il12b^"-/-"))) +
  scale_fill_manual(values = c("#440154FF", viridis(3)[2], "#CC9900"), labels = c("WT", "APPPS1", bquote(APPPS1.il12b^"-/-"))) +
  theme_classic() +
  ggtitle("Neuronal cell types proportions") + 
  ylab("percent (%)") +
  theme(axis.text.x = element_text(size =15, angle = 45, color = "black", vjust = 0.5, family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size =15, color = "black", family = "Helvetica"),
        axis.title.y = element_text(size =15, family = "Helvetica"),
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 15, family = "Helvetica"),
        legend.key.size = unit(1, "cm")) +
  ylim(0, 35) +
  stat_compare_means(aes(group = sample), label = "p.format", method = "anova", show.legend = FALSE, size = 6, family = "Helvetica") +
  guides(color = FALSE) + scale_x_discrete(breaks=c("Dentate_Gyrus", "CA1", "CA2_3", "subiculum", "Unidentified_Neurons", "Inhibitory_Neurons"),
                                           labels=c("Dentate Gyrus", "CA1", "CA2/3", "subiculum", "unidentified Neurons", "Inhibitory Neurons"))


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Neurons_cell_type_proportion_cell_type_boxplot_p_format.png", 
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

celltypes <- plot.data2$cell_type %>% unique()

for(i in celltypes){
  print(i)
  res.aov <- aov(proportion ~ sample, data = plot.data2[plot.data2$cell_type %in% i, ])
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

write.csv(df, file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Neurons_cell_type_proportion_cell_type_raw_data.csv")




##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Neurons_cell_type_proportion_session_info.txt")
