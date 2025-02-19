library(viridis)
library(dplyr)
library(SingleCellExperiment)
library(slingshot, quietly = FALSE)
library(dplyr)
library(Seurat)
library(SeuratDisk)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)


load(file = "/nfs/team292/sk27/tmp/20200210_figures_for_meeting_report_cleaned_Seurat_object_res0.8.Robj")

data9set_cleaned.SO$cell_type <- data9set_cleaned.SO@active.ident
table(data9set_cleaned.SO$cell_type)


# ------------------------------- #
# Supple Fig 6A
# ------------------------------- #

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

plot.data3 <- plot.data2[!plot.data2$cell_type %in% c("Unidentified_Neurons"), ]

g1 <- ggplot(plot.data3, aes(x=cell_type, y=proportion2, fill = sample, color=sample)) +
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
  guides(color = FALSE) + scale_x_discrete(breaks=c("Dentate_Gyrus", "CA1", "CA2_3", "subiculum",  "Inhibitory_Neurons"),
                                           labels=c("Dentate Gyrus", "CA1", "CA2/3", "subiculum",  "Inhibitory Neurons"))

g1
ggsave(filename = "~/MDC/Final_revisions/Supple_Fig6_A.pdf", 
       plot = g1, 
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)



f10 <- FeaturePlot(data9set_cleaned.SO, features = "Pvt1", order = TRUE)

ggsave(filename = "~/MDC/Final_revisions/Supple_Fig7_G.pdf", 
       plot = f10, 
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)


f10 <- FeaturePlot(data9set_cleaned.SO, features = "Neat1", order = TRUE)

ggsave(filename = "~/MDC/Final_revisions/Supple_Fig7_F.pdf", 
       plot = f10, 
       scale = 1, width = 7, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)




#################################
# Neuron proportion box plot 
#################################
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

table(plot.data2$cell_type)

plot.data2 <- plot.data2[plot.data2$cell_type %in% c("CA1", "CA2_3", "Dentate_Gyrus", "subiculum", "Inhibitory_Neurons"), ]

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
  guides(color = FALSE) + scale_x_discrete(breaks=c("Dentate_Gyrus", "CA1", "CA2_3", "subiculum",  "Inhibitory_Neurons"),
                                           labels=c("Dentate Gyrus", "CA1", "CA2/3", "subiculum",  "Inhibitory Neurons"))
g1

data.df <- ggplot_build(g1)$data

data_1 <- data.df[[1]] %>% as.data.frame()
data_1$outliers <- NULL
write.csv(data_1, file = "~/MDC/Final_revisions/Extend_Fig6a_boxplot_1st.csv")
data_2 <- data.df[[2]] %>% as.data.frame()
write.csv(data_2, file = "~/MDC/Final_revisions/Extend_Fig6a_boxplot_2nd.csv")
data_3 <- data.df[[3]] %>% as.data.frame()
write.csv(data_3, file = "~/MDC/Final_revisions/Extend_Fig6a_boxplot_3rd.csv")


g2 <- ggplot(plot.data2, aes(x=cell_type, y=proportion2, fill = sample, color=sample)) +
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
  guides(color = FALSE) + scale_x_discrete(breaks=c("Dentate_Gyrus", "CA1", "CA2_3", "subiculum",  "Inhibitory_Neurons"),
                                           labels=c("Dentate Gyrus", "CA1", "CA2/3", "subiculum",  "Inhibitory Neurons"))







