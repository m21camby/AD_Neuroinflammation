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
library(gridExtra)



load(file = "/nfs/team292/sk27/tmp/20200210_figures_for_meeting_report_cleaned_Seurat_object_res0.8.Robj")



data9set_cleaned.SO$cell_type <- data9set_cleaned.SO@active.ident
table(data9set_cleaned.SO$cell_type)

#data9set_cleaned.SO$cell_type <- ifelse(data9set_cleaned.SO$cell_type %in% c("Dentate Gyrus", "CA1 Neurons", "CA2/CA3 Neuron", "Subiculum", "Neurons"), "Excitatory Neurons",
#                                               ifelse(data9set_cleaned.SO$cell_type %in% c("Astrocytes"), "Astrocytes",
#                                        ifelse(data9set_cleaned.SO$cell_type %in% c("Inhibitory Interneurons"), "Inhibitory Interneurons",
#                                                      ifelse(data9set_cleaned.SO$cell_type %in% c("Microglia"), "Microglia",
#                                                             ifelse(data9set_cleaned.SO$cell_type %in% c("Oligo"), "Oligodendrocytes",
#                                                                    ifelse(data9set_cleaned.SO$cell_type %in% c("OPC"), "OPC", "Rest"
#                                                                    ))))))

meta.df <- data9set_cleaned.SO@meta.data %>% mutate(cell_type = case_when(seurat_clusters %in% c(0, 10) ~ "Dentate Gyrus",
                                                                          seurat_clusters %in% c(2, 13, 18, 40) ~ "CA1",
                                                                          seurat_clusters %in% c(9, 21, 27, 32) ~ "CA2/3",
                                                                          seurat_clusters %in% c(11, 14, 15, 17, 19, 20) ~ "subiculum",
                                                                          seurat_clusters %in% c(24, 29) ~ "Unidentified\nNeurons",
                                                                          seurat_clusters %in% c(7, 16, 22, 30) ~ "Inhibitory Neurons",
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
data9set_cleaned.SO$cell_type <- factor(data9set_cleaned.SO$cell_type, levels = c("Dentate Gyrus", "CA1", "CA2/3", "subiculum",  "Unidentified\nNeurons", 
                                                                                  "Inhibitory Neurons", "MOL", "MFOL", "OPC", "NFOL", "Microglia","Astrocytes",
                                                                                  "Vascular", "VLMC", "Choroid", "Fibroblast", "Cajal", "Pericyte", "Macrophage"))



# --------------------- #
# Fig 2.B
# --------------------- #

d1 <- DimPlot(data9set_cleaned.SO, reduction = "umap",  group.by = "cell_type",
              label = FALSE, label.size = 6, 
              cols = c("#CC0000", "#FF0000", "#990000", "#660000", "#FF3366", 
                       "#FF6600",
                       "#006666", "#669999","#99CCCC","#9966FF", 
                       "#99CC00","#FFCC66", "#000033", "#000033","#CC9900", "#3399FF", "#FF9900", "#3399FF","#000033")) + 
  theme(legend.position = "none", 
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) + labs(title = NULL)  
d1

ggsave(filename = "~/MDC/Final_revisions/Fig2_B.pdf", 
       plot = d1, 
       scale = 1, width = 10, height = 10, units = "in", device = cairo_pdf,
       dpi = 300)



# --------------------- #
# Fig 2.C
# --------------------- #
data9set_cleaned.SO$cell_type <- ifelse(data9set_cleaned.SO$cell_type %in% c("Dentate Gyrus", "CA1 Neurons", "CA2/CA3 Neuron", "Subiculum", "Neurons"), "Excitatory Neurons",
                                        ifelse(data9set_cleaned.SO$cell_type %in% c("Inhibitory Interneurons"), "Inhibitory Interneurons",
                                               ifelse(data9set_cleaned.SO$cell_type %in% c("Astrocytes"), "Astrocytes",
                                                      ifelse(data9set_cleaned.SO$cell_type %in% c("Microglia"), "Microglia",
                                                             ifelse(data9set_cleaned.SO$cell_type %in% c("Oligo"), "Oligodendrocytes",
                                                                    ifelse(data9set_cleaned.SO$cell_type %in% c("OPC"), "OPC", "Rest"
                                                                           
                                                                    ))))))



plot.data <- data9set_cleaned.SO@meta.data %>%
  dplyr::select(sample, gemgroup, cell_type = cell_type) %>%
  dplyr::mutate(cell_type = cell_type) %>%
  dplyr::group_by(gemgroup, cell_type, sample) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::group_by(gemgroup) %>%
  dplyr::mutate(gemgroup_total = sum(count)) %>%
  dplyr::mutate(cell_prop = count / gemgroup_total)

plot.data <- plot.data[which(plot.data$cell_type %in% c("Excitatory Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", "Inhibitory Interneurons", "OPC")),]


plot.data$cell_type <- factor(plot.data$cell_type, 
                              levels = c("Excitatory Neurons","Inhibitory Interneurons", "Astrocytes", "Microglia","Oligodendrocytes",   "OPC"))

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
  ylim(0, 0.6) + guides(color = FALSE)
#+
#  stat_compare_means(aes(group = sample), label = "p.signif", method = "anova", show.legend = FALSE, size = 6, family = "Helvetica") +
#  guides(color = FALSE)
g3

ggsave(filename = "~/MDC/Final_revisions/Fig2_C.pdf",
       plot = g3,
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)


# --------------------- #
# Fig 2. box plot table
# --------------------- #
data9set_cleaned.SO$broad_cell_type <- ifelse(data9set_cleaned.SO$cell_type %in% c("Dentate Gyrus", "CA1 Neurons", 'CA2/CA3 Neuron',
                                                                                   ' Subiculum', 'Neurons'), "Excitatory Neurons", 
                                              ifelse(data9set_cleaned.SO$cell_type %in% c('Inhibitory Interneurons'), "Inhibitory Neurons",
                                                     ifelse(data9set_cleaned.SO$cell_type %in% c('Oligo'), "Oligodendrocytes",
                                                            ifelse(data9set_cleaned.SO$cell_type %in% c('Microglia'), "Microglia", 
                                                                   ifelse(data9set_cleaned.SO$cell_type %in% c('OPC'), "OPC", 
                                                                          ifelse(data9set_cleaned.SO$cell_type %in% c('Astrocytes'), "Astrocytes", 
                                                                                 "Rest"))))))

table(data9set_cleaned.SO$broad_cell_type)

plot.data <- data9set_cleaned.SO@meta.data %>%
  dplyr::select(sample, gemgroup, broad_cell_type = broad_cell_type) %>%
  dplyr::mutate(broad_cell_type = broad_cell_type) %>%
  dplyr::group_by(gemgroup, broad_cell_type, sample) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::group_by(gemgroup) %>%
  dplyr::mutate(gemgroup_total = sum(count)) %>%
  dplyr::mutate(cell_prop = count / gemgroup_total)

plot.data <- plot.data[which(plot.data$broad_cell_type %in% c("Excitatory Neurons", "Astrocytes", "Oligodendrocytes", "Microglia", "Inhibitory Neurons", "OPC", "Rest")),]
plot.data
g1 <- ggplot(plot.data, aes(x=broad_cell_type, y=cell_prop, fill = sample, color=sample)) +
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
g1

g2 <- ggplot(plot.data, aes(x=broad_cell_type, y=cell_prop, fill = sample, color=sample)) + 
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
g2

data.df <- ggplot_build(g2)$data

data_1 <- data.df[[1]] %>% as.data.frame()
data_1$outliers <- NULL
write.csv(data_1, file = "~/MDC/Final_revisions/Fig2c_boxplot_1st.csv")
data_2 <- data.df[[2]] %>% as.data.frame()
write.csv(data_2, file = "~/MDC/Final_revisions/Fig2c_boxplot_2nd.csv")
data_3 <- data.df[[3]] %>% as.data.frame()
write.csv(data_3, file = "~/MDC/Final_revisions/Fig2c_boxplot_3rd.csv")















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

# Again with big cell types

data9set_cleaned.SO$cell_types <- ifelse(data9set_cleaned.SO$cell_type %in% c("CA1", "CA2_3", "Dentate_Gyrus","Unidentified_Neurons", "subiculum"), "Excitatory Neurons",
                                         ifelse(data9set_cleaned.SO$cell_type %in% "Inhibitory_Neurons", "Inhibitory Neurons",
                                                ifelse(data9set_cleaned.SO$cell_type %in% c("Microglia"), "Microglia",
                                                       ifelse(data9set_cleaned.SO$cell_type %in% c("Macrophage"), "Macrophages",
                                                              ifelse(data9set_cleaned.SO$cell_type %in% c("MOL", "MFOL", "NFOL"), "Oligodendrocytes",
                                                                     ifelse(data9set_cleaned.SO$cell_type %in% c("OPC"), "OPC",
                                                                            ifelse(data9set_cleaned.SO$cell_type %in% c("Astrocytes"), "Asctrocytes",
                                                                                   ifelse(data9set_cleaned.SO$cell_type %in% c("Cajal"), "Cajal Retzius",
                                                                                          ifelse(data9set_cleaned.SO$cell_type %in% c("Fibroblast"), "Fibroblast",
                                                                                                 ifelse(data9set_cleaned.SO$cell_type %in% c("Vascular", "VLMC", "Pericyte"), "Vascular cells", "Choroid Plexus"))))))))))



data9set_cleaned.SO$cell_types <- factor(data9set_cleaned.SO$cell_types, levels = c("Excitatory Neurons", "Inhibitory Neurons", "Cajal Retzius", "Choroid Plexus", "Asctrocytes", "Microglia", "Macrophages", "Oligodendrocytes", "OPC", "Fibroblast", "Vascular cells"))



# meta data
data9set_cleaned.meta <- data9set_cleaned.SO@meta.data %>% as.data.frame

v1 <- VlnPlot(data9set_cleaned.SO, features = "nCount_RNA", group.by = "cell_types", pt.size = 0, cols = c("#CC0000", "#FF6600", "#FF9900","#CC9900", "#FFCC66", "#99CC00","#003300", "#99CCCC", "#009999", "#3399FF","#000033")) + 
  ggtitle("nUMIs per cell by cell type") + ylim(c(0, 31000)) +
  theme(legend.position = "none", axis.title.x = element_blank()) + 
  annotate("text", x = 1, y = 30500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Excitatory Neurons", ]$nCount_RNA)) +
  annotate("text", x = 2, y = 30500, label = round(median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Inhibitory Neurons", ]$nCount_RNA))) +
  annotate("text", x = 3, y = 30500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Cajal Retzius", ]$nCount_RNA)) +
  annotate("text", x = 4, y = 30500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Choroid Plexus", ]$nCount_RNA)) +
  annotate("text", x = 5, y = 30500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Asctrocytes", ]$nCount_RNA)) +
  annotate("text", x = 6, y = 30500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Microglia", ]$nCount_RNA)) +
  annotate("text", x = 7, y = 30500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Macrophages", ]$nCount_RNA)) +
  annotate("text", x = 8, y = 30500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Oligodendrocytes", ]$nCount_RNA)) +
  annotate("text", x = 9, y = 30500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "OPC", ]$nCount_RNA)) +
  annotate("text", x = 10, y = 30500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Fibroblast", ]$nCount_RNA)) +
  annotate("text", x = 11, y = 30500, label = round(median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Vascular cells", ]$nCount_RNA)))

ggsave(filename = "~/MDC/Final_revisions/Supple_Fig1_E.pdf", 
       plot = v1, 
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)

v1 <- VlnPlot(data9set_cleaned.SO, features = "nFeature_RNA", group.by = "cell_types", pt.size = 0, cols = c("#CC0000", "#FF6600", "#FF9900","#CC9900", "#FFCC66", "#99CC00","#003300", "#99CCCC", "#009999", "#3399FF","#000033")) + 
  ggtitle("nGenes per cell by cell type") + ylim(c(0, 8000)) +
  theme(legend.position = "none", axis.title.x = element_blank()) + 
  annotate("text", x = 1, y = 7500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Excitatory Neurons", ]$nFeature_RNA)) +
  annotate("text", x = 2, y = 7500, label = round(median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Inhibitory Neurons", ]$nFeature_RNA))) +
  annotate("text", x = 3, y = 7500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Cajal Retzius", ]$nFeature_RNA)) +
  annotate("text", x = 4, y = 7500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Choroid Plexus", ]$nFeature_RNA)) +
  annotate("text", x = 5, y = 7500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Asctrocytes", ]$nFeature_RNA)) +
  annotate("text", x = 6, y = 7500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Microglia", ]$nFeature_RNA)) +
  annotate("text", x = 7, y = 7500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Macrophages", ]$nFeature_RNA)) +
  annotate("text", x = 8, y = 7500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Oligodendrocytes", ]$nFeature_RNA)) +
  annotate("text", x = 9, y = 7500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "OPC", ]$nFeature_RNA)) +
  annotate("text", x = 10, y = 7500, label = median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Fibroblast", ]$nFeature_RNA)) +
  annotate("text", x = 11, y = 7500, label = round(median(data9set_cleaned.meta[data9set_cleaned.meta$cell_types %in% "Vascular cells", ]$nFeature_RNA)))


ggsave(filename = "~/MDC/Final_revisions/Supple_Fig1_D.pdf", 
       plot = v1, 
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)





# AS Fig3

data9set_cleaned_AS.SO <- subset(data9set_cleaned.SO, subset = seurat_clusters %in% c(4))

f1 <- FeaturePlot(data9set_AS.SO, features = c("Gfap"), cols = c("#FFCC00",'#003366'))
f2 <- FeaturePlot(data9set_cleaned_AS.SO, features = c("Gfap"), cols = c("#FFCC00",'#003366'))
f2

ggsave(filename = "~/MDC/Final_revisions/Supple_Fig3_G.pdf", 
       plot = f2, 
       scale = 1, width = 7, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)
