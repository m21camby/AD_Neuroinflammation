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




# --------------------- #
# oligo Fig4A
# --------------------- #

data9set_OL_OPC.SO <- subset(data9set_cleaned.SO, subset = seurat_clusters %in% c(1, 5, 6, 12, 38))


f1 <- FeaturePlot(data9set_OL_OPC.SO, features = c("Pdgfra")) + 
  scale_x_continuous(limits = c(10,25)) + 
  scale_y_continuous(limits = c(-20, 10)) + 
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", limits=c(0, 5), midpoint = 2.5, breaks=c(0,2.5,5)) + 
  ggtitle("OPC (Pdgfra)") + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

f2 <- FeaturePlot(data9set_OL_OPC.SO, features = c("Fyn")) + 
  scale_x_continuous(limits = c(10,25)) + 
  scale_y_continuous(limits = c(-20, 10)) + 
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", limits=c(0, 5), midpoint = 2.5, breaks=c(0,2.5,5)) + 
  ggtitle("COPs (Fyn)") + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

f3 <- FeaturePlot(data9set_OL_OPC.SO, features = c("Tcf7l2")) + 
  scale_x_continuous(limits = c(10,25)) + 
  scale_y_continuous(limits = c(-20, 10)) + 
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", limits=c(0, 5), midpoint = 2.5, breaks=c(0,2.5,5)) + 
  ggtitle("NFOL (Tcf7l2)") + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

f4 <- FeaturePlot(data9set_OL_OPC.SO, features = c("Opalin")) + 
  scale_x_continuous(limits = c(10,25)) + 
  scale_y_continuous(limits = c(-20, 10)) + 
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", limits=c(0, 5), midpoint = 2.5, breaks=c(0,2.5,5)) + 
  ggtitle("MFOL (Opalin)") + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

f5 <- FeaturePlot(data9set_OL_OPC.SO, features = c("Apod")) + 
  scale_x_continuous(limits = c(10,25)) + 
  scale_y_continuous(limits = c(-20, 10)) + 
  scale_color_gradient2(low="#003366", high='#FFCC00', mid="#66FFFF", limits=c(0, 5), midpoint = 2.5, breaks=c(0,2.5,5)) + 
  ggtitle("MOL (Apod)") + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank())

g1 <- arrangeGrob(f1, f2, f3, f4, f5, ncol = 5, widths= c(1,1,1,1,1.35))
g1
ggsave(filename = "~/MDC/Final_revisions/Fig4_A.pdf",
       plot = g1,
       scale = 1, width = 12, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)


# --------------------- #
# Fig 4B
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



data9set_sub.SO <- subset(data9set_cleaned.SO, subset = cell_type %in% c("Oligo", "OPC"))

data9set_sub.meta <- data9set_sub.SO@meta.data %>% as.data.frame

data9set_sub.meta$cell_barcode <- rownames(data9set_sub.meta)

data9set_sub.meta <- data9set_sub.meta %>% mutate(cell_type_precise = case_when(seurat_clusters %in% 6 ~ "MOL",
                                                                                seurat_clusters %in% c(1,5) ~ "MFOL",
                                                                                seurat_clusters %in% 12 ~ "OPC",
                                                                                seurat_clusters %in% 38 ~ "NFOL"))

rownames(data9set_sub.meta) <- data9set_sub.meta$cell_barcode

group_name <- data9set_sub.SO@meta.data$seurat_clusters %>% as.character

group_name <- factor(group_name, levels = c(6, 1, 5, 38, 12))

space <- readRDS("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20201027_IL12b_downstream_analysis_scorpius_space.rda")
traj2 <- readRDS("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20201027_IL12b_downstream_analysis_scorpius_trajectory_rev.rda")


palette <- 
  setNames(
    c("#999999", "#009999", "#006666", "#99CCCC", "red"),
    c(1,5,6,12,38)
  )


d1 <- draw_trajectory_plot(
  space,
  progression_group = group_name,
  path = traj2$path,
  point_size = 2,
  point_alpha = 0.8,
  path_size = 1, 
  contour = FALSE,
  progression_group_palette = palette) + 
  ggtitle("Oligo differentiation trajectory") + 
  scale_color_manual(name = "Cell type", labels = c("MOL(cluster 6)", "MFOL(cluster 1)", "MFOL(cluster 5)",  "NFOL(cluster 38)", "OPC(cluster 12)"),
                     values = c("#006666", "#669999", "#009999", "#9966FF",  "#99CCCC")) + 
  theme(axis.text = element_text(size = 12, color = "black", family = "helvetica"),
        axis.title = element_text(size = 12, color = "black", family = "helvetica"),
        plot.title = element_text(hjust = 0.5, size = 15, color = "black", family = "helvetica"),
        legend.title = element_text(size = 12, color = "black", family = "helvetica"),
        legend.text = element_text(size = 11, color = "black", family = "helvetica"))



# --------------------- #
# Fig 4D
# --------------------- #

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
ggsave(filename = "~/MDC/Final_revisions/Fig2_Oligo_Cell_Percentage_dot_plot_ANOVA_box_dot_plot_ANOVA_Final_anova.pdf",
       plot = g1,
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)

# without stat
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
  #stat_compare_means(aes(group = sample), label = "p.signif", method = "anova", show.legend = FALSE, size = 6, family = "Helvetica") +
  guides(color = FALSE)
g1
ggsave(filename = "~/MDC/Final_revisions/Fig4_D.pdf",
       plot = g1,
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)



# --------------------- #
# check box plot 
# --------------------- #
table(plot.data_OL$cell_type)


data.df <- ggplot_build(g1)$data

data_1 <- data.df[[1]] %>% as.data.frame()
data_1$outliers <- NULL
write.csv(data_1, file = "~/MDC/Final_revisions/Fig4d_boxplot_1st.csv")
data_2 <- data.df[[2]] %>% as.data.frame()
write.csv(data_2, file = "~/MDC/Final_revisions/Fig4d_boxplot_2nd.csv")
data_3 <- data.df[[3]] %>% as.data.frame()
write.csv(data_3, file = "~/MDC/Final_revisions/Fig4d_boxplot_3rd.csv")






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

g1

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

