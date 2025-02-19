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




# Fig 3 


data9set_cleaned_Ctrl.SO <- subset(data9set_cleaned.SO, subset = sample %in% "Ctrl")
data9set_cleaned_AD.SO <- subset(data9set_cleaned.SO, subset = sample %in% "AD")
data9set_cleaned_ADp40KO.SO <- subset(data9set_cleaned.SO, subset = sample %in% "ADp40KO")



f1 <- FeaturePlot(data9set_cleaned_Ctrl.SO, features = "Il12rb1", cols = c("#CCCCCC", "#0000FF"), order = TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
f1 <- f1 + ggtitle("WT") + labs(y = "Il12rb1") + scale_color_continuous(high = "blue", low = "gray", limit = c(0,3)) + 
  theme(legend.position = "none", 
        axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title.x = element_blank())
f2 <- FeaturePlot(data9set_cleaned_AD.SO, features = "Il12rb1", cols = c("#CCCCCC", "#0000FF"), sort.cell =  TRUE, order = TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
f2 <- f2 + ggtitle("AD") + scale_color_continuous(high = "blue", low = "gray", limit = c(0,3)) + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank())

f3 <- FeaturePlot(data9set_cleaned_ADp40KO.SO, features = "Il12rb1", cols = c("#CCCCCC", "#0000FF"), order = TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
f3 <- f3 + ggtitle("AD.p40KO") + scale_color_continuous(high = "blue", low = "gray", limit = c(0,3)) + 
  theme(legend.position = "none",axis.line=element_blank(),  axis.ticks=element_blank(), axis.text=element_blank(), axis.title = element_blank())


g1 <- arrangeGrob(f1, f3, f2, ncol = 3, widths= c(1.1,1.05,1.2))


ggsave(filename = "~/MDC/Final_revisions/Fig3_A.pdf", 
       plot = g1, 
       scale = 1, width = 15, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

f1 <- FeaturePlot(data9set_cleaned_Ctrl.SO, features = "Il12rb2", cols = c("#CCCCCC", "#0000FF"), order = TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
f1 <- f1 + ggtitle("WT") + labs(y = "Il12rb2") + scale_color_continuous(high = "blue", low = "gray", limit = c(0,3)) + 
  theme(legend.position = "none", 
        axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title.x = element_blank())


f2 <- FeaturePlot(data9set_cleaned_AD.SO, features = "Il12rb2", cols = c("#CCCCCC", "#0000FF"), sort.cell =  TRUE, order = TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
f2 <- f2 + ggtitle("AD") + scale_color_continuous(high = "blue", low = "gray", limit = c(0,3)) + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank())

f3 <- FeaturePlot(data9set_cleaned_ADp40KO.SO, features = "Il12rb2", cols = c("#CCCCCC", "#0000FF"), order = TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
f3 <- f3 + ggtitle("AD.p40KO") + scale_color_continuous(high = "blue", low = "gray", limit = c(0,3)) + 
  theme(legend.position = "none",axis.line=element_blank(),  axis.ticks=element_blank(), axis.text=element_blank(), axis.title = element_blank())


g2 <- arrangeGrob(f1, f3, f2, ncol = 3, widths= c(1.1,1.05,1.2))


ggsave(filename = "~/MDC/Final_revisions/Fig3_B.pdf", 
       plot = g2, 
       scale = 1, width = 15, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)



# Il23r

f1 <- FeaturePlot(data9set_cleaned_Ctrl.SO, features = "Il23r", cols = c("#CCCCCC", "#0000FF"), order = TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
f1 <- f1 + ggtitle("WT") + labs(y = "Il23r") + scale_color_continuous(high = "blue", low = "gray", limit = c(0,3)) +
  theme(legend.position = "none", 
        axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title.x = element_blank())


f2 <- FeaturePlot(data9set_cleaned_AD.SO, features = "Il23r", cols = c("#CCCCCC", "#0000FF"), sort.cell =  TRUE, order = TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
f2 <- f2 + ggtitle("AD") + scale_color_continuous(high = "blue", low = "gray", limit = c(0,3)) + 
  theme(axis.line=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(),
        axis.title = element_blank())

f3 <- FeaturePlot(data9set_cleaned_ADp40KO.SO, features = "Il23r", cols = c("#CCCCCC", "#0000FF"), order = TRUE, pt.size = 0.5, min.cutoff = 0, max.cutoff = 5)
f3 <- f3 + ggtitle("AD.p40KO") + scale_color_continuous(high = "blue", low = "gray", limit = c(0,3)) + 
  theme(legend.position = "none",axis.line=element_blank(),  axis.ticks=element_blank(), axis.text=element_blank(), axis.title = element_blank())

g5 <- arrangeGrob(f1, f3, f2, ncol = 3, widths= c(1.1,1.05,1.2))

ggsave(filename = "~/MDC/Final_revisions/Fig3_C.pdf", 
       plot = g5, 
       scale = 1, width = 15, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)





meta.df <- meta.df %>% mutate(large_cell_type = case_when(cell_type %in% c("Dentate Gyrus", "CA1", "CA2/3","subiculum","Unidentified\nNeurons") ~ "Excitatory Neurons",
                                                          cell_type %in% c("Inhibitory Neurons") ~ "Inhibitory Neurons",
                                                          cell_type %in% c("MOL", "MFOL","NFOL") ~ "Oligodendrocytes",
                                                          cell_type %in% c("Astrocytes") ~ "Astrocytes",
                                                          cell_type %in% c("Microglia") ~ "Microglia",
                                                          cell_type %in% c("OPC") ~ "OPC",
                                                          cell_type %in% c("Vascular", "VLMC","Choroid","Fibroblast","Cajal","Pericyte","Macrophage") ~ "Rest"))


meta.df$large_cell_type <- factor(meta.df$large_cell_type, levels = c("Excitatory Neurons", "Inhibitory Neurons", "Oligodendrocytes", "OPC", "Microglia", "Astrocytes", "Rest"))




meta.df$Il12rb1 <- data9set_cleaned.SO@assays$RNA@data[rownames(data9set_cleaned.SO@assays$RNA@data) %in% "Il12rb1", ]

meta.df2 <- meta.df[!meta.df$large_cell_type %in% c("Rest"), ]

v1 <- ggplot(meta.df2, aes(x=large_cell_type, y=Il12rb1, color=large_cell_type))  + 
  geom_violin(trim=FALSE, width = 1) + scale_color_manual(values=c("#CC0000", "#FF6600", "#99CCCC", "#009999", "#99CC00","#FFCC66", "purple")) + 
  theme_classic() + 
  geom_jitter(position=position_jitter(0.2), size = 1.5, alpha = 1) + scale_y_continuous(expand = c(0,0),limits=c(0, 4)) +  
  ggtitle("Il12rb1") + ylab("expression level") + 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12, family = "helvetica", color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, size = 12, family = "helvetica", color = "black"),
        axis.text.y = element_text(size = 12, family = "helvetica", color = "black"),
        plot.title = element_text(family = "helvetica", color = "black", hjust = 0.5, size = 20)) + coord_cartesian(ylim = c(-0.1,4))

ggsave(filename = "~/MDC/Final_revisions/Fig3_A_violin.pdf", 
       plot = v1, 
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)



meta.df$Il12rb2 <- data9set_cleaned.SO@assays$RNA@data[rownames(data9set_cleaned.SO@assays$RNA@data) %in% "Il12rb2", ]

meta.df2 <- meta.df[!meta.df$large_cell_type %in% c("Rest"), ]

v2 <- ggplot(meta.df2, aes(x=large_cell_type, y=Il12rb2, color=large_cell_type))  + 
  geom_violin(trim=FALSE, width = 1) + scale_color_manual(values=c("#CC0000", "#FF6600", "#99CCCC", "#009999", "#99CC00","#FFCC66", "purple")) + 
  theme_classic() + 
  geom_jitter(position=position_jitter(0.2), size = 1.5, alpha = 1) + scale_y_continuous(expand = c(0,0),limits=c(0, 4)) +  
  ggtitle("Il12rb2") + ylab("expression level") + 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12, family = "helvetica", color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, size = 12, family = "helvetica", color = "black"),
        axis.text.y = element_text(size = 12, family = "helvetica", color = "black"),
        plot.title = element_text(family = "helvetica", color = "black", hjust = 0.5, size = 20)) + coord_cartesian(ylim = c(-0.1,4))
v2

ggsave(filename = "~/MDC/Final_revisions/Fig3_B_violin.pdf", 
       plot = v2, 
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)

meta.df$Il23r <- data9set_cleaned.SO@assays$RNA@data[rownames(data9set_cleaned.SO@assays$RNA@data) %in% "Il23r", ]

meta.df2 <- meta.df[!meta.df$large_cell_type %in% c("Rest"), ]

v3 <- ggplot(meta.df2, aes(x=large_cell_type, y=Il23r, color=large_cell_type))  + 
  geom_violin(trim=FALSE, width = 1) + scale_color_manual(values=c("#CC0000", "#FF6600", "#99CCCC", "#009999", "#99CC00","#FFCC66", "purple")) + 
  theme_classic() + 
  geom_jitter(position=position_jitter(0.2), size = 1.5, alpha = 1) + scale_y_continuous(expand = c(0,0),limits=c(0, 4)) +  
  ggtitle("Il23r") + ylab("expression level") + 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12, family = "helvetica", color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, size = 12, family = "helvetica", color = "black"),
        axis.text.y = element_text(size = 12, family = "helvetica", color = "black"),
        plot.title = element_text(family = "helvetica", color = "black", hjust = 0.5, size = 20)) + coord_cartesian(ylim = c(-0.1,4))
v3
ggsave(filename = "~/MDC/Final_revisions/Fig3_C_violin.pdf", 
       plot = v3, 
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)


