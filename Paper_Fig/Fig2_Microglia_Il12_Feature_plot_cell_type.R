
library(gridExtra)
library(dplyr)
library(ggplot2)
library(tidyr)
library(Seurat, lib.loc = "/data/rajewsky/home/skim/R/")
library(ggrepel)
library(cowplot)


load(file="/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200210_figures_for_meeting_report_cleaned_Seurat_object_res0.8.Robj")

# New assigned cell type

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



tmp <- data9set_cleaned.SO@reductions$umap@cell.embeddings

meta.df$UMAP1 <- tmp[, 1]
meta.df$UMAP2 <- tmp[, 2]

meta.df <- meta.df %>% mutate(large_cell_type = case_when(cell_type %in% c("Dentate_Gyrus", "CA1", "CA2_3","subiculum","Unidentified_Neurons") ~ "Excitatory Neurons",
                                                          cell_type %in% c("Inhibitory_Neurons") ~ "Inhibitory Neurons",
                                                          cell_type %in% c("MOL", "MFOL","NFOL") ~ "Oligodendrocytes",
                                                          cell_type %in% c("Astrocytes") ~ "Astrocytes",
                                                          cell_type %in% c("Microglia") ~ "Microglia",
                                                          cell_type %in% c("OPC") ~ "OPC",
                                                          cell_type %in% c("Vascular", "VLMC","Choroid","Fibroblast","Cajal","Pericyte","Macrophage") ~ "Rest"))


meta.df$large_cell_type <- factor(meta.df$large_cell_type, levels = c("Excitatory Neurons", "Inhibitory Neurons", "Oligodendrocytes", "OPC", "Microglia", "Astrocytes", "Rest"))


g1 <- ggplot(meta.df, aes(x = UMAP1, y = UMAP2, color = large_cell_type)) + 
  geom_point() + 
  theme_cowplot() + theme(legend.position = "None") + 
  scale_color_manual(values=c("#CC0000", "#FF6600", "#99CCCC", "#009999", "#99CC00","#FFCC66","purple")) + 
  annotate(geom="text", x=-10, y=0, label="Excitatory Neurons",
             color="black", size = 6, fontface = 2) +
  annotate(geom="text", x=-10, y=-17, label="Inhibitory Neurons",
           color="black", size = 6, fontface = 2) +
  annotate(geom="text", x=17, y= 2, label="Oligodendrocytes",
           color="black", size = 6, fontface = 2) +
  annotate(geom="text", x=21, y=-18, label="OPC",
           color="black", size = 6, fontface = 2) +
  annotate(geom="text", x= 1, y=23, label="Microglia",
           color="black", size = 6, fontface = 2) +
  annotate(geom="text", x=7, y=-25, label="Astrocytes",
           color="black", size = 6, fontface = 2)  +
  annotate(geom="text", x=14, y=14, label="Rest",
           color="black", size = 6, fontface = 2)



ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_Il12_Feature_plot_cell_type_UMAP.pdf",
       plot = g1,
       scale = 1, width = 10, height = 10, units = "in", device = cairo_pdf,
       dpi = 300)




