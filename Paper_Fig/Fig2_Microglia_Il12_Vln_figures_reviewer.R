

library(Seurat, lib.loc = "/data/rajewsky/home/skim/R/")
library(ggplot2)
library(dplyr)
library(tidyverse)
source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

# In this plot, I remove zero counts and only retain expressing cell and show by violin plot

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


data9set_cleaned.SO$cell_type <- meta.df$cell_type


# give levels to large cell type
data9set_cleaned.SO$cell_type <- factor(data9set_cleaned.SO$cell_type, levels = c("Dentate_Gyrus", "CA1", "CA2_3", "subiculum",  "Unidentified_Neurons", 
                                                                                  "Inhibitory_Neurons", "MOL", "MFOL", "OPC", "NFOL", "Microglia","Astrocytes",
                                                                                  "Vascular", "VLMC", "Choroid", "Fibroblast", "Cajal", "Pericyte", "Macrophage"))



meta.df <- meta.df %>% mutate(large_cell_type = case_when(cell_type %in% c("Dentate_Gyrus", "CA1", "CA2_3","subiculum","Unidentified_Neurons") ~ "Excitatory Neurons",
                                                          cell_type %in% c("Inhibitory_Neurons") ~ "Inhibitory Neurons",
                                                          cell_type %in% c("MOL", "MFOL","NFOL") ~ "Oligodendrocytes",
                                                          cell_type %in% c("Astrocytes") ~ "Astrocytes",
                                                          cell_type %in% c("Microglia") ~ "Microglia",
                                                          cell_type %in% c("OPC") ~ "OPC",
                                                          cell_type %in% c("Vascular", "VLMC","Choroid","Fibroblast","Cajal","Pericyte","Macrophage") ~ "Rest"))


meta.df$large_cell_type <- factor(meta.df$large_cell_type, levels = c("Excitatory Neurons", "Inhibitory Neurons", "Oligodendrocytes", "OPC", "Microglia", "Astrocytes", "Rest"))


#######################
# Il12b
#######################

#Il12b2 <- Il12b[Il12b$Il12b != 0, ]
meta.df$Il12b <- data9set_cleaned.SO@assays$RNA@data[rownames(data9set_cleaned.SO@assays$RNA@data) %in% "Il12b", ]

v1 <- ggplot(meta.df, aes(x=large_cell_type, y=Il12b, color=large_cell_type))  + 
  geom_violin(trim=FALSE, width = 1) + scale_color_manual(values=c("#CC0000", "#FF6600", "#99CCCC", "#009999", "#99CC00","#FFCC66", "purple")) + 
  theme_classic() + 
  geom_jitter(position=position_jitter(0.2), size = 1.5, alpha = 0.5) + scale_y_continuous(expand = c(0,0),limits=c(0, 4)) +  
  ggtitle("Il12b") + ylab("expression level") + 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12, family = "helvetica", color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, size = 12, family = "helvetica", color = "black"),
        axis.text.y = element_text(size = 12, family = "helvetica", color = "black"),
        plot.title = element_text(family = "helvetica", color = "black", hjust = 0.5, size = 20)) + coord_cartesian(ylim = c(-0.1,4))
v1

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_Il12_Vln_figures_reviewer_Il12b.pdf", 
       plot = v1, 
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)

#######################
# Il12a
#######################
meta.df$Il12a <- data9set_cleaned.SO@assays$RNA@data[rownames(data9set_cleaned.SO@assays$RNA@data) %in% "Il12a", ]

v1 <- ggplot(meta.df, aes(x=large_cell_type, y=Il12a, color=large_cell_type))  + 
  geom_violin(trim=FALSE, width = 1) + scale_color_manual(values=c("#CC0000", "#FF6600", "#99CCCC", "#009999", "#99CC00","#FFCC66", "purple")) + 
  theme_classic() + 
  geom_jitter(position=position_jitter(0.2), size = 1.5, alpha = 0.5) + scale_y_continuous(expand = c(0,0),limits=c(0, 4)) +  
  ggtitle("Il12a") + ylab("expression level") + 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12, family = "helvetica", color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, size = 12, family = "helvetica", color = "black"),
        axis.text.y = element_text(size = 12, family = "helvetica", color = "black"),
        plot.title = element_text(family = "helvetica", color = "black", hjust = 0.5, size = 20)) + coord_cartesian(ylim = c(-0.1,4))
v1

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_Il12_Vln_figures_reviewer_Il12a.pdf", 
       plot = v1, 
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)

#######################
# Il23a
#######################
meta.df$Il23a <- data9set_cleaned.SO@assays$RNA@data[rownames(data9set_cleaned.SO@assays$RNA@data) %in% "Il23a", ]

v1 <- ggplot(meta.df, aes(x=large_cell_type, y=Il23a, color=large_cell_type))  + 
  geom_violin(trim=FALSE, width = 1) + scale_color_manual(values=c("#CC0000", "#FF6600", "#99CCCC", "#009999", "#99CC00","#FFCC66", "purple")) + 
  theme_classic() + 
  geom_jitter(position=position_jitter(0.2), size = 1.5, alpha = 0.5) + scale_y_continuous(expand = c(0,0),limits=c(0, 4)) +  
  ggtitle("Il23a") + ylab("expression level") + 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12, family = "helvetica", color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, size = 12, family = "helvetica", color = "black"),
        axis.text.y = element_text(size = 12, family = "helvetica", color = "black"),
        plot.title = element_text(family = "helvetica", color = "black", hjust = 0.5, size = 20)) + coord_cartesian(ylim = c(-0.1,4))
v1

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_Il12_Vln_figures_reviewer_Il23a.pdf", 
       plot = v1, 
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)

#######################
# Il12rb1
#######################

meta.df$Il12rb1 <- data9set_cleaned.SO@assays$RNA@data[rownames(data9set_cleaned.SO@assays$RNA@data) %in% "Il12rb1", ]

v1 <- ggplot(meta.df, aes(x=large_cell_type, y=Il12rb1, color=large_cell_type))  + 
  geom_violin(trim=FALSE, width = 1) + scale_color_manual(values=c("#CC0000", "#FF6600", "#99CCCC", "#009999", "#99CC00","#FFCC66", "purple")) + 
  theme_classic() + 
  geom_jitter(position=position_jitter(0.2), size = 1.5, alpha = 0.5) + scale_y_continuous(expand = c(0,0),limits=c(0, 4)) +  
  ggtitle("Il12rb1") + ylab("expression level") + 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12, family = "helvetica", color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, size = 12, family = "helvetica", color = "black"),
        axis.text.y = element_text(size = 12, family = "helvetica", color = "black"),
        plot.title = element_text(family = "helvetica", color = "black", hjust = 0.5, size = 20)) + coord_cartesian(ylim = c(-0.1,4))
v1

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_Il12_Vln_figures_reviewer_Il12rb1.pdf", 
       plot = v1, 
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)

#######################
# Il12rb2
#######################

meta.df$Il12rb2 <- data9set_cleaned.SO@assays$RNA@data[rownames(data9set_cleaned.SO@assays$RNA@data) %in% "Il12rb2", ]

v1 <- ggplot(meta.df, aes(x=large_cell_type, y=Il12rb2, color=large_cell_type))  + 
  geom_violin(trim=FALSE, width = 1) + scale_color_manual(values=c("#CC0000", "#FF6600", "#99CCCC", "#009999", "#99CC00","#FFCC66", "purple")) + 
  theme_classic() + 
  geom_jitter(position=position_jitter(0.2), size = 1.5, alpha = 0.5) + scale_y_continuous(expand = c(0,0),limits=c(0, 4)) +  
  ggtitle("Il12rb2") + ylab("expression level") + 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12, family = "helvetica", color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, size = 12, family = "helvetica", color = "black"),
        axis.text.y = element_text(size = 12, family = "helvetica", color = "black"),
        plot.title = element_text(family = "helvetica", color = "black", hjust = 0.5, size = 20)) + coord_cartesian(ylim = c(-0.1,4))
v1

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_Il12_Vln_figures_reviewer_Il12rb2.pdf", 
       plot = v1, 
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)

#######################
# Il23r
#######################

meta.df$Il23r <- data9set_cleaned.SO@assays$RNA@data[rownames(data9set_cleaned.SO@assays$RNA@data) %in% "Il23r", ]

v1 <- ggplot(meta.df, aes(x=large_cell_type, y=Il23r, color=large_cell_type))  + 
  geom_violin(trim=FALSE, width = 1) + scale_color_manual(values=c("#CC0000", "#FF6600", "#99CCCC", "#009999", "#99CC00","#FFCC66", "purple")) + 
  theme_classic() + 
  geom_jitter(position=position_jitter(0.2), size = 1.5, alpha = 0.5) + scale_y_continuous(expand = c(0,0),limits=c(0, 4)) +  
  ggtitle("Il23r") + ylab("expression level") + 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12, family = "helvetica", color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, size = 12, family = "helvetica", color = "black"),
        axis.text.y = element_text(size = 12, family = "helvetica", color = "black"),
        plot.title = element_text(family = "helvetica", color = "black", hjust = 0.5, size = 20)) + coord_cartesian(ylim = c(-0.1,4)) 
v1

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_Il12_Vln_figures_reviewer_Il23r.pdf", 
       plot = v1, 
       scale = 1, width = 9, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)


##########################################
# Separate by genotype and all molecules
##########################################

meta.df2 <- meta.df %>% gather(molecule,  value, Il12b:Il23r)

meta.df2 <- meta.df2 %>% mutate(genotype = case_when(sample %in% c("Ctrl") ~ "WT", 
                                                     sample %in% c("AD") ~ "APPPS1",
                                                     sample %in% c("ADp40KO") ~ "APPPS1.il12b-/-"))


meta.df2$genotype <- factor(meta.df2$genotype, levels = c("WT","APPPS1","APPPS1.il12b-/-"))

v2 <- ggplot(meta.df2, aes(x=large_cell_type, y=value, color=large_cell_type))  + 
  geom_violin(trim=FALSE, width = 1) + scale_color_manual(values=c("#CC0000", "#FF6600", "#99CCCC", "#009999", "#99CC00","#FFCC66", "purple")) + 
  theme_classic() + 
  geom_jitter(position=position_jitter(0.2), size = 1.5, alpha = 0.5) + 
  #scale_y_continuous(expand = c(0,0),limits=c(0, 4)) + 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12, family = "helvetica", color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, size = 12, family = "helvetica", color = "black"),
        axis.text.y = element_blank(),
        #axis.text.y = element_text(size = 12, family = "helvetica", color = "black"),
        plot.title = element_text(family = "helvetica", color = "black", hjust = 0.5)) + coord_cartesian(ylim = c(-0.1,4)) +
  facet_grid(molecule ~  genotype)
v2

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_Il12_Vln_figures_reviewer_separate_genotype.pdf", 
       plot = v2, 
       scale = 1, width = 16, height = 12, units = "in", device = cairo_pdf,
       dpi = 300)


################################
# avg expression of cell per cell type
################################

# Il12b
plot.data <- meta.df %>%
  dplyr::select(sample, Il12b, cell_type = large_cell_type) %>%
  dplyr::group_by(cell_type, sample) %>% 
  dplyr::summarise(Il12b = mean(Il12b)) 

# add Il12a
plot.data2 <- meta.df %>%
  dplyr::select(sample, Il12a, cell_type = large_cell_type) %>%
  dplyr::group_by(cell_type, sample) %>% 
  dplyr::summarise(Il12a = mean(Il12a)) 
plot.data$Il12a <- plot.data2$Il12a

# add Il12rb1
plot.data2 <- meta.df %>%
  dplyr::select(sample, Il12rb1, cell_type = large_cell_type) %>%
  dplyr::group_by(cell_type, sample) %>% 
  dplyr::summarise(Il12rb1 = mean(Il12rb1)) 
plot.data$Il12rb1 <- plot.data2$Il12rb1

# add Il12rb2
plot.data2 <- meta.df %>%
  dplyr::select(sample, Il12rb2, cell_type = large_cell_type) %>%
  dplyr::group_by(cell_type, sample) %>% 
  dplyr::summarise(Il12rb2 = mean(Il12rb2)) 
plot.data$Il12rb2 <- plot.data2$Il12rb2

# add Il23a
plot.data2 <- meta.df %>%
  dplyr::select(sample, Il23a, cell_type = large_cell_type) %>%
  dplyr::group_by(cell_type, sample) %>% 
  dplyr::summarise(Il23a = mean(Il23a)) 
plot.data$Il23a <- plot.data2$Il23a

# add Il23r
plot.data2 <- meta.df %>%
  dplyr::select(sample, Il23r, cell_type = large_cell_type) %>%
  dplyr::group_by(cell_type, sample) %>% 
  dplyr::summarise(Il23r = mean(Il23r)) 
plot.data$Il23r <- plot.data2$Il23r

# preprocessing data frame for figure
expression.final <- plot.data %>% gather(molecule, expression, Il12b:Il23r)

################################
# percent of cell per cell type
################################


plot.data <- meta.df %>%
  dplyr::select(sample, Il12b, cell_type = large_cell_type) %>%
  dplyr::group_by(cell_type, sample) %>% 
  dplyr::summarise(count = n()) 

# add Il12b
plot.data2 <- meta.df %>%
  dplyr::select(sample, Il12b, cell_type = large_cell_type) %>%
  dplyr::group_by(cell_type, sample) %>% 
  dplyr::summarise(Il12b = sum(Il12b > 0, na.rm = TRUE)) 

plot.data$Il12b <- plot.data2$Il12b

# add Il12a
plot.data2 <- meta.df %>%
  dplyr::select(sample, Il12a, cell_type = large_cell_type) %>%
  dplyr::group_by(cell_type, sample) %>% 
  dplyr::summarise(Il12a = sum(Il12a > 0, na.rm = TRUE)) 

plot.data$Il12a <- plot.data2$Il12a

# add Il12rb1
plot.data2 <- meta.df %>%
  dplyr::select(sample, Il12rb1, cell_type = large_cell_type) %>%
  dplyr::group_by(cell_type, sample) %>% 
  dplyr::summarise(Il12rb1 = sum(Il12rb1 > 0, na.rm = TRUE)) 

plot.data$Il12rb1 <- plot.data2$Il12rb1

# add Il12rb2
plot.data2 <- meta.df %>%
  dplyr::select(sample, Il12rb2, cell_type = large_cell_type) %>%
  dplyr::group_by(cell_type, sample) %>% 
  dplyr::summarise(Il12rb2 = sum(Il12rb2 > 0, na.rm = TRUE)) 

plot.data$Il12rb2 <- plot.data2$Il12rb2

# add Il23a
plot.data2 <- meta.df %>%
  dplyr::select(sample, Il23a, cell_type = large_cell_type) %>%
  dplyr::group_by(cell_type, sample) %>% 
  dplyr::summarise(Il23a = sum(Il23a > 0, na.rm = TRUE)) 

plot.data$Il23a <- plot.data2$Il23a

# add Il23r
plot.data2 <- meta.df %>%
  dplyr::select(sample, Il23r, cell_type = large_cell_type) %>%
  dplyr::group_by(cell_type, sample) %>% 
  dplyr::summarise(Il23r = sum(Il23r > 0, na.rm = TRUE)) 

plot.data$Il23r <- plot.data2$Il23r

# preprocessing data frame for figure
plot.data3 <- apply(plot.data[, c(4:9)], 2, function(x) x/plot.data$count) %>% as.data.frame()
plot.final <- cbind(plot.data[, c(1,2)], plot.data3)
plot.final2 <- plot.final %>% gather(molecule, percent, Il12b:Il23r)

# add average expression 
plot.final2$expression <- expression.final$expression

plot.final2$molecule <- factor(plot.final2$molecule, levels = rev(c("Il12a", "Il12b", "Il12rb1", "Il12rb2","Il23a", "Il23r")))

g1 <- ggplot(plot.final2, 
       aes(x=sample, y=molecule,
           colour=expression,
           size=ifelse(percent==0, NA, percent))) +
  geom_point() + facet_grid(~  cell_type) + theme_cowplot() + 
  scale_x_discrete(position = "bottom", labels=c("Ctrl" = "WT", "AD" = "APPPS1",
                                              "ADp40KO" = "APPPS1.il12b-/-")) + 
  scale_color_viridis(limits = c(0,0.1)) +
  theme(axis.title.x = element_blank(), 
        axis.title.y =  element_blank(),
        axis.text.y = element_text(size = 11, family = "Helvetica", color = "black"),
        axis.text.x = element_text(angle=90, vjust = 0.5, size = 11, family = "Helvetica", color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_text(angle=0, family = "Helvetica", color = "black"),
        strip.text.y.left = element_text(angle=0, size = 11, family = "Helvetica", color = "black"),
        #strip.placement = "outside",
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(), 
        #panel.spacing = unit(c(0,0,0, -2), "lines"),
        #panel.spacing.x = unit(-2, "cm"),
        #plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  labs(size="percent",col="expression")

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_Il12_Vln_figures_reviewer_dot_plot.pdf", 
       plot = g1, 
       scale = 1, width = 14, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

##################
# control dot plot
##################
meta.df$Gapdh <- data9set_cleaned.SO@assays$RNA@data[rownames(data9set_cleaned.SO@assays$RNA@data) %in% "Gapdh", ]
meta.df$Actb <- data9set_cleaned.SO@assays$RNA@data[rownames(data9set_cleaned.SO@assays$RNA@data) %in% "Actb", ]
meta.df$Pgk1 <- data9set_cleaned.SO@assays$RNA@data[rownames(data9set_cleaned.SO@assays$RNA@data) %in% "Pgk1", ]

################################
# avg expression of cell per cell type
################################

# Gapdh
plot.data <- meta.df %>%
  dplyr::select(sample, Gapdh, cell_type = large_cell_type) %>%
  dplyr::group_by(cell_type, sample) %>% 
  dplyr::summarise(Gapdh = mean(Gapdh)) 

# add Actb
plot.data2 <- meta.df %>%
  dplyr::select(sample, Actb, cell_type = large_cell_type) %>%
  dplyr::group_by(cell_type, sample) %>% 
  dplyr::summarise(Actb = mean(Actb)) 
plot.data$Actb <- plot.data2$Actb

# add Pgk1
plot.data2 <- meta.df %>%
  dplyr::select(sample, Pgk1, cell_type = large_cell_type) %>%
  dplyr::group_by(cell_type, sample) %>% 
  dplyr::summarise(Pgk1 = mean(Pgk1)) 
plot.data$Pgk1 <- plot.data2$Pgk1

# preprocessing data frame for figure
expression.final <- plot.data %>% gather(molecule, expression, Gapdh:Pgk1)

################################
# percent of cell per cell type
################################


plot.data <- meta.df %>%
  dplyr::select(sample, Gapdh, cell_type = large_cell_type) %>%
  dplyr::group_by(cell_type, sample) %>% 
  dplyr::summarise(count = n()) 

# add Gapdh
plot.data2 <- meta.df %>%
  dplyr::select(sample, Gapdh, cell_type = large_cell_type) %>%
  dplyr::group_by(cell_type, sample) %>% 
  dplyr::summarise(Gapdh = sum(Gapdh > 0, na.rm = TRUE)) 

plot.data$Gapdh <- plot.data2$Gapdh

# add Actb
plot.data2 <- meta.df %>%
  dplyr::select(sample, Actb, cell_type = large_cell_type) %>%
  dplyr::group_by(cell_type, sample) %>% 
  dplyr::summarise(Actb = sum(Actb > 0, na.rm = TRUE)) 

plot.data$Actb <- plot.data2$Actb

# add Pgk1
plot.data2 <- meta.df %>%
  dplyr::select(sample, Pgk1, cell_type = large_cell_type) %>%
  dplyr::group_by(cell_type, sample) %>% 
  dplyr::summarise(Pgk1 = sum(Pgk1 > 0, na.rm = TRUE)) 

plot.data$Pgk1 <- plot.data2$Pgk1

# preprocessing data frame for figure
plot.data3 <- apply(plot.data[, c(4:6)], 2, function(x) x/plot.data$count) %>% as.data.frame()
plot.final <- cbind(plot.data[, c(1,2)], plot.data3)
plot.final2 <- plot.final %>% gather(molecule, percent, Gapdh:Pgk1)

# add average expression 
plot.final2$expression <- expression.final$expression

plot.final2$molecule <- factor(plot.final2$molecule, levels = rev(c("Gapdh", "Actb", "Pgk1")))

g1 <- ggplot(plot.final2, 
             aes(x=sample, y=molecule,
                 colour=expression,
                 size=ifelse(percent==0, NA, percent))) +
  geom_point() + facet_grid(~  cell_type) + theme_cowplot() + 
  scale_x_discrete(position = "bottom", labels=c("Ctrl" = "WT", "AD" = "APPPS1",
                                                 "ADp40KO" = "APPPS1.il12b-/-")) + 
  scale_color_viridis(limits = c(0,1.1)) +
  theme(axis.title.x = element_blank(), 
        axis.title.y =  element_blank(),
        axis.text.y = element_text(size = 11, family = "Helvetica", color = "black"),
        axis.text.x = element_text(angle=90, vjust = 0.5, size = 11, family = "Helvetica", color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_text(angle=0, family = "Helvetica", color = "black"),
        strip.text.y.left = element_text(angle=0, size = 11, family = "Helvetica", color = "black"),
        #strip.placement = "outside",
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(), 
        #panel.spacing = unit(c(0,0,0, -2), "lines"),
        #panel.spacing.x = unit(-2, "cm"),
        #plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  labs(size="percent",col="expression")
g1
ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_Il12_Vln_figures_reviewer_dot_plot_ctrl.pdf", 
       plot = g1, 
       scale = 1, width = 14, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)


