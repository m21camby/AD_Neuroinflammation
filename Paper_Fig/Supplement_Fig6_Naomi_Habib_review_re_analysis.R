
Naomi.meta <- Naomi.SO@meta.data



DimPlot(object = Naomi.SO, reduction = "umap", label = TRUE)

############################
# Check marker and cell type
############################

# 23 MG
# 7, 13, 18 IN
# 17 Oligo
# 14, 26 OPC
# 2, 8, 19, 22, 27 AS
# 0,1,21, 6, 11, 25, 4,5,20, 10,24, 9, 3, 12, 16 EX
# 22, 15 Rest


# IN
FeaturePlot(Naomi.SO, features = "Gad2", order = TRUE)
# Oligo
FeaturePlot(Naomi.SO, features = "Mog", order = TRUE)
FeaturePlot(Naomi.SO, features = "Mbp", order = FALSE)
# MG
FeaturePlot(Naomi.SO, features = "Ly86", order = TRUE)
# MG
FeaturePlot(Naomi.SO, features = "Hexb", order = TRUE)
# OPC
FeaturePlot(Naomi.SO, features = "Pdgfra", order = TRUE)
# AS
FeaturePlot(Naomi.SO, features = "Gfap", order = TRUE)
FeaturePlot(Naomi.SO, features = "Slc1a3", order = TRUE)
FeaturePlot(Naomi.SO, features = "Rmst", order = TRUE)

# Rest
FeaturePlot(Naomi.SO, features = "Flt1", order = TRUE)
FeaturePlot(Naomi.SO, features = "Cped1", order = TRUE)

Naomi.meta <- Naomi.meta %>% mutate(cell_type = case_when(seurat_clusters %in% c(23) ~ "MG",
                                                            seurat_clusters %in% c(7, 13, 18) ~ "IN",
                                                            seurat_clusters %in% c(17) ~ "Oligo",
                                                            seurat_clusters %in% c(14, 26) ~ "OPC",
                                                            seurat_clusters %in% c(2, 8, 19, 22, 27) ~ "AS",
                                                            seurat_clusters %in% c(0,1,21, 6, 11, 25, 4,5,20, 10,24, 9, 3, 12, 16) ~ "EX",
                                                            seurat_clusters %in% c(22, 15) ~ "Rest"))

table(Naomi.meta$cell_type)

tmp <- Naomi.SO@reductions$umap@cell.embeddings

Naomi.meta$UMAP1 <- tmp[, 1]
Naomi.meta$UMAP2 <- tmp[, 2]

Naomi.meta$cell_type <- factor(Naomi.meta$cell_type, levels = c("EX", "IN", "Oligo", "OPC", "MG", "AS", "Rest"))

g1 <- ggplot(Naomi.meta, aes(x = UMAP1, y = UMAP2, color = cell_type)) +
  geom_point() +
  theme_cowplot() + theme(legend.position = "None") +
  scale_color_manual(values=c("#CC0000", "#FF6600", "#99CCCC", "#009999", "#99CC00","#FFCC66","purple")) +
  annotate(geom="text", x=-7, y=0, label="Excitatory Neurons",
           color="black", size = 5, fontface = 1) +
  annotate(geom="text", x=7, y=13, label="Inhibitory Neurons",
           color="black", size = 5, fontface = 1) +
  annotate(geom="text", x=19, y= 12, label="Oligodendrocytes",
           color="black", size = 5, fontface = 1) +
  annotate(geom="text", x=12, y=20, label="OPC",
           color="black", size = 5, fontface = 1) +
  annotate(geom="text", x= 12, y=5.5, label="Microglia",
           color="black", size = 5, fontface = 1) +
  annotate(geom="text", x=18, y=-3.5, label="Astrocytes",
           color="black", size = 5, fontface = 1)  +
  annotate(geom="text", x=4, y=25, label="Rest",
           color="black", size = 5, fontface = 1) +
  annotate(geom="text", x=9, y=-20, label="55,317 cells",
           color="black", size = 5, fontface = 1) +
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
g1

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig6_Naomi_Habib_review_re_analysis_UMAP.pdf",
       plot = g1,
       scale = 1, width = 10, height = 10, units = "in", device = cairo_pdf,
       dpi = 300)

# Il12b is not exist
# IL12A
f1 <- FeaturePlot(Naomi.SO, features = "Il12a",order = TRUE, pt.size = 0.5)  +
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig6_Naomi_Habib_review_re_analysis_UMAP_Il12a.pdf",
       plot = f1,
       scale = 1, width = 8, height = 7, units = "in", device = cairo_pdf,
       dpi = 300)

# Il12rb1
f1 <- FeaturePlot(Naomi.SO, features = "Il12rb1",order = TRUE, pt.size = 0.5)  +
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig6_Naomi_Habib_review_re_analysis_UMAP_Il12rb1.pdf",
       plot = f1,
       scale = 1, width = 8, height = 7, units = "in", device = cairo_pdf,
       dpi = 300)

# Il12rb2
f1 <- FeaturePlot(Naomi.SO, features = "Il12rb2",order = TRUE, pt.size = 0.5)  +
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig6_Naomi_Habib_review_re_analysis_UMAP_Il12rb2.pdf",
       plot = f1,
       scale = 1, width = 8, height = 7, units = "in", device = cairo_pdf,
       dpi = 300)

# Il23a is not exist
# Il23r
f1 <- FeaturePlot(Naomi.SO, features = "Il23r",order = TRUE, pt.size = 0.5)  +
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig6_Naomi_Habib_review_re_analysis_UMAP_Il23r.pdf",
       plot = f1,
       scale = 1, width = 8, height = 7, units = "in", device = cairo_pdf,
       dpi = 300)

# Gad2
f1 <- FeaturePlot(Naomi.SO, features = "Gad2",order = TRUE, pt.size = 0.5)  +
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig6_Naomi_Habib_review_re_analysis_UMAP_Gad2_IN.pdf",
       plot = f1,
       scale = 1, width = 8, height = 7, units = "in", device = cairo_pdf,
       dpi = 300)

# Mog
f1 <- FeaturePlot(Naomi.SO, features = "Mog",order = FALSE, pt.size = 0.5)  +
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig6_Naomi_Habib_review_re_analysis_UMAP_Mog_Oligo.pdf",
       plot = f1,
       scale = 1, width = 8, height = 7, units = "in", device = cairo_pdf,
       dpi = 300)

# Pdgfra
f1 <- FeaturePlot(Naomi.SO, features = "Pdgfra",order = TRUE, pt.size = 0.5)  +
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig6_Naomi_Habib_review_re_analysis_UMAP_Pdgfra_OPC.pdf",
       plot = f1,
       scale = 1, width = 8, height = 7, units = "in", device = cairo_pdf,
       dpi = 300)

# Ly86
f1 <- FeaturePlot(Naomi.SO, features = "Ly86",order = TRUE, pt.size = 0.5)  +
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig6_Naomi_Habib_review_re_analysis_UMAP_Ly86_MG.pdf",
       plot = f1,
       scale = 1, width = 8, height = 7, units = "in", device = cairo_pdf,
       dpi = 300)

# Slc1a3 AS
f1 <- FeaturePlot(Naomi.SO, features = "Slc1a3",order = FALSE, pt.size = 0.5)  +
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig6_Naomi_Habib_review_re_analysis_UMAP_Slc1a3_AS.pdf",
       plot = f1,
       scale = 1, width = 8, height = 7, units = "in", device = cairo_pdf,
       dpi = 300)

# Kcnip4 
f1 <- FeaturePlot(Naomi.SO, features = "Kcnip4",order = FALSE, pt.size = 0.5)  +
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig6_Naomi_Habib_review_re_analysis_UMAP_Kcnip4_EX.pdf",
       plot = f1,
       scale = 1, width = 8, height = 7, units = "in", device = cairo_pdf,
       dpi = 300)


####################
# dot plot
####################

Naomi.meta$Il12a <- Naomi.SO@assays$RNA@data[rownames(Naomi.SO@assays$RNA@data) %in% "Il12a", ]
Naomi.meta$Il12rb1 <- Naomi.SO@assays$RNA@data[rownames(Naomi.SO@assays$RNA@data) %in% "Il12rb1", ]
Naomi.meta$Il12rb2 <- Naomi.SO@assays$RNA@data[rownames(Naomi.SO@assays$RNA@data) %in% "Il12rb2", ]
Naomi.meta$Il23r <- Naomi.SO@assays$RNA@data[rownames(Naomi.SO@assays$RNA@data) %in% "Il23r", ]



################################
# avg expression of cell per cell type
################################

# IL12a
plot.data <- Naomi.meta %>%
  dplyr::select(exp, Il12a, cell_type = cell_type) %>%
  dplyr::group_by(cell_type, exp) %>%
  dplyr::summarise(Il12a = mean(Il12a))

# add IL12rb1
plot.data2 <- Naomi.meta %>%
  dplyr::select(exp, Il12rb1, cell_type = cell_type) %>%
  dplyr::group_by(cell_type, exp) %>%
  dplyr::summarise(Il12rb1 = mean(Il12rb1))
plot.data$Il12rb1 <- plot.data2$Il12rb1

# add IL12rb2
plot.data2 <- Naomi.meta %>%
  dplyr::select(exp, Il12rb2, cell_type = cell_type) %>%
  dplyr::group_by(cell_type, exp) %>%
  dplyr::summarise(Il12rb2 = mean(Il12rb2))
plot.data$Il12rb2 <- plot.data2$Il12rb2

# add Il23r
plot.data2 <- Naomi.meta %>%
  dplyr::select(exp, Il23r, cell_type = cell_type) %>%
  dplyr::group_by(cell_type, exp) %>%
  dplyr::summarise(Il23r = mean(Il23r))
plot.data$Il23r <- plot.data2$Il23r

expression.final <- plot.data %>% tidyr::gather(molecule, expression, Il12a:Il23r)

################################
# percent of cell per cell type
################################


plot.data <- Naomi.meta %>%
  dplyr::select(exp, Il12a, cell_type = cell_type) %>%
  dplyr::group_by(cell_type, exp) %>%
  dplyr::summarise(count = n())

# add Il12a
plot.data2 <- Naomi.meta %>%
  dplyr::select(exp, Il12a, cell_type = cell_type) %>%
  dplyr::group_by(cell_type, exp) %>%
  dplyr::summarise(Il12a = sum(Il12a > 0, na.rm = TRUE))
plot.data$Il12a <- plot.data2$Il12a

# add Il12rb1
plot.data2 <- Naomi.meta %>%
  dplyr::select(exp, Il12rb1, cell_type = cell_type) %>%
  dplyr::group_by(cell_type, exp) %>%
  dplyr::summarise(Il12rb1 = sum(Il12rb1 > 0, na.rm = TRUE))
plot.data$Il12rb1 <- plot.data2$Il12rb1

# add Il12rb2
plot.data2 <- Naomi.meta %>%
  dplyr::select(exp, Il12rb2, cell_type = cell_type) %>%
  dplyr::group_by(cell_type, exp) %>%
  dplyr::summarise(Il12rb2 = sum(Il12rb2 > 0, na.rm = TRUE))
plot.data$Il12rb2 <- plot.data2$Il12rb2

# add Il23r
plot.data2 <- Naomi.meta %>%
  dplyr::select(exp, Il23r, cell_type = cell_type) %>%
  dplyr::group_by(cell_type, exp) %>%
  dplyr::summarise(Il23r = sum(Il23r > 0, na.rm = TRUE))
plot.data$Il23r <- plot.data2$Il23r

# preprocessing data frame for figure
plot.data3 <- apply(plot.data[, c(4:7)], 2, function(x) x/plot.data$count) %>% as.data.frame()
plot.final <- cbind(plot.data[, c(1,2)], plot.data3)
plot.final2 <- plot.final %>% tidyr::gather(molecule, percent, Il12a:Il23r)

# add average expression 
plot.final2$expression <- expression.final$expression

plot.final2$molecule <- factor(plot.final2$molecule, levels = rev(c("Il12a",  "Il12rb1", "Il12rb2","Il23r")))

plot.final2$exp <- factor(plot.final2$exp, levels = c("WT","AD", "Untreated"))

g1 <- ggplot(plot.final2,
             aes(x=exp, y=molecule,
                 colour=expression,
                 size=ifelse(percent==0, NA, percent))) +
  geom_point() + facet_grid(~  cell_type) + theme_cowplot() +
  scale_x_discrete(position = "bottom", labels=c("WT" = "WT", "AD" = "AD", "Untreated" = "Untreated")) +
  scale_color_viridis(limits = c(0,0.04)) +
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

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig6_Naomi_Habib_review_re_analysis_UMAP_dot_plot.pdf",
       plot = g1,
       scale = 1, width = 14, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)
