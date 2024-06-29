#! /bin/env RScript
# written by SJK at 8. Mar. 2020

.libPaths(c("/data/murphy/shared_libs/R", "/usr/local/lib/R/site-librar","/usr/local/lib/R/library","/usr/lib64/R/library","/usr/share/R/library", .libPaths() ))

library(Seurat)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)
library(tidyr)

source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

load("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200117_Clustering_DE_by_MAST_Seurat_object.Robj")
data9set.SO@meta.data$sample <- ifelse((data9set.SO$gemgroup == 1 | data9set.SO$gemgroup == 4 | data9set.SO$gemgroup == 7), "Ctrl", 
                                       ifelse((data9set.SO$gemgroup == 2 | data9set.SO$gemgroup == 5 | data9set.SO$gemgroup == 8), "ADp40KO", "AD"))
data9set.SO@meta.data$sample <- factor(data9set.SO@meta.data$sample, levels = c("Ctrl", "AD", "ADp40KO"))

data9set_cleaned.SO  <- subset(data9set.SO, subset = seurat_clusters %in% c(seq(0,24), seq(26,30), 32, seq(35,40)))

data9set_cleaned.SO <- RenameIdents(data9set_cleaned.SO, '0' = "Dentate Gyrus", '1' = "Oligo", '2' = "CA1 Neurons", '3' = "Microglia", 
                                    '4' = "Astrocytes", '5' = "Oligo", '6' = "Oligo", '7' = "Inhibitory Interneurons", '8' = "Microglia", 
                                    '9' = "CA2/CA3 Neurons", '10' = "Dentate Gyrus", '11' = "Subiculum", '12' = "OPC", 
                                    '13' = "CA1 Neurons", '14' = "Subiculum", '15' = "Subiculum", 
                                    '16' = "Inhibitory Interneurons", '17' = "Subiculum", '18' = "CA1 Neurons",
                                    '19' = "Subiculum", '20' = "Subiculum", '21' = "CA2/CA3 Neurons", '22' = "Inhibitory Interneurons",
                                    '23' = "Fibroblast", '24' = "Neurons", '26' = "Choroid Plexus", 
                                    '27' = "CA2/CA3 Neurons", '28' = "Cajal Retzius", '29' = "Neurons", '30' = "Inhibitory Interneurons",
                                    '32' = "Neurons", 
                                    '35' = "Pericytes", '36' = "Vascular Endothelial",'37' = "Macrophage", '38' = "OPC",
                                    '39' = "VLMC", '40' = "Neurons")

data9set_cleaned.SO@active.ident <- factor(data9set_cleaned.SO@active.ident, 
                                           levels = c("Dentate Gyrus", "CA1 Neurons", "CA2/CA3 Neurons", "Subiculum",
                                                      "Inhibitory Interneurons", "Neurons", "Cajal Retzius", "Choroid Plexus", "Astrocytes",
                                                      "Microglia","Macrophage", "Oligo", "OPC", "Fibroblast", 
                                                      "Vascular Endothelial", "VLMC", "Pericytes"))



data9set_cleaned.SO$cell_type <- data9set_cleaned.SO@active.ident


d1 <- DimPlot(data9set_cleaned.SO, reduction = "umap", label = TRUE, label.size = 5) + theme(legend.position = "none")

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Marker_genes_Heatmap_region_UMAP.png", 
       plot = d1, 
       scale = 1, width = 10, height = 10, units = "in", device = "png",
       dpi = 300)


data9set_cleaned.SO.subset <- subset(data9set_cleaned.SO, downsample = 100)

###################################
# 2 heatmap by cell type region
###################################
data9set_scaled.df3 <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200308_Scaled_data_10000_Cell_type_by_region.csv", row.names = 1, check.names = FALSE)


# 3-2. extract meta data
data9set_cleaned.SO.data <- data9set_cleaned.SO@meta.data
data9set_cleaned.SO.data$Cell <- rownames(data9set_cleaned.SO.data)
data9set_cleaned.SO.data.data_sub <- data9set_cleaned.SO.data[,c(9,10)]

# extract marker genes
GM.df <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200210_FindMarkers_each_cluster_data9set_marker_top10.csv")
GM.df <- as.data.frame(GM.df)
GM.df2 <- GM.df[c(seq(1,5), seq(1,5) + 100, seq(1,5) + 20, seq(1,5) + 130, seq(1,5) + 180, 
                  seq(1,5) + 90, seq(1,5) + 210, seq(1,5) + 270, seq(1,5) + 110, 
                  seq(1,5) + 140, seq(1,5) + 150,seq(1,5) + 170, seq(1,5) + 190, seq(1,5) + 200,
                  seq(1,5) + 70, seq(1,5) + 160, seq(1,5) + 220, seq(1,5) + 300,
                  seq(1,5) + 240, seq(1,5) + 290, seq(1,5) + 320, seq(1,5) + 400, seq(1,5) + 280, 
                  seq(1,5) + 40, seq(1,5) + 30, seq(1,5) + 80, seq(1,5) + 370,
                  seq(1,5) + 10, seq(1,5) + 50, seq(1,5) + 60, seq(1,5) + 120, seq(1,5) + 380,
                  seq(1,5) + 230, seq(1,5) + 360, seq(1,5) + 390, seq(1,5) + 350, seq(1,5) + 280, seq(1,5) + 260), ]

GM.df2$gene <- as.character(GM.df2$gene)
# remove Cck, Atp6v0c, Cst3, Tmsb4x, Gm42418 which are not exist in scale.data
# 1,2,6,7,8, 11,12,16,17,21, 27,28,31,32,36,  41,42,46,51,57, 71,72,76,81,Gad2, 91,92,102,104,106,
# 111:115, 186:190, 116:120, 121,122,126,129,Hexb, 131:135, 136,137,138,141,146, 151,152,153,157,158,
# 161:165, 166:170, 171:175, 176:180

GM.df3 <- GM.df2[c(1,2,6,7,8, 11,12,16,17,21, 27,28,33,34,36,  41,42,46,52,57, 71,72,76,81),] 
GM.df4 <- GM.df2[c(91,102,104,106,107,
         111:115, 186:190, 116:120, 121,122,126,129),] 
GM.df5 <- GM.df2[c(131:135, 136,137,138,146,150, 151,152,155,157,158,
         161:165, 166:170, 171:175, 176:180),]


genes <- c(GM.df3$gene, "Gad2", GM.df4$gene, "Hexb", GM.df5$gene)

data9set_scaled.df3_sub <- data9set_scaled.df3[genes, ]


data9set_scaled.df3_sub$gene <- rownames(data9set_scaled.df3_sub)

tail(colnames(data9set_scaled.df3_sub))

data9set_scaled.df3_sub.long <- gather(data = data9set_scaled.df3_sub, z_score, key = Cell, "AATTTCCAGCATGCAG-1":"TTTACGTCAATTTCGG-9")


data9set_scaled.df3_sub.long2 <- left_join(data9set_scaled.df3_sub.long, data9set_cleaned.SO.data.data_sub, by = "Cell")
data9set_scaled.df3_sub.long2$gene <- factor(data9set_scaled.df3_sub.long2$gene, levels = rev(genes))

# 3-5. plot
g1 <- ggplot(data = data9set_scaled.df3_sub.long2, mapping = aes(x = cell_type, y = gene, fill = z_score)) +
  geom_tile() +  theme_classic() + theme(axis.title.y = element_blank(),
                                         axis.title.x = element_blank(), 
                                         axis.text.x = element_text(angle = 90, size = 12, color = "black", vjust = 0.1), 
                                         axis.text.y = element_text(size = 11, color = "black"),
                                         axis.line = element_line(color = "white"),
                                         axis.ticks = element_line(color = "white"),
                                         legend.title = element_blank(),
                                         legend.text = element_text(size = 11, color = "black")) + 
  scale_fill_viridis_c(limits=c(-3, 3), na.value = "#FDE725FF", labels = c("low", " ", " "," ", " "," ","max"))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Marker_genes_Heatmap_region_cell_type.png", 
       plot = g1, 
       scale = 1, width = 10, height = 16, units = "in", device = "png",
       dpi = 300)

# save doHeatmap
d1 <- DoHeatmap(data9set_cleaned.SO.subset, features = genes, angle = 90, group.bar.height = 0) + 
  scale_fill_viridis() + 
  theme(legend.position = "none")

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Marker_genes_doHeatmap_region_cell_type.png", 
       plot = d1, 
       scale = 1, width = 10, height = 16, units = "in", device = "png",
       dpi = 300)

############################
# 3. Heatmap by cell type
############################

# 3-1. load data
data9set_scaled.df2 <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200308_Scaled_data_10000_Cell_type.csv", row.names = 1, check.names = FALSE)
load(file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200221_figure1_for_meeting_report_Seurat_object.Robj")
data9set_cleaned.SO$cell_type <- data9set_cleaned.SO@active.ident


# 3-2. extract meta data
data9set_cleaned.SO_meta.data <- data9set_cleaned.SO@meta.data
data9set_cleaned.SO_meta.data$Cell <- rownames(data9set_cleaned.SO_meta.data)
data9set_cleaned.SO_meta.data_sub <- data9set_cleaned.SO_meta.data[,c(9,10)]

data9set_cleaned.SO.subset <- subset(data9set_cleaned.SO, downsample = 100)

# Identify markers for each cell type (don't need for plot but just for double-check)
#data9set_cleaned_all_markers <- FindAllMarkers(data9set_cleaned.SO)
#write.csv(data9set_cleaned_all_markers, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Marker_genes_Heatmap_Cluster_markers.csv")
#data9set_cleaned_all_markers <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Marker_genes_Heatmap_Cluster_markers.csv", row.names = 1)
#top5 <- data9set_cleaned_all_markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
#top20 <- data9set_cleaned_all_markers %>% group_by(cluster) %>% top_n(20, avg_logFC)

# 3-3. genes to use in heatmap
Ex <- c("Kcnip4", "Slit3", "Grin2a", "Tafa1", "Hs6st3")
In <- c("Grip1", "Erbb4", "Nxph1", "Kcnmb2", "Kcnc2")
Cr <- c("Reln", "Dach1","Ndnf", "Thsd7b", "Cdh4")
Cp <- c("Htr2c", "Otx2os1", "Ttr", "Enpp2", "Spag16")
As <- c("Slc1a2", "Gpc5", "Wdr17", "Slc1a3", "Rgs20")
Mg <- c("Inpp5d", "Tgfbr1", "Ly86", "Fyb", "Hexb")
Ma <- c("F13a1", "Mrc1", "Rbpj", "Stab1", "P2rx7")
Ol <- c("St18", "Prr5l", "Mbp", "Plp1","Mog")
Op <- c("Lhfpl3", "Sox6", "Ptprz1","Vcan", "Pdgfra")
Fb <- c("Slc7a11", "Ranbp3l", "Slc6a20a", "Cped1", "Atp1a2")
Vc <- c("Flt1", "Ebf1","Slco1a4", "Ly6a", "Atp13a5")

genes <- c(Ex, In, Cr, Cp, As, Mg, Ma, Ol, Op, Fb, Vc)

# 3-4. modify data frame
data9set_scaled.df2_sub <- data9set_scaled.df2[genes, ]
data9set_scaled.df2_sub$gene <- rownames(data9set_scaled.df2_sub)
data9set_scaled.df2_sub.long <- gather(data = data9set_scaled.df2_sub, z_score, key = Cell, "ACACGCGCACGTATAC-1":"TTTCACACAAACGGCA-9")
data9set_scaled.df2_sub.long2 <- left_join(data9set_scaled.df2_sub.long, data9set_cleaned.SO_meta.data_sub, by = "Cell")
data9set_scaled.df2_sub.long2$gene <- factor(data9set_scaled.df2_sub.long2$gene, levels = rev(genes))

# 3-5. plot
g1 <- ggplot(data = data9set_scaled.df2_sub.long2, mapping = aes(x = cell_type, y = gene, fill = z_score)) +
  geom_tile() +  theme_classic() + theme(axis.title.y = element_blank(),
                                         axis.title.x = element_blank(), 
                                         axis.text.x = element_text(angle = 90, size = 12, color = "black"), 
                                         axis.text.y = element_text(size = 11, color = "black"),
                                         axis.line = element_line(color = "white"),
                                         axis.ticks = element_line(color = "white"),
                                         legend.title = element_blank(),
                                         legend.text = element_text(size = 11, color = "black")) + 
  scale_fill_viridis_c(limits=c(-2, 2), na.value = "#FDE725FF", labels = c("low", " ", " "," ", "max"))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Marker_genes_Heatmap_cell_type.png", 
       plot = g1, 
       scale = 1, width = 6, height = 10, units = "in", device = "png",
       dpi = 300)

# save doHeatmap
d1 <- DoHeatmap(data9set_cleaned.SO.subset, features = genes, angle = 90, group.bar.height = 0) + 
  scale_fill_viridis() + 
  theme(legend.position = "none")

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Marker_genes_doHeatmap_cell_type.png", 
       plot = d1, 
       scale = 1, width = 6, height = 10, units = "in", device = "png",
       dpi = 300)

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Marker_genes_Heatmap_session_info.txt")

















