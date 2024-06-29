#! /bin/env RScript
# written by SJK at 3. May. 2020

# This file is for heatmap of below genes request from Shirin
# Excitatory Neurons: Grin1, Rbfox3, Slc17a7
# Inhibitory Neurons: Grip1, Erbb4, Gad2
# Microglia : Hexb, Tgfbr1, Inpp5d
# Macrophages: F13a1, Mrc1, Stab1
# Choroid Plexus: Ttr, Htr2c, Spag16
# Cajal Rutzius: ReIn, Dach1, Thsd7b
# Astrocytes: Slc1a2, Clu, Gja1
# Oligodendrocyte: Mbp, Plp, Mog
# OPC: Vcan, Lhfpl3, Pdgfra
# Fibroblast: Ranbpl3, Flt1, Cped1
# Vascular cells: Flt1, Cldn5, Ly6

.libPaths(c("/home/skim/R/usr_lib", .libPaths()))


library(Seurat)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)
library(tidyr)

source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

############################
# Heatmap by cell type
############################

# load data
data9set_scaled.df2 <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200308_Scaled_data_10000_Cell_type.csv", row.names = 1, check.names = FALSE)
load(file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200221_figure1_for_meeting_report_Seurat_object.Robj")
data9set_cleaned.SO$cell_type <- data9set_cleaned.SO@active.ident


# extract meta data
data9set_cleaned.SO_meta.data <- data9set_cleaned.SO@meta.data
data9set_cleaned.SO_meta.data$Cell <- rownames(data9set_cleaned.SO_meta.data)
data9set_cleaned.SO_meta.data_sub <- data9set_cleaned.SO_meta.data[,c(9,10)]

data9set_cleaned.SO.subset <- subset(data9set_cleaned.SO, downsample = 100)


# genes to use in heatmap
Ex <- c("Grin2a", "Tafa1", "Kcnip4") # modified Grin1 - > Grin2a, Rbfox3 -> Tafa1, Slc17a7 -> Kcnip4
In <- c("Grip1", "Erbb4", "Gad2")
Mg <- c("Inpp5d", "Tgfbr1", "Hexb")
Ma <- c("F13a1", "Mrc1", "Stab1")
Cp <- c("Htr2c", "Ttr", "Spag16")
Cr <- c("Reln", "Dach1","Thsd7b") # modified ReIn -> Reln
As <- c("Slc1a2", "Slc1a3", "Gfap") # modified Clu-> Gfap, Gja1 -> Slc1a3
Ol <- c("Mbp", "Plp1","Mog")
Op <- c("Lhfpl3","Vcan", "Pdgfra")
Fb <- c("Slc7a11", "Ranbp3l","Cped1") # modified Flt1 -> Slc7a11
Vc <- c("Flt1", "Cldn5","Ly6a")

genes <- c(Ex, In, Cr, Cp, As, Mg, Ma, Ol, Op, Fb, Vc)


# modify data frame
data9set_scaled.df2_sub <- data9set_scaled.df2[genes, ]
data9set_scaled.df2_sub$gene <- rownames(data9set_scaled.df2_sub)
data9set_scaled.df2_sub.long <- gather(data = data9set_scaled.df2_sub, z_score, key = Cell, "ACACGCGCACGTATAC-1":"TTTCACACAAACGGCA-9")
data9set_scaled.df2_sub.long2 <- left_join(data9set_scaled.df2_sub.long, data9set_cleaned.SO_meta.data_sub, by = "Cell")
data9set_scaled.df2_sub.long2$gene <- factor(data9set_scaled.df2_sub.long2$gene, levels = rev(genes))

# plot
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

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Marker_genes_Heatmap_cell_type_Modified.png", 
       plot = g1, 
       scale = 1, width = 6, height = 10, units = "in", device = "png",
       dpi = 300)


##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Marker_genes_Heatmap_Modified_session_info.txt")

















