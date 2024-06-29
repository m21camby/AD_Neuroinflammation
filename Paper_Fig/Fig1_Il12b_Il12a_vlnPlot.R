##############################################

# This goes to Fig2_Microglia_Il12_Vln_figures

##############################################


library(Seurat)
library(gridExtra)
library(dplyr)
library(gridExtra)
library(ggplot2)

load(file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200221_figure1_for_meeting_report_Seurat_object.Robj")

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


##########
# EX
##########
data9set_sub.SO <- subset(data9set_cleaned.SO, subset = cell_type %in% c("Dentate_Gyrus","subiculum","CA1", "CA2_3","Unidentified_Neurons"))
data9set_sub.SO <- subset(data9set_sub.SO, subset = sample %in% c("AD"))
EX_Il12b <- data9set_sub.SO@assays$RNA@data[rownames(data9set_sub.SO@assays$RNA@data) %in% "Il12b", ] %>% as.data.frame
colnames(EX_Il12b) <- "Il12b"
EX_Il12b$celltype <- "EX"

###########
# IN
###########
data9set_sub.SO <- subset(data9set_cleaned.SO, subset = cell_type %in% c("Inhibitory_Neurons"))
data9set_sub.SO <- subset(data9set_sub.SO, subset = sample %in% c("AD"))
IN_Il12b <- data9set_sub.SO@assays$RNA@data[rownames(data9set_sub.SO@assays$RNA@data) %in% "Il12b", ] %>% as.data.frame
colnames(IN_Il12b) <- "Il12b"
IN_Il12b$celltype <- "IN"


##############
# OL
##############
data9set_sub.SO <- subset(data9set_cleaned.SO, subset = cell_type %in% c("MOL", "MFOL"))
data9set_sub.SO <- subset(data9set_sub.SO, subset = sample %in% c("AD"))
OL_Il12b <- data9set_sub.SO@assays$RNA@data[rownames(data9set_sub.SO@assays$RNA@data) %in% "Il12b", ] %>% as.data.frame
colnames(OL_Il12b) <- "Il12b"
OL_Il12b$celltype <- "OL"

###############
# OPC
###############
data9set_sub.SO <- subset(data9set_cleaned.SO, subset = cell_type %in% c("OPC", "NFOL"))
data9set_sub.SO <- subset(data9set_sub.SO, subset = sample %in% c("AD"))
OPC_Il12b <- data9set_sub.SO@assays$RNA@data[rownames(data9set_sub.SO@assays$RNA@data) %in% "Il12b", ] %>% as.data.frame
colnames(OPC_Il12b) <- "Il12b"
OPC_Il12b$celltype <- "OPC"

###############
# Microglia
###############
data9set_sub.SO <- subset(data9set_cleaned.SO, subset = cell_type %in% c("Microglia"))
data9set_sub.SO <- subset(data9set_sub.SO, subset = sample %in% c("AD"))
MG_Il12b <- data9set_sub.SO@assays$RNA@data[rownames(data9set_sub.SO@assays$RNA@data) %in% "Il12b", ] %>% as.data.frame
colnames(MG_Il12b) <- "Il12b"
MG_Il12b$celltype <- "MG"

################
# Astrocytes
################
data9set_sub.SO <- subset(data9set_cleaned.SO, subset = cell_type %in% c("Astrocytes"))
data9set_sub.SO <- subset(data9set_sub.SO, subset = sample %in% c("AD"))
AS_Il12b <- data9set_sub.SO@assays$RNA@data[rownames(data9set_sub.SO@assays$RNA@data) %in% "Il12b", ] %>% as.data.frame
colnames(AS_Il12b) <- "Il12b"
AS_Il12b$celltype <- "AS"

Il12b <- rbind(EX_Il12b, IN_Il12b, OL_Il12b, OPC_Il12b, MG_Il12b, AS_Il12b)
Il12b2 <- Il12b[Il12b$Il12b != 0, ]

ggplot(Il12b2, aes(x=celltype, y=Il12b, color=celltype)) + 
  geom_violin(trim=FALSE) + 
  geom_jitter(position=position_jitter(0.2), size = 0.3, alpha = 0.3) + ggtitle("Il12b normalized counts excluding zero counts") 

##############################################

# This goes to Fig2_Microglia_Il12_Vln_figures

##############################################