#! /bin/env RScript
# written by SJK at 14 Jan. 2021

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(gridExtra)
source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")


# Analysis

# #Naomi <- read.table(file = "/data/rajewsky/home/skim/Microglia_Heppner/published_DAA_Habib/GSE143758_Admouse_Hippocampus_7m_AllNuclei_UMIcounts.txt.gz", row.names = 1, header = TRUE)
# Naomi <- data.table::fread(file = "/data/rajewsky/home/skim/Microglia_Heppner/published_DAA_Habib/GSE143758_Admouse_Hippocampus_7m_AllNuclei_UMIcounts.txt.gz",showProgress = TRUE, nThread = 16)
# 
# Naomi.SO <- CreateSeuratObject(counts =  Naomi,  min.cells = 3, min.features = 200, project = "DAA")
# 
# Naomi.SO[["percent.mt"]] <- PercentageFeatureSet(object = Naomi.SO, pattern = "^mt-")
# Naomi.SO <- subset(x = Naomi.SO, subset = nCount_RNA < 6000 & nCount_RNA > 200 & percent.mt < 5)
# Naomi.df <- Naomi.SO@meta.data
# 
# CB_WT <- rownames(Naomi.df)[grepl("WT|Wt", rownames(Naomi.df))]
# CB_AD <- rownames(Naomi.df)[grepl("AD|Ad", rownames(Naomi.df))]
# CB_Unknwon <- rownames(Naomi.df)[grepl("Untreated", rownames(Naomi.df))]
# Naomi.df$ID <- rownames(Naomi.df)
# Naomi.df <- Naomi.df %>% mutate(exp = case_when(ID %in% CB_WT ~ "WT",
#                                                 ID %in% CB_AD ~ "AD",
#                                                 ID %in% CB_Unknwon ~ "Untreated"))
# Naomi.SO$exp <- Naomi.df$exp
# Naomi.df2 <- Naomi.SO@meta.data
# Naomi.SO <- NormalizeData(object = Naomi.SO, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
# Naomi.SO <- FindVariableFeatures(object = Naomi.SO, nfeatures = 2000, selection.method = "vst", verbose = FALSE)
# Naomi.SO <- ScaleData(object = Naomi.SO, vars.to.regress = c("percent.mt", "nCount_RNA"))
# 
# # PCA
# Naomi.SO <- RunPCA(Naomi.SO, features = VariableFeatures(object = Naomi.SO), verbose = FALSE, ndims.print = "None")
# #ElbowPlot(Naomi.SO, ndims = 50)
# Naomi.SO <- FindNeighbors(object = Naomi.SO, reduction = "pca", dims = 1:20)
# Naomi.SO <- FindClusters(object = Naomi.SO, resolution = 1, verbose = FALSE)
# Naomi.SO <- RunUMAP(object = Naomi.SO, dims = 1:20, n.neighbors = 20, min.dist = 0.35, n.epochs = 500, spread = 2)
# d5 <- DimPlot(object = Naomi.SO, reduction = "umap")
# save(Naomi.SO, file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig6_Naomi_Habib_confirm.Robj")

load(file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig6_Naomi_Habib_confirm.Robj")

f1 <- FeaturePlot(Naomi.SO, features = "Il12rb1", order = TRUE, pt.size = 1, min.cutoff = 0, max.cutoff = 5)
f1 <- f1 + theme(legend.position = "none", 
                 axis.line=element_blank(), 
                 axis.ticks=element_blank(), 
                 axis.text=element_blank(),
                 axis.title = element_blank())

f2 <- FeaturePlot(Naomi.SO, features = "Il12rb2", order = TRUE, pt.size = 1, min.cutoff = 0, max.cutoff = 5)
f2 <- f2 + theme(legend.position = "none", 
                 axis.line=element_blank(), 
                 axis.ticks=element_blank(), 
                 axis.text=element_blank(),
                 axis.title = element_blank())

g1 <- arrangeGrob(f1, f2, ncol = 2)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig6_Naomi_Habib_confirm_Il12rb1_Il12rb2.pdf",
       plot = g1,
       scale = 1, width = 10, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)


f1 <- FeaturePlot(Naomi.SO, features = "Il12a", order = TRUE, pt.size = 1, min.cutoff = 0, max.cutoff = 5)
f1 <- f1 + theme(legend.position = "none", 
                 axis.line=element_blank(), 
                 axis.ticks=element_blank(), 
                 axis.text=element_blank(),
                 axis.title = element_blank())


f2 <- FeaturePlot(Naomi.SO, features = "Il23r", order = TRUE, pt.size = 1, min.cutoff = 0, max.cutoff = 5)
f2 <- f2 + theme(legend.position = "none", 
                 axis.line=element_blank(), 
                 axis.ticks=element_blank(), 
                 axis.text=element_blank(),
                 axis.title = element_blank())

g1 <- arrangeGrob(f1, f2, ncol = 2)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig6_Naomi_Habib_confirm_Il12a_Il23r.pdf",
       plot = g1,
       scale = 1, width = 10, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)


f1 <- FeaturePlot(Naomi.SO, features = "C4b", order = FALSE, pt.size = 1, min.cutoff = 0, max.cutoff = 5)
f1 <- f1 + theme(legend.position = "none", 
                 axis.line=element_blank(), 
                 axis.ticks=element_blank(), 
                 axis.text=element_blank(),
                 axis.title = element_blank())

f2 <- FeaturePlot(Naomi.SO, features = "Kcnip4", order = FALSE, pt.size = 1, min.cutoff = 0, max.cutoff = 5)
f2 <- f2 + theme(legend.position = "none", 
                 axis.line=element_blank(), 
                 axis.ticks=element_blank(), 
                 axis.text=element_blank(),
                 axis.title = element_blank())

g1 <- arrangeGrob(f1, f2, ncol = 2)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig6_Naomi_Habib_confirm_C4b_Kcnip4.pdf",
       plot = g1,
       scale = 1, width = 10, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)


f1 <- FeaturePlot(Naomi.SO, features = "Mbp", order = FALSE, pt.size = 1, min.cutoff = 0, max.cutoff = 5)
f1 <- f1 + theme(legend.position = "none", 
                 axis.line=element_blank(), 
                 axis.ticks=element_blank(), 
                 axis.text=element_blank(),
                 axis.title = element_blank())

f2 <- FeaturePlot(Naomi.SO, features = "Snap25", order = FALSE, pt.size = 1, min.cutoff = 0, max.cutoff = 5)
f2 <- f2 + theme(legend.position = "none", 
                 axis.line=element_blank(), 
                 axis.ticks=element_blank(), 
                 axis.text=element_blank(),
                 axis.title = element_blank())

g1 <- arrangeGrob(f1, f2, ncol = 2)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig6_Naomi_Habib_confirm_Mbp_Snap25.pdf",
       plot = g1,
       scale = 1, width = 10, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig6_Naomi_Habib_confirm_session_info.txt")

