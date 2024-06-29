#! /bin/env RScript
# written by SJK at 08. Nov. 2020

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)

# If you have extra time, would you mind creating featureplots for the following genes:
# Il6, Il11, Il27a, Gm45928, Ctf2 /Gm494, Ctf1, Csf3/Csfg, Clcf1/Bsf3/Nnt1

load(file="/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200210_figures_for_meeting_report_cleaned_Seurat_object_res0.8.Robj")

# Il6
f1 <- FeaturePlot(data9set_cleaned.SO, features = "Il6", order = TRUE, pt.size = 1)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig3_Featureplots_genes_Il6.png",
       plot = f1,
       scale = 1, width = 7, height = 6, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig3_Featureplots_genes_Il6.pdf",
       plot = f1,
       scale = 1, width = 7, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)

# Il11
f1 <- FeaturePlot(data9set_cleaned.SO, features = "Il11", order = TRUE, pt.size = 1)
f1
ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig3_Featureplots_genes_Il11.png",
       plot = f1,
       scale = 1, width = 7, height = 6, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig3_Featureplots_genes_Il11.pdf",
       plot = f1,
       scale = 1, width = 7, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)

# no Il27a

# Gm45928
f1 <- FeaturePlot(data9set_cleaned.SO, features = "Gm45928", order = TRUE, pt.size = 1)
f1
ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig3_Featureplots_genes_Gm45928.png",
       plot = f1,
       scale = 1, width = 7, height = 6, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig3_Featureplots_genes_Gm45928.pdf",
       plot = f1,
       scale = 1, width = 7, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)

# Ctf2
f1 <- FeaturePlot(data9set_cleaned.SO, features = "Ctf2", order = TRUE, pt.size = 1)
f1
ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig3_Featureplots_genes_Ctf2.png",
       plot = f1,
       scale = 1, width = 7, height = 6, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig3_Featureplots_genes_Ctf2.pdf",
       plot = f1,
       scale = 1, width = 7, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)

# no Gm494

# Ctf1
f1 <- FeaturePlot(data9set_cleaned.SO, features = "Ctf1", order = TRUE, pt.size = 1)
f1
ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig3_Featureplots_genes_Ctf1.png",
       plot = f1,
       scale = 1, width = 7, height = 6, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig3_Featureplots_genes_Ctf1.pdf",
       plot = f1,
       scale = 1, width = 7, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)

# Csf3
f1 <- FeaturePlot(data9set_cleaned.SO, features = "Csf3", order = TRUE, pt.size = 1)
f1
ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig3_Featureplots_genes_Csf3.png",
       plot = f1,
       scale = 1, width = 7, height = 6, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig3_Featureplots_genes_Csf3.pdf",
       plot = f1,
       scale = 1, width = 7, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)

# no Csfg

# no Clcf1

# no Bsf3

# no Nnt1
