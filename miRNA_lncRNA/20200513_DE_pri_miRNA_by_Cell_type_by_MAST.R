#! /bin/env RScript
# written by SJK at 13. May. 2020
# For DE by each cell type including pri-miRNA

.libPaths(c("/home/skim/R/usr_lib", .libPaths()))

library(Seurat)
library(dplyr)
library(gridExtra)
library(tidyverse)
library(MAST)

source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

load(file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/miRNA_lncRNA/R_Scripts/20200505_pri_miRNA_dgCMatrix_Create_Seurat_with_intergenic_pri_miRNA.obj")

data9set.meta.df <- data9set.SO@meta.data

data9set.SO$cell_type <- data9set.SO@active.ident



#subset.SO <- subset(data9set.SO, subset = cell_type %in% "CA2/CA3 Neuron")
#DE_AD_Ctrl <- FindMarkers(subset.SO, ident.1 =  "AD", ident.2 = "Ctrl", group.by = "sample", test.use = "MAST", logfc.threshold = 0, min.pct = 0)
#DE_AD_ADp40KO <- FindMarkers(subset.SO, ident.1 =  "AD", ident.2 = "ADp40KO", group.by = "sample", test.use = "MAST", logfc.threshold = 0, min.pct = 0)
#DE_ADp40KO_Ctrl <- FindMarkers(subset.SO, ident.1 =  "ADp40KO", ident.2 = "Ctrl", group.by = "sample", test.use = "MAST", logfc.threshold = 0, min.pct = 0)
#write.csv(DE_ADp40KO_Ctrl, paste0("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/miRNA_lncRNA/R_Scripts/DE_csv_files/20200513_pri_miRNA_Cell_type_DE_by_MAST_CA2_CA3_ADp40KO_Ctrl.csv"))
#unique(data9set.SO@meta.data$cell_type)[9:16]
#gsub("/", "_", unique(data9set.SO@meta.data$cell_type)[9])



 
# 3. run each cluster DE and save csv 
for(i in unique(data9set.SO@meta.data$cell_type)[10:16]){
  print("DE start")
  print(i)
  # 3-1. subset Seurat object
  subset.SO <- subset(data9set.SO, subset = cell_type %in% i)
  
  # 3-2. AD vs Ctrl
  DE_AD_Ctrl <- FindMarkers(subset.SO, ident.1 =  "AD", ident.2 = "Ctrl", group.by = "sample", test.use = "MAST", logfc.threshold = 0, min.pct = 0)
  write.csv(DE_AD_Ctrl, paste0("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/miRNA_lncRNA/R_Scripts/DE_csv_files/20200513_pri_miRNA_Cell_type_DE_by_MAST_", i,"AD_Ctrl.csv"))

  # 3-3. AD vs ADp40KO
  DE_AD_ADp40KO <- FindMarkers(subset.SO, ident.1 =  "AD", ident.2 = "ADp40KO", group.by = "sample", test.use = "MAST", logfc.threshold = 0, min.pct = 0)
  write.csv(DE_AD_ADp40KO, paste0("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/miRNA_lncRNA/R_Scripts/DE_csv_files/20200513_pri_miRNA_Cell_type_DE_by_MAST_",i,"AD_ADp40KO.csv"))

  # 3-4. ADp40KO vs Ctrl 
  DE_ADp40KO_Ctrl <- FindMarkers(subset.SO, ident.1 =  "ADp40KO", ident.2 = "Ctrl", group.by = "sample", test.use = "MAST", logfc.threshold = 0, min.pct = 0)
  write.csv(DE_ADp40KO_Ctrl, paste0("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/miRNA_lncRNA/R_Scripts/DE_csv_files/20200513_pri_miRNA_Cell_type_DE_by_MAST_",i,"ADp40KO_Ctrl.csv"))

}

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/miRNA_lncRNA/R_Scripts/20200513_DE_pri_miRNA_by_Cell_type_by_MAST_session_info.txt")



