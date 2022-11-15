#! /bin/env RScript
# written by SJK at 25. 07. 2020
# This script is for edgeR glmQLFit with CDR for DE analysis
# As I previously analyzed DE by MAST, very few genes were statistically significant.
# I tested edgeR and it showed improved detection of statistically significant genes.

library(Seurat)
library(dplyr)
library(tidyverse)
library(irr)
library(fgsea)
library(tidyr)
library(scran)
library(edgeR)


for(i in unique(data9set_cleaned.SO@meta.data$cell_type)){
  print("DE start")
  print(i)
  # 3-1. subset Seurat object
  subset.SO <- subset(data9set_cleaned.SO, subset = cell_type %in% i)

  # convert Seurat object to sce object
  sce <- as.SingleCellExperiment(subset.SO)

  # convert to dge by edgeR
  dge <- scran::convertTo(sce, type = "edgeR")

  dge$samples$group <- dge$samples$sample

  dge <- calcNormFactors(dge)

  meta <- data.frame(id = rownames(subset.SO@meta.data),
                     group = subset.SO$sample,
                     stringsAsFactors = FALSE)

  cdr <- scale(colMeans(subset.SO@assays$RNA@counts))

  design <- model.matrix(~ cdr + meta$group)

  dge <- estimateDisp(dge, design = design)

  dge$samples$group <- relevel(dge$samples$group, ref="Ctrl")

  fit <- glmQLFit(dge, design = design)


  # 3-2. AD vs Ctrl
  qlf <- glmQLFTest(fit, contrast = c(0,0,1,0))
  tt <- topTags(qlf, n = Inf)
  resNoFilt <- topTags(tt, n=nrow(tt$table))
  resNoFilt <- resNoFilt %>% as.data.frame
  resNoFilt$gene <- rownames(resNoFilt)
  resNoFilt_final <- resNoFilt %>% as.data.frame
  # save csv
  write.csv(resNoFilt_final, paste0("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/", i, "_AD_Ctrl.csv"))

  # 3-3. AD vs ADp40KO
  qlf <- glmQLFTest(fit, contrast = c(0,0,1,-1))
  tt <- topTags(qlf, n = Inf)
  resNoFilt <- topTags(tt, n=nrow(tt$table))
  resNoFilt <- resNoFilt %>% as.data.frame
  # save csv
  write.csv(resNoFilt_final, paste0("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/", i, "_AD_ADp40KO.csv"))

  # 3-4. ADp40KO vs Ctrl 
  qlf <- glmQLFTest(fit, contrast = c(0,0,0,1))
  tt <- topTags(qlf, n = Inf)
  resNoFilt <- topTags(tt, n=nrow(tt$table))
  resNoFilt <- resNoFilt %>% as.data.frame
  resNoFilt$gene <- rownames(resNoFilt)
  resNoFilt_final <- resNoFilt %>% as.data.frame
  # save csv
  write.csv(resNoFilt_final, paste0("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/", i, "_ADp40KO_Ctrl.csv"))

}


