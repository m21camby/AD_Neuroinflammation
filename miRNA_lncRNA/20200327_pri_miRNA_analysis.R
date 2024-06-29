
SP064_023_pri_miRNA <- Read10X(data.dir = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/miRNA_lncRNA/SP064_023_pri_miRNA/outs/filtered_feature_bc_matrix/")
SP064_023_pri_miRNA.SO <- CreateSeuratObject(counts = SP064_023_pri_miRNA, project = "pbmc3k", min.cells = 3)

SP064_023_pri_miRNA_rowSums.df <- as.data.frame(Matrix::rowSums(GetAssayData(SP064_023_pri_miRNA.SO, slot = "counts")))

SP064_023_pri_miRNA_rowSums.df$gene <- rownames(SP064_023_pri_miRNA_rowSums.df)

sum(SP064_023_pri_miRNA_rowSums.df$`Matrix::rowSums(GetAssayData(SP064_023_pri_miRNA.SO, slot = "counts"))`)


SP064_023_pri_miRNA.df <- as.data.frame(SP064_023_pri_miRNA.SO[["RNA"]]@counts["let-7i", ])

SP064_023_pri_miRNA.df$cell_barcode <- paste0(rownames(SP064_023_pri_miRNA.df), "-2")

colnames(SP064_023_pri_miRNA.df) <- c("let7i", "cell_barcode")

    