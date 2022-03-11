# Written by SJK 08.08 2019

# This R script is for DroNc-seq AD project
# This file aims for GSEA Analysis of scRNA-seq data


library(biomaRt)
library(fgsea)
library(dplyr)

scRNA_GSEA <- function(DE.df, 
                      pathway = "Reactome",
                      metric = "logFC"){

if(pathway == "Reactome"){
  pathways.hallmark <- gmtPathways("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/ext_files/ReactomePathways.gmt")
}
if(pathway == "Hall"){
  pathways.hallmark <- gmtPathways("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/ext_files/h.all.v7.0.symbols.gmt")
}
if(pathway == "GO"){
  pathways.hallmark <- gmtPathways("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/GSEA_gmt_files/c5.go.v7.2.symbols.gmt")
}
if(pathway == "Immune"){
  pathways.hallmark <- gmtPathways("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/GSEA_gmt_files/c7.all.v7.2.symbols.gmt")
}
  
#################
# human genes
# error temporary
#################
#mart <- useDataset("mmusculus_gene_ensembl", mart=useMart("ensembl"))
#bm <- getBM(attributes=c("ensembl_gene_id", "hsapiens_homolog_associated_gene_name"),  
#            filters = 'mgi_symbol',
#            values = rownames(DE.df),
#            mart=mart)

##############################################################
# due to error on bioconductor I use below code for temporary
##############################################################
bm <- read.table("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200113_CCI_test_GENE_BM.df")


###############################################################


# only distinct gene names
bm <- bm %>% distinct() %>% as_tibble() %>% na_if("") %>% na.omit() %>% as.data.frame()
# for Data Frame
bm.df <- data.frame(ensembl_gene_id= bm$ensembl_gene_id, hsapiens_homolog_associated_gene_name = bm$hsapiens_homolog_associated_gene_name, stringsAsFactors = FALSE)

#names(bm.df) <- gsub("\\s"," ",names(bm.df))
#names(bm.df) <- names(bm.df) %>% stringr::str_replace_all("\\s","_")

##################
# mouse genes
# error temporary
##################

#bm2 <- getBM(attributes=c("ensembl_gene_id", "mgi_symbol"),  
#             filters = 'mgi_symbol',
#             values = rownames(DE.df),
#             mart=mart)

##############################################################
# due to error on bioconductor I use below code for temporary
##############################################################
bm2 <- read.table("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200113_CCI_test_GENE_BM2.df")


#############################################



bm2.df <- as.data.frame(bm2)

total <- merge(bm.df,bm2.df, by = "ensembl_gene_id")

total <- data.frame(Gene = total$mgi_symbol, HGene = total$hsapiens_homolog_associated_gene_name)

# only unique gene
total_uni.df <- distinct(total, HGene, .keep_all = TRUE)

# merge 
DE.df$Gene <- rownames(DE.df)
total_HUMAN <- merge(DE.df, total_uni.df, by = "Gene")

# Creating rnk file
if(metric == "logFC"){
  total_HUMAN$metric= total_HUMAN$avg_logFC
}
if(metric == "padj"){
  total_HUMAN$fcsign <- sign(total_HUMAN$avg_logFC)
  total_HUMAN$logP = -log10(total_HUMAN$p_val_adj)
  total_HUMAN$metric= total_HUMAN$logP/total_HUMAN$fcsign
}
if(metric == "pvalue"){
  total_HUMAN$fcsign <- sign(total_HUMAN$avg_logFC)
  total_HUMAN$logP = -log10(total_HUMAN$p_val)
  total_HUMAN$metric= total_HUMAN$logP/total_HUMAN$fcsign
}

metric <- total_HUMAN[,c("HGene", "metric")]
colnames(metric) <- c("Gene", "metric")

metric <- metric %>% distinct()

ranks <- tibble::deframe(metric)

ranks <<- ranks
fgseaRes <- fgsea(pathways.hallmark, ranks, minSize=15, maxSize=500, nperm=10000)

return(fgseaRes)
}

scRNA_GSEA_edgeR <- function(DE.df, 
                       pathway = "Reactome",
                       metric = "logFC"){
  
  if(pathway == "Reactome"){
    pathways.hallmark <- gmtPathways("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/ext_files/ReactomePathways.gmt")
  }
  if(pathway == "Hall"){
    pathways.hallmark <- gmtPathways("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/ext_files/h.all.v7.0.symbols.gmt")
  }
  if(pathway == "GO"){
    pathways.hallmark <- gmtPathways("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/GSEA_gmt_files/c5.go.v7.2.symbols.gmt")
  }
  if(pathway == "Immune"){
    pathways.hallmark <- gmtPathways("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/GSEA_gmt_files/c7.all.v7.2.symbols.gmt")
  }
  
  #################
  # human genes
  # error temporary
  #################
  #mart <- useDataset("mmusculus_gene_ensembl", mart=useMart("ensembl"))
  #bm <- getBM(attributes=c("ensembl_gene_id", "hsapiens_homolog_associated_gene_name"),  
  #            filters = 'mgi_symbol',
  #            values = rownames(DE.df),
  #            mart=mart)
  
  ##############################################################
  # due to error on bioconductor I use below code for temporary
  ##############################################################
  bm <- read.table("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200113_CCI_test_GENE_BM.df")
  
  
  ###############################################################
  
  
  # only distinct gene names
  bm <- bm %>% distinct() %>% as_tibble() %>% na_if("") %>% na.omit() %>% as.data.frame()
  # for Data Frame
  bm.df <- data.frame(ensembl_gene_id= bm$ensembl_gene_id, hsapiens_homolog_associated_gene_name = bm$hsapiens_homolog_associated_gene_name, stringsAsFactors = FALSE)
  
  #names(bm.df) <- gsub("\\s"," ",names(bm.df))
  #names(bm.df) <- names(bm.df) %>% stringr::str_replace_all("\\s","_")
  
  ##################
  # mouse genes
  # error temporary
  ##################
  
  #bm2 <- getBM(attributes=c("ensembl_gene_id", "mgi_symbol"),  
  #             filters = 'mgi_symbol',
  #             values = rownames(DE.df),
  #             mart=mart)
  
  ##############################################################
  # due to error on bioconductor I use below code for temporary
  ##############################################################
  bm2 <- read.table("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200113_CCI_test_GENE_BM2.df")
  
  
  #############################################
  
  
  
  bm2.df <- as.data.frame(bm2)
  
  total <- merge(bm.df,bm2.df, by = "ensembl_gene_id")
  
  total <- data.frame(Gene = total$mgi_symbol, HGene = total$hsapiens_homolog_associated_gene_name)
  
  # only unique gene
  total_uni.df <- distinct(total, HGene, .keep_all = TRUE)
  
  # merge 
  DE.df$Gene <- rownames(DE.df)
  total_HUMAN <- merge(DE.df, total_uni.df, by = "Gene")
  
  # Creating rnk file
  if(metric == "logFC"){
    total_HUMAN$metric= total_HUMAN$logFC
  }
  if(metric == "padj"){
    total_HUMAN$fcsign <- sign(total_HUMAN$logFC)
    total_HUMAN$logP = -log10(total_HUMAN$FDR)
    total_HUMAN$metric= total_HUMAN$logP/total_HUMAN$fcsign
  }
  if(metric == "pvalue"){
    total_HUMAN$fcsign <- sign(total_HUMAN$avg_logFC)
    total_HUMAN$logP = -log10(total_HUMAN$p_val)
    total_HUMAN$metric= total_HUMAN$logP/total_HUMAN$fcsign
  }
  
  metric <- total_HUMAN[,c("HGene", "metric")]
  colnames(metric) <- c("Gene", "metric")
  
  metric <- metric %>% distinct()
  
  ranks <- tibble::deframe(metric)
  
  ranks <<- ranks
  fgseaRes <- fgsea(pathways.hallmark, ranks, minSize=15, maxSize=500, nperm=10000)
  
  return(fgseaRes)
}


