#! /bin/env RScript
# written by SJK at 25. Mar. 2020

.libPaths(c("/home/skim/R/usr_lib", .libPaths()))

#library(Seurat)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyr)
library(topGO)

source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

#######################
# Term explanation
#######################

# The top 10 ranked GO terms of BP and MF according to enriched genes in cluster. 
# GO, gene ontology; BP, biological process; MF, molecular function; CC, cellular component.
# In here, I list the top 10 significant GO terms identified by the elim method
# GO terms which the number of Annotated is more than 500 are removed
# The elim method was design to be more conservative then the classic method

# Term explanation
# Annotated : number of genes in our list which are annotated with the GO-term.
# Signnificant: it is the number of genes enriched in a GO term. (number of significantly DE genes belonging to your input which are annotated with the GO-term.)
# Expected : show an estimate of the number of genes a node of size Annotated would have if the significant genes were to be randomly selected from the gene universe (Under random chance, number of genes that would be expected to be significantly DE and annotated with that term)
# gene ratio: it is the percentage of total enriched genes of cluster in the given GO term
# p-value: it is calcluated from Fisher elim algorithm

# In here, I get gene universe from DE analysis of AD vs Ctrl in Microglia. (genes from pct.1 > 0.01 OR pct.2 > 0.01)
# Cluster significant genes are avg_logFC > 1 & p_val_adj < 0.01 from cluster markers calcualted from FindAllmarkers function

####################
# cluster function
####################
cluster_GO <- function(genesOfInterest, geneUniverse){
  
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse
  
  onts = c( "MF", "BP", "CC" )
  tab <- as.list(onts)
  names(tab) <- onts
  
  for(i in 1:3){
    
    GOdata <- new("topGOdata",
                  description = "GOanalysis",
                  ontology = onts[i],
                  allGenes = geneList,
                  annot = annFUN.org,
                  mapping = "org.Mm.eg.db",
                  ID = "SYMBOL",
                  nodeSize = 20)
    
    res.result1 <- runTest(GOdata, statistic = "fisher", algorithm = "elim")
    res.result2 <- runTest(GOdata, statistic = "fisher", algorithm = "classic")
    
    
    tab[[i]] <- cbind(data.frame(category = onts[i]), GenTable(GOdata, Fisher.elim = res.result1,
                                                               Fisher.classic = res.result2,
                                                               orderBy = "Fisher.elim" , topNodes = 30))
    
  }
  
  topGOResults <- plyr::rbind.fill(tab)
  topGOResults.df <- as.data.frame(topGOResults)
  topGOResults.df$gene_ratio <- topGOResults.df$Significant / topGOResults.df$Annotated
  
  # modification appropriate for plot
  topGOResults.df$Fisher.elim <- as.numeric(topGOResults.df$Fisher.elim)
  topGOResults.df$Fisher.classic <- as.numeric(topGOResults.df$Fisher.classic)
  topGOResults.df$Term <- factor(topGOResults.df$Term, levels = rev(unique(topGOResults.df$Term)))
  
  return(topGOResults.df)
}

scRNA_TopGO_plot2 <- function(topGOResults.df){
  
  topGOResults.df$Fisher.elim <- as.numeric(topGOResults.df$Fisher.elim)
  topGOResults.df$Fisher.elim <- -log10(topGOResults.df$Fisher.elim)
  topGOResults.df$Term <- factor(topGOResults.df$Term, levels = rev(unique(topGOResults.df$Term)))
  
  ggplot(topGOResults.df, aes(x=gene_ratio,
                              y=Term,
                              colour=Fisher.elim,
                              size=Significant)) +
    geom_point() +
    expand_limits(x=0) +
    labs(x="gene ratio", y="GO term", colour="-log10(p-value)", size="Significant") +
    theme_minimal() + theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))
  
}
######################
# GO figure function
######################

scRNA_TopGO_plot2 <- function(topGOResults.df, title = "title"){
  
  topGOResults.df$Fisher.elim <- as.numeric(topGOResults.df$Fisher.elim)
  topGOResults.df$Fisher.elim <- -log10(topGOResults.df$Fisher.elim)
  topGOResults.df$Term <- factor(topGOResults.df$Term, levels = rev(unique(topGOResults.df$Term)))
  
  ggplot(topGOResults.df, aes(x=gene_ratio,
                              y=Term,
                              colour=Fisher.elim,
                              size=Significant)) +
    geom_point() +
    expand_limits(x=0) +
    labs(x="gene ratio", y="GO term", colour="-log10(p-value)", size="Significant") + 
    ggtitle(title) + 
    theme_minimal() + theme(axis.text.x = element_text(size = 10, color = "black"), axis.text.y = element_text(size = 10, color = "black"), axis.title = element_text(size = 12, color = "black"), title = element_text(hjust = 0.5)) + guides(
      size = guide_legend(order = 1),
      fill = guide_legend(order = 0)
    )
  
}

# load cluster marker genes
markers <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200210_FindMarkers_each_cluster_data9set_marker.csv", row.names = 1)

# extract Oligo marker genes
MG_markers <- markers[markers$cluster %in% c("3", "8"), ]

###################
# cluster 3
###################
# gene of interest
# In here, I extract genes padj < 0.01 & logFC > 1
genesOfInterest <- as.character(MG_markers[which(MG_markers$cluster %in% "3" &MG_markers$p_val_adj < 0.01 & MG_markers$avg_logFC > 1), ]$gene) 

# background genes
data9set_MG_DE_AD_Ctrl_markers <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/Cell_type/20200221_Cell_type_DE_by_MAST_MicrogliaAD_Ctrl.csv", row.names = 1)
geneUniverse <- rownames(data9set_MG_DE_AD_Ctrl_markers[data9set_MG_DE_AD_Ctrl_markers$pct.1 > 0.01 | data9set_MG_DE_AD_Ctrl_markers$pct.2 > 0.01, ])

# run GO analysis
topGOResults.df <- cluster_GO(genesOfInterest, geneUniverse)

# MF subset
topGOResults_MF.df1 <- topGOResults.df[which(topGOResults.df$category == "MF" & topGOResults.df$Annotated < 500), ][c(1:10), ]
g1 <- scRNA_TopGO_plot2(topGOResults_MF.df1, title = "Cluster 3 Molecular function")

# BP subset
topGOResults_BP.df1 <- topGOResults.df[which(topGOResults.df$category == "BP" & topGOResults.df$Annotated < 500), ][c(1:10), ]
g2 <- scRNA_TopGO_plot2(topGOResults_BP.df1, title = "Cluster 3 Biological process")


g3 <- arrangeGrob(g1, g2, ncol = 2)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_cluster3_GO_figure.png",
       plot = g3,
       scale = 1, width = 15, height = 5, units = "in", device = "png",
       dpi = 300)

###################
# cluster 8
###################
# gene of interest
# In here, I extract genes padj < 0.01 & logFC > 1
genesOfInterest <- as.character(MG_markers[which(MG_markers$cluster %in% "8" & MG_markers$p_val_adj < 0.01 & MG_markers$avg_logFC > 1), ]$gene) 

# run GO analysis
topGOResults.df <- cluster_GO(genesOfInterest, geneUniverse)

# MF subset
topGOResults_MF.df2 <- topGOResults.df[which(topGOResults.df$category == "MF" & topGOResults.df$Annotated < 500), ][c(1:10), ]
g4 <- scRNA_TopGO_plot2(topGOResults_MF.df2, title = "Cluster 8 Molecular function")

# BP subset
topGOResults_BP.df2 <- topGOResults.df[which(topGOResults.df$category == "BP" & topGOResults.df$Annotated < 500), ][c(1:10), ]
g5 <- scRNA_TopGO_plot2(topGOResults_BP.df2, title = "Cluster 8 Biological process")


g6 <- arrangeGrob(g4, g5, ncol = 2)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_cluster8_GO_figure.png",
       plot = g6,
       scale = 1, width = 15, height = 5, units = "in", device = "png",
       dpi = 300)


##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_cluster_GO_figures_session_info.txt")


