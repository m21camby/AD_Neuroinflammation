#! /bin/env RScript
# written by SJK at 25. Mar. 2020

.libPaths(c("/home/skim/R/usr_lib", .libPaths()))

library(Seurat)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyr)
library(topGO)
library(ggrepel)

source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

####################################################################
# This analysis is for comparison 
# with DAM cluster from Ido Amit's paper (figureS1 (D) & (E))
# Although some markers show similar expression, in here we double checked Genes and GO and compare with his data 
####################################################################

# load DE data between cluster 8 and 3
data9set.SO_MG.marker8vs3 <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200215_MG_2nd_analysis_Cluster8_vs_3_marker.csv", row.names = 1)
data9set.SO_MG.marker8vs3_padj <- data9set.SO_MG.marker8vs3[which(data9set.SO_MG.marker8vs3$p_val_adj < 0.01), ]
data9set.SO_MG.marker8vs3_padj$gene <- rownames(data9set.SO_MG.marker8vs3_padj)

# some padj value is zero. So I set as lowest value 
data9set.SO_MG.marker8vs3_padj$p_val_adj <- ifelse(data9set.SO_MG.marker8vs3_padj$p_val_adj == 0, 1e-320, data9set.SO_MG.marker8vs3_padj$p_val_adj)

# DAM plus gene and minus genes list
DAM_plus <- c("Cst7", "Lpl", "Apoe", "Clec7a", "Ank", "Axl", "Spp1", "Itgax", "Igf1", "Csf1", "Cybb", "Fam20c", "Gpnmb", "Gm11428")
DAM_minus <- c("P2ry12", "Tmem119", "Cx3cr1", "Serinc3")

# subset of DAM plus DF and mins DF
data9set.SO_MG.marker8vs3_padj_plus.df <- data9set.SO_MG.marker8vs3_padj[data9set.SO_MG.marker8vs3_padj$gene %in% DAM_plus,]
# padj > 10
data9set.SO_MG.marker8vs3_padj_minus.df <- data9set.SO_MG.marker8vs3_padj[data9set.SO_MG.marker8vs3_padj$gene %in% DAM_minus,]

# volcano plot
g1 <- ggplot(data9set.SO_MG.marker8vs3_padj) + geom_point(aes(x = avg_logFC, y = -log10(p_val_adj)), size = 0.5) + 
  geom_point(data = data9set.SO_MG.marker8vs3_padj_plus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), color = "red") +
  geom_point(data = data9set.SO_MG.marker8vs3_padj_minus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), color = "blue") +
  ggtitle("cluster8 versus cluster3 in Microglia") + 
  theme_classic() + 
  xlab("log2 fold change") + 
  ylab("-log10(adjusted p-value)") +
  scale_y_continuous(expand = c(0,0), limits = c(0,350)) + 
  theme(legend.position = "none",
        plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12, color = "black")) +
  geom_text_repel(data = data9set.SO_MG.marker8vs3_padj_plus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), label = data9set.SO_MG.marker8vs3_padj_plus.df$gene, force = 30, nudge_y = -0.2, size = 4) + 
  geom_text_repel(data = data9set.SO_MG.marker8vs3_padj_minus.df, aes(x = avg_logFC, y = -log10(p_val_adj)), label = data9set.SO_MG.marker8vs3_padj_minus.df$gene, force = 10, size = 4)

################################
# scatter plot
################################
# Scatterplot showing the average molecules count (log2 scale) of sample compared with other sample (y axis).
load(file="/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200210_figures_for_meeting_report_cleaned_Seurat_object_res0.8.Robj")

# split by sample
data9set_cleaned_MG_c3.SO <- subset(data9set_cleaned.SO, subset = seurat_clusters %in% 3)
data9set_cleaned_MG_c8.SO <- subset(data9set_cleaned.SO, subset = seurat_clusters %in% 8)

# extract UMI counts and calculate average per each gene
MG_c3_count.df <- data.frame(Matrix::rowMeans(GetAssayData(data9set_cleaned_MG_c3.SO, slot = "counts")))
MG_c8_count.df <- data.frame(Matrix::rowMeans(GetAssayData(data9set_cleaned_MG_c8.SO, slot = "counts")))

MG_count.df <- cbind(MG_c3_count.df, MG_c8_count.df)

# remove genes where all zero in samples
MG_count.df <- MG_count.df[apply(MG_count.df == 0, 1, sum) != 3, ]
colnames(MG_count.df) <- c("cluster3", "cluster8")
# add pseudocount 0.001 and logarithm 
MG_count.df <- log(MG_count.df + 0.001)

MG_count.df$gene <- rownames(MG_count.df)
MG_count.df <- as.data.frame(MG_count.df)

# extract genes from volcano plots
MG_count_plus.df <- MG_count.df[MG_count.df$gene %in% DAM_plus, ]
MG_count_minus.df <- MG_count.df[MG_count.df$gene %in% DAM_minus, ]


g2 <- ggplot(MG_count.df, aes(x = cluster3, y = cluster8)) + geom_point(size = 0.5) + 
  geom_point(data = MG_count_plus.df,  aes(x = cluster3, y = cluster8), color = "red", size = 1) +
  geom_point(data = MG_count_minus.df,  aes(x = cluster3, y = cluster8), color = "blue", size = 1) +
  ggtitle("cluster8 versus cluster3 in Microglia") + 
  theme_classic() + 
  xlab("log (average UMI counts) cluster3") + 
  ylab("log (average UMI counts) cluster8") + 
  theme(plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12, color = "black")) +
  geom_abline(slope = 1, color="red", linetype="dashed", size=.5) + 
  geom_text_repel(data = MG_count_plus.df, aes(x = cluster3, y = cluster8), label = MG_count_plus.df$gene, force = 20, nudge_y = .8) + 
  geom_text_repel(data = MG_count_minus.df, aes(x = cluster3, y = cluster8), label = MG_count_minus.df$gene, force = 10, nudge_y = -1, nudge_x = 0.5)

g3 <- arrangeGrob(g1, g2, ncol = 2)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_DAM_comparison_volcano_scatter.png",
       plot = g3,
       scale = 1, width = 11, height = 5, units = "in", device = "png",
       dpi = 300)


###########################
# GO for DE genes
###########################

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
# ggplot function
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

###################
# cluster 8
###################
# gene of interest
# In here, I extract genes padj < 0.01 & logFC > 1
genesOfInterest <- as.character(data9set.SO_MG.marker8vs3_padj[which(data9set.SO_MG.marker8vs3_padj$p_val_adj < 0.01 & data9set.SO_MG.marker8vs3_padj$avg_logFC > 0.5), ]$gene) 

# background genes
geneUniverse <- rownames(data9set.SO_MG.marker8vs3)

# run GO analysis
topGOResults.df <- cluster_GO(genesOfInterest, geneUniverse)

# MF subset
topGOResults_MF.df1 <- topGOResults.df[which(topGOResults.df$category == "MF" & topGOResults.df$Annotated < 500), ][c(1:10), ]
g4 <- scRNA_TopGO_plot2(topGOResults_MF.df1, title = "Cluster 8 Molecular function")

# BP subset
topGOResults_BP.df1 <- topGOResults.df[which(topGOResults.df$category == "BP" & topGOResults.df$Annotated < 500), ][c(1:10), ]
g5 <- scRNA_TopGO_plot2(topGOResults_BP.df1, title = "Cluster 8 Biological process")


g6 <- arrangeGrob(g4, g5, ncol = 2)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_DAM_comparison_cluster8_GO_figure.png",
       plot = g6,
       scale = 1, width = 15, height = 5, units = "in", device = "png",
       dpi = 300)

###################
# cluster 3
###################
# gene of interest
# In here, I extract genes padj < 0.01 & logFC > 1
genesOfInterest <- as.character(data9set.SO_MG.marker8vs3_padj[which(data9set.SO_MG.marker8vs3_padj$p_val_adj < 0.01 & data9set.SO_MG.marker8vs3_padj$avg_logFC < -0.5), ]$gene) 

# background genes
geneUniverse <- rownames(data9set.SO_MG.marker8vs3)

# run GO analysis
topGOResults.df <- cluster_GO(genesOfInterest, geneUniverse)

# MF subset
topGOResults_MF.df2 <- topGOResults.df[which(topGOResults.df$category == "MF" & topGOResults.df$Annotated < 500), ][c(1:10), ]
g7 <- scRNA_TopGO_plot2(topGOResults_MF.df2, title = "Cluster 3 Molecular function")

# BP subset
topGOResults_BP.df2 <- topGOResults.df[which(topGOResults.df$category == "BP" & topGOResults.df$Annotated < 500), ][c(1:10), ]
g8 <- scRNA_TopGO_plot2(topGOResults_BP.df2, title = "Cluster 3 Biological process")


g9 <- arrangeGrob(g7, g8, ncol = 2)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_DAM_comparison_cluster3_GO_figure.png",
       plot = g9,
       scale = 1, width = 15, height = 5, units = "in", device = "png",
       dpi = 300)

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_DAM_comparison_session_info.txt")



