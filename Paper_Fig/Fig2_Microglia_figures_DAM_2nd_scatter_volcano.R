#! /bin/env RScript
# written by SJK at 22. Oct. 2020

library(gridExtra)
library(dplyr)
library(ggplot2)
library(tidyr)
library(Seurat)
library(ggrepel)
library(readxl)
library(tidyr)
library(scran)
library(edgeR)
library(topGO)

source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

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
    
    # export GO ID from BP
    if(i == 1){
      allGO_MF <- genesInTerm(GOdata)
    }
    if(i == 2){
      allGO_BP <- genesInTerm(GOdata)
    }
    
    if(i == 3){
      allGO_CC <- genesInTerm(GOdata)
    }
    
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
  
  for(i in 1:30){
    ID <- topGOResults.df$GO.ID[i]
    #print(ID)
    #print(paste(genesOfInterest[genesOfInterest %in% allGO[ID][[1]]], collapse=', ' ))
    topGOResults.df[i, "genes"] <- paste(genesOfInterest[genesOfInterest %in% allGO_MF[ID][[1]]], collapse=', ' )
  }
  for(i in 31:60){
    ID <- topGOResults.df$GO.ID[i]
    #print(ID)
    #print(paste(genesOfInterest[genesOfInterest %in% allGO[ID][[1]]], collapse=', ' ))
    topGOResults.df[i, "genes"] <- paste(genesOfInterest[genesOfInterest %in% allGO_BP[ID][[1]]], collapse=', ' )
  }
  for(i in 61:90){
    ID <- topGOResults.df$GO.ID[i]
    #print(ID)
    #print(paste(genesOfInterest[genesOfInterest %in% allGO[ID][[1]]], collapse=', ' ))
    topGOResults.df[i, "genes"] <- paste(genesOfInterest[genesOfInterest %in% allGO_CC[ID][[1]]], collapse=', ' )
  }
  
  return(topGOResults.df)
}

seurat_mouse_genes <- read.table(file = "/data/rajewsky/home/skim/Microglia_Heppner/Bulk_transcriptome/R_Scripts/20200312_Bulk_Correlation_mouse_genes_bio_type_biomaRt.txt")

load(file="/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200210_figures_for_meeting_report_cleaned_Seurat_object_res0.8.Robj")

data9set_cleaned_MG.SO <- subset(data9set_cleaned.SO, subset = seurat_clusters %in% c(3, 8))

###########################
# Cluster 8 AD AD vs Cluster 3 WT

data9set_cleaned_3.SO <- subset(data9set_cleaned.SO, subset = seurat_clusters %in% c(3))
data9set_cleaned_3_WT.SO <- subset(data9set_cleaned_3.SO, subset = sample %in% c("Ctrl"))

data9set_cleaned_8.SO <- subset(data9set_cleaned.SO, subset = seurat_clusters %in% c(8))
data9set_cleaned_8_AD.SO <- subset(data9set_cleaned_8.SO, subset = sample %in% c("AD"))

data9set_cleaned_3_WT_8_AD.SO <- merge(data9set_cleaned_3_WT.SO, data9set_cleaned_8_AD.SO)

##############################


data9set_cleaned_3.SO <- subset(data9set_cleaned.SO, subset = seurat_clusters %in% c(3))
cluster3.df <- rowMeans(GetAssayData(data9set_cleaned_3.SO, slot = "data")) %>% as.data.frame
colnames(cluster3.df) <- "cluster3"

data9set_cleaned_8.SO <- subset(data9set_cleaned.SO, subset = seurat_clusters %in% c(8))
cluster8.df <- rowMeans(GetAssayData(data9set_cleaned_8.SO, slot = "data")) %>% as.data.frame
colnames(cluster8.df) <- "cluster8"

clusters.df <- cbind(cluster3.df, cluster8.df)
clusters.df$gene <- rownames(clusters.df)

data9set_cleaned_20.SO <- subset(data9set_cleaned.SO, subset = seurat_clusters %in% c(20))
FeaturePlot(data9set_cleaned_20.SO, features = "Erbb4", split.by = "sample")
c20DE <- FindMarkers(data9set_cleaned_20.SO, ident.1 = "AD", ident.2 = "Ctrl", group.by = "sample", test.use = "MAST")
c20DE$gene <- rownames(c20DE)
c20DE2 <- FindMarkers(data9set_cleaned_20.SO, ident.1 = "AD", ident.2 = "ADp40KO", group.by = "sample", test.use = "MAST")
c20DE2$gene <- rownames(c20DE2)


#######################
# edgeR cluster 3 vs 8
#######################

# sce <- as.SingleCellExperiment(data9set_cleaned_MG.SO)
# dge <- scran::convertTo(sce, type = "edgeR")
# 
# dge$samples$seurat_clusters <- dge$samples$seurat_clusters %>% as.character
# 
# dge$samples$group <- dge$samples$seurat_clusters
# table(dge$samples$group)
# 
# dge <- calcNormFactors(dge)
# 
# data9set_cleaned_MG.SO$seurat_clusters <- data9set_cleaned_MG.SO$seurat_clusters %>% as.character
# 
# meta <- data.frame(id = rownames(data9set_cleaned_MG.SO@meta.data),
#                    group = data9set_cleaned_MG.SO$seurat_clusters,
#                    stringsAsFactors = FALSE)
# 
# cdr <- scale(colMeans(data9set_cleaned_MG.SO@assays$RNA@counts))
# 
# design <- model.matrix(~ cdr + meta$group)
# 
# #design2 <- model.matrix(~0+meta$group)
# 
# dge <- estimateDisp(dge, design = design)
# 
# dge$samples$group <- dge$samples$group %>% as.factor
# 
# 
# dge$samples$group <- relevel(dge$samples$group, ref="3")
# 
# fit <- glmQLFit(dge, design = design)
# 
# qlf <- glmQLFTest(fit)
# tt <- topTags(qlf, n = Inf)
# 
# resNoFilt <- topTags(tt, n=nrow(tt$table))
# 
# resNoFilt <- resNoFilt %>% as.data.frame
# 
# resNoFilt$gene <- rownames(resNoFilt)
# 
# resNoFilt_final <- resNoFilt %>% as.data.frame
# 
# write.csv(resNoFilt_final, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cluster/Fig2_Microglia_figures_DAM_edgeR_DE_test_3_vs_8.csv")

########################################
# edgeR cluster 3 (WT) vs 8 (AD) only
########################################

# sce <- as.SingleCellExperiment(data9set_cleaned_3_WT_8_AD.SO)
# dge <- scran::convertTo(sce, type = "edgeR")
# 
# dge$samples$seurat_clusters <- dge$samples$seurat_clusters %>% as.character
# 
# dge$samples$group <- dge$samples$seurat_clusters
# table(dge$samples$group)
# 
# dge <- calcNormFactors(dge)
# 
# data9set_cleaned_3_WT_8_AD.SO$seurat_clusters <- data9set_cleaned_3_WT_8_AD.SO$seurat_clusters %>% as.character
# 
# meta <- data.frame(id = rownames(data9set_cleaned_3_WT_8_AD.SO@meta.data),
#                    group = data9set_cleaned_3_WT_8_AD.SO$seurat_clusters,
#                    stringsAsFactors = FALSE)
# 
# cdr <- scale(colMeans(data9set_cleaned_3_WT_8_AD.SO@assays$RNA@counts))
# 
# design <- model.matrix(~ cdr + meta$group)
# 
# #design2 <- model.matrix(~0+meta$group)
# 
# dge <- estimateDisp(dge, design = design)
# 
# dge$samples$group <- dge$samples$group %>% as.factor
# 
# 
# dge$samples$group <- relevel(dge$samples$group, ref="3")
# 
# fit <- glmQLFit(dge, design = design)
# 
# qlf <- glmQLFTest(fit)
# tt <- topTags(qlf, n = Inf)
# 
# resNoFilt <- topTags(tt, n=nrow(tt$table))
# 
# resNoFilt <- resNoFilt %>% as.data.frame
# 
# resNoFilt$gene <- rownames(resNoFilt)
# 
# resNoFilt_final <- resNoFilt %>% as.data.frame
# 
# write.csv(resNoFilt_final, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cluster/Fig2_Microglia_figures_DAM_edgeR_DE_test_3_WT_vs_8_AD.csv")





resNoFilt_final <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cluster/Fig2_Microglia_figures_DAM_edgeR_DE_test_3_vs_8.csv", row.names = 1)

resNoFilt_final_up_0.5 <- resNoFilt_final[which(resNoFilt_final$logFC > 0.5 & resNoFilt_final$FDR < 0.01), ]
resNoFilt_final_up_0.5_sub <- resNoFilt_final_up_0.5[,c("logFC", "gene", "FDR")]

resNoFilt_final_down_0.5 <- resNoFilt_final[which(resNoFilt_final$logFC < -0.5 & resNoFilt_final$FDR < 0.01), ]
resNoFilt_final_down_0.5_sub <- resNoFilt_final_down_0.5[,c("logFC", "gene", "FDR")]

resNoFilt_final_0.5_sub <- rbind(resNoFilt_final_up_0.5_sub, resNoFilt_final_down_0.5_sub)

resNoFilt_final_DAM <- resNoFilt_final[resNoFilt_final$gene %in% DAM_markers3$`Gene name`, ]
resNoFilt_final_DAM_sub <- resNoFilt_final_DAM[,c("logFC", "gene", "FDR")]

resNoFilt_final_sub <- resNoFilt_final[,c("logFC", "gene", "FDR")]

##########################
# DAM from Ido Amit paper
##########################

#####################################
# Table 2
# Cluster 3 vs 1 for Figure 1
#####################################

DAM_markers <- read_xlsx("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DAM_markers/mmc2.xlsx")
DAM_markers_UP <- DAM_markers[DAM_markers$`up/down` == 1, ]
DAM_markers_DOWN <- DAM_markers[DAM_markers$`up/down` == -1, ]

DAM_markers_UP$logFC <- log2(DAM_markers_UP$`Microglia3  (average UMI count)`/DAM_markers_UP$`Microglia1 (average UMI count)`)
colnames(DAM_markers_UP) <- c("gene", "-log10(p-value)", "MG1_avg_UMI", "MG2_avg_UMI", "MG3_avg_UMI", "up/down", "logFC_DAM")
DAM_markers_UP <- DAM_markers_UP[!is.infinite(DAM_markers_UP$logFC_DAM), ]


DAM_markers_DOWN$logFC <- log2(DAM_markers_DOWN$`Microglia3  (average UMI count)`/DAM_markers_DOWN$`Microglia1 (average UMI count)`)
colnames(DAM_markers_DOWN) <- c("gene", "-log10(p-value)", "MG1_avg_UMI", "MG2_avg_UMI", "MG3_avg_UMI", "up/down", "logFC_DAM")

DAM_All_markers <- rbind(DAM_markers_UP, DAM_markers_DOWN)

# DAM from Ido data
DAM_UP.df <- right_join(resNoFilt_final_sub, DAM_markers_UP, by = "gene")

DAM_DOWN.df <- right_join(resNoFilt_final_sub, DAM_markers_DOWN, by = "gene")

ggplot(DAM_UP.df, aes(x = logFC_DAM, y = logFC)) + geom_point()

ggplot(DAM_DOWN.df, aes(x = logFC_DAM, y = logFC)) + geom_point()

DAM_All.df <- rbind(DAM_UP.df, DAM_DOWN.df)

DAM_All.df <- left_join(DAM_All.df, clusters.df, by = "gene")

DAM_All.df$norm_exp <- ifelse(DAM_All.df$`up/down` == "1", DAM_All.df$cluster8, DAM_All.df$cluster3)

DAM_All.df_sig <- DAM_All.df[DAM_All.df$FDR < 0.05, ]

DAM_All.df_1 <- DAM_All.df[which(DAM_All.df$logFC_DAM > 6.5 & DAM_All.df$logFC < 0.2), ]
DAM_All.df_1 <- DAM_All.df_1[DAM_All.df_1$gene != "Atp6v0d2", ]
DAM_All.df_2 <- DAM_All.df[which(DAM_All.df$logFC_DAM > 6 & DAM_All.df$logFC > 1.2), ]
DAM_All.df_3 <- DAM_All.df[which(DAM_All.df$logFC_DAM > 2 & DAM_All.df$logFC > 2), ]
DAM_All.df_4 <- DAM_All.df[which(DAM_All.df$logFC_DAM > 6.5 & DAM_All.df$logFC > 0.2 & DAM_All.df$logFC < 1), ]
DAM_All.df_5 <- DAM_All.df[which(DAM_All.df$logFC_DAM < -1.5 & DAM_All.df$logFC < -0.7 ), ]
DAM_All.df_5 <- DAM_All.df_5[DAM_All.df_5$gene != "P2ry12", ]
DAM_All.df_6 <- DAM_All.df[DAM_All.df$gene == "P2ry12", ]


ggplot(DAM_All.df, aes(x = logFC_DAM, y = logFC, color = norm_exp)) + 
  theme_bw() + 
  geom_point(alpha = 0.7, size = 2) + 
  #geom_point(data = DAM_All.df_sig, aes(x = logFC_DAM, y = logFC), alpha = 0.7, shape = 21, size = 2, stroke = 0.3, color = "black") + 
  xlim(c(-4, 9)) + ylim(c(-4, 9)) +
  ggtitle("S.table2 417 genes from  490 genes of DAM vs homeostatic from Ido Amit (77 genes are not detected in our data)") +
  scale_color_gradient(low="darkblue", high="red", limits = c(0,3), oob = scales::squish) + 
  geom_text_repel(data = DAM_All.df_1, aes(label = gene), family = "helvetica", color = "black", force = 3, nudge_x = 0.5,nudge_y = -1) + 
  geom_text_repel(data = DAM_All.df_2, aes(label = gene), family = "helvetica", color = "black", force = 3, nudge_y = 0.5) + 
  geom_text_repel(data = DAM_All.df_3, aes(label = gene), family = "helvetica", color = "black", force = 3, nudge_y = 0.5) + 
  geom_text_repel(data = DAM_All.df_4, aes(label = gene), family = "helvetica", color = "black", force = 1, nudge_y = 0.2) + 
  geom_text_repel(data = DAM_All.df_5, aes(label = gene), family = "helvetica", color = "black", force = 2, nudge_y = -0.5)

# GO analysis

genes <- DAM_All.df[is.na(DAM_All.df$logFC), ]$gene
table(seurat_mouse_genes[seurat_mouse_genes$mgi_symbol %in% DAM_All.df$gene, ]$gene_biotype)

geneUniverse <- rownames(data9set_cleaned_MG.SO@assays$RNA@counts[rowSums(data9set_cleaned_MG.SO@assays$RNA@counts) > 0, ]) %>% as.character
geneUniverse <- c(geneUniverse, DAM_All.df$gene)
topGOResults.df <- cluster_GO(genes, geneUniverse)


# For supplementary figure

g1 <- ggplot(DAM_All.df, aes(x = logFC_DAM, y = logFC, color = norm_exp)) + 
  theme_bw() + 
  theme(axis.title = element_text(size = 15, color = "black", family = "helvetica"),
        axis.text = element_text(size = 15, color = "black", family = "helvetica"),
        plot.title = element_text(size = 15, color = "black", family = "helvetica", hjust = 0.5),
        legend.title = element_text(size = 15, color = "black", family = "helvetica"),
        legend.text = element_text(size = 15, color = "black", family = "helvetica")) + 
  geom_point(alpha = 0.7, size = 2) + 
  #geom_point(data = DAM_All.df_sig, aes(x = logFC_DAM, y = logFC), alpha = 0.7, shape = 21, size = 2, stroke = 0.3, color = "black") + 
  xlim(c(-4, 9)) + ylim(c(-4, 9)) +
  xlab("logFC (DAM vs Homeostatic Microglia)") + ylab("logFC (Cluster 8 vs Cluster 3)") + 
  ggtitle("418 genes from 490 genes of Keren-Shaul et al., S.Table2") +
  scale_color_gradient(name = "Nor. Exp", low="darkblue", high="red", limits = c(0,3), oob = scales::squish) + 
  geom_text_repel(data = DAM_All.df_1, aes(label = gene), family = "helvetica", color = "black", force = 3, nudge_x = 0.1, nudge_y = -1) + 
  geom_text_repel(data = DAM_All.df_2, aes(label = gene), family = "helvetica", color = "black", force = 3, nudge_y = 0.5) + 
  geom_text_repel(data = DAM_All.df_3, aes(label = gene), family = "helvetica", color = "black", force = 3, nudge_y = 0.5) + 
  geom_text_repel(data = DAM_All.df_4, aes(label = gene), family = "helvetica", color = "black", force = 1, nudge_y = 0.2) + 
  geom_text_repel(data = DAM_All.df_5, aes(label = gene), family = "helvetica", color = "black", force = 2, nudge_y = -0.5) + 
  geom_text_repel(data = DAM_All.df_6, aes(label = gene), family = "helvetica", color = "black", force = 2, nudge_x = -0.6, nudge_y = 0.4)


                      

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_DAM_2nd_scatter_volcano_DAM_scatter.png",
       plot = g1,
       scale = 1, width = 9, height = 6, units = "in", device = "png",
       dpi = 300)





##################
# our DAM markers
##################

cor(DAM_All.df[!is.na(DAM_All.df$logFC), ]$logFC, DAM_All.df[!is.na(DAM_All.df$logFC), ]$logFC_DAM, method = c("spearman"))
our_to_DAM.df <- left_join(resNoFilt_final_0.5_sub, DAM_All_markers, by = "gene")
our_to_DAM.df <- left_join(our_to_DAM.df, clusters.df, by = "gene")
our_to_DAM.df$norm_exp <- ifelse(our_to_DAM.df$`up/down` == "1", our_to_DAM.df$cluster8, our_to_DAM.df$cluster3)
our_to_DAM.df1 <- our_to_DAM.df[which(our_to_DAM.df$logFC_DAM > 6), ] 
our_to_DAM.df2 <- our_to_DAM.df[which(our_to_DAM.df$logFC_DAM < -2), ] 
our_to_DAM.df2 <- our_to_DAM.df2[our_to_DAM.df2$gene != "P2ry12", ]
our_to_DAM.df3 <- our_to_DAM.df[our_to_DAM.df$gene == "P2ry12", ]

ggplot(our_to_DAM.df, aes(x = logFC_DAM, y = logFC, color = norm_exp)) + 
  geom_point(alpha = 0.7, size = 2) + 
  theme_bw() +
  xlim(c(-4, 8.5)) + ylim(c(-4, 8.5)) +
  ggtitle("our 303 logFC > |0.5| genes DAM vs homeostatic from Ido Amit (239 genes are not detected in DAM data)") + 
  scale_color_gradient(low="darkblue", high="red", limits = c(0,3), oob = scales::squish) + 
  geom_text_repel(data = our_to_DAM.df1, aes(label = gene), family = "helvetica", color = "black", force = 3, nudge_y = 0.5) + 
  geom_text_repel(data = our_to_DAM.df2, aes(label = gene), family = "helvetica", color = "black", force = 5)


cor(our_to_DAM.df[!is.na(our_to_DAM.df$logFC_DAM), ]$logFC, our_to_DAM.df[!is.na(our_to_DAM.df$logFC_DAM), ]$logFC_DAM, method = c("spearman"))


genes <- our_to_DAM.df[is.na(our_to_DAM.df$logFC_DAM), ]$gene
table(seurat_mouse_genes[seurat_mouse_genes$mgi_symbol %in% genes, ]$gene_biotype)
geneUniverse <- rownames(data9set_cleaned_MG.SO@assays$RNA@counts[rowSums(data9set_cleaned_MG.SO@assays$RNA@counts) > 0, ]) %>% as.character
topGOResults.df2 <- cluster_GO(genes, geneUniverse)

our_to_DAM_lnc.df <- our_to_DAM.df[our_to_DAM.df$gene %in% seurat_mouse_genes[seurat_mouse_genes$gene_biotype %in% "lncRNA", ]$mgi_symbol , ]

write.csv(our_to_DAM_lnc.df, file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_DAM_2nd_scatter_volcano_DAM_scatter_our_lncRNA.csv")

g1 <- ggplot(our_to_DAM.df, aes(x = logFC_DAM, y = logFC, color = norm_exp)) + 
  theme_bw() + 
  theme(axis.title = element_text(size = 15, color = "black", family = "helvetica"),
        axis.text = element_text(size = 15, color = "black", family = "helvetica"),
        plot.title = element_text(size = 15, color = "black", family = "helvetica", hjust = 0.5),
        legend.title = element_text(size = 15, color = "black", family = "helvetica"),
        legend.text = element_text(size = 15, color = "black", family = "helvetica")) + 
  geom_point(alpha = 0.7, size = 2) + 
  #geom_point(data = DAM_All.df_sig, aes(x = logFC_DAM, y = logFC), alpha = 0.7, shape = 21, size = 2, stroke = 0.3, color = "black") + 
  xlim(c(-4, 9)) + ylim(c(-4, 9)) +
  xlab("logFC (DAM vs Homeostatic Microglia)") + ylab("logFC (Cluster 8 vs Cluster 3)") + 
  ggtitle("68 genes from 303 genes of logFC > |0.5| Cluster 8 vs Cluster 3") +
  scale_color_gradient(name = "Nor. Exp", low="darkblue", high="red", limits = c(0,3), oob = scales::squish) + 
  geom_text_repel(data = our_to_DAM.df1, aes(label = gene), family = "helvetica", color = "black", force = 3, nudge_x = 0.5, nudge_y = 0.5) + 
  geom_text_repel(data = our_to_DAM.df2, aes(label = gene), family = "helvetica", color = "black", force = 5, nudge_x = -0.4, nudge_y = -0.2) + 
  geom_text_repel(data = our_to_DAM.df3, aes(label = gene), family = "helvetica", color = "black", force = 5, nudge_x = 0.3, nudge_y = 0.3)






#####################################
# Scatter plot
######################################
# From below
#data9set_cleaned_3_WT.SO <- subset(data9set_cleaned_3.SO, subset = sample %in% c("Ctrl"))
#data9set_cleaned_8_AD.SO <- subset(data9set_cleaned_8.SO, subset = sample %in% c("AD"))

# extract UMI counts and calculate average per each gene
MG_c3_count.df <- data.frame(Matrix::rowMeans(GetAssayData(data9set_cleaned_3_WT.SO, slot = "counts")))
MG_c8_count.df <- data.frame(Matrix::rowMeans(GetAssayData(data9set_cleaned_8_AD.SO, slot = "counts")))

MG_c3_count_cpm.df <- apply(MG_c3_count.df, 2, function(x) (x/sum(x)) * 1000000)
MG_c8_count_cpm.df <- apply(MG_c8_count.df, 2, function(x) (x/sum(x)) * 1000000)

MG_c3_count_cpm.df <- data.frame(Matrix::rowMeans(MG_c3_count_cpm.df))
MG_c8_count_cpm.df <- data.frame(Matrix::rowMeans(MG_c8_count_cpm.df))

MG_count_cpm.df <- cbind(MG_c3_count_cpm.df, MG_c8_count_cpm.df)
MG_count_cpm.df <- as.data.frame(MG_count_cpm.df)

# remove genes where all zero in samples
MG_count_cpm.df <- MG_count_cpm.df[rowSums(MG_count_cpm.df) > 0, ]
colnames(MG_count_cpm.df) <- c("cluster3", "cluster8")

# add pseudocount 0.001 and logarithm 
MG_count_cpm.df <- log2(MG_count_cpm.df + 1)
MG_count_cpm.df$gene <- rownames(MG_count_cpm.df)

MG_count_cpm_plus.df <- MG_count_cpm.df[MG_count_cpm.df$gene %in% DAM_markers_UP$gene, ]
MG_count_cpm_minus.df <- MG_count_cpm.df[MG_count_cpm.df$gene %in% DAM_markers_DOWN$gene, ]

MG_count_plus_label1.df <- MG_count_cpm.df[MG_count_cpm.df$gene %in% c("Myo1e", "Mamdc2","Axl", "Cst7", "Lpl", "Apoe", "Clec7a", "Ank",  "Spp1", "Igf1", "Csf1", "Gm11428"), ]
MG_count_plus_label2.df <- MG_count_cpm.df[MG_count_cpm.df$gene %in% c("Il12b"), ]
MG_count_plus_label3.df <- MG_count_cpm.df[MG_count_cpm.df$gene %in% c("Fam20c", "Itgax", "Cybb"), ]

MG_count_mius_label1.df <- MG_count_cpm.df[MG_count_cpm.df$gene %in% c("Malat1","Tmem119", "P2ry12", "Serinc3", "Cx3cr1"), ]

g1 <- ggplot(MG_count_cpm.df, aes(x = cluster3, y = cluster8)) + geom_point(size = 1, color = "#666666", alpha = 0.2) + 
  geom_point(data = MG_count_cpm_plus.df,  aes(x = cluster3, y = cluster8), color = "red", size = 1.4, alpha = 0.5) +
  geom_point(data = MG_count_cpm_minus.df,  aes(x = cluster3, y = cluster8), color = "blue", size = 1.4, alpha = 0.5) +
  ggtitle("cluster3 WT vs cluster8 APPPS1 in Microglia") + 
  theme_classic() + xlab("log2 (average CPM + 1) cluster3 WT") +
  ylab("log2 (average CPM + 1) cluster8 APPPS1") +
  theme(plot.title = element_text(size = 15, hjust = 0.5, family = "helvetica"),
        axis.title.x = element_text(size = 15, family = "helvetica", color = "#E69F00"),
        axis.title.y = element_text(size = 15, family = "helvetica", color = "#56B4E9"),
        axis.text = element_text(size = 15, color = "black", family = "helvetica")) +
  geom_abline(slope = 1, color="#FFFFCC", linetype="dashed", size=.8) + 
  geom_text_repel(data = MG_count_plus_label1.df, aes(x = cluster3, y = cluster8), label = MG_count_plus_label1.df$gene, force = 20, nudge_y = 1) +
  geom_text_repel(data = MG_count_plus_label2.df, aes(x = cluster3, y = cluster8), label = MG_count_plus_label2.df$gene, nudge_x = -0.5) +
  geom_text_repel(data = MG_count_plus_label3.df, aes(x = cluster3, y = cluster8), label = MG_count_plus_label3.df$gene, force = 10, nudge_y = 0.5, nudge_x = 0.5) + 
  geom_text_repel(data = MG_count_mius_label1.df, aes(x = cluster3, y = cluster8), label = MG_count_mius_label1.df$gene, force = 10, nudge_y = -0.7, nudge_x = 0.12) 
  

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_DAM_2nd_scatter_Cluster8_3.png",
       plot = g1,
       scale = 1, width = 7.5, height = 6.5, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_DAM_2nd_scatter_Cluster8_3.svg",
       plot = g1,
       scale = 1, width = 7.5, height = 6.5, units = "in", device = "svg",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_figures_DAM_2nd_scatter_Cluster8_3.pdf",
       plot = g1,
       scale = 1, width = 7.5, height = 6.5, units = "in", device = cairo_pdf,
       dpi = 300)




#####################################
# Not used

# Table 3
# DAM vs homeostatic for Figure 2
#####################################
DAM_markers2 <- read_xlsx("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DAM_markers/mmc3.xlsx")

# Remove NaN
DAM_markers2 <- DAM_markers2[complete.cases(DAM_markers2), ]
DAM_markers2 <- DAM_markers2[which(!DAM_markers2$`Fold-change (DAM to homeostatic microglia)` %in% "NaN" ), ]
DAM_markers2 <- DAM_markers2[which(!DAM_markers2$`DAM FDR p-value` %in% "NaN" ), ]

# -log10(FDR) > 2
DAM_markers2 <- DAM_markers2[DAM_markers2$`DAM FDR p-value` > 2, ]
# homeostatic > 0
DAM_markers2 <- DAM_markers2[DAM_markers2$`Average UMI per cell (homeostatic microglia)` > 0, ]

# foldchange -> logfoldchange 
DAM_markers2$`Fold-change (DAM to homeostatic microglia)` <- DAM_markers2$`Fold-change (DAM to homeostatic microglia)` %>% as.numeric
# only up 
DAM_markers2 <- DAM_markers2[DAM_markers2$`Fold-change (DAM to homeostatic microglia)` > 0, ]

DAM_markers2$logFC <- log(DAM_markers2$`Fold-change (DAM to homeostatic microglia)`)

# log foldchange > 0.5
DAM_markers3 <- DAM_markers2[DAM_markers2$logFC > 0.5, ]
colnames(DAM_markers3) <- c("gene", "Avg_UMI_homeo", "Fold_change", "-log10(p-value)", "FDR_DAM", "logFC_DAM")


# merge two results
# DAM from Ido data
all.df <- full_join(resNoFilt_final_DAM_sub, DAM_markers3, by = "gene")
ggplot(all.df, aes(x = logFC_DAM, y = logFC)) + geom_point() + xlim(c(-0.2, 2.5)) + ylim(c(-0.2, 2.5))



# DAM from our data
DAM_markers2 <- read_xlsx("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DAM_markers/mmc3.xlsx")

# Remove NaN
DAM_markers2 <- DAM_markers2[complete.cases(DAM_markers2), ]
DAM_markers2 <- DAM_markers2[which(!DAM_markers2$`Fold-change (DAM to homeostatic microglia)` %in% "NaN" ), ]
DAM_markers2 <- DAM_markers2[which(!DAM_markers2$`DAM FDR p-value` %in% "NaN" ), ]
# foldchange -> logfoldchange 
DAM_markers2$`Fold-change (DAM to homeostatic microglia)` <- DAM_markers2$`Fold-change (DAM to homeostatic microglia)` %>% as.numeric

DAM_markers2$logFC <- log(DAM_markers2$`Fold-change (DAM to homeostatic microglia)`)

colnames(DAM_markers2) <- c("gene", "Avg_UMI_homeo", "Fold_change", "-log10(p-value)", "FDR_DAM", "logFC_DAM")

all.df2 <- left_join(resNoFilt_final_up_0.5_sub, DAM_markers2, by = "gene")
ggplot(all.df2, aes(x = logFC_DAM, y = logFC)) + geom_point() + xlim(c(-5, 2.5)) + ylim(c(0, 2.5))


##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig4_Volcano_Plots_for_All_cell_type_session_info.txt")


