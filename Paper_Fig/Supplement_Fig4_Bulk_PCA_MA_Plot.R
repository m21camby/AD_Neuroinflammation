
#! /bin/env RScript
# written by SJK at 18. Dec. 2020

library(biomaRt)
library(DESeq2)
library(ggcorrplot)
library(corrplot)
library("Hmisc")
library(ggplot2)
library(plotly)
library(gridExtra)
library(grid)
library(ggrepel)
library(dplyr)
library(viridis)

source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/Remove_duplicated_genes.R")
source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

DGEm <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/Bulk_transcriptome/Bulk_by_ENSEMBL/Final_DF.csv", row.names = 1, check.names = FALSE, sep = "\t")


# remove genes where zero count is more than 2 out of 9 samples
DGEm <- DGEm[apply(DGEm == 0, 1, sum) < 7, ]
colnames(DGEm) <- c("ADp40KO_2", "WT_3", "ADp40KO_1", "WT_2", "WT_1", "AD_2", "AD_1", "ADp40KO_3", "AD_3")


##############################################
# extract gene information
##############################################
# mart <- useDataset("mmusculus_gene_ensembl", mart = useMart("ensembl"))
# 
# listAttributes(mart)
# 
# mouse_genes <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
#                   filters = 'ensembl_gene_id',
#                   values = rownames(DGEm),
#                   mart = mart)
# 
# 
# write.table(mouse_genes, file = "/data/rajewsky/home/skim/Microglia_Heppner/Bulk_transcriptome/R_Scripts/Supple_Fig3_Bulk_figures_mouse_genes.txt")
# mouse_genes <- read.table(file = "/data/rajewsky/home/skim/Microglia_Heppner/Bulk_transcriptome/R_Scripts/Supple_Fig3_Bulk_figures_mouse_genes.txt", stringsAsFactors = FALSE)
# 
# # remove duplicated ENSEMBL and mgi_symbol information
# mouse_genes_final.df <- Remove_duplicated_no_id_version(mouse_genes)
# 
# DGEm$ensembl_gene_id<- rownames(DGEm)
# DGEm <- left_join(DGEm, mouse_genes_final.df, by = "ensembl_gene_id")
# # remove NA
# DGEm$mgi_symbol <- ifelse(is.na(DGEm$mgi_symbol), DGEm$ensembl_gene_id, DGEm$mgi_symbol)
# 
# 
# rownames(DGEm) <- DGEm$mgi_symbol
# DGEm$mgi_symbol <- NULL
# DGEm$ensembl_gene_id <- NULL
# 
# write.csv(DGEm, "/data/rajewsky/home/skim/Microglia_Heppner/Bulk_transcriptome/R_Scripts/Supple_Fig3_Bulk_figures_Modified_Final_DF.csv")
##############################################

DGEm <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/Bulk_transcriptome/R_Scripts/Supple_Fig3_Bulk_figures_Modified_Final_DF.csv", row.names = 1)

sample_info <- data.frame(condition = c("ADp40KO","WT","ADp40KO","WT","WT","AD","AD","ADp40KO","AD"), row.names = names(DGEm))
dds <- DESeqDataSetFromMatrix(DGEm, colData = sample_info, design = ~ condition)
vsd <- vst(dds, blind = FALSE)

#plotPCA(vsd)
pc <- plotPCA(vsd, returnData =TRUE)
pc$group <- factor(pc$group, levels = c("WT", "AD", "ADp40KO"))

g1 <- ggplot(pc ,aes(x=PC1,y=PC2, color=group)) + geom_point(size = 5) +
  xlab("PC1: 54% variance") + 
  ylab("PC2: 31% variance") +
  #  geom_text(aes(label=name),hjust=.5, vjust=2, size = 3) +
  theme(axis.text = element_text(face = "bold", size = 15, color = "black", family = "helvetica"),
        axis.title = element_text(face ="bold", size = 15, color = "black", family = "helvetica"),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(),
        panel.grid.major = element_line(color="white"),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        legend.key = element_rect(fill = "transparent", colour = "transparent")) + 
  scale_color_manual(values = c("#440154FF", viridis(3)[2], "#CC9900"), labels = c("WT", "APPPS1", "APPPS1.Il12b.-/-"))
g1

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig4_Bulk_figures_PCA.pdf",
       plot = g1,
       scale = 1, width = 7, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)



###############################################
# DESeq2
###############################################
dds$condition <- relevel(dds$condition, ref = "WT")
dds.ds <- DESeq(dds)
dds.ds <- estimateSizeFactors(dds.ds)
resultsNames(dds.ds)

# Calculate WT normalized count
sizeFactors(dds.ds)
# WT_1: 1.1135411
# WT_2: 0.9270373
# WT_3: 0.8843946
DGEm_WT <- DGEm[,c("WT_1", "WT_2", "WT_3")]

DGEm_WT$WT_1_Norm <- DGEm_WT$WT_1 / 1.1135411
DGEm_WT$WT_2_Norm <- DGEm_WT$WT_2 / 0.9270373
DGEm_WT$WT_3_Norm <- DGEm_WT$WT_3 / 0.8843946



MA_mRNA_Plot <- function(DF, Min_exp = 1, y_axis = "exp"){
  
  mc1 <- subset(DF, log2FoldChange > 1 & baseMean > Min_exp & padj < .05)
  mc2 <- subset(DF, log2FoldChange < -1 & baseMean > Min_exp & padj < .05)
  
  ggplot(DF, aes(x = WT_norm, y = log2FoldChange)) +
    geom_point(size = 2, color = "darkgray", alpha = 0.2) +
    geom_point(data = mc1, color = "#990000", size = 2, alpha = 0.2) +
    geom_point(data = mc2, color = "#003366", size = 2, alpha = 0.2) +
    xlab("Normalized Counts (WT)") + ylab(paste0("Fold Change (log2)", y_axis)) +
    scale_x_continuous(trans='log10') +
    geom_hline(yintercept = 0, size = 0.3) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "darkred", size = 0.1) + geom_hline(yintercept = -1, linetype = "dashed", color = "darkred", size = 0.1) +
    theme(axis.title = element_text(face = "bold", size = 12, color = "black", family = "helvetica"), 
          axis.text = element_text(face = "bold", size = 12, color = "black",  family = "helvetica"), 
          plot.title = element_text(hjust = 0.5, color = "black",  family = "helvetica"),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))
}


#################################################
# AD vs WT
#################################################
res <- results(dds.ds, name="condition_AD_vs_WT", test="Wald")
res.df <- as.data.frame(res)
res.df$gene <- rownames(res.df)
res.df$WT_norm <- rowMeans(DGEm_WT[,c(4:6)])

res.df1 <- res.df[res.df$gene %in% c("Il12b"), ] 
res.df2 <-res.df[res.df$gene %in% c("Il12rb1"), ] 
res.df3 <- res.df[which(res.df$gene %in% c("Il12rb2")), ]
res.df4 <- res.df[which(res.df$gene %in% c("Il23a")), ]
res.df5 <- res.df[which(res.df$gene %in% c("Il12a")), ]


g2 <- MA_mRNA_Plot(res.df, Min_exp = 1, y_axis = "") + ggtitle("APPPS1 vs WT") + coord_cartesian(y = c(-10,10)) + 
  geom_point(data = res.df1, color = "#990000", size = 2, alpha = 1) +
  geom_text_repel(data = res.df1, aes(label = gene), family = "helvetica", color = "black", nudge_y = 1) + 
  geom_point(data = res.df2, color = "black", size = 2, alpha = 1) +
  geom_text_repel(data = res.df2, aes(label = gene), family = "helvetica", color = "black", nudge_y = -0.5, nudge_x = 0.5) + 
  geom_point(data = res.df3, color = "black", size = 2, alpha = 1) +
  geom_text_repel(data = res.df3, aes(label = gene), family = "helvetica", color = "black", nudge_y = -0.5, nudge_x =  -0.5) +
  geom_point(data = res.df4, color = "black", size = 2, alpha = 1) +
  geom_text_repel(data = res.df4, aes(label = gene), family = "helvetica", color = "black", nudge_y = -1, nudge_x =  -0)  +
  geom_point(data = res.df5, color = "black", size = 2, alpha = 1) +
  geom_text_repel(data = res.df5, aes(label = gene), family = "helvetica", color = "black", nudge_y = 1, nudge_x =  0) 


g2

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig4_Bulk_PCA_MA_Plot_AD_WT.pdf",
       plot = g2,
       scale = 1, width = 6, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

#################################################
# ADp40KO vs WT
#################################################
res <- results(dds.ds, name="condition_ADp40KO_vs_WT", test="Wald")
res.df <- as.data.frame(res)
res.df$gene <- rownames(res.df)
res.df$WT_norm <- rowMeans(DGEm_WT[,c(4:6)])

res.df1 <- res.df[res.df$gene %in% c("Il12b"), ] 
res.df2 <-res.df[res.df$gene %in% c("Il12rb1"), ] 
res.df3 <- res.df[which(res.df$gene %in% c("Il12rb2")), ]
res.df4 <- res.df[which(res.df$gene %in% c("Il23a")), ]
res.df5 <- res.df[which(res.df$gene %in% c("Il12a")), ]



g2 <- MA_mRNA_Plot(res.df, Min_exp = 1, y_axis = "") + ggtitle("APPPS1.Il12b.-/- vs WT") + coord_cartesian(y = c(-10,10)) + 
  geom_point(data = res.df1, color = "black", size = 2, alpha = 1) +
  geom_text_repel(data = res.df1, aes(label = gene), family = "helvetica", color = "black", nudge_y = 1.5) + 
  geom_point(data = res.df2, color = "black", size = 2, alpha = 1) +
  geom_text_repel(data = res.df2, aes(label = gene), family = "helvetica", color = "black", nudge_y = -0.5, nudge_x = 0.5) + 
  geom_point(data = res.df3, color = "black", size = 2, alpha = 1) +
  geom_text_repel(data = res.df3, aes(label = gene), family = "helvetica", color = "black", nudge_y = 0.7, nudge_x =  -0.3) +
  geom_point(data = res.df4, color = "black", size = 2, alpha = 1) +
  geom_text_repel(data = res.df4, aes(label = gene), family = "helvetica", color = "black", nudge_y = 1, nudge_x =  -0)  +
  geom_point(data = res.df5, color = "black", size = 2, alpha = 1) +
  geom_text_repel(data = res.df5, aes(label = gene), family = "helvetica", color = "black", nudge_y = -0.7, nudge_x = -0.5) 


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig4_Bulk_PCA_MA_Plot_ADp40KO_WT.pdf",
       plot = g2,
       scale = 1, width = 6, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

#################################################
# ADp40KO vs AD
#################################################
res <- results(dds.ds, contrast=c("condition", "ADp40KO", "AD"), test="Wald")
res.df <- as.data.frame(res)
res.df$gene <- rownames(res.df)

res.df$WT_norm <- rowMeans(DGEm_WT[,c(4:6)])

res.df1 <- res.df[res.df$gene %in% c("Il12b"), ] 
res.df2 <-res.df[res.df$gene %in% c("Il12rb1"), ] 
res.df3 <- res.df[which(res.df$gene %in% c("Il12rb2")), ]
res.df4 <- res.df[which(res.df$gene %in% c("Il23a")), ]
res.df5 <- res.df[which(res.df$gene %in% c("Il12a")), ]


g2 <- MA_mRNA_Plot(res.df, Min_exp = 1, y_axis = "") + ggtitle("APPPS1.Il12b.-/- vs APPPS1") + coord_cartesian(y = c(-10,10)) + 
  geom_point(data = res.df1, color = "#003366", size = 2, alpha = 1) +
  geom_text_repel(data = res.df1, aes(label = gene), family = "helvetica", color = "black", nudge_y = 1.5) + 
  geom_point(data = res.df2, color = "black", size = 2, alpha = 1) +
  geom_text_repel(data = res.df2, aes(label = gene), family = "helvetica", color = "black", nudge_y = -0.5, nudge_x = 0.5) + 
  geom_point(data = res.df3, color = "black", size = 2, alpha = 1) +
  geom_text_repel(data = res.df3, aes(label = gene), family = "helvetica", color = "black", nudge_y = 1, nudge_x =  -0.2) +
  geom_point(data = res.df4, color = "black", size = 2, alpha = 1) +
  geom_text_repel(data = res.df4, aes(label = gene), family = "helvetica", color = "black", nudge_y = 1, nudge_x =  -0)  +
  geom_point(data = res.df5, color = "black", size = 2, alpha = 1) +
  geom_text_repel(data = res.df5, aes(label = gene), family = "helvetica", color = "black", nudge_y = -0.7, nudge_x = -0.5) 


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig4_Bulk_PCA_MA_Plot_ADp40KO_AD.pdf",
       plot = g2,
       scale = 1, width = 6, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)


##########################################
# correlation with scRNA-seq
##########################################

gene_length.df <- read.table("/data/rajewsky/home/skim/Microglia_Heppner/Bulk_transcriptome/GTFtools_0.6.9/Mus_musculus.GRCm38.98_add_chr_no_filtered.isoform.length.txt", row.names = 1, header = TRUE)

mouse_genes <- read.table(file = "/data/rajewsky/home/skim/Microglia_Heppner/Bulk_transcriptome/R_Scripts/Supple_Fig3_Bulk_figures_mouse_genes.txt", stringsAsFactors = FALSE)
colnames(mouse_genes) <- c("gene", "mgi_symbol")

gene_length_gene.df <- left_join(gene_length.df, mouse_genes, by = "gene")

gene_length_gene.df <- gene_length_gene.df[!is.na(gene_length_gene.df$mgi_symbol), ]

# calculate median for each gene
gene_length_median.df <- gene_length_gene.df %>% group_by(mgi_symbol) %>% summarise(Median=median(length))

# Remove empty rows
gene_length_median.df <- gene_length_median.df[-1, ]

#########################
# DGEm and length combine
#########################
DGEm$mgi_symbol <- rownames(DGEm)

DGEm_length <- inner_join(DGEm, gene_length_median.df, by = c("mgi_symbol"))
rownames(DGEm_length) <- DGEm_length$mgi_symbol
DGEm_length$mgi_symbol <- NULL

DGEm_length_info.df <- data.frame(length = DGEm_length$Median, row.names = rownames(DGEm_length))
DGEm_length$Median <-NULL

tpm <- function(counts, lengths) {
  rate <- counts / lengths # RPK
  rate.df <- rate / sum(rate) * 1e6
  return(rate.df)
}

DGEm_length.df <- apply(DGEm_length, 2, function(x) tpm(x, DGEm_length_info.df$length)) %>% as.data.frame

DGEm_length.df <- DGEm_length.df + 1

# log2 
DGEm_length_log2.df <- as.data.frame(lapply(DGEm_length.df, log2), row.names = rownames(DGEm_length.df))

############################
# scRNA-seq cpm
############################
cpm.df <- read.table(file = "/data/rajewsky/home/skim/Microglia_Heppner/Bulk_transcriptome/R_Scripts/20200312_Bulk_Correlation_scRNA_data_cpm.txt")
#cpm <- apply(data9set_1.df,2, function(x) (x/sum(x))*1000000)

# 1. remove if non-expressed genes
cpm.df <- cpm.df[apply(cpm.df == 0, 1, sum) != 9, ]

# 3. rearrange column order
cpm.df <- cpm.df[, c(1,4,7,3,6,9,2,5,8)]

# 4. rename column name
colnames(cpm.df) <- c("WT_1", "WT_2", "WT_3", "APPPS1_1", "APPPS1_2", "APPPS1_3", "APPPS1.p40KO_1", "APPPS1.p40KO_2", "APPPS1.p40KO_3")

# 5. add pseudocount 1
cpm.df <- cpm.df + 1

# log2 
cpm_log2.df <- as.data.frame(lapply(cpm.df, log2), row.names = rownames(cpm.df))


###################################
# merged bulk and snRNA-seq
###################################

cpm_genotype.df <- data.frame(WT_snRNA = rowMeans(cpm_log2.df[, colnames(cpm_log2.df) %in% c("WT_1", "WT_2","WT_3")]),
                              APPPS1_snRNA = rowMeans(cpm_log2.df[, colnames(cpm_log2.df) %in% c("APPPS1_1", "APPPS1_2","APPPS1_3")]),
                              APPPS1.il12b_KO_snRNA = rowMeans(cpm_log2.df[, colnames(cpm_log2.df) %in% c("APPPS1.p40KO_1", "APPPS1.p40KO_2","APPPS1.p40KO_3")]))


tpm_genotype.df <- data.frame(WT_bulk = rowMeans(DGEm_length_log2.df[, colnames(DGEm_length_log2.df) %in% c("WT_1", "WT_2","WT_3")]),
                              APPPS1_bulk = rowMeans(DGEm_length_log2.df[, colnames(DGEm_length_log2.df) %in% c("AD_1", "AD_2","AD_3")]),
                              APPPS1.il12b_KO_bulk = rowMeans(DGEm_length_log2.df[, colnames(DGEm_length_log2.df) %in% c("ADp40KO_1", "ADp40KO_2","ADp40KO_3")]))


cpm_genotype.df$gene <- rownames(cpm_genotype.df)
tpm_genotype.df$gene <- rownames(tpm_genotype.df)

genotype.df <- inner_join(cpm_genotype.df, tpm_genotype.df, by = "gene")
rownames(genotype.df) <- genotype.df$gene
genotype.df$gene <- NULL


#################################
# correlation test
#################################

cor_res <- rcorr(as.matrix(genotype.df))

png(height=1400, width=1400, units="px", res=300, file="/data/rajewsky/home/skim/Microglia_Heppner/Bulk_transcriptome/R_Scripts/Supple_Fig3_Bulk_figures_corr.png", type = "cairo")

corrplot(cor_res$r, type = "upper", 
         tl.col = "black", tl.srt = 90, tl.cex = 0.7, number.cex=.6, number.digits = 3, 
         addCoef.col = "grey", order="original")

dev.off()

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/Bulk_transcriptome/R_Scripts/Supple_Fig3_Bulk_figures_session_info.txt")





