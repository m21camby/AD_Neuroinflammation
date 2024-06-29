#! /bin/env RScript
# written by SJK at 24. Mar. 2020

.libPaths(c("/home/skim/R/usr_lib", .libPaths()))

library(Seurat)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)
library(tidyr)
library(grid)
library(ggrepel)
source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")


seurat_mouse_genes <- read.table(file = "/data/rajewsky/home/skim/Microglia_Heppner/Bulk_transcriptome/R_Scripts/20200312_Bulk_Correlation_mouse_genes_bio_type_biomaRt.txt")

seurat_mouse_genes$gene_biotype_global <- ifelse(seurat_mouse_genes$gene_biotype %in% c("IG_V_gene", "IG_C_gene"), "IG_gene",
                                                 ifelse(seurat_mouse_genes$gene_biotype %in% c("pseudogene",
                                                                                               "unprocessed_pseudogene",
                                                                                               "processed_pseudogene",
                                                                                               "transcribed_unprocessed_pseudogene",
                                                                                               "transcribed_processed_pseudogene",
                                                                                               "unitary_pseudogene", 
                                                                                               "transcribed_unitary_pseudogene",
                                                                                               "polymorphic_pseudogene", 
                                                                                               "TR_V_pseudogene"), "pseudogene",
                                                        ifelse(seurat_mouse_genes$gene_biotype %in% c("TR_V_gene", "TR_C_gene"), "TR_gene",
                                                               ifelse(seurat_mouse_genes$gene_biotype %in% c("protein_coding"), "protein_coding",
                                                                      ifelse(seurat_mouse_genes$gene_biotype %in% c("TEC"), "TEC",
                                                                             ifelse(seurat_mouse_genes$gene_biotype %in% c("lncRNA"), "lncRNA", "ncRNA"))))))

seurat_mouse_genes_stat.df <- as.data.frame(table(seurat_mouse_genes$gene_biotype_global))

seurat_mouse_genes_stat.df$Var1 <- factor(seurat_mouse_genes_stat.df$Var1, levels = c("IG_gene", "TR_gene", "ncRNA", "TEC", "pseudogene", "lncRNA", "protein_coding"))




data9set_cleaned_sum.df_summary <- read.table("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200430_lncRNAs_2nd_analysis_data9set_cleaned_sum.df_summary", sep = "\t")


data9set_cleaned_sum.df_summary$gene_biotype_global <- factor(data9set_cleaned_sum.df_summary$gene_biotype_global, levels = c("IG_gene", "TR_gene", "ncRNA", "TEC", "pseudogene", "lncRNA", "protein_coding"))


g1 <- ggplot(seurat_mouse_genes_stat.df, aes(x = Var1, y = Freq, fill = Var1)) +
  geom_col(stat="identity") +  
  theme_classic() + coord_flip() + scale_y_log10(expand = c(0,0)) + 
  scale_fill_manual(values=c("darkred","#6633CC","#56B4E9","#666666","darkgreen","#E69F00", "#003366")) + 
  ggtitle("Total number of detected genes by biotypes") + 
  xlab("RNA biotypes") + ylab("number of genes") +
  theme(plot.title = element_text(size = 15),
        axis.text = element_text(size = 12, color = "black"), 
        axis.title = element_text(size = 15, color = "black"),
        legend.position = "none") +
  geom_text(aes(x = Var1, y = Freq, label = Freq),
            position = position_dodge(width = 0), hjust = 1.5, size = 5, color = "white")


g2 <- ggplot(data9set_cleaned_sum.df_summary, aes(fill=gene_biotype_global, y=percent, x="gene_biotype")) + 
  geom_bar(position="fill", stat="identity") + 
  theme_classic() + 
  geom_text(aes(y=0.95, label="lncRNA: 13.4%"), vjust=1.6, color="white", size=5) + 
  geom_text(aes(y=0.5, label="protein coding:"), vjust=1.6, color="white", size=5) + 
  geom_text(aes(y=0.42, label="85.7%"), vjust=1.6, color="white", size=5) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,1.2), breaks = c(0,1)) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.title = element_blank(), 
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12)) + 
  annotate(geom="text", x=1, y=1.08, label="rest: 0.9%", color="black", size = 5) + 
  scale_fill_manual(values=c("darkred","#6633CC","#56B4E9","#666666","darkgreen","#E69F00", "#003366")) +
  geom_segment(aes(x = 1, y = 0.999, xend = 1, yend = 1.05),color='black',size=0.5) 

blank <- grid.rect(gp=gpar(col="white"))


g12 <- arrangeGrob(g1, blank, g2, ncol = 3, widths = c(0.9, 0.15, 0.55))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig7_lncRNA_biotypes_plot.png",
       plot = g12,
       scale = 1, width = 11, height = 4, units = "in", device = "png",
       dpi = 300)


##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig7_lncRNA_biotypes_session_info.txt")


