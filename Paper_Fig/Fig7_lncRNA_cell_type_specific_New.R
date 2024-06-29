#! /bin/env RScript
# written by SJK at 10. Feb. 2021
# This file is for specific lncRNA figure
# adding a supplementary figure of an overview of lncRNA. from Shirin's email

library(Seurat)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)
library(tidyr)
library(grid)
library(ggrepel)
source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")


##############################
# 1. load biotype information
##############################

seurat_mouse_genes <- read.table(file = "/data/rajewsky/home/skim/Microglia_Heppner/Bulk_transcriptome/R_Scripts/20200312_Bulk_Correlation_mouse_genes_bio_type_biomaRt.txt")

# narrow down biotype
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


#################
# 2. bulk data load
#################
DGEm <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/Bulk_transcriptome/R_Scripts/20200311_Bulk_Analysis_Modified_Final_DF.csv", row.names =1, sep = ",", check.names = FALSE)

DGEm$gene <- rownames(DGEm)
DGEm2 <- data.frame(WT = DGEm$`WT-1` + DGEm$`WT-2` + DGEm$`WT-3`, AD = DGEm$`AD-1` + DGEm$`AD-2`+ DGEm$`AD-3`, ADp40KO = DGEm$`ADp40KO-1` + DGEm$`ADp40KO-2` + DGEm$`ADp40KO-3`, row.names = rownames(DGEm))
DGEm2$mgi_symbol <- rownames(DGEm2)


# only lncRNAs
data_lncRNA_sum_final.df <- DGEm2[DGEm2$mgi_symbol %in% seurat_mouse_genes[seurat_mouse_genes$gene_biotype_global %in% "lncRNA", ]$mgi_symbol , ]

# 1st high expressed (WT is higher in 15 normalized counts)
data_lncRNA_cleaned_sum_final.df <- data_lncRNA_sum_final.df[data_lncRNA_sum_final.df$WT > 15, ]


#################
# 3. snRNA avg data load
#################
avg.df <- read.table(file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200430_lncRNAs_2nd_analysis_avg.txt", sep = "\t")

avg.df$gene <- rownames(avg.df)

# Higly expressed in bulk
avg_lncRNA.df <- avg.df[avg.df$gene %in% data_lncRNA_cleaned_sum_final.df$mgi_symbol, ] 
avg_lncRNA.df$gene <- NULL
avg_lncRNA.df <- avg_lncRNA.df[rowSums(avg_lncRNA.df) !=0, ]

###########################
# 4. calculate z-score (high)
###########################
avg_lncRNA_z_score.df <- as.data.frame(t(apply(avg_lncRNA.df, 1, function(x) (x - mean(x)) / sd(x))))

avg_lncRNA_top10_z_score.df <- rbind(avg_lncRNA_z_score.df[rev(order(avg_lncRNA_z_score.df$Excitatory_Neuron))[1:10],],
                                     avg_lncRNA_z_score.df[rev(order(avg_lncRNA_z_score.df$Inhibitory_Neuron))[1:10],],
                                     avg_lncRNA_z_score.df[rev(order(avg_lncRNA_z_score.df$Oligodendrocytes))[1:10],],
                                     avg_lncRNA_z_score.df[rev(order(avg_lncRNA_z_score.df$OPC))[1:10],],
                                     avg_lncRNA_z_score.df[rev(order(avg_lncRNA_z_score.df$Microglia))[1:10],],
                                     avg_lncRNA_z_score.df[rev(order(avg_lncRNA_z_score.df$Astrocytes))[1:10],])


avg_lncRNA_top10_z_score.df$gene <- rownames(avg_lncRNA_top10_z_score.df)
avg_lncRNA_top10_z_score.df <- avg_lncRNA_top10_z_score.df[,c(7,1:6)]


# for ggplot modification
avg_lncRNA_top10_z_score.df2 <- gather(avg_lncRNA_top10_z_score.df, Cell_type, z_score, Excitatory_Neuron:Astrocytes)

avg_lncRNA_top10_z_score.df2$gene <- factor(avg_lncRNA_top10_z_score.df2$gene, levels = rev(rownames(avg_lncRNA_top10_z_score.df)))

avg_lncRNA_top10_z_score.df2$Cell_type <- factor(avg_lncRNA_top10_z_score.df2$Cell_type, levels = 
                                                   c("Excitatory_Neuron", "Inhibitory_Neuron",
                                                     "Oligodendrocytes","OPC",
                                                     "Microglia", "Astrocytes"))


g1 <- ggplot(data = avg_lncRNA_top10_z_score.df2, mapping = aes(x = Cell_type, y = gene, fill = z_score)) +
  geom_tile()+  ggtitle("Cell-type specific lncRNAs") + 
  theme_classic() + theme(axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5),
                          axis.title.x = element_blank(), 
                          axis.text.x = element_text(angle = 90, vjust = 0, size = 12, color = "black"), 
                          axis.text.y = element_text(size = 8, color = "black"),
                          axis.line = element_line(color = "white"),
                          axis.ticks = element_line(color = "white"),
                          legend.title = element_blank(),
                          legend.text = element_text(size = 11, color = "black")) + 
  scale_fill_viridis_c(limits=c(-2, 2), na.value = "#FDE725FF", labels = c("-2", " ", " "," ", "2"))


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig7_lncRNA_cell_type_specific_plot_New.pdf",
       plot = g1,
       scale = 1, width = 5, height = 7.5, units = "in", device = cairo_pdf,
       dpi = 300)


#####################
# 5. Interesting AD lncRNA
#####################

avg_lncRNA_z_score.df2 <- avg_lncRNA_z_score.df[rownames(avg_lncRNA_z_score.df) %in% c("Neat1", "Pvt1", "Rmst", "Sox2ot", "Snhg1"), ]

avg_lncRNA_z_score.df2$gene <- rownames(avg_lncRNA_z_score.df2)
avg_lncRNA_z_score.df2 <- avg_lncRNA_z_score.df2[,c(7,1:6)]

# for ggplot modification
avg_lncRNA_z_score.df3 <- gather(avg_lncRNA_z_score.df2, Cell_type, z_score, Excitatory_Neuron:Astrocytes)

avg_lncRNA_z_score.df3$gene <- factor(avg_lncRNA_z_score.df3$gene, levels = rev(c("Snhg1", "Pvt1", "Neat1","Sox2ot", "Rmst")))

avg_lncRNA_z_score.df3$Cell_type <- factor(avg_lncRNA_z_score.df3$Cell_type, levels = 
                                                   c("Excitatory_Neuron", "Inhibitory_Neuron",
                                                     "Oligodendrocytes","OPC",
                                                     "Microglia", "Astrocytes"))

g2 <- ggplot(data = avg_lncRNA_z_score.df3, mapping = aes(x = Cell_type, y = gene, fill = z_score)) +
  geom_tile()+  ggtitle("Cell-type specific AD-related lncRNAs") + 
  theme_classic() + theme(axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5),
                          axis.title.x = element_blank(), 
                          axis.text.x = element_text(angle = 90, vjust = 0, size = 12, color = "black"), 
                          axis.text.y = element_text(size = 8, color = "black"),
                          axis.line = element_line(color = "white"),
                          axis.ticks = element_line(color = "white"),
                          legend.title = element_blank(),
                          legend.text = element_text(size = 11, color = "black")) + 
  scale_fill_viridis_c(limits=c(-2, 2), na.value = "#FDE725FF", labels = c("-2", " ", " "," ", "2"))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig7_lncRNA_cell_type_specific_AD_plot.pdf",
       plot = g2,
       scale = 1, width = 5, height = 2.5, units = "in", device = cairo_pdf,
       dpi = 300)


##################
# 6. All (60 specific + AD-related)
##################

avg_lnc_all.df <- rbind(avg_lncRNA_top10_z_score.df, avg_lncRNA_z_score.df2)

avg_lnc_all.df2 <- gather(avg_lnc_all.df, Cell_type, z_score, Excitatory_Neuron:Astrocytes)

avg_lnc_all.df2$gene <- factor(avg_lnc_all.df2$gene, levels = c("Rmst", avg_lnc_all.df2[c(51:60), ]$gene, 
                                                                avg_lnc_all.df2[c(41:50), ]$gene, "Pvt1", 
                                                                avg_lnc_all.df2[c(31:40), ]$gene, "Sox2ot",
                                                                avg_lnc_all.df2[c(21:30), ]$gene, "Neat1",
                                                                avg_lnc_all.df2[c(11:20), ]$gene, "Snhg1",
                                                                avg_lnc_all.df2[c(1:10), ]$gene))
                                 
avg_lnc_all.df2$Cell_type <- factor(avg_lnc_all.df2$Cell_type, levels = 
                                             c("Excitatory_Neuron", "Inhibitory_Neuron",
                                               "Oligodendrocytes","OPC",
                                               "Microglia", "Astrocytes"))


g3 <- ggplot(data = avg_lnc_all.df2, mapping = aes(x = Cell_type, y = gene, fill = z_score)) +
  geom_tile()+  ggtitle("Cell-type specific lncRNAs") + 
  theme_classic() + theme(axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5),
                          axis.title.x = element_blank(), 
                          axis.text.x = element_text(angle = 90, vjust = 0, size = 12, color = "black"), 
                          axis.text.y = element_text(size = 8, color = "black"),
                          axis.line = element_line(color = "white"),
                          axis.ticks = element_line(color = "white"),
                          legend.title = element_blank(),
                          legend.text = element_text(size = 11, color = "black")) + 
  scale_fill_viridis_c(limits=c(-2, 2), na.value = "#FDE725FF", labels = c("-2", " ", " "," ", "2"))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig7_lncRNA_cell_type_specific_plot_New_All.pdf",
       plot = g3,
       scale = 1, width = 5, height = 8, units = "in", device = cairo_pdf,
       dpi = 300)


#####################
# 7. All + Rmst only
#####################
avg_lnc_all.df <- rbind(avg_lncRNA_top10_z_score.df, avg_lncRNA_z_score.df2[avg_lncRNA_z_score.df2$gene %in% "Rmst", ])

avg_lnc_all.df2 <- gather(avg_lnc_all.df, Cell_type, z_score, Excitatory_Neuron:Astrocytes)

avg_lnc_all.df2$gene <- factor(avg_lnc_all.df2$gene, levels = c("Rmst", avg_lnc_all.df2[c(51:60), ]$gene, 
                                                                avg_lnc_all.df2[c(41:50), ]$gene,
                                                                avg_lnc_all.df2[c(31:40), ]$gene,
                                                                avg_lnc_all.df2[c(21:30), ]$gene,
                                                                avg_lnc_all.df2[c(11:20), ]$gene,
                                                                avg_lnc_all.df2[c(1:10), ]$gene))

avg_lnc_all.df2$Cell_type <- factor(avg_lnc_all.df2$Cell_type, levels = 
                                      c("Excitatory_Neuron", "Inhibitory_Neuron",
                                        "Oligodendrocytes","OPC",
                                        "Microglia", "Astrocytes"))


g4 <- ggplot(data = avg_lnc_all.df2, mapping = aes(x = Cell_type, y = gene, fill = z_score)) +
  geom_tile()+  ggtitle("Cell-type specific lncRNAs") + 
  theme_classic() + theme(axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5),
                          axis.title.x = element_blank(), 
                          axis.text.x = element_text(angle = 90, vjust = 0, size = 12, color = "black"), 
                          axis.text.y = element_text(size = 8, color = "black"),
                          axis.line = element_line(color = "white"),
                          axis.ticks = element_line(color = "white"),
                          legend.title = element_blank(),
                          legend.text = element_text(size = 11, color = "black")) + 
  scale_fill_viridis_c(limits=c(-2, 2), na.value = "#FDE725FF", labels = c("-2", " ", " "," ", "2"))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig7_lncRNA_cell_type_specific_plot_New_Rmst_add.pdf",
       plot = g4,
       scale = 1, width = 5, height = 7.5, units = "in", device = cairo_pdf,
       dpi = 300)

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig7_lncRNA_cell_type_specific_New_session_info.txt")

