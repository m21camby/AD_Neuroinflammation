#! /bin/env RScript
# written by SJK at 31. Oct. 2020
# 14. Dec 2020
# I add ADp40KO vs WT in GO from last Friday meeting

library(gridExtra)
library(dplyr)
library(ggplot2)
library(tidyr)
library(Seurat)
library(ggrepel)
library(readxl)

source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

Vln_plot <- function(df = DF){
  ggplot(df, aes(x=pathway, y=logFC)) + 
    geom_violin(trim=FALSE, adjust = 3) + coord_flip() + geom_jitter(data = df[df$FDR < 0.01, ], aes(color = logFC > 0),  position=position_jitter(0.1), alpha = 0.3) + 
    theme(axis.title.y = element_blank()) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) + theme_minimal() + 
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 10, family = "helvetica", color = "black"),
          axis.text = element_text(size = 8, family = "helvetica", color = "black"),
          panel.grid.minor = element_blank())
}

###################
# MG Data loading
###################
MG_AD_WT_UP <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/GO_csv_files/edgeR_GO/cell_type_up_0.25/Microglia_AD_Ctrl.csv_GO.csv", row.names = 1)
MG_AD_WT_UP$group <- "APPPS1 vs WT Up"

MG_AD_WT_DOWN <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/GO_csv_files/edgeR_GO/cell_type_down_0.25/Microglia_AD_Ctrl.csv_GO.csv", row.names = 1)
MG_AD_WT_DOWN$group <- "WT vs APPPS1 Up"

MG_ADp40KO_WT_UP <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/GO_csv_files/edgeR_GO/cell_type_up_0.25/Microglia_ADp40KO_Ctrl.csv_GO.csv", row.names = 1)
MG_ADp40KO_WT_UP$group <- "APPPS1.Il12b-/- vs WT Up"

MG_ADp40KO_WT_DOWN <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/GO_csv_files/edgeR_GO/cell_type_down_0.25/Microglia_ADp40KO_Ctrl.csv_GO.csv", row.names = 1)
MG_ADp40KO_WT_DOWN$group <- "WT vs APPPS1.Il12b Up"


# merge
MG_AD_WT <- rbind(MG_AD_WT_UP, MG_AD_WT_DOWN, MG_ADp40KO_WT_UP, MG_ADp40KO_WT_DOWN)

g1 <- ggplot(MG_AD_WT[MG_AD_WT$Term %in% c("amyloid-beta binding", "regulation of neuroinflammatory response", 
                                           "response to lipoprotein particle", "inflammatory response", 
                                           "heparin binding", "cellular response to fatty acid",
                                           "positive regulation of ERK1 and ERK2 cas...", "amyloid-beta clearance",
                                           "SH3 domain binding", "positive regulation of MAPK cascade", 
                                           "ubiquitin protein ligase binding", "ubiquitin protein ligase activity"), ], 
             aes(x=group, y=Term,
                 colour=Fisher.elim,
                 size=gene_ratio)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="group", y="GO term", colour="p-value", size="gene_ratio") +
  #ggtitle("Microglia APPPS1 vs WT") + 
  theme_cowplot() + 
  theme(axis.text.x =  element_text(size = 10, angle = 45, vjust = 0.5, family = "helvetica", color = "black"),
        axis.text.y = element_text(size = 10, family = "helvetica", color = "black"),
        axis.title.y = element_text(size = 10, family = "helvetica", color = "black"),
        axis.title.x =  element_blank(),
        plot.title = element_text(size = 10, family = "helvetica", color = "black", hjust = 0.5)) + 
  scale_color_gradient(low="#003366", high="#FF9900", trans = 'reverse') + 
  guides(color = guide_colourbar(order = 1), 
         size = guide_legend(order = 2)) + guides(color = guide_colorbar(reverse=TRUE))

g1
##########################
# Vln plot
##########################

MG_AD_WT_DE <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/Microglia_AD_Ctrl.csv", row.names = 1)

genes <- MG_AD_WT[MG_AD_WT$Term %in% c("positive regulation of cell migration"), ]$genes %>% as.character
genes

MG_AD_WT_DE1 <- MG_AD_WT_DE[MG_AD_WT_DE$gene %in% c("Apbb2", "Apoe", "Ldlrad3", "Ldlr", "Gsap", "Fcgr2b", "Scarb1", "Trem2"), ]
MG_AD_WT_DE1$pathway <- "amyloid-beta binding"

MG_AD_WT_DE2 <- MG_AD_WT_DE[MG_AD_WT_DE$gene %in% c("Igf1", "Cst7", "Ldlr", "Trem2"), ]
MG_AD_WT_DE2$pathway <- "regulation of neuroinflammatory response"

MG_AD_WT_DE3 <- MG_AD_WT_DE[MG_AD_WT_DE$gene %in% c("Cd9", "Apoe", "Lpl", "Ldlr"), ]
MG_AD_WT_DE3$pathway <- "response to lipoprotein particle"

MG_AD_WT_DE4 <- MG_AD_WT_DE[MG_AD_WT_DE$gene %in% c("Igf1","Hdac9","Hif1a","Ccl4","Apoe","Otulin","Axl","Csf1","Lpl","Ccl3","Abcd2","Acer3","Cst7","Ldlr","Abhd12","Cd180","Ly86","Ccl6","Nfe2l2","Fcgr2b","Cdk19","Plgrkt","C3ar1","Trem2"), ]
MG_AD_WT_DE4$pathway <- "inflammatory response"

MG_AD_WT_DE5 <- MG_AD_WT_DE[MG_AD_WT_DE$gene %in% c("Nrp1", "Serpine2", "Apoe", "Lpl", "Ptprc", "Fgf1", "Aplp2", "Nrp2"), ]
MG_AD_WT_DE5$pathway <- "heparin binding"

MG_AD_WT_DE6 <- MG_AD_WT_DE[MG_AD_WT_DE$gene %in% c("Acaca", "Lpl", "Gnas", "Ldlr"), ]
MG_AD_WT_DE6$pathway <- "cellular response to fatty acid"

MG_AD_WT_DE7 <- MG_AD_WT_DE[MG_AD_WT_DE$gene %in% c("Nrp1", "Flt1", "Ccl4", "Apoe", "Ccl3", "Ptprc", "Gcnt2", "Ccl6", "Npnt", "Trem2"), ]
MG_AD_WT_DE7$pathway <- "positive regulation of ERK1 and ERK2 cas..."

MG_AD_WT_DE8 <- MG_AD_WT_DE[MG_AD_WT_DE$gene %in% c("Apoe", "Ldlr", "Igf1r", "Trem2"), ]
MG_AD_WT_DE8$pathway <- "amyloid-beta clearance"

MG_AD_WT_DE9 <- MG_AD_WT_DE[MG_AD_WT_DE$gene %in% c("Elmo1", "Dock4", "Prkce", "Lyn", "Afap1l1", "Ncf1", "Khdrbs3"), ]
MG_AD_WT_DE9$pathway <- "SH3 domain binding"

MG_AD_WT_DE10 <- MG_AD_WT_DE[MG_AD_WT_DE$gene %in% c("Prkca", "Bank1", "Fgd2", "Camk2d", "Pten", "Prkce", "Prkcd", "Rock2", "Syk", "Csf1r", "Arrb1", "Pxn", "Gab1", "Ncf1", "Pik3r5"), ]
MG_AD_WT_DE10$pathway <- "positive regulation of MAPK cascade"

MG_AD_WT_DE11 <- MG_AD_WT_DE[MG_AD_WT_DE$gene %in% c("Lyn", "Syk", "Arrb1", "Pxn", "Traf5", "Casp8", "Sipa1l1", "Hspa8"), ]
MG_AD_WT_DE11$pathway <- "ubiquitin protein ligase binding"

MG_AD_WT_DE12 <- MG_AD_WT_DE[MG_AD_WT_DE$gene %in% c("Mycbp2", "Asb2", "Rnf216", "Rnf130", "Znrf1"), ]
MG_AD_WT_DE12$pathway <- "ubiquitin protein ligase activity"


MG_AD_WT_DE_All <- rbind(MG_AD_WT_DE1, MG_AD_WT_DE2, MG_AD_WT_DE3,
                         MG_AD_WT_DE4, MG_AD_WT_DE5, MG_AD_WT_DE6,
                         MG_AD_WT_DE7, MG_AD_WT_DE8, MG_AD_WT_DE9,
                         MG_AD_WT_DE10, MG_AD_WT_DE11, MG_AD_WT_DE12)

MG_AD_WT_DE_All$pathway <- factor(MG_AD_WT_DE_All$pathway, levels = c("amyloid-beta binding", "amyloid-beta clearance", "cellular response to fatty acid",
                                                                          "heparin binding", "inflammatory response", "positive regulation of ERK1 and ERK2 cas...",
                                                                          "regulation of neuroinflammatory response","response to lipoprotein particle","positive regulation of MAPK cascade", 
                                                                          "SH3 domain binding", "ubiquitin protein ligase activity", "ubiquitin protein ligase binding"))

g2 <- Vln_plot(MG_AD_WT_DE_All) + scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"), name = "logFC", labels = c("minus", "plus"))  + theme(axis.text.y = element_blank())

blankPlot <- ggplot() + geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

g3 <- grid.arrange(g1, blankPlot, g2, blankPlot,  ncol = 2, widths = c(1.6, 0.9), heights = c(0.0, 3.5, 1.4), layout_matrix = cbind(c(1,1,1), c(2,3,4)))
g3
g3 <- arrangeGrob(g1, blankPlot, g2, blankPlot,  ncol = 2, widths = c(1.6, 0.9), heights = c(0, 3.5, 0.8), layout_matrix = cbind(c(1,1,1), c(2,3,4)))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_GO_vlnplot_Final.pdf",
       plot = g3,
       scale = 1, width = 10, height = 5.5, units = "in", device = cairo_pdf,
       dpi = 300)

MG_AD_WT_csv <- MG_AD_WT[MG_AD_WT$Term %in% c("amyloid-beta binding", "regulation of neuroinflammatory response", 
                                              "response to lipoprotein particle", "inflammatory response", 
                                              "heparin binding", "cellular response to fatty acid",
                                              "positive regulation of ERK1 and ERK2 cas...", "amyloid-beta clearance",
                                              "SH3 domain binding", "positive regulation of MAPK cascade", 
                                              "ubiquitin protein ligase binding", "ubiquitin protein ligase activity"), ]

write.csv(MG_AD_WT_csv, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_GO_vlnplot_Final.csv")

#######################
# Shirin asked GO
#######################

g1 <- ggplot(MG_AD_WT[MG_AD_WT$Term %in% c("amyloid-beta binding", "regulation of neuroinflammatory response", 
                                           "response to lipoprotein particle", "inflammatory response", 
                                           "positive regulation of cell migration"), ], 
             aes(x=group, y=Term,
                 colour=-log10(Fisher.elim),
                 size=gene_ratio)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="group", y="GO term", colour="-log10(p-value)", size="gene_ratio") +
  ggtitle("Microglia APPPS1 vs WT") + 
  theme_minimal() + 
  theme(axis.text.x =  element_text(size = 10, angle = 45, vjust = 0.5, family = "helvetica", color = "black"),
        axis.text.y = element_text(size = 10, family = "helvetica", color = "black"),
        axis.title.y = element_text(size = 10, family = "helvetica", color = "black"),
        axis.title.x =  element_blank(),
        plot.title = element_text(size = 10, family = "helvetica", color = "black", hjust = 0.5)) + 
  scale_color_gradient(low="#003366", high="#FF9900") + 
  guides(color = guide_colourbar(order = 1), 
         size = guide_legend(order = 2))


MG_AD_WT_DE1 <- MG_AD_WT_DE[MG_AD_WT_DE$gene %in% c("Apbb2", "Apoe", "Ldlrad3", "Ldlr", "Gsap", "Fcgr2b", "Scarb1", "Trem2"), ]
MG_AD_WT_DE1$pathway <- "amyloid-beta binding"

MG_AD_WT_DE2 <- MG_AD_WT_DE[MG_AD_WT_DE$gene %in% c("Igf1", "Cst7", "Ldlr", "Trem2"), ]
MG_AD_WT_DE2$pathway <- "regulation of neuroinflammatory response"

MG_AD_WT_DE3 <- MG_AD_WT_DE[MG_AD_WT_DE$gene %in% c("Cd9", "Apoe", "Lpl", "Ldlr"), ]
MG_AD_WT_DE3$pathway <- "response to lipoprotein particle"

MG_AD_WT_DE4 <- MG_AD_WT_DE[MG_AD_WT_DE$gene %in% c("Igf1","Hdac9","Hif1a","Ccl4","Apoe","Otulin","Axl","Csf1","Lpl","Ccl3","Abcd2","Acer3","Cst7","Ldlr","Abhd12","Cd180","Ly86","Ccl6","Nfe2l2","Fcgr2b","Cdk19","Plgrkt","C3ar1","Trem2"), ]
MG_AD_WT_DE4$pathway <- "inflammatory response"

genes1 <- c("Igf1","Nrp1","Myo5a","Hdac9","Hif1a","Flt1","Ccl4","Tgfbr2","Akt3","Csf1","Bcl2","Ptprc","Gcnt2","Thy1","Myo1f","Anxa3","Cd274","Atp8a1","Fgf1","Gab2","Igf1r","Nfe2l2","Il6st","P2rx4","C3ar1","Trem2","Prkca","Dock4","Numb","Fer","Camk2d","Pld1","Prkce","Wnk1","Lyn","Rreb1","Rock2","Dock8","Prex1","Csf1r","Cass4","Sparc","Gab1","Fermt3","Smo")

MG_AD_WT_DE5 <- MG_AD_WT_DE[MG_AD_WT_DE$gene %in% c(genes1), ]
MG_AD_WT_DE5$pathway <- "positive regulation of cell migration"

MG_AD_WT_DE_All <- rbind(MG_AD_WT_DE1, MG_AD_WT_DE2, MG_AD_WT_DE3,
                         MG_AD_WT_DE4, MG_AD_WT_DE5)


g2 <- Vln_plot(MG_AD_WT_DE_All) + ylim(c(-1.5,1.5)) + 
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"), name = "logFC", labels = c("minus", "plus"))  + 
  theme(axis.text.y = element_blank())

blankPlot <- ggplot() + geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

# Shirin
# g3 <- arrangeGrob(g1, blankPlot, g2, blankPlot,  ncol = 2, widths = c(1.6, 0.9), heights = c(0.2, 3.5, 0.6), layout_matrix = cbind(c(1,1,1), c(2,3,4)))

g3 <- arrangeGrob(g1, blankPlot, g2, blankPlot,  ncol = 2, widths = c(1.6, 0.9), heights = c(0, 3.5, 0.7), layout_matrix = cbind(c(1,1,1), c(2,3,4)))

g3 <- grid.arrange(g1, blankPlot, g2, blankPlot,  ncol = 2, widths = c(1.6, 0.9), heights = c(0.1, 3.5, 0.6), layout_matrix = cbind(c(1,1,1), c(2,3,4)))


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_GO_vlnplot_Shirin.png",
       plot = g3,
       scale = 1, width = 10, height = 4, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_GO_vlnplot_Final.pdf",
       plot = g3,
       scale = 1, width = 10, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig2_Microglia_GO_vlnplot_Shirin.svg",
       plot = g3,
       scale = 1, width = 10, height = 4, units = "in", device = "svg",
       dpi = 300)


