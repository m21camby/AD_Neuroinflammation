#! /bin/env RScript
# written by SJK at 25. Oct. 2020

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
    geom_violin(trim=FALSE) + coord_flip() + geom_jitter(data = df[df$FDR < 0.01, ], aes(color = logFC > 0),  position=position_jitter(0.1), alpha = 0.3) + theme(axis.title.y = element_blank()) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) + theme_minimal() + 
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 10, family = "helvetica", color = "black"),
          axis.text = element_text(size = 8, family = "helvetica", color = "black"),
          panel.grid.minor = element_blank())
}


###################
# MFOL Data loading
###################

MFOL_AD_ADp40KO_UP <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/GO_csv_files/edgeR_GO/cell_type_up_0.25/MFOL_AD_ADp40KO.csv_GO.csv", row.names = 1)
MFOL_AD_ADp40KO_UP$group <- "APPPS1 Up-GO"

MFOL_AD_ADp40KO_DOWN <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/GO_csv_files/edgeR_GO/cell_type_down_0.25/MFOL_AD_ADp40KO.csv_GO.csv", row.names = 1)

MFOL_AD_ADp40KO_DOWN$group <- "APPPS1.il12b.-/- Up-GO"

MFOL_AD_ADp40KO <- rbind(MFOL_AD_ADp40KO_UP, MFOL_AD_ADp40KO_DOWN)

MFOL_AD_ADp40KO[MFOL_AD_ADp40KO$Term %in% "cell projection", ]$genes

MFOL_AD_ADp40KO_DE <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/MFOL_AD_ADp40KO.csv", row.names = 1)

####################
# GO plot
####################

g1 <- ggplot(MFOL_AD_ADp40KO[MFOL_AD_ADp40KO$Term %in% c("cell morphogenesis", "programmed cell death", "cell death"), ], 
       aes(x=group, y=Term,
           colour=-log10(Fisher.elim),
           size=gene_ratio)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="group", y="GO term", colour="-log10(p-value)", size="gene_ratio") +
  ggtitle("MFOL APPPS1 vs APPPS1.il12b.-/-") + 
  theme_minimal() + 
  theme(axis.text.x =  element_text(size = 10, angle = 45, vjust = 0.5, family = "helvetica", color = "black"),
        axis.text.y = element_text(size = 10, family = "helvetica", color = "black"),
        axis.title.y = element_text(size = 10, family = "helvetica", color = "black"),
        axis.title.x =  element_blank(),
        plot.title = element_text(size = 10, family = "helvetica", color = "black")) + 
  scale_color_gradient(low="#003366", high="#FF9900") + 
  guides(color = guide_colourbar(order = 1), 
         size = guide_legend(order = 2))


##########################
# Vln plot
##########################

MFOL_AD_ADp40KO_DE1 <- MFOL_AD_ADp40KO_DE[MFOL_AD_ADp40KO_DE$gene %in% c("Shtn1", "Cdh20", "Gas7", "Ephb1", "Dcc", "Syne1"), ]
MFOL_AD_ADp40KO_DE1$pathway <- "cell morphogenesis"

MFOL_AD_ADp40KO_DE2 <- MFOL_AD_ADp40KO_DE[MFOL_AD_ADp40KO_DE$gene %in% c("Elmo1", "Oxr1", "Hif3a", "Gli2", "Eya4", "Cst3"), ]
MFOL_AD_ADp40KO_DE2$pathway <- "programmed cell death"

MFOL_AD_ADp40KO_DE3 <- MFOL_AD_ADp40KO_DE[MFOL_AD_ADp40KO_DE$gene %in% c("Elmo1", "Oxr1", "Hif3a", "Gli2", "Eya4", "Cst3"), ]
MFOL_AD_ADp40KO_DE3$pathway <- "cell death"

MFOL_AD_ADp40KO_DE_all <- rbind(MFOL_AD_ADp40KO_DE1, MFOL_AD_ADp40KO_DE2, MFOL_AD_ADp40KO_DE3)

MFOL_AD_ADp40KO_DE_all$pathway <- factor(MFOL_AD_ADp40KO_DE_all$pathway, levels = rev(c("cell morphogenesis", "programmed cell death", "cell death")))


g2 <- Vln_plot(MFOL_AD_ADp40KO_DE_all) + scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"), name = "logFC", labels = c("minus", "plus")) + theme(axis.text.y = element_blank())

blankPlot <- ggplot() + geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

g3 <- arrangeGrob(g1, blankPlot, g2, blankPlot,  ncol = 2, widths = c(1.5, 1), heights = c(0.2, 3.5, 1), layout_matrix = cbind(c(1,1,1), c(2,3,4)))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_GO_vlnplot_MFOL.png",
       plot = g3,
       scale = 1, width = 9, height = 5, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_GO_vlnplot_MFOL.pdf",
       plot = g3,
       scale = 1, width = 9, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_GO_vlnplot_MFOL.svg",
       plot = g3,
       scale = 1, width = 9, height = 5, units = "in", device = "svg",
       dpi = 300)

###################
# MOL Data loading
###################
MOL_AD_ADp40KO_UP <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/GO_csv_files/edgeR_GO/cell_type_up_0.25/MOL_AD_ADp40KO.csv_GO.csv", row.names = 1)
MOL_AD_ADp40KO_UP$group <- "APPPS1 Up-GO"

MOL_AD_ADp40KO_DOWN <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/GO_csv_files/edgeR_GO/cell_type_down_0.25/MOL_AD_ADp40KO.csv_GO.csv", row.names = 1)

MOL_AD_ADp40KO_DOWN$group <- "APPPS1.il12b.-/- Up-GO"

MOL_AD_ADp40KO <- rbind(MOL_AD_ADp40KO_UP, MOL_AD_ADp40KO_DOWN)

MOL_AD_ADp40KO[MOL_AD_ADp40KO$Term %in% "negative regulation of cellular process", ]$genes

MOL_AD_ADp40KO_DE <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/MOL_AD_ADp40KO.csv", row.names = 1)

####################
# GO plot
####################
g1 <- ggplot(MOL_AD_ADp40KO[MOL_AD_ADp40KO$Term %in% c("cell projection", "cell part", "cell morphogenesis"), ], 
             aes(x=group, y=Term,
                 colour=-log10(Fisher.elim),
                 size=gene_ratio)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="group", y="GO term", colour="-log10(p-value)", size="gene_ratio") +
  ggtitle("MOL APPPS1 vs APPPS1.il12b.-/-") + 
  theme_minimal() + 
  theme(axis.text.x =  element_text(size = 10, angle = 45, vjust = 0.5, family = "helvetica", color = "black"),
        axis.text.y = element_text(size = 10, family = "helvetica", color = "black"),
        axis.title.y = element_text(size = 10, family = "helvetica", color = "black"),
        axis.title.x =  element_blank(),
        plot.title = element_text(size = 10, family = "helvetica", color = "black")) + 
  scale_color_gradient(low="#003366", high="#FF9900") + 
  guides(color = guide_colourbar(order = 1), 
         size = guide_legend(order = 2))

##########################
# Vln plot
##########################

MOL_AD_ADp40KO_DE1 <- MOL_AD_ADp40KO_DE[MOL_AD_ADp40KO_DE$gene %in% c("Shtn1", "Dcc", "Samd4", "Gas7", "Cadm2"), ]
MOL_AD_ADp40KO_DE1$pathway <- "cell morphogenesis"

MOL_AD_ADp40KO_DE2 <- MOL_AD_ADp40KO_DE[MOL_AD_ADp40KO_DE$gene %in% c("Rnf220", "St18", "Shtn1", "Rcan2", "Ano4", "Pttg1", "Dcc", "Fth1", "Zdhhc14", "Samd4", "Piezo2", "Gas7", "Cadm2", "Cdh20"), ]
MOL_AD_ADp40KO_DE2$pathway <- "cell part"

MOL_AD_ADp40KO_DE3 <- MOL_AD_ADp40KO_DE[MOL_AD_ADp40KO_DE$gene %in% c("Shtn1", "Dcc", "Samd4", "Gas7", "Cadm2"), ]
MOL_AD_ADp40KO_DE3$pathway <- "cell projection"


MOL_AD_ADp40KO_DE_all <- rbind(MOL_AD_ADp40KO_DE1, MOL_AD_ADp40KO_DE2, MOL_AD_ADp40KO_DE3)

MOL_AD_ADp40KO_DE_all$pathway <- factor(MOL_AD_ADp40KO_DE_all$pathway, levels = rev(c("cell projection", "cell part", "cell morphogenesis")))


g2 <- Vln_plot(MOL_AD_ADp40KO_DE_all) + scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"), name = "logFC", labels = c("minus", "plus")) + theme(axis.text.y = element_blank())

blankPlot <- ggplot() + geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

g3 <- arrangeGrob(g1, blankPlot, g2, blankPlot,  ncol = 2, widths = c(1.2, 1), heights = c(0.2, 3.5, 1), layout_matrix = cbind(c(1,1,1), c(2,3,4)))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_GO_vlnplot_MOL.png",
       plot = g3,
       scale = 1, width = 10, height = 4.5, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_GO_vlnplot_MOL.pdf",
       plot = g3,
       scale = 1, width = 10, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_GO_vlnplot_MOL.svg",
       plot = g3,
       scale = 1, width = 10, height = 4.5, units = "in", device = "svg",
       dpi = 300)

##################
# For my Thesis
##################

g1 <- ggplot(MFOL_AD_ADp40KO[MFOL_AD_ADp40KO$Term %in% c("actin filament binding","regulation of ERK1 and ERK2 cascade","cell morphogenesis", "programmed cell death", "cell death", "cellular response to oxidative stress","dephosphorylation"), ], 
             aes(x=group, y=Term,
                 colour=-log10(Fisher.elim),
                 size=gene_ratio)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="group", y="GO term", colour="-log10(p-value)", size="gene_ratio") +
  ggtitle("MFOL APPPS1 vs APPPS1.il12b.-/-") + 
  theme_minimal() + 
  theme(axis.text.x =  element_text(size = 10, angle = 45, vjust = 0.5, family = "helvetica", color = "black"),
        axis.text.y = element_text(size = 10, family = "helvetica", color = "black"),
        axis.title.y = element_text(size = 10, family = "helvetica", color = "black"),
        axis.title.x =  element_blank(),
        plot.title = element_text(size = 10, family = "helvetica", color = "black")) + 
  scale_color_gradient(low="#003366", high="#FF9900") + 
  guides(color = guide_colourbar(order = 1), 
         size = guide_legend(order = 2))


##########################
# Vln plot
##########################

MFOL_AD_ADp40KO_DE1 <- MFOL_AD_ADp40KO_DE[MFOL_AD_ADp40KO_DE$gene %in% c("Shtn1", "Gas7", "Syne1"), ]
MFOL_AD_ADp40KO_DE1$pathway <- "actin filament binding"

MFOL_AD_ADp40KO_DE2 <- MFOL_AD_ADp40KO_DE[MFOL_AD_ADp40KO_DE$gene %in% c("Shtn1", "Cdh20", "Gas7", "Ephb1", "Dcc", "Syne1"), ]
MFOL_AD_ADp40KO_DE2$pathway <- "cell morphogenesis"

MFOL_AD_ADp40KO_DE3 <- MFOL_AD_ADp40KO_DE[MFOL_AD_ADp40KO_DE$gene %in% c("Ephb1", "Dcc"), ]
MFOL_AD_ADp40KO_DE3$pathway <- "regulation of ERK1 and ERK2 cascade"

MFOL_AD_ADp40KO_DE4 <- MFOL_AD_ADp40KO_DE[MFOL_AD_ADp40KO_DE$gene %in% c("Elmo1", "Oxr1", "Hif3a", "Gli2", "Eya4", "Cst3"), ]
MFOL_AD_ADp40KO_DE4$pathway <- "programmed cell death"

MFOL_AD_ADp40KO_DE5 <- MFOL_AD_ADp40KO_DE[MFOL_AD_ADp40KO_DE$gene %in% c("Elmo1", "Oxr1", "Hif3a", "Gli2", "Eya4", "Cst3"), ]
MFOL_AD_ADp40KO_DE5$pathway <- "cell death"

MFOL_AD_ADp40KO_DE6 <- MFOL_AD_ADp40KO_DE[MFOL_AD_ADp40KO_DE$gene %in% c("Slc8a1", "Ptprk", "Oxr1"), ]
MFOL_AD_ADp40KO_DE6$pathway <- "cellular response to oxidative stress"

MFOL_AD_ADp40KO_DE7 <- MFOL_AD_ADp40KO_DE[MFOL_AD_ADp40KO_DE$gene %in% c("Pcdh11x", "Ptprk", "Eya4"), ]
MFOL_AD_ADp40KO_DE7$pathway <- "dephosphorylation"


MFOL_AD_ADp40KO_DE_all <- rbind(MFOL_AD_ADp40KO_DE1, MFOL_AD_ADp40KO_DE2, MFOL_AD_ADp40KO_DE3,
                                MFOL_AD_ADp40KO_DE4, MFOL_AD_ADp40KO_DE5, MFOL_AD_ADp40KO_DE6,
                                MFOL_AD_ADp40KO_DE7)

MFOL_AD_ADp40KO_DE_all$pathway <- factor(MFOL_AD_ADp40KO_DE_all$pathway, levels = rev(c("regulation of ERK1 and ERK2 cascade","cell morphogenesis","actin filament binding","programmed cell death","dephosphorylation", "cellular response to oxidative stress", "cell death")))

g2 <- Vln_plot(MFOL_AD_ADp40KO_DE_all) + scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"), name = "logFC", labels = c("minus", "plus")) + theme(axis.text.y = element_blank())

blankPlot <- ggplot() + geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

g3 <- arrangeGrob(g1, blankPlot, g2, blankPlot,  ncol = 2, widths = c(1.5, 1), heights = c(0.2, 3.5, 0.9), layout_matrix = cbind(c(1,1,1), c(2,3,4)))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_GO_vlnplot_MFOL_Thesis.png",
       plot = g3,
       scale = 1, width = 10, height = 4.5, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_GO_vlnplot_MFOL_Thesis.pdf",
       plot = g3,
       scale = 1, width = 10, height = 4.5, units = "in", device = cairo_pdf,
       dpi = 300)

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_GO_vlnplot_session_info.txt")



