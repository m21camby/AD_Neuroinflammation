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
library(cowplot)

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


MFOL_UP <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/GO_csv_files/edgeR_GO/cell_type_up_0.25/MFOL_AD_ADp40KO.csv_GO.csv", row.names = 1, stringsAsFactors = F)
MFOL_UP <- MFOL_UP[MFOL_UP$Fisher.elim < 0.05, ]
MFOL_UP <- MFOL_UP[MFOL_UP$Significant > 1, ]
MFOL_UP$cell_type <- "MFOL"
MFOL_UP$exp <- "APPPS1 UP"
# negative regulation of cell death

MFOL_DOWN <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/GO_csv_files/edgeR_GO/cell_type_down_0.25/MFOL_AD_ADp40KO.csv_GO.csv", row.names = 1, stringsAsFactors = F)
MFOL_DOWN <- MFOL_DOWN[MFOL_DOWN$Fisher.elim < 0.05, ]
MFOL_DOWN <- MFOL_DOWN[MFOL_DOWN$Significant > 1, ]
MFOL_DOWN$cell_type <- "MFOL"
MFOL_DOWN$exp <- "APPPS1.Il12b-/- UP"
# synaptic signaling

MOL_UP <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/GO_csv_files/edgeR_GO/cell_type_up_0.25/MOL_AD_ADp40KO.csv_GO.csv", row.names = 1, stringsAsFactors = F)
MOL_UP <- MOL_UP[MOL_UP$Fisher.elim < 0.05, ]
MOL_UP <- MOL_UP[MOL_UP$Significant > 1, ]
MOL_UP$cell_type <- "MOL"
MOL_UP$exp <- "APPPS1 UP"
# negative regulation of cell differentiat...
# negative regulation of programmed cell d...
# negative regulation of cellular process

MOL_DOWN <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/GO_csv_files/edgeR_GO/cell_type_down_0.25/MOL_AD_ADp40KO.csv_GO.csv", row.names = 1, stringsAsFactors = F)
MOL_DOWN <- MOL_DOWN[MOL_DOWN$Fisher.elim < 0.05, ]
MOL_DOWN <- MOL_DOWN[MOL_DOWN$Significant > 1, ]
MOL_DOWN$cell_type <- "MOL"
MOL_DOWN$exp <- "APPPS1.Il12b-/- UP"
# negative regulation of cell proliferatio...


OPC_UP <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/GO_csv_files/edgeR_GO/cell_type_up_0.25/OPC_AD_ADp40KO.csv_GO.csv", row.names = 1, stringsAsFactors = F)
OPC_UP <- OPC_UP[OPC_UP$Fisher.elim < 0.05, ]
OPC_UP <- OPC_UP[OPC_UP$Significant > 1, ]
OPC_UP$cell_type <- "OPC"
OPC_UP$exp <- "APPPS1 UP"
# regulation of programmed cell death
# negative regulation of neuron differenti...
# developmental maturation
# oligodendrocyte differentiation
# myelination

OPC_DOWN <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/GO_csv_files/edgeR_GO/cell_type_down_0.25/OPC_AD_ADp40KO.csv_GO.csv", row.names = 1, stringsAsFactors = F)
OPC_DOWN <- OPC_DOWN[OPC_DOWN$Fisher.elim < 0.05, ]
OPC_DOWN <- OPC_DOWN[OPC_DOWN$Significant > 1, ]
OPC_DOWN$cell_type <- "OPC"
OPC_DOWN$exp <- "APPPS1.Il12b-/- UP"
# negative regulation of cell proliferatio...
# cell growth
# DNA-binding transcription activator acti...


MFOL_Terms <- c("programmed cell death", "negative regulation of cell death", "cellular response to oxidative stress", "dephosphorylation",
                "actin filament binding", "cell morphogenesis", "regulation of ERK1 and ERK2 cascade", "synaptic signaling")


MOL_Terms <- c("negative regulation of cellular process", "negative regulation of programmed cell d...",
               "negative regulation of cell proliferatio...")

OPC_Terms <- c("regulation of programmed cell death", "negative regulation of neuron differenti...", "developmental maturation", "myelination",
               "cell growth")

All_Terms <- c(MFOL_Terms, MOL_Terms, OPC_Terms) %>% unique

All_GO <- rbind(MOL_UP, MOL_DOWN, MFOL_UP, MFOL_DOWN, OPC_UP, OPC_DOWN)

All_GO <- All_GO[All_GO$Term %in% All_Terms, ]

# element rename
All_GO[All_GO == "negative regulation of neuron differenti..."] <- "negative regulation of neuron differentiation"
All_GO[All_GO == "negative regulation of programmed cell d..."] <- "negative regulation of programmed cell death"
All_GO[All_GO == "negative regulation of cell proliferatio..."] <- "negative regulation of cell proliferation"



# order change
All_GO$cell_type <- factor(All_GO$cell_type, levels = c("OPC", "MFOL", "MOL"))

All_GO$Term <- factor(All_GO$Term, levels = c("developmental maturation", "myelination","negative regulation of neuron differentiation","regulation of programmed cell death",
                                              "cell growth", "negative regulation of cell proliferation",
                                              "cellular response to oxidative stress", "dephosphorylation", "negative regulation of cellular process", "programmed cell death", "negative regulation of cell death",
                                              "synaptic signaling", "regulation of ERK1 and ERK2 cascade", "actin filament binding", "cell morphogenesis",
                                              "negative regulation of programmed cell death"))

####################
# GO plot
####################
g1 <- ggplot(All_GO, 
       aes(x=exp, y=Term,
           colour=Fisher.elim,
           size=gene_ratio)) +
  geom_point() +
  expand_limits(x=0) + facet_wrap(~cell_type) + 
  labs(x="group", y="GO term", colour="p-value", size="gene_ratio") +
  theme_cowplot() + 
  theme(axis.title.y = element_blank(),
        axis.title.x =  element_blank(),
        axis.text.x = element_text(size = 12, angle = 90, vjust=0.5, family = "helvetica", color = "black"),
        axis.text.y = element_text(size = 10, family = "helvetica", color = "black"),
        strip.text.y = element_text(size = 10, angle = 0, hjust = 0, vjust= 0, family = "helvetica", color = "black"),
        strip.background = element_rect(color = "white", fill="white"),
        panel.spacing = unit(c(-1), "lines")) + 
  scale_color_gradient(low="#003366", high="#FF9900", trans = 'reverse') + 
  guides(color = guide_colourbar(order = 1), 
         size = guide_legend(order = 2)) + guides(color = guide_colorbar(reverse=TRUE))


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_GO_dotplot_New.pdf",
       plot = g1,
       scale = 1, width = 8, height = 8, units = "in", device = cairo_pdf,
       dpi = 300)

write.csv(All_GO, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_GO_dotplot_New.csv")


##########################
# Vln plot
##########################

# 1. OPC
# Data loading
OPC_AD_ADp40KO_DE <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/OPC_AD_ADp40KO.csv", row.names = 1)

OPC_AD_ADp40KO_DE1 <- OPC_AD_ADp40KO_DE[OPC_AD_ADp40KO_DE$gene %in% c("Sema5a", "Cntn4", "Olig2"), ]
OPC_AD_ADp40KO_DE1$pathway <- "negative regulation of neuron differentiation"

OPC_AD_ADp40KO_DE2 <- OPC_AD_ADp40KO_DE[OPC_AD_ADp40KO_DE$gene %in% c("Plcb1", "Kcnma1", "Sez6l"), ]
OPC_AD_ADp40KO_DE2$pathway <- "developmental maturation"

OPC_AD_ADp40KO_DE3 <- OPC_AD_ADp40KO_DE[OPC_AD_ADp40KO_DE$gene %in% c("Clu", "Olig2"), ]
OPC_AD_ADp40KO_DE3$pathway <- "myelination"

OPC_AD_ADp40KO_DE4 <- OPC_AD_ADp40KO_DE[OPC_AD_ADp40KO_DE$gene %in% c("Il12b", "Kcnma1", "Sema5a", "Clu", "Sycp2"), ]
OPC_AD_ADp40KO_DE4$pathway <- "regulation of programmed cell death"

OPC_AD_ADp40KO_DE5 <- OPC_AD_ADp40KO_DE[OPC_AD_ADp40KO_DE$gene %in% c("Tmem108", "Pttg1"), ]
OPC_AD_ADp40KO_DE5$pathway <- "cell growth"

OPC_AD_ADp40KO_DE6 <- OPC_AD_ADp40KO_DE[OPC_AD_ADp40KO_DE$gene %in% c("Pttg1", "St18"), ]
OPC_AD_ADp40KO_DE6$pathway <- "negative regulation of cell proliferation"

# all
OPC_AD_ADp40KO_DE_all <- rbind(OPC_AD_ADp40KO_DE1, OPC_AD_ADp40KO_DE2, OPC_AD_ADp40KO_DE3,
                               OPC_AD_ADp40KO_DE4, OPC_AD_ADp40KO_DE5, OPC_AD_ADp40KO_DE6)

OPC_AD_ADp40KO_DE_all$cell_type <- "OPC"

# MFOL
# Data loading
MFOL_AD_ADp40KO_DE <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/MFOL_AD_ADp40KO.csv", row.names = 1)


MFOL_AD_ADp40KO_DE1 <- MFOL_AD_ADp40KO_DE[MFOL_AD_ADp40KO_DE$gene %in% c("Shtn1", "Gas7", "Syne1"), ]
MFOL_AD_ADp40KO_DE1$pathway <- "actin filament binding"

MFOL_AD_ADp40KO_DE2 <- MFOL_AD_ADp40KO_DE[MFOL_AD_ADp40KO_DE$gene %in% c("Shtn1", "Cdh20", "Gas7", "Ephb1", "Dcc", "Syne1"), ]
MFOL_AD_ADp40KO_DE2$pathway <- "cell morphogenesis"

MFOL_AD_ADp40KO_DE3 <- MFOL_AD_ADp40KO_DE[MFOL_AD_ADp40KO_DE$gene %in% c("Ephb1", "Dcc"), ]
MFOL_AD_ADp40KO_DE3$pathway <- "regulation of ERK1 and ERK2 cascade"

MFOL_AD_ADp40KO_DE4 <- MFOL_AD_ADp40KO_DE[MFOL_AD_ADp40KO_DE$gene %in% c("Elmo1", "Oxr1", "Hif3a", "Gli2", "Eya4", "Cst3"), ]
MFOL_AD_ADp40KO_DE4$pathway <- "programmed cell death"

MFOL_AD_ADp40KO_DE5 <- MFOL_AD_ADp40KO_DE[MFOL_AD_ADp40KO_DE$gene %in% c(	"Oxr1", "Gli2", "Eya4", "Cst3"), ]
MFOL_AD_ADp40KO_DE5$pathway <- "negative regulation of cell death"

MFOL_AD_ADp40KO_DE6 <- MFOL_AD_ADp40KO_DE[MFOL_AD_ADp40KO_DE$gene %in% c("Slc8a1", "Ptprk", "Oxr1"), ]
MFOL_AD_ADp40KO_DE6$pathway <- "cellular response to oxidative stress"

MFOL_AD_ADp40KO_DE7 <- MFOL_AD_ADp40KO_DE[MFOL_AD_ADp40KO_DE$gene %in% c("Pcdh11x", "Ptprk", "Eya4"), ]
MFOL_AD_ADp40KO_DE7$pathway <- "dephosphorylation"

MFOL_AD_ADp40KO_DE8 <- MFOL_AD_ADp40KO_DE[MFOL_AD_ADp40KO_DE$gene %in% c("Otud7a", "Spock1", "Slc8a1", "Ptgds", "Pcdh11x", "Ptprk", "Oxr1", "Gli2", "Eya4", "Cst3"), ]
MFOL_AD_ADp40KO_DE8$pathway <- "negative regulation of cellular process"

MFOL_AD_ADp40KO_DE9 <- MFOL_AD_ADp40KO_DE[MFOL_AD_ADp40KO_DE$gene %in% c("Anks1b", "Ephb1", "Dcc", "Syne1"), ]
MFOL_AD_ADp40KO_DE9$pathway <- "synaptic signaling"

MFOL_AD_ADp40KO_DE_all <- rbind(MFOL_AD_ADp40KO_DE1, MFOL_AD_ADp40KO_DE2, MFOL_AD_ADp40KO_DE3,
                                MFOL_AD_ADp40KO_DE4, MFOL_AD_ADp40KO_DE5, MFOL_AD_ADp40KO_DE6,
                                MFOL_AD_ADp40KO_DE7, MFOL_AD_ADp40KO_DE8, MFOL_AD_ADp40KO_DE9)

MFOL_AD_ADp40KO_DE_all$cell_type <- "MFOL"

# MOL
# Data loading
MOL_AD_ADp40KO_DE <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/MOL_AD_ADp40KO.csv", row.names = 1)

MOL_AD_ADp40KO_DE1 <- MOL_AD_ADp40KO_DE[MOL_AD_ADp40KO_DE$gene %in% c("Sox10", "Gli2", "Sgk3"), ]
MOL_AD_ADp40KO_DE1$pathway <- "negative regulation of programmed cell death"

MOL_AD_ADp40KO_DE1$cell_type <- "MOL"


All_DE <- rbind(OPC_AD_ADp40KO_DE_all, MFOL_AD_ADp40KO_DE_all, MOL_AD_ADp40KO_DE1)

All_DE$cell_type <- factor(All_DE$cell_type, levels = c("OPC", "MFOL", "MOL"))

All_DE$pathway <- factor(All_DE$pathway, levels = c("developmental maturation", "myelination","negative regulation of neuron differentiation","regulation of programmed cell death",
                                              "cell growth", "negative regulation of cell proliferation",
                                              "cellular response to oxidative stress", "dephosphorylation", "negative regulation of cellular process", "programmed cell death", "negative regulation of cell death",
                                              "synaptic signaling", "regulation of ERK1 and ERK2 cascade", "actin filament binding", "cell morphogenesis",
                                              "negative regulation of programmed cell death"))


All_DE$cell_type <- factor(All_DE$cell_type, levels = rev(c("OPC", "MFOL", "MOL")))

g2 <- ggplot(All_DE, aes(x=pathway, y=logFC)) + 
  geom_violin(trim=FALSE) + 
  coord_flip() + 
  geom_jitter(data = All_DE, aes(color = logFC > 0),  position=position_jitter(0.1), alpha = 0.3) + 
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"), name = "logFC", labels = c("minus", "plus")) + 
  theme_minimal() + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10, family = "helvetica", color = "black"),
        axis.text = element_text(size = 8, family = "helvetica", color = "black"),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(c(-0.2), "lines")) + facet_grid(cell_type ~., scales = "free", space = "free")
g2


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_GO_vlnplot_New.pdf",
       plot = g2,
       scale = 1, width = 6, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)


#################################
# Combined
#################################

g2 <- ggplot(All_DE, aes(x=pathway, y=logFC)) + 
  geom_violin(trim=FALSE) + 
  coord_flip() + 
  geom_jitter(data = All_DE, aes(color = logFC > 0),  position=position_jitter(0.1), alpha = 0.3) + 
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"), name = "logFC", labels = c("minus", "plus")) + 
  theme_minimal() + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(), 
        axis.title.x = element_text(size = 10, family = "helvetica", color = "black"),
        axis.text.x = element_text(size = 8, family = "helvetica", color = "black"),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(c(-0.2), "lines")) + facet_grid(cell_type ~., scales = "free", space = "free")
g2

blankPlot <- ggplot() + geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

g3 <- grid.arrange(g1, blankPlot, g2, blankPlot,  ncol = 2, widths = c(1.6, 0.9), heights = c(0.3, 3.4, 1.4), layout_matrix = cbind(c(1,1,1), c(2,3,4)))
g3
g3 <- arrangeGrob(g1, blankPlot, g2, blankPlot,  ncol = 2, widths = c(1.6, 0.9), heights = c(0.2, 4, 0.8), layout_matrix = cbind(c(1,1,1), c(2,3,4)))


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig3_Oligo_OPC_GO_dotplot_vlnplot_New.pdf",
       plot = g3,
       scale = 1, width = 12, height = 8, units = "in", device = cairo_pdf,
       dpi = 300)









