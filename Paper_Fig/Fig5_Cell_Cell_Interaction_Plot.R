#! /bin/env RScript
# written by SJK at 21. Oct. 2020

##################################
# fig reference
# https://www.nature.com/articles/s41467-020-17834-w
# https://www.nature.com/articles/s41587-020-0602-4
##################################

##################################
# Text: 
# Dot plot of the predicted interactions between neuronal cell and glial cells
# in WT(W), APPPS1(A) and APPPS1.il12b.-/-(K).
# P values are indicated by the circle sizes, as shown in the scale on the right (permutation test). 
# The means of the average expression level of interacting molecule 1 in cluster 1 and interacting molecule 2 in cluster 2 are indicated by the colour. 
# Assays were carried out at the mRNA level but were extrapolated to protein interactions.

# we used the accumulated ligand/receptor interaction database CellPhoneDB (www.cellphonedb.org) to identify alterations of the molecular interactions

# Dot plots of ligand-receptors in neuronal and glial cells Significantly altered expression 
# To map the complex interactions of neuronal and glial cells, 
# we inferred all possible intercellular communications by the expression of ligand–receptor pairs in both cell populations with CellPhoneDB

library(Seurat)
library(ggplot2)
library(dplyr)
library(DT)
library(gridExtra)
library(tidyverse)
library(irr)
library(viridis)
library(grid)

source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")


###################################
# loading cellphoneDB results file
###################################

sigmeanvalue_Ctrl.df <- read.table("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/out/cellphoneDB_Ctrl_cell_type/significant_means.txt", sep = "\t", header = TRUE, check.names = FALSE)
sigmeanvalue_AD.df <- read.table("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/out/cellphoneDB_AD_cell_type/significant_means.txt", sep = "\t", header = TRUE, check.names = FALSE)
sigmeanvalue_ADp40KO.df <- read.table("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/out/cellphoneDB_ADp40KO_cell_type/significant_means.txt", sep = "\t", header = TRUE, check.names = FALSE)

pvalue_Ctrl.df <- read.table("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/out/cellphoneDB_Ctrl_cell_type/pvalues.txt", sep = "\t", header = TRUE, check.names = FALSE)
pvalue_AD.df <- read.table("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/out/cellphoneDB_AD_cell_type/pvalues.txt", sep = "\t", header = TRUE, check.names = FALSE)
pvalue_ADp40KO.df <- read.table("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/out/cellphoneDB_ADp40KO_cell_type/pvalues.txt", sep = "\t", header = TRUE, check.names = FALSE)

###################################
# function to extract interesting receptor-ligand and cell types
###################################

cellphone_Cell_pathway <- function(pathway = "BDNF_SORT1", cell_type = c("MFOL", "MOL")){
  Results_Ctrl <- sigmeanvalue_Ctrl.df[sigmeanvalue_Ctrl.df$interacting_pair	== pathway, c(13:373)] %>% t %>% as.data.frame
  colnames(Results_Ctrl) <- "WT_means"
  Results_Ctrl$interacting_cells <- rownames(Results_Ctrl)
  
  Results_AD <- sigmeanvalue_AD.df[sigmeanvalue_AD.df$interacting_pair	== pathway, c(13:373)] %>% t %>% as.data.frame
  colnames(Results_AD) <- "AD_means"
  Results_AD$interacting_cells <- rownames(Results_AD)
  
  Results_ADp40KO <- sigmeanvalue_ADp40KO.df[sigmeanvalue_ADp40KO.df$interacting_pair	== pathway, c(13:373)] %>% t %>% as.data.frame
  colnames(Results_ADp40KO) <- "ADp40KO_means"
  Results_ADp40KO$interacting_cells <- rownames(Results_ADp40KO)
  
  Results <- inner_join(Results_Ctrl, Results_AD, by = "interacting_cells")
  
  Results <- inner_join(Results, Results_ADp40KO, by = "interacting_cells")
  
  # subset MOL, MFOL
  Results <- Results[grepl(paste(cell_type, collapse="|"), Results$interacting_cells), ]
  
  # remove non informtive cell types 
  Results <- Results[!grepl(paste(c("Cajal", "Choroid", "Fibroblast", "Pericyte", "VLMC", "Vascular", "Macrophage"), collapse="|"), Results$interacting_cells), ]
  
  rownames(Results) <- Results$interacting_cells
  Results$interacting_cells <- NULL
  
  Results[is.na(Results)] <- 0
  
  Results$interacting_cells <- rownames(Results)
  
  Results <- Results %>% gather(exp, means, "WT_means":"ADp40KO_means")
  
  P_Results_Ctrl <- pvalue_Ctrl.df[pvalue_Ctrl.df$interacting_pair	== pathway, c(12:372)] %>% t %>% as.data.frame
  colnames(P_Results_Ctrl) <- "WT_pvalue"
  P_Results_Ctrl$interacting_cells <- rownames(P_Results_Ctrl)
  
  P_Results_AD <- pvalue_AD.df[pvalue_AD.df$interacting_pair	== pathway, c(12:372)] %>% t %>% as.data.frame
  colnames(P_Results_AD) <- "AD_pvalue"
  P_Results_AD$interacting_cells <- rownames(P_Results_AD)
  
  P_Results_ADp40KO <- pvalue_ADp40KO.df[pvalue_ADp40KO.df$interacting_pair	== pathway, c(12:372)] %>% t %>% as.data.frame
  colnames(P_Results_ADp40KO) <- "ADp40KO_pvalue"
  P_Results_ADp40KO$interacting_cells <- rownames(P_Results_ADp40KO)
  
  P_Results <- inner_join(P_Results_Ctrl, P_Results_AD, by = "interacting_cells")
  
  P_Results <- inner_join(P_Results, P_Results_ADp40KO, by = "interacting_cells")
  
  rownames(P_Results) <- P_Results$interacting_cells
  P_Results$interacting_cells <- NULL
  
  P_Results[is.na(P_Results)] <- 0
  
  P_Results <- P_Results[rownames(P_Results) %in% Results$interacting_cells, ]
  
  P_Results$interacting_cells <- rownames(P_Results)
  
  P_Results <- P_Results %>% gather(exp, pvalue, "WT_pvalue":"ADp40KO_pvalue")
  
  Results_final <- cbind(Results, P_Results)
  Results_final$exp <- str_sub(Results_final$exp, end = -7)
  
  Results_final <- Results_final[, c(1,2,3,6)]
  Results_final$logpvalue <- -log10(Results_final$pvalue)
  
  Results_final[sapply(Results_final, is.infinite)] <- 3
  
  return(Results_final)
  
}

###################################
# function to modify for plot of extract interesting receptor-ligand and cell types
###################################

cellphone_Cell_pathway_modification_plot <- function(Pathway = "CADM1_CADM1",
                                                     mouse_Pathway = "Cadm1_Cadm1",
                                                     interesting_cell_type = "Microglia"){
  
  cell_types <- c("Dentate_Gyrus", "CA1", "CA2_3", "subiculum", "Inhibitory_Neurons", "MOL", "MFOL", "NFOL", "OPC", "Microglia", "Astrocytes")
  
  
  PATHWAY_final <- cellphone_Cell_pathway(pathway = Pathway, cell_type = cell_types)
  PATHWAY_final$ligand_receptor <- mouse_Pathway
  PATHWAY_final$sender <- sapply(strsplit(PATHWAY_final$interacting_cells,"\\|"), `[`, 1)
  PATHWAY_final$receiver <- sapply(strsplit(PATHWAY_final$interacting_cells,"\\|"), `[`, 2)
  PATHWAY_final <- PATHWAY_final[!grepl("Unidenti", PATHWAY_final$interacting_cells), ]
  
  PATHWAY_final <- PATHWAY_final[grepl(interesting_cell_type, PATHWAY_final$interacting_cells), ]
  PATHWAY_final_s <- PATHWAY_final[PATHWAY_final$sender %in% interesting_cell_type, ]
  PATHWAY_final_r <- PATHWAY_final[PATHWAY_final$receiver %in% interesting_cell_type, ]
  
  # Remove Both Microglia
  PATHWAY_final_s <- PATHWAY_final_s[PATHWAY_final_s$interacting_cells != paste0(interesting_cell_type, "|",interesting_cell_type), ]
  PATHWAY_final_r <- PATHWAY_final_r[PATHWAY_final_r$interacting_cells != paste0(interesting_cell_type, "|",interesting_cell_type), ]
  
  
  # Name change
  PATHWAY_final_s[PATHWAY_final_s == "CA2_3"] <- "CA2/3"
  PATHWAY_final_s[PATHWAY_final_s == "Dentate_Gyrus"] <- "Dentate Gyrus"
  PATHWAY_final_s[PATHWAY_final_s == "Inhibitory_Neurons"] <- "Inhibitory Neurons"
  
  PATHWAY_final_r[PATHWAY_final_r == "CA2_3"] <- "CA2/3"
  PATHWAY_final_r[PATHWAY_final_r == "Dentate_Gyrus"] <- "Dentate Gyrus"
  PATHWAY_final_r[PATHWAY_final_r == "Inhibitory_Neurons"] <- "Inhibitory Neurons"
  
  # Remove Microglia for labeling and left only one
  if(interesting_cell_type == "Inhibitory_Neurons"){
    PATHWAY_final_r$receiver <- ifelse(PATHWAY_final_r$sender %in% "MOL", "Inhibitory Neurons", " ")
  }else if(interesting_cell_type == "CA2_3"){
    PATHWAY_final_r$receiver <- ifelse(PATHWAY_final_r$sender %in% "Inhibitory Neurons", "CA2/3", " ")
  } else {
    PATHWAY_final_r$receiver <- ifelse(PATHWAY_final_r$sender %in% "Inhibitory Neurons", interesting_cell_type, " ")
    }
  
  
  
  # modified cell_types
  cell_types2 <- c("Dentate Gyrus", "CA1", "CA2/3", "subiculum", "Inhibitory Neurons", "MOL", "MFOL", "NFOL", "OPC", "Microglia", "Astrocytes")
  
  
  rest_cell_types <- cell_types2[!cell_types2 %in% interesting_cell_type]
  
  # give levels
  PATHWAY_final_s$receiver <- factor(PATHWAY_final_s$receiver, 
                                     levels = rev(rest_cell_types))
  
  PATHWAY_final_r$sender <- factor(PATHWAY_final_r$sender, 
                                   levels = rest_cell_types)
  
  
  PATHWAY_final_s$exp <- factor(PATHWAY_final_s$exp, levels = c("WT", "AD", "ADp40KO"))
  PATHWAY_final_r$exp <- factor(PATHWAY_final_r$exp, levels = c("WT", "AD", "ADp40KO"))
  
  return(list(PATHWAY_final_s, PATHWAY_final_r))
  
}

###################################
# function to extract upper plot legend
###################################

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

###################################
# function to uppper plot
###################################

cellphonedb_plot1 <- function(exp.df, range = c(0,15)){
  
  # Input as receiver cell type
  
  ggplot(exp.df, aes(x=exp,
                     y=receiver,
                     colour=means,
                     size=logpvalue)) +
    geom_point() +
    expand_limits(x=0) + labs(x="cluster", y="Interacting cell types", colour="mean exp", size="-log10(p value)") +
    facet_grid(sender ~ ligand_receptor, scales="free", space = "free", switch = "y") + 
    theme_minimal() + 
    theme(axis.title.x = element_blank(), 
          axis.title.y =  element_blank(),
          axis.text.y = element_text(size = 11, family = "Helvetica", color = "black"),
          axis.text.x = element_text(size = 11, family = "Helvetica", color = "black"),
          strip.background = element_blank(),
          strip.text.x = element_text(angle=45, size = 12, family = "Helvetica", color = "black"),
          strip.text.y.left = element_text(angle=0, size = 11, family = "Helvetica", color = "black"),
          strip.placement = "outside",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          #panel.spacing = unit(c(0, 0, 0, 0), "lines"),
          panel.spacing.x = unit(-2, "cm"),
          plot.margin = margin(1, 1, 1, 1, "cm")
    ) + scale_color_viridis(limits = range) + 
    scale_x_discrete(position = "top", labels=c("WT" = "W", "AD" = "A",
                                                "ADp40KO" = "K"), 
                     expand = expansion(mult = c(0.2, 1))) 
  
  
}

###################################
# function to lower plot
###################################

cellphonedb_plot2 <- function(exp.df, range = c(0,15)){
  
  # Input as sender cell type
  
  ggplot(exp.df, aes(x=exp,
                     y=receiver,
                     colour=means,
                     size=logpvalue)) +
    geom_point() +
    expand_limits(x=0) + labs(x="cluster", y="Interacting cell types", colour="mean exp", size="-log10(p value)") +
    facet_grid(sender ~ ligand_receptor, scales="free", space = "free", switch = "y") + 
    theme_minimal() + 
    theme(axis.title.x = element_blank(), 
          axis.title.y =  element_blank(),
          axis.text.y = element_text(size = 11, family = "Helvetica", color = "black"),
          axis.text.x = element_text(size = 11, family = "Helvetica", color = "black"),
          strip.background = element_blank(),
          strip.text.x = element_text(angle=45, family = "Helvetica", color = "black"),
          strip.text.y.left = element_text(angle=0, size = 11, family = "Helvetica", color = "black"),
          strip.placement = "outside",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          #panel.spacing = unit(c(0,0,0, -2), "lines"),
          panel.spacing.x = unit(-2, "cm"),
          plot.margin = margin(1, 1, 1, 1, "cm")
    ) + scale_color_viridis(limits = range) + 
    scale_x_discrete(position = "top", labels=c("WT" = "W", "AD" = "A",
                                                "ADp40KO" = "K"), 
                     expand = expansion(mult = c(0.2, 1))) 

}


###################
# 1. Microglia
###################
# s is sender and r is receiver

MG_CADM1 <- cellphone_Cell_pathway_modification_plot(Pathway = "CADM1_CADM1",
                                                 mouse_Pathway = "Cadm1_Cadm1",
                                                 interesting_cell_type = "Microglia")
MG_CADM1_s <- MG_CADM1[1] %>% as.data.frame
MG_CADM1_r <- MG_CADM1[2] %>% as.data.frame

MG_PTN <- cellphone_Cell_pathway_modification_plot(Pathway = "PTN_PTPRZ1",
                                                     mouse_Pathway = "Ptn_Ptprz1",
                                                     interesting_cell_type = "Microglia")
MG_PTN_s <- MG_PTN[1] %>% as.data.frame
MG_PTN_r <- MG_PTN[2] %>% as.data.frame

MG_r <- rbind(MG_CADM1_r, MG_PTN_r)
MG_s <- rbind(MG_CADM1_s, MG_PTN_s)

################
# Figure
################

c1 <- cellphonedb_plot1(MG_r, range = c(0,30))
c2 <- cellphonedb_plot2(MG_s, range = c(0,15))

# extract legend
legend <- get_legend(c1)

c1 <- c1 + theme(legend.position="none") + theme(plot.margin = margin(0, 0, 0, 0, "cm"),
                                                 legend.margin=margin(0,0,0,0),
                                                 legend.box.margin=margin(-10,-10,-10,-10)) 

c2 <- c2 + theme(legend.position="none") + theme(plot.margin = margin(0, 0, 0, 0, "cm"), 
                                                 strip.text.x = element_blank(),
                                                 axis.text.x = element_blank())
# create blank plot
blankPlot <- ggplot() + geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

g1 <- arrangeGrob(c1, c2, blankPlot, legend, blankPlot, blankPlot, ncol = 2, 
             widths = c(2.7, .5), heights = c(3.5, 3.5, 2.5, 2.5), 
             layout_matrix = cbind(c(1,1,2,2), c(3,4,5,6)))


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_MG.png",
       plot = g1,
       scale = 1, width = 6.3, height = 7, units = "in", device = "png",
       dpi = 300)


###################
# 2. Astrocytes
###################

AS_FGF1 <- cellphone_Cell_pathway_modification_plot(Pathway = "FGF1_FGFR2",
                                                     mouse_Pathway = "Fgf1_Fgfr2",
                                                     interesting_cell_type = "Astrocytes")
AS_FGF1_s <- AS_FGF1[1] %>% as.data.frame
AS_FGF1_r <- AS_FGF1[2] %>% as.data.frame

AS_NRG2 <- cellphone_Cell_pathway_modification_plot(Pathway = "NRG2_ERBB4",
                                                    mouse_Pathway = "Nrg2_Erbb4",
                                                    interesting_cell_type = "Astrocytes")
AS_NRG2_s <- AS_NRG2[1] %>% as.data.frame
AS_NRG2_r <- AS_NRG2[2] %>% as.data.frame


AS_r <- rbind(AS_FGF1_r, AS_NRG2_r)
AS_s <- rbind(AS_FGF1_s, AS_NRG2_s)

# Figure
c1 <- cellphonedb_plot1(AS_r, range = c(0,40))
c2 <- cellphonedb_plot2(AS_s, range = c(0,40))

c1 <- c1 + guides(color = guide_colourbar(order = 1),
                  size = guide_legend(order = 2))

legend <- get_legend(c1)

c1 <- c1 + theme(legend.position="none") + 
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10)) 


c2 <- c2 + theme(legend.position="none") + 
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        strip.text.x = element_blank(),
        axis.text.x = element_blank()) 

blankPlot <- ggplot() + geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

g1 <- arrangeGrob(c1, c2, blankPlot, legend, blankPlot, blankPlot, ncol = 2, 
                  widths = c(2.7, .5), heights = c(3.5, 3.5, 2.5, 2.5), 
                  layout_matrix = cbind(c(1,1,2,2), c(3,4,5,6)))


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_AS.png",
       plot = g1,
       scale = 1, width = 6.3, height = 7, units = "in", device = "png",
       dpi = 300)


###################
# 3. MOL
###################
MOL_ERBB3 <- cellphone_Cell_pathway_modification_plot(Pathway = "ERBB3_NRG1",
                                                     mouse_Pathway = "Erbb3_Nrg1",
                                                     interesting_cell_type = "MOL")
MOL_ERBB3_s <- MOL_ERBB3[1] %>% as.data.frame
MOL_ERBB3_r <- MOL_ERBB3[2] %>% as.data.frame


MOL_NRG2 <- cellphone_Cell_pathway_modification_plot(Pathway = "NRG2_ERBB4",
                                                    mouse_Pathway = "Nrg2_Erbb4",
                                                    interesting_cell_type = "MOL")
MOL_NRG2_s <- MOL_NRG2[1] %>% as.data.frame
MOL_NRG2_r <- MOL_NRG2[2] %>% as.data.frame

MOL_NRG1 <- cellphone_Cell_pathway_modification_plot(Pathway = "NRG1_ERBB4",
                                                     mouse_Pathway = "Nrg1_Erbb4",
                                                     interesting_cell_type = "MOL")
MOL_NRG1_s <- MOL_NRG1[1] %>% as.data.frame
MOL_NRG1_r <- MOL_NRG1[2] %>% as.data.frame

MOL_r <- rbind(MOL_ERBB3_r, MOL_NRG2_r, MOL_NRG1_r)
MOL_s <- rbind(MOL_ERBB3_s, MOL_NRG2_s, MOL_NRG1_s)


# Figure
c1 <- cellphonedb_plot1(MOL_r, range = c(0,18))
c2 <- cellphonedb_plot2(MOL_s, range = c(0,18))

c1 <- c1 + guides(color = guide_colourbar(order = 1),
                  size = guide_legend(order = 2))

legend <- get_legend(c1)

c1 <- c1 + theme(legend.position="none") + 
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10)) 


c2 <- c2 + theme(legend.position="none") + 
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        strip.text.x = element_blank(),
        axis.text.x = element_blank()) 

blankPlot <- ggplot() + geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

g1 <- arrangeGrob(c1, c2, blankPlot, legend, blankPlot, blankPlot, ncol = 2, 
                  widths = c(2.7, .5), heights = c(3.5, 3.5, 2.5, 2.5), 
                  layout_matrix = cbind(c(1,1,2,2), c(3,4,5,6)))


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_MOL.png",
       plot = g1,
       scale = 1, width = 7, height = 7, units = "in", device = "png",
       dpi = 300)

# New Figure (new scale)
c1 <- cellphonedb_plot1(MOL_r, range = c(0,20))
c2 <- cellphonedb_plot2(MOL_s, range = c(0,20))

c1 <- c1 + guides(color = guide_colourbar(order = 1),
                  size = guide_legend(order = 2))

legend <- get_legend(c1)

c1 <- c1 + theme(legend.position="none") + 
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10)) 


c2 <- c2 + theme(legend.position="none") + 
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        strip.text.x = element_blank(),
        axis.text.x = element_blank()) 

blankPlot <- ggplot() + geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

g1 <- arrangeGrob(c1, c2, blankPlot, legend, blankPlot, blankPlot, ncol = 2, 
                  widths = c(2.7, .5), heights = c(3.5, 3.5, 2.5, 2.5), 
                  layout_matrix = cbind(c(1,1,2,2), c(3,4,5,6)))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_MOL_New_scale.pdf",
       plot = g1,
       scale = 1, width = 7, height = 7, units = "in", device = cairo_pdf,
       dpi = 300)

###################
# 4. MFOL
###################

MFOL_FGF1 <- cellphone_Cell_pathway_modification_plot(Pathway = "FGF1_FGFR2",
                                                      mouse_Pathway = "Fgf1_Fgfr2",
                                                      interesting_cell_type = "MFOL")
MFOL_FGF1_s <- MFOL_FGF1[1] %>% as.data.frame
MFOL_FGF1_r <- MFOL_FGF1[2] %>% as.data.frame

MFOL_ERBB4 <- cellphone_Cell_pathway_modification_plot(Pathway = "ERBB4_NRG4",
                                                      mouse_Pathway = "Erbb4_Nrg4",
                                                      interesting_cell_type = "MFOL")
MFOL_ERBB4_s <- MFOL_ERBB4[1] %>% as.data.frame
MFOL_ERBB4_r <- MFOL_ERBB4[2] %>% as.data.frame

MFOL_r <- rbind(MFOL_FGF1_r, MFOL_ERBB4_r)
MFOL_s <- rbind(MFOL_FGF1_s, MFOL_ERBB4_s)

# Figure
c1 <- cellphonedb_plot1(MFOL_r, range = c(0,11))
c2 <- cellphonedb_plot2(MFOL_s, range = c(0,11))

c1 <- c1 + guides(color = guide_colourbar(order = 1),
                  size = guide_legend(order = 2))

legend <- get_legend(c1)

c1 <- c1 + theme(legend.position="none") + 
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10)) 


c2 <- c2 + theme(legend.position="none") + 
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        strip.text.x = element_blank(),
        axis.text.x = element_blank()) 

blankPlot <- ggplot() + geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

g1 <- arrangeGrob(c1, c2, blankPlot, legend, blankPlot, blankPlot, ncol = 2, 
                  widths = c(2.7, .5), heights = c(3.5, 3.5, 2.5, 2.5), 
                  layout_matrix = cbind(c(1,1,2,2), c(3,4,5,6)))


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_MFOL.png",
       plot = g1,
       scale = 1, width = 6.3, height = 7, units = "in", device = "png",
       dpi = 300)

# New Figure with new scale
c1 <- cellphonedb_plot1(MFOL_r, range = c(0,20))
c2 <- cellphonedb_plot2(MFOL_s, range = c(0,20))

c1 <- c1 + guides(color = guide_colourbar(order = 1),
                  size = guide_legend(order = 2))

legend <- get_legend(c1)

c1 <- c1 + theme(legend.position="none") + 
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10)) 


c2 <- c2 + theme(legend.position="none") + 
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        strip.text.x = element_blank(),
        axis.text.x = element_blank()) 

blankPlot <- ggplot() + geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

g1 <- arrangeGrob(c1, c2, blankPlot, legend, blankPlot, blankPlot, ncol = 2, 
                  widths = c(2.7, .5), heights = c(3.5, 3.5, 2.5, 2.5), 
                  layout_matrix = cbind(c(1,1,2,2), c(3,4,5,6)))


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_MFOL_New_scale.pdf",
       plot = g1,
       scale = 1, width = 6.3, height = 7, units = "in", device = cairo_pdf,
       dpi = 300)


###################
# 5. NFOL
###################

NFOL_FGFR2 <- cellphone_Cell_pathway_modification_plot(Pathway = "FGFR2_EPHA4",
                                                      mouse_Pathway = "Fgfr2_Epha4",
                                                      interesting_cell_type = "NFOL")
NFOL_FGFR2_s <- NFOL_FGFR2[1] %>% as.data.frame
NFOL_FGFR2_r <- NFOL_FGFR2[2] %>% as.data.frame

NFOL_FGFR1 <- cellphone_Cell_pathway_modification_plot(Pathway = "FGFR1_FGFR2",
                                                       mouse_Pathway = "Fgfr1_Fgfr2",
                                                       interesting_cell_type = "NFOL")
NFOL_FGFR1_s <- NFOL_FGFR1[1] %>% as.data.frame
NFOL_FGFR1_r <- NFOL_FGFR1[2] %>% as.data.frame

NFOL_r <- rbind(NFOL_FGFR2_r, NFOL_FGFR1_r)
NFOL_s <- rbind(NFOL_FGFR2_s, NFOL_FGFR1_s)

# Figure
c1 <- cellphonedb_plot1(NFOL_r, range = c(0,11))
c2 <- cellphonedb_plot2(NFOL_s, range = c(0,11))

c1 <- c1 + guides(color = guide_colourbar(order = 1),
                  size = guide_legend(order = 2))

legend <- get_legend(c1)

c1 <- c1 + theme(legend.position="none") + 
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10)) 


c2 <- c2 + theme(legend.position="none") + 
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        strip.text.x = element_blank(),
        axis.text.x = element_blank()) 

blankPlot <- ggplot() + geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

g1 <- arrangeGrob(c1, c2, blankPlot, legend, blankPlot, blankPlot, ncol = 2, 
                  widths = c(2.7, .5), heights = c(3.5, 3.5, 2.5, 2.5), 
                  layout_matrix = cbind(c(1,1,2,2), c(3,4,5,6)))


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_NFOL.png",
       plot = g1,
       scale = 1, width = 6.3, height = 7, units = "in", device = "png",
       dpi = 300)


###################
# 6. OPC
###################

# OPC and NFOL was switched. 
# Old one is NFOL with OPC name
# while new scale is change to NFOL

OPC_ERBB3 <- cellphone_Cell_pathway_modification_plot(Pathway = "ERBB3_NRG1",
                                                       mouse_Pathway = "Erbb3_Nrg1",
                                                       interesting_cell_type = "OPC")
OPC_ERBB3_s <- OPC_ERBB3[1] %>% as.data.frame
OPC_ERBB3_r <- OPC_ERBB3[2] %>% as.data.frame

OPC_ERBB4 <- cellphone_Cell_pathway_modification_plot(Pathway = "ERBB4_NRG4",
                                                       mouse_Pathway = "Erbb4_Nrg4",
                                                       interesting_cell_type = "OPC")
OPC_ERBB4_s <- OPC_ERBB4[1] %>% as.data.frame
OPC_ERBB4_r <- OPC_ERBB4[2] %>% as.data.frame

OPC_r <- rbind(OPC_ERBB3_r, OPC_ERBB4_r)
OPC_s <- rbind(OPC_ERBB3_s, OPC_ERBB4_s)

# Figure
c1 <- cellphonedb_plot1(OPC_r, range = c(0,19))
c2 <- cellphonedb_plot2(OPC_s, range = c(0,19))

c1 <- c1 + guides(color = guide_colourbar(order = 1),
                  size = guide_legend(order = 2))

legend <- get_legend(c1)

c1 <- c1 + theme(legend.position="none") + 
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10)) 


c2 <- c2 + theme(legend.position="none") + 
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        strip.text.x = element_blank(),
        axis.text.x = element_blank()) 

blankPlot <- ggplot() + geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

g1 <- arrangeGrob(c1, c2, blankPlot, legend, blankPlot, blankPlot, ncol = 2, 
                  widths = c(2.7, .5), heights = c(3.5, 3.5, 2.5, 2.5), 
                  layout_matrix = cbind(c(1,1,2,2), c(3,4,5,6)))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_OPC.png",
       plot = g1,
       scale = 1, width = 6.3, height = 7, units = "in", device = "png",
       dpi = 300)

# New Figure New scale
c1 <- cellphonedb_plot1(OPC_r, range = c(0,20))
c2 <- cellphonedb_plot2(OPC_s, range = c(0,20))

c1 <- c1 + guides(color = guide_colourbar(order = 1),
                  size = guide_legend(order = 2))

legend <- get_legend(c1)

c1 <- c1 + theme(legend.position="none") + 
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10)) 


c2 <- c2 + theme(legend.position="none") + 
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        strip.text.x = element_blank(),
        axis.text.x = element_blank()) 

blankPlot <- ggplot() + geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

g1 <- arrangeGrob(c1, c2, blankPlot, legend, blankPlot, blankPlot, ncol = 2, 
                  widths = c(2.7, .5), heights = c(3.5, 3.5, 2.5, 2.5), 
                  layout_matrix = cbind(c(1,1,2,2), c(3,4,5,6)))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_NFOL_New_scale.pdf",
       plot = g1,
       scale = 1, width = 6.3, height = 7, units = "in", device = cairo_pdf,
       dpi = 300)

###################
# 6. Inhibitory
###################

IN_GRIN2B <- cellphone_Cell_pathway_modification_plot(Pathway = "VEGFA_GRIN2B",
                                                      mouse_Pathway = "Vegfa_Grin2b",
                                                      interesting_cell_type = "Inhibitory_Neurons")
IN_GRIN2B_s <- IN_GRIN2B[1] %>% as.data.frame
IN_GRIN2B_r <- IN_GRIN2B[2] %>% as.data.frame

IN_EPHB2 <- cellphone_Cell_pathway_modification_plot(Pathway = "VEGFA_EPHB2",
                                                      mouse_Pathway = "Vegfa_Ephb2",
                                                      interesting_cell_type = "Inhibitory_Neurons")
IN_EPHB2_s <- IN_EPHB2[1] %>% as.data.frame
IN_EPHB2_r <- IN_EPHB2[2] %>% as.data.frame

IN_NRP2 <- cellphone_Cell_pathway_modification_plot(Pathway = "NRP2_VEGFA",
                                                     mouse_Pathway = "Nrp2_Vegfa",
                                                     interesting_cell_type = "Inhibitory_Neurons")
IN_NRP2_s <- IN_NRP2[1] %>% as.data.frame
IN_NRP2_r <- IN_NRP2[2] %>% as.data.frame

IN_NRP1 <- cellphone_Cell_pathway_modification_plot(Pathway = "NRP1_VEGFA",
                                                    mouse_Pathway = "Nrp1_Vegfa",
                                                    interesting_cell_type = "Inhibitory_Neurons")
IN_NRP1_s <- IN_NRP1[1] %>% as.data.frame
IN_NRP1_r <- IN_NRP1[2] %>% as.data.frame

IN_r <- rbind(IN_GRIN2B_r, IN_EPHB2_r, IN_NRP1_r, IN_NRP2_r)
IN_s <- rbind(IN_GRIN2B_s, IN_EPHB2_s, IN_NRP1_s, IN_NRP2_s)

# Figure
c1 <- cellphonedb_plot1(IN_r, range = c(0,12))
c2 <- cellphonedb_plot2(IN_s, range = c(0,12))

c1 <- c1 + guides(color = guide_colourbar(order = 1),
                  size = guide_legend(order = 2))

legend <- get_legend(c1)

c1 <- c1 + theme(legend.position="none") + 
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10)) 


c2 <- c2 + theme(legend.position="none") + 
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        strip.text.x = element_blank(),
        axis.text.x = element_blank()) 

blankPlot <- ggplot() + geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

g1 <- arrangeGrob(c1, c2, blankPlot, legend, blankPlot, blankPlot, ncol = 2, 
                  widths = c(2.7, .5), heights = c(3.5, 3.5, 2.5, 2.5), 
                  layout_matrix = cbind(c(1,1,2,2), c(3,4,5,6)))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_IN.png",
       plot = g1,
       scale = 1, width = 9, height = 7, units = "in", device = "png",
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_IN.pdf",
       plot = g1,
       scale = 1, width = 9, height = 7, units = "in", device = cairo_pdf,
       dpi = 300)

###################
# 7. CA2/3
###################

CA23_ERBB4 <- cellphone_Cell_pathway_modification_plot(Pathway = "ERBB4_NRG4",
                                                      mouse_Pathway = "Erbb4_Nrg4",
                                                      interesting_cell_type = "CA2_3")
CA23_ERBB4_s <- CA23_ERBB4[1] %>% as.data.frame
CA23_ERBB4_r <- CA23_ERBB4[2] %>% as.data.frame

CA2_3_NRG3 <- cellphone_Cell_pathway_modification_plot(Pathway = "NRG3_ERBB4",
                                                     mouse_Pathway = "Nrg3_Erbb4",
                                                     interesting_cell_type = "CA2_3")
CA2_3_NRG3_s <- CA2_3_NRG3[1] %>% as.data.frame
CA2_3_NRG3_r <- CA2_3_NRG3[2] %>% as.data.frame


CA23_r <- rbind(CA23_ERBB4_r, CA2_3_NRG3_r)
CA23_s <- rbind(CA23_ERBB4_s, CA2_3_NRG3_s)

# Figure
c1 <- cellphonedb_plot1(CA23_r, range = c(0,70))
c2 <- cellphonedb_plot2(CA23_s, range = c(0,70))

c1 <- c1 + guides(color = guide_colourbar(order = 1),
                  size = guide_legend(order = 2))

legend <- get_legend(c1)

c1 <- c1 + theme(legend.position="none") + 
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10)) 


c2 <- c2 + theme(legend.position="none") + 
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        strip.text.x = element_blank(),
        axis.text.x = element_blank()) 

blankPlot <- ggplot() + geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

g1 <- arrangeGrob(c1, c2, blankPlot, legend, blankPlot, blankPlot, ncol = 2, 
                  widths = c(2.7, .5), heights = c(3.5, 3.5, 2.5, 2.5), 
                  layout_matrix = cbind(c(1,1,2,2), c(3,4,5,6)))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_CA23.png",
       plot = g1,
       scale = 1, width = 6.3, height = 7, units = "in", device = "png",
       dpi = 300)

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_Cell_Cell_Interaction_Plot_session_info.txt")


