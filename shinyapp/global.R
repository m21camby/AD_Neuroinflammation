.libPaths( c( "/shinyserver/libs/AD_Neuroinflammation/prod/site-library/" , .libPaths() ) )
#.libPaths( c( "/shinyserver/libs/AD_Neuroinflammation/prod/site-library" , .libPaths() ) )
# Change to /shinyserver/apps/AD_Neuroinflammation later when Dan change

lib.path=.libPaths() #"/shinyserver/libs/AD_Neuroinflammation/prod/site-library/")

library(ggplot2, lib.loc = lib.path)
library(gridExtra, lib.loc = lib.path)
library(data.table, lib.loc = lib.path)
library(dplyr, lib.loc = lib.path)
library(stringr, lib.loc = lib.path)
library(fst, lib.loc = lib.path)
library(DT, lib.loc = lib.path)
library(tidyr, lib.loc = lib.path)
library(tibble, lib.loc = lib.path)
library(gghighlight, lib.loc = lib.path)
library(viridis, lib.loc = lib.path)
library(viridisLite, lib.loc = lib.path)
library(scales, lib.loc = lib.path)
library(RColorBrewer, lib.loc = lib.path)

######################
# 1. Data loading
######################

# QC meta data loading

df <- readRDS('/shinyserver/apps/AD_Neuroinflammation/UMAP_final.rda')

# DGE matrix selected genes
DGE <- readRDS('/shinyserver/apps/AD_Neuroinflammation/DGE_selected.rda')

# local density DataFrame
g1.df <- readRDS("/shinyserver/apps/AD_Neuroinflammation/Density1.rda")
g2.df <- readRDS("/shinyserver/apps/AD_Neuroinflammation/Density2.rda")
g3.df <- readRDS("/shinyserver/apps/AD_Neuroinflammation/Density3.rda")

# Merge meta data and gene DataFrame
gene.df <- data.frame(df, DGE)

# Trajectory data
traj.df <- readRDS("/shinyserver/apps/AD_Neuroinflammation/traj.rda")

expression.df <- readRDS("/shinyserver/apps/AD_Neuroinflammation/oligo_genes.rda")


######################
# 2. Features
######################

# QC features
QC_Features <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "Density_AD_WT", "Density_AD_ADp40KO", "Density_ADp40KO_WT")

# cell type features
celltype_Features <- c("Astrocytes", "subiculum", "MFOL", "OPC", "MOL", "Unidentified_Neurons",
                    "CA1", "Dentate_Gyrus", "Inhibitory_Neurons", "Microglia",
                    "NFOL", "CA2_3", "Pericyte", "Cajal", "VLMC", "Fibroblast", "Vascular", "Macrophage",
                    "Choroid")

# gene features
genes_Features <- colnames(DGE)


# negative regulation of oligodendrocyte differentiation
gene1 = c("Bmp4","Ctnnb1","Daam2","Dusp10","Nf1","Notch1","Sirt2","Tmem98")
# postive regulation of oligodendrocyte differentiation
gene2 = c("Aspa","Enpp2","Ptn", "Nkx2-2","Ptpra","Ptprz1","Qk","Tenm4", "Zfp365")
# negative regulation of myelination
gene3 = c("Jam2","Mtmr2","Pten","Tmem98", "Epha4")
# positive regulation of myelination
gene4 = c("Dicer1", "Mag","Myrf","Nrg1","Pard3","Sox10","Tppp","Trf")

# OPC marker genes
gene5 = c("Cdo1", "Pdgfra", "Rlbp1")
# NFOL marker genes
gene6 = c("Fyn", "Enpp6", "Tmem163")
# MOL marker genes
gene7 = c("Mbp", "Ndrg1", "Opalin", "Mal", "Mog", "Ppp1r14a", "Apod", "Trf", "Cryab")
# Housekeeping genes
gene8 = c("Gapdh","Hsp90ab1", "Pgk1", "Actb")

Oligo_Genes <- c(gene1, gene2, gene3, gene4, gene5, gene6, gene7, gene8) %>% unique

#############################
# 3. Plots functions
#############################

# 3-1. QC UMAP plots
QC1 <- function(input_feature = "nCount_RNA"){

  if(input_feature %in% c("nCount_RNA", "nFeature_RNA")){
    ggplot(data = df, aes_string(x = "UMAP1", y = "UMAP2", color = input_feature)) +
      theme_classic() + scale_color_viridis() +
      geom_point(size = 0.1)
  }else if(input_feature == "percent.mt"){
    ggplot(data = df, aes_string(x = "UMAP1", y = "UMAP2", color = input_feature)) +
      theme_classic() + scale_color_viridis(limits = c(0, 2), oob = scales::squish) +
      geom_point(size = 0.1)
  }else if(input_feature == "Density_AD_WT"){
    ggplot(g1.df, aes(x=UMAP_1,y=UMAP_2,colour=log2ratio)) +
      geom_point(size=0.2) +
      scale_colour_gradient2(low="blue", mid='gray',high="red") +
      theme_void() +
      theme(legend.position = c(0.1, 0.8),
            legend.text = element_text(size=12)) +
      guides(colour=guide_colorbar(title=paste0("WT",' vs ',"APPPS1"),
                                   title.position='right',
                                   title.hjust=.5,
                                   title.theme=element_text(angle=90,size=8),
                                   barwidth=.5))
  }else if(input_feature == "Density_AD_ADp40KO"){
    ggplot(g2.df, aes(x=UMAP_1,y=UMAP_2,colour=log2ratio)) +
      geom_point(size=0.2) +
      scale_colour_gradient2(low="blue", mid='gray',high="red") +
      theme_void() +
      theme(legend.position = c(0.1, 0.8),
            legend.text = element_text(size=12)) +
      guides(colour=guide_colorbar(title=paste0("APPPS1",' vs ',"APPPS1.il12b-/-"),
                                   title.position='right',
                                   title.hjust=.5,
                                   title.theme=element_text(angle=90,size=8),
                                   barwidth=.5))
  }else{
    ggplot(g3.df, aes(x=UMAP_1,y=UMAP_2,colour=log2ratio)) +
      geom_point(size=0.2) +
      scale_colour_gradient2(low="blue", mid='gray',high="red") +
      theme_void() +
      theme(legend.position = c(0.1, 0.8),
            legend.text = element_text(size=12)) +
      guides(colour=guide_colorbar(title=paste0("WT",' vs ',"APPPS1.il12b-/-"),
                                   title.position='right',
                                   title.hjust=.5,
                                   title.theme=element_text(angle=90,size=8),
                                   barwidth=.5))
  }
}

# 3-2. QC Violin plots
nb.cols <- 19
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

Violin_plot <- function(input_feature = "nCount_RNA"){
  if(input_feature %in% c("nCount_RNA", "nFeature_RNA", "percent.mt")){
    ggplot(df, aes_string(x="cell_type", y = input_feature, color = "cell_type")) +
      geom_violin() + theme_classic() +
      geom_jitter(shape=16, size = 0.05, alpha = 0.1, position=position_jitter(0.1)) +
      scale_color_manual(values = mycolors) +
      theme(legend.position = "none",
            axis.text.x =  element_text(angle = 90, vjust = 0.5),
            axis.title.x = element_blank())

  }else{
    ggplot() + theme_void()
  }
}

# 3-3. UMAP cell type plots
UMAP_celltype <- function(input_feature = "MOL"){
  ggplot(data = df, aes(x = UMAP1, y = UMAP2)) + theme_classic() +
  geom_point(size = 0.1, color = "darkgray") + scale_color_viridis_d() +
    theme(legend.position = "none") +
  geom_point(data = df[df$cell_type %in% input_feature, ], aes(x = UMAP1, y = UMAP2, color = "red"), size = 0.1)
}

# 3-4. UMAP gene plots
UMAP_genes <- function(input_feature = "Il12b"){
  g1 <- ggplot(data = gene.df[gene.df$sample %in% "Ctrl", ], aes_string(x = "UMAP1", y = "UMAP2", color = input_feature, order = input_feature)) +
    theme_classic() + ggtitle("WT") +
    geom_point(size = 0.5) + scale_color_viridis() + theme(plot.title = element_text(hjust = 0.5)) +
    gghighlight(input_feature > 1)

  g2 <- ggplot(data = gene.df[gene.df$sample %in% "AD", ], aes_string(x = "UMAP1", y = "UMAP2", color = input_feature, order = input_feature)) +
    theme_classic() + ggtitle("APPPS1") +
    geom_point(size = 0.5) + scale_color_viridis() + theme(plot.title = element_text(hjust = 0.5)) +
    gghighlight(input_feature > 1)

  g3 <- ggplot(data = gene.df[gene.df$sample %in% "ADp40KO", ], aes_string(x = "UMAP1", y = "UMAP2", color = input_feature, order = input_feature)) +
    theme_classic() + ggtitle("APPPS1.Il12b-/-") +
    geom_point(size = 0.5) + scale_color_viridis() + theme(plot.title = element_text(hjust = 0.5)) +
    gghighlight(input_feature > 1)

  g4 <- ggplot(data = gene.df, aes_string(x = "UMAP1", y = "UMAP2", color = input_feature, order = input_feature)) +
    theme_classic() + ggtitle("All samples") +
    geom_point(size = 0.5) + scale_color_viridis() + theme(plot.title = element_text(hjust = 0.5)) +
    gghighlight(input_feature > 1)

  grid.arrange(g1, g2, g3, g4, ncol = 2)
}

DE_genes <- function(input_feature = "Il12b"){

  ggplot(gene.df, aes_string(x= "sample", y = input_feature, color = "sample")) +
    theme_classic() +
    geom_violin(trim=TRUE) + scale_color_viridis_d() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 10),
          axis.title.x = element_blank()) + facet_wrap(~ cell_type)
}


pseudograph <- function(input_feature = "Pdgfra"){
  gene.df <- expression.df[rownames(expression.df) %in% input_feature, ] %>% t
  traj_gene.df <- cbind(traj.df, gene.df)
  colnames(traj_gene.df) <- c("pseudotime", "orig.ident", "sample", "cell_type", "cell_barcode", "cell_type_precise", input_feature)

  ggplot(traj_gene.df, aes_string(x = "pseudotime", y = input_feature, color = "sample")) +
    geom_point(size = 0.1, alpha = 0.1) + ylab("gene expression") +
    stat_smooth(method = "gam", formula = y ~ s(x), size = 1, se = FALSE) +
    scale_color_manual(labels = c("WT", "APPPS1", "APPPS1.il12b.-/-"), values = c(viridis(3)[1], viridis(3)[2], "#CC9966")) +
    theme_classic() +
    theme(axis.text = element_text(size = 12, color = "black", family = "helvetica"),
          axis.title = element_text(size = 12, color = "black", family = "helvetica"),
          #legend.position = "none",
          legend.title = element_blank(),
          legend.text = element_text(size = 11, color = "black", family = "helvetica"))
}








