#! /bin/env RScript
# written by SJK at 4. Oct. 2020
# This is for entropy analysis 
# As Shirin asked for entropy analysis by genotype, I created by whole plot


library(viridis)
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)

source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

load(file="/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200210_figures_for_meeting_report_cleaned_Seurat_object_res0.8.Robj")


meta.df <- data9set_cleaned.SO@meta.data %>% mutate(cell_type = case_when(seurat_clusters %in% c(0, 10) ~ "Dentate_Gyrus",
                                                                          seurat_clusters %in% c(2, 13, 18, 40) ~ "CA1",
                                                                          seurat_clusters %in% c(9, 21, 27, 32) ~ "CA2_3",
                                                                          seurat_clusters %in% c(11, 14, 15, 17, 19, 20) ~ "subiculum",
                                                                          seurat_clusters %in% c(24, 29) ~ "Unidentified_Neurons",
                                                                          seurat_clusters %in% c(7, 16, 22, 30) ~ "Inhibitory_Neurons",
                                                                          seurat_clusters %in% c(1,5, 6) ~ "OL", 
                                                                          seurat_clusters %in% c(12, 38) ~ "OPC",
                                                                          seurat_clusters %in% c(3, 8) ~ "Microglia",
                                                                          seurat_clusters %in% c(4) ~ "Astrocytes",
                                                                          seurat_clusters %in% c(36) ~ "Vascular",
                                                                          seurat_clusters %in% c(39) ~ "VLMC",
                                                                          seurat_clusters %in% c(26) ~ "Choroid",
                                                                          seurat_clusters %in% c(23) ~ "Fibroblast",
                                                                          seurat_clusters %in% c(28) ~ "Cajal",
                                                                          seurat_clusters %in% c(35) ~ "Pericyte",
                                                                          seurat_clusters %in% c(37) ~ "Macrophage"))


# Neuronal large cell type

meta.df2 <- meta.df %>% mutate(large_Neuronal_cell_type = case_when(cell_type %in% c("Dentate_Gyrus", "CA1","CA2_3", "subiculum", "Unidentified_Neurons") ~ "EX",
                                                                    cell_type %in% c("Inhibitory_Neurons") ~ "IN",
                                                                    cell_type %in% c("OL") ~ "OL", 
                                                                    cell_type %in% c("OPC") ~ "OPC",
                                                                    cell_type %in% c("Microglia") ~ "Microglia",
                                                                    cell_type %in% c("Astrocytes") ~ "Astrocytes",
                                                                    cell_type %in% c("Vascular", "VLMC") ~ "Vascular",
                                                                    cell_type %in% c("Choroid") ~ "Choroid",
                                                                    cell_type %in% c("Fibroblast") ~ "Fibroblast",
                                                                    cell_type %in% c("Cajal") ~ "Cajal",
                                                                    cell_type %in% c("Pericyte") ~ "Pericyte",
                                                                    cell_type %in% c("Macrophage") ~ "Macrophage"))

data9set_cleaned.SO$cell_type <- meta.df$cell_type

# give levels to large cell type
data9set_cleaned.SO$cell_type <- factor(data9set_cleaned.SO$cell_type, levels = c("Dentate_Gyrus", "CA1", "CA2_3", "subiculum",  "Unidentified_Neurons", 
                                                                                  "Inhibitory_Neurons", "OL", "OPC", "Microglia","Astrocytes",
                                                                                  "Vascular", "VLMC", "Choroid", "Fibroblast", "Cajal", "Pericyte", "Macrophage"))

data9set_cleaned.SO$large_cell_type <- meta.df2$large_Neuronal_cell_type

data9set_cleaned.SO$large_cell_type <- factor(data9set_cleaned.SO$large_cell_type, levels = c("EX", "IN", "OL", "OPC", "Microglia","Astrocytes",
                                                                                  "Vascular", "Choroid", "Fibroblast", "Cajal", "Pericyte", "Macrophage"))



get_entropy <- function(sobj,
                        groups,
                        k.param=30,
                        reduction.type='pca',
                        dims.use=1:30, sobj_meta.data = data9set.meta) {
  
  require(RANN)
  
  data.use <- Embeddings(sobj,
                         reduction.type=reduction.type, 
                         dims.use=dims.use)
  message('constructing knn network')
  my.knn <- nn2(data=data.use,
                k=k.param,
                searchtype='standard',
                eps=0)
  
  get_ent <- function(neighbors, md, group, qref) {
    q <- table(md[neighbors, group])/length(neighbors)
    sum(q*log(q/qref),na.rm=TRUE)
  }
  
  message('calculcating entropies')
  
  entropies <- list()
  for (group in groups) {
    qref <- table(sobj_meta.data[,group])/dim(sobj_meta.data)[1]
    entropies[[paste0(group,'_entropy')]] <- apply(my.knn$nn.idx, 1,
                                                   get_ent, sobj_meta.data, group, qref)
  }
  data.frame(entropies, row.names= rownames(sobj_meta.data))
  #data.frame(entropies)
}

set.seed(1)

data9set.meta <- data9set_cleaned.SO@meta.data %>% as.data.frame

data9set.suffle_meta <- data9set.meta
data9set.suffle_meta$gemgroup = sample(data9set.suffle_meta$gemgroup, replace=FALSE)
data9set.suffle_meta$sample = sample(data9set.suffle_meta$sample, replace=FALSE)

entropy.dt <- get_entropy(data9set_cleaned.SO, groups = c('sample','gemgroup'), sobj_meta.data = data9set.meta)

# Shuffled data
entropy_suffled.dt <- get_entropy(data9set_cleaned.SO, groups = c('sample','gemgroup'), sobj_meta.data = data9set.suffle_meta)


entropy2.dt <- cbind(data9set.meta, entropy.dt)

whole2.dt <- rbind(data.frame(entropy = entropy2.dt[,11], test = "genotype", cell_type = entropy2.dt[,9]),
                   data.frame(entropy = entropy2.dt[,12], test = "orig.ident", cell_type = entropy2.dt[,9]))



whole2.dt$cell_type <- factor(whole2.dt$cell_type, levels = c("Dentate_Gyrus", "CA1", "CA2_3", "subiculum",  "Unidentified_Neurons", 
                                                              "Inhibitory_Neurons", "Cajal", "Choroid",
                                                              "Astrocytes", "Microglia","Macrophage", "OL", "OPC", 
                                                              "Fibroblast",  "Pericyte", "Vascular", "VLMC"))


whole3.dt <- rbind(data.frame(entropy = entropy2.dt[,11], test = "genotype", cell_type = entropy2.dt[,10]),
                   data.frame(entropy = entropy2.dt[,12], test = "orig.ident", cell_type = entropy2.dt[,10]))

whole3.dt$cell_type <- factor(whole3.dt$cell_type, levels = c("EX", "IN", "Cajal", "Choroid",
                                                              "Astrocytes", "Microglia","Macrophage", "OL", "OPC", 
                                                              "Fibroblast",  "Pericyte", "Vascular"))



whole.dt <- rbind(data.frame(entropy = entropy.dt[,1], test = "genotype", orig = "original"),
                  data.frame(entropy = entropy.dt[,2], test = "orig.ident", orig = "original"),
                  data.frame(entropy = entropy_suffled.dt[,1], test = "genotype", orig = "randomized"),
                  data.frame(entropy = entropy_suffled.dt[,2], test = "orig.ident", orig = "randomized"))

g0 <- ggplot(whole.dt, aes(x = test, y = entropy, fill = orig)) + geom_boxplot(coef = 6) + 
  theme_classic() +
  ylab("relative entropy") + ylim(c(0,2)) + 
  theme(axis.text.x = element_text(size =15, angle = 45, color = "black", vjust = 1, hjust = 1, family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size =15, color = "black", family = "Helvetica"),
        axis.title.y = element_text(size =15, family = "Helvetica"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15, family = "Helvetica"),
        legend.key.size = unit(1, "cm")) + scale_fill_manual(values=c("#E69F00", "#56B4E9"))  
g0
ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig1_Entropy_by_genotype_all_shuffled.pdf", 
       plot = g0, 
       scale = 1, width = 6.5, height = 4, units = "in", device = cairo_pdf,
       dpi = 300)

cols = c("#CC0000", "#CC0000", "#CC0000", "#CC0000", "#CC0000", 
         "#FF6600", "#FF9900","#CC9900", "#FFCC66", "#99CC00","#003300", "#99CCCC", "#009999", "#3399FF","#3399FF", "#000033","#000033")


g1 <- ggplot(whole2.dt, aes(x = cell_type, y = entropy, fill = cell_type)) + geom_boxplot(coef = 6) + theme_classic() +
  ylab("relative entropy") + ylim(c(0,1.5)) + 
  theme(axis.text.x = element_text(size =12, angle = 50, color = "black", vjust = 1, hjust = 1, family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size =15, color = "black", family = "Helvetica"),
        axis.title.y = element_text(size =15, family = "Helvetica"),
        legend.position = "none",
        strip.text = element_text(size=15, family = "Helvetica")) + facet_wrap(~test, dir = "v") + 
  scale_fill_manual(values=cols) + 
  scale_x_discrete(labels=c("Dentate_Gyrus" = "Dentate Gyrus", "CA1", "CA2_3" = "CA2/3", "subiculum",  
                            "Unidentified_Neurons" = "Unidentified Neurons", 
                            "Inhibitory_Neurons" = "Inhibitory Neurons", "Cajal" = "Cajal Retzius", 
                            "Choroid" = "Choroid Plexus",
                            "Astrocytes", "Microglia","Macrophage", 
                            "OL" = "Oligodendrocyte", "OPC", 
                            "Fibroblast",  "Pericyte", "Vascular", "VLMC"))
                                                                                                                         

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig1_Entropy_by_genotype_all.png", 
       plot = g1, 
       scale = 1, width = 8, height = 8, units = "in", device = "png",
       dpi = 300)


cols = c("#CC0000", 
         "#FF6600", "#FF9900","#CC9900", "#FFCC66", "#99CC00","#003300", "#99CCCC", "#009999", "#3399FF","#3399FF", "#000033")


g2 <- ggplot(whole3.dt, aes(x = cell_type, y = entropy, fill = cell_type)) + geom_boxplot(coef = 6) + theme_classic() +
  ylab("relative entropy") + ylim(c(0,1.5)) + 
  theme(axis.text.x = element_text(size =12, angle = 50, color = "black", vjust = 1, hjust = 1, family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size =15, color = "black", family = "Helvetica"),
        axis.title.y = element_text(size =15, family = "Helvetica"),
        legend.position = "none",
        strip.text = element_text(size=15, family = "Helvetica")) + facet_wrap(~test, dir = "v") + 
  scale_fill_manual(values=cols) + 
  scale_x_discrete(labels=c("EX" = "Excitatory Neurons", "IN" = "Inhibitory Neurons", "Cajal" = "Cajal Retzius", 
                            "Choroid" = "Choroid Plexus",
                            "Astrocytes", "Microglia","Macrophage", 
                            "OL" = "Oligodendrocyte", "OPC", 
                            "Fibroblast",  "Pericyte", "Vascular"))


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig1_Entropy_by_genotype_all_large_cell_type.pdf", 
       plot = g2, 
       scale = 1, width = 8, height = 8, units = "in", device = cairo_pdf,
       dpi = 300)



###################
# By genotype
###################

# Ctrl

data9set_sub.SO <- subset(data9set_cleaned.SO, subset = sample %in% c("Ctrl"))

data9set_ctrl.meta <- data9set_sub.SO@meta.data %>% as.data.frame

entropy_ctrl.dt <- get_entropy(data9set_sub.SO, groups = c('gemgroup'), sobj_meta.data = data9set_ctrl.meta)

entropy2_ctrl.dt <- cbind(data9set_ctrl.meta, entropy_ctrl.dt)

whole2_ctrl.dt <- data.frame(entropy = entropy2_ctrl.dt[,10], test = "orig.ident", cell_type = entropy2_ctrl.dt[,9])

whole2_ctrl.dt$exp <- "WT"


ggplot(whole2_ctrl.dt, aes(x = cell_type, y = entropy, fill = cell_type)) + geom_boxplot(coef = 6) + theme_classic() +
  ylab("relative entropy") + ylim(c(0,1.5)) + 
  theme(axis.text.x = element_text(size =12, angle = 50, color = "black", vjust = 1, hjust = 1, family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size =15, color = "black", family = "Helvetica"),
        axis.title.y = element_text(size =15, family = "Helvetica"),
        legend.position = "none",
        strip.text = element_text(size=15, family = "Helvetica")) + 
  scale_fill_manual(values=cols) + 
  scale_x_discrete(labels=c("Dentate_Gyrus" = "Dentate Gyrus", "CA1", "CA2_3" = "CA2/3", "subiculum",  
                            "Unidentified_Neurons" = "Unidentified Neurons", 
                            "Inhibitory_Neurons" = "Inhibitory Neurons", "Cajal" = "Cajal Retzius", 
                            "Choroid" = "Choroid Plexus",
                            "Astrocytes", "Microglia","Macrophage", 
                            "OL" = "Oligodendrocyte", "OPC", 
                            "Fibroblast",  "Pericyte", "Vascular", "VLMC"))

# AD

data9set_sub_AD.SO <- subset(data9set_cleaned.SO, subset = sample %in% c("AD"))

data9set_AD.meta <- data9set_sub_AD.SO@meta.data %>% as.data.frame

entropy_AD.dt <- get_entropy(data9set_sub_AD.SO, groups = c('gemgroup'), sobj_meta.data = data9set_AD.meta)

entropy2_AD.dt <- cbind(data9set_AD.meta, entropy_AD.dt)

whole2_AD.dt <- data.frame(entropy = entropy2_AD.dt[,10], test = "orig.ident", cell_type = entropy2_AD.dt[,9])

whole2_AD.dt$exp <- "APPPS1"

ggplot(whole2_AD.dt, aes(x = cell_type, y = entropy, fill = cell_type)) + geom_boxplot(coef = 6) + theme_classic() +
  ylab("relative entropy") + ylim(c(0,1.5)) + 
  theme(axis.text.x = element_text(size =12, angle = 50, color = "black", vjust = 1, hjust = 1, family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size =15, color = "black", family = "Helvetica"),
        axis.title.y = element_text(size =15, family = "Helvetica"),
        legend.position = "none",
        strip.text = element_text(size=15, family = "Helvetica")) + 
  scale_fill_manual(values=cols) + 
  scale_x_discrete(labels=c("Dentate_Gyrus" = "Dentate Gyrus", "CA1", "CA2_3" = "CA2/3", "subiculum",  
                            "Unidentified_Neurons" = "Unidentified Neurons", 
                            "Inhibitory_Neurons" = "Inhibitory Neurons", "Cajal" = "Cajal Retzius", 
                            "Choroid" = "Choroid Plexus",
                            "Astrocytes", "Microglia","Macrophage", 
                            "OL" = "Oligodendrocyte", "OPC", 
                            "Fibroblast",  "Pericyte", "Vascular", "VLMC"))

# ADp40KO

data9set_sub_ADp40KO.SO <- subset(data9set_cleaned.SO, subset = sample %in% c("ADp40KO"))

data9set_ADp40KO.meta <- data9set_sub_ADp40KO.SO@meta.data %>% as.data.frame

entropy_ADp40KO.dt <- get_entropy(data9set_sub_ADp40KO.SO, groups = c('gemgroup'), sobj_meta.data = data9set_ADp40KO.meta)

entropy2_ADp40KO.dt <- cbind(data9set_ADp40KO.meta, entropy_ADp40KO.dt)

whole2_ADp40KO.dt <- data.frame(entropy = entropy2_ADp40KO.dt[,10], test = "orig.ident", cell_type = entropy2_ADp40KO.dt[,9])

whole2_ADp40KO.dt$exp <- "APPPS1.il12b-/-"

ggplot(whole2_ADp40KO.dt, aes(x = cell_type, y = entropy, fill = cell_type)) + geom_boxplot(coef = 6) + theme_classic() +
  ylab("relative entropy") + ylim(c(0,1.5)) + 
  theme(axis.text.x = element_text(size =12, angle = 50, color = "black", vjust = 1, hjust = 1, family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size =15, color = "black", family = "Helvetica"),
        axis.title.y = element_text(size =15, family = "Helvetica"),
        legend.position = "none",
        strip.text = element_text(size=15, family = "Helvetica")) + 
  scale_fill_manual(values=cols) + 
  scale_x_discrete(labels=c("Dentate_Gyrus" = "Dentate Gyrus", "CA1", "CA2_3" = "CA2/3", "subiculum",  
                            "Unidentified_Neurons" = "Unidentified Neurons", 
                            "Inhibitory_Neurons" = "Inhibitory Neurons", "Cajal" = "Cajal Retzius", 
                            "Choroid" = "Choroid Plexus",
                            "Astrocytes", "Microglia","Macrophage", 
                            "OL" = "Oligodendrocyte", "OPC", 
                            "Fibroblast",  "Pericyte", "Vascular", "VLMC"))



whole2_all.dt <- rbind(whole2_ctrl.dt, whole2_AD.dt, whole2_ADp40KO.dt)


whole2_all.dt$cell_type <- factor(whole2_all.dt$cell_type, levels = c("Dentate_Gyrus", "CA1", "CA2_3", "subiculum",  "Unidentified_Neurons", 
                                                              "Inhibitory_Neurons", "Cajal", "Choroid",
                                                              "Astrocytes", "Microglia","Macrophage", "OL", "OPC", 
                                                              "Fibroblast",  "Pericyte", "Vascular", "VLMC"))

whole2_all.dt$exp <- factor(whole2_all.dt$exp, levels = c("WT", "APPPS1", "APPPS1.il12b-/-"))

g2 <- ggplot(whole2_all.dt, aes(x = cell_type, y = entropy, fill = exp)) + geom_boxplot(coef = 6) + theme_classic() +
  ylab("relative entropy") + ylim(c(0,1.5)) + 
  theme(axis.text.x = element_text(size =12, angle = 50, color = "black", vjust = 1, hjust = 1, family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size =15, color = "black", family = "Helvetica"),
        axis.title.y = element_text(size =15, family = "Helvetica"),
        legend.position = c(0.9,0.9),
        legend.title = element_blank(),
        strip.text = element_text(size=15, family = "Helvetica")) + 
  scale_fill_manual(values= c("#440154FF", viridis(3)[2], "#CC9900")) + 
  scale_x_discrete(labels=c("Dentate_Gyrus" = "Dentate Gyrus", "CA1", "CA2_3" = "CA2/3", "subiculum",  
                            "Unidentified_Neurons" = "Unidentified Neurons", 
                            "Inhibitory_Neurons" = "Inhibitory Neurons", "Cajal" = "Cajal Retzius", 
                            "Choroid" = "Choroid Plexus",
                            "Astrocytes", "Microglia","Macrophage", 
                            "OL" = "Oligodendrocyte", "OPC", 
                            "Fibroblast",  "Pericyte", "Vascular", "VLMC"))

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig1_Entropy_by_genotype.png", 
       plot = g2, 
       scale = 1, width = 10, height = 5, units = "in", device = "png",
       dpi = 300)


##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Supplement_Fig1_Entropy_by_genotype_session_info.txt")
