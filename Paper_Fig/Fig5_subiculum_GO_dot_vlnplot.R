



SB_UP <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/GO_csv_files/edgeR_GO/cell_type_up_0.25/subiculum_AD_ADp40KO.csv_GO.csv", row.names = 1, stringsAsFactors = F)
SB_UP <- SB_UP[SB_UP$Fisher.elim < 0.05, ]
SB_UP <- SB_UP[SB_UP$Significant > 1, ]
SB_UP$cell_type <- "subiculum"
SB_UP$exp <- "APPPS1 UP"


SB_DOWN <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/GO_csv_files/edgeR_GO/cell_type_down_0.25/subiculum_AD_ADp40KO.csv_GO.csv", row.names = 1, stringsAsFactors = F)
SB_DOWN <- SB_DOWN[SB_DOWN$Fisher.elim < 0.05, ]
SB_DOWN <- SB_DOWN[SB_DOWN$Significant > 1, ]
SB_DOWN$cell_type <- "subiculum"
SB_DOWN$exp <- "APPPS1.Il12b-/- UP"


UP_spec <- setdiff(SB_UP$Term, SB_DOWN$Term)
SB_UP_spec <- SB_UP[SB_UP$Term %in% UP_spec, ]
SB_UP_terms <- c("calcium ion regulated exocytosis", "release of sequestered calcium ion into ...",
                 "potassium ion transport", "voltage-gated cation channel activity", "transmitter-gated ion channel activity",
                 "regulation of synaptic transmission, GAB...",
                 "G-protein alpha-subunit binding", "GTPase regulator activity",
                 "synapse assembly", "regulation of postsynapse organization",
                 "positive regulation of Wnt signaling pat...")

# calcium ion regulated exocytosis
# release of sequestered calcium ion into
# potassium ion transport
# cation channel activity
# voltage-gated cation channel activity
# transmitter-gated ion channel activity
# regulation of synaptic transmission, GAB

# G-protein alpha-subunit binding
# Rho GTPase binding
# GTPase regulator activity

# synapse organization
# synapse assembly
# regulation of postsynapse organization

# negative regulation of cellular response...
# positive regulation of Wnt signaling pat

DOWN_spec <- setdiff(SB_DOWN$Term, SB_UP$Term)
SB_DOWN_spec <- SB_DOWN[SB_DOWN$Term %in% DOWN_spec, ]

SB_DOWN_terms <- c("SH2 domain binding", 
                   "potassium channel regulator activity","cellular response to calcium ion", 
                   "ion channel inhibitor activity", "calcium-dependent phospholipid binding",
                   "dendrite development", "positive regulation of dendrite extensio...",
                   "negative regulation of canonical Wnt sig...", 
                   "axon ensheathment",
                   "glutamate receptor signaling pathway",
                   "positive regulation of MAPK cascade", "regulation of MAP kinase activity",
                   "phospholipid metabolic process",  "cholesterol biosynthetic process")

# SH2 domain binding
# cellular response to calcium ion
# calcium-dependent cell-cell adhesion via...
# ion channel inhibitor activity
# potassium channel regulator activity
# calcium-dependent phospholipid bindin...

# dendrite development
# positive regulation of dendrite extensio...

# negative regulation of canonical Wnt sig...	
# axon ensheathment

# glutamate receptor signaling pathwa
# positive regulation of MAPK cascade
# 	regulation of MAP kinase activity

# phospholipid metabolic process
# phospholipid binding
# cholesterol biosynthetic process


SB_UP_Sub <- SB_UP[SB_UP$Term %in% SB_UP_terms, ]
SB_DOWN_Sub <- SB_DOWN[SB_DOWN$Term %in% SB_DOWN_terms, ]

SB_Sub <- rbind(SB_UP_Sub, SB_DOWN_Sub)

SB_Sub[SB_Sub == "regulation of synaptic transmission, GAB..."] <- "regulation of synaptic transmission, GABAergic"
SB_Sub[SB_Sub == "positive regulation of Wnt signaling pat..."] <- "positive regulation of Wnt signaling pathway"
SB_Sub[SB_Sub == "release of sequestered calcium ion into ..."] <- "release of sequestered calcium ion into cytosol"
SB_Sub[SB_Sub == "negative regulation of canonical Wnt sig..."] <- "negative regulation of canonical Wnt signaling"
SB_Sub[SB_Sub == "positive regulation of dendrite extensio..."] <- "positive regulation of dendrite extension"

SB_UP_terms[SB_UP_terms == "release of sequestered calcium ion into ..."] <- "release of sequestered calcium ion into cytosol"
SB_UP_terms[SB_UP_terms == "regulation of synaptic transmission, GAB..."] <- "regulation of synaptic transmission, GABAergic"
SB_UP_terms[SB_UP_terms == "positive regulation of Wnt signaling pat..."] <- "positive regulation of Wnt signaling pathway"

SB_DOWN_terms[SB_DOWN_terms == "negative regulation of canonical Wnt sig..."] <- "negative regulation of canonical Wnt signaling"
SB_DOWN_terms[SB_DOWN_terms == "positive regulation of dendrite extensio..."] <- "positive regulation of dendrite extension"


# order
SB_Sub$Term <- factor(SB_Sub$Term, levels = c(SB_UP_terms, SB_DOWN_terms))


g1 <- ggplot(SB_Sub, 
             aes(x=exp, y=Term,
                 colour=Fisher.elim,
                 size=gene_ratio)) +
  geom_point() +
  expand_limits(x=0) + 
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
g1

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_subiculum_GO_dotplot.pdf",
       plot = g1,
       scale = 1, width = 6, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)


write.csv(SB_Sub, file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_subiculum_GO.csv")

# Violin plot

# Data loading
subi_AD_ADp40KO_DE <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/DE_csv_files/edgeR/cell_type/subiculum_AD_ADp40KO.csv", row.names = 1)

# 1st
gene <- SB_Sub[1, ]$genes
gene <- strsplit(gene, ",")  
gene <- unlist(gene, use.names=FALSE)
gene <- trimws(gene)

subi_AD_ADp40KO_DE1 <- subi_AD_ADp40KO_DE[subi_AD_ADp40KO_DE$gene %in% c(gene), ]
subi_AD_ADp40KO_DE1$pathway <- SB_Sub[1, ]$Term

for(i in c(2:25)){
  gene <- SB_Sub[i, ]$genes
  gene <- strsplit(gene, ",")  
  gene <- unlist(gene, use.names=FALSE)
  gene <- trimws(gene)
  
  subi_AD_ADp40KO_DE2 <- subi_AD_ADp40KO_DE[subi_AD_ADp40KO_DE$gene %in% c(gene), ]
  subi_AD_ADp40KO_DE2$pathway <- SB_Sub[i, ]$Term
  
  subi_AD_ADp40KO_DE1 <- rbind(subi_AD_ADp40KO_DE1, subi_AD_ADp40KO_DE2)
}

g2 <- ggplot(subi_AD_ADp40KO_DE1, aes(x=pathway, y=logFC)) + 
  geom_violin(trim=FALSE) + 
  coord_flip() + 
  geom_jitter(data = subi_AD_ADp40KO_DE1, aes(color = logFC > 0),  position=position_jitter(0.1), alpha = 0.3) + 
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"), name = "logFC", labels = c("minus", "plus")) + 
  theme_minimal() + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10, family = "helvetica", color = "black"),
        axis.text = element_text(size = 8, family = "helvetica", color = "black"),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(c(-0.2), "lines"))
g2


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_subiculum_GO_vlnplot.pdf",
       plot = g2,
       scale = 1, width = 6, height = 5, units = "in", device = cairo_pdf,
       dpi = 300)


# Combined plot
g2 <- ggplot(subi_AD_ADp40KO_DE1, aes(x=pathway, y=logFC)) + 
  geom_violin(trim=FALSE) + 
  coord_flip() + 
  geom_jitter(data = subi_AD_ADp40KO_DE1, aes(color = logFC > 0),  position=position_jitter(0.1), alpha = 0.3) + 
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"), name = "logFC", labels = c("minus", "plus")) + 
  theme_minimal() + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(), 
        axis.title.x = element_text(size = 10, family = "helvetica", color = "black"),
        axis.text.x = element_text(size = 8, family = "helvetica", color = "black"),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(c(-0.2), "lines"))

blankPlot <- ggplot() + geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

g3 <- grid.arrange(g1, blankPlot, g2, blankPlot,  ncol = 2, widths = c(1.6, 0.9), heights = c(0.3, 3.4, 1.4), layout_matrix = cbind(c(1,1,1), c(2,3,4)))
g3
g3 <- arrangeGrob(g1, blankPlot, g2, blankPlot,  ncol = 2, widths = c(1.6, 0.9), heights = c(0.0, 4, 0.7), layout_matrix = cbind(c(1,1,1), c(2,3,4)))


ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_subiculum_GO_dot_vln_plot.pdf",
       plot = g3,
       scale = 1, width = 12, height = 8, units = "in", device = cairo_pdf,
       dpi = 300)
######################################
# Shirin Asked 
######################################
# Also included APPPS1 vs WT & APPPS1.Il12b-/- vs WT

SB_UP <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/GO_csv_files/edgeR_GO/cell_type_up_0.25/subiculum_AD_ADp40KO.csv_GO.csv", row.names = 1, stringsAsFactors = F)
SB_UP <- SB_UP[SB_UP$Fisher.elim < 0.05, ]
SB_UP <- SB_UP[SB_UP$Significant > 1, ]
SB_UP$cell_type <- "subiculum"
SB_UP$exp <- "APPPS1 vs APPPS1.Il12b-/- UP"


SB_DOWN <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/GO_csv_files/edgeR_GO/cell_type_down_0.25/subiculum_AD_ADp40KO.csv_GO.csv", row.names = 1, stringsAsFactors = F)
SB_DOWN <- SB_DOWN[SB_DOWN$Fisher.elim < 0.05, ]
SB_DOWN <- SB_DOWN[SB_DOWN$Significant > 1, ]
SB_DOWN$cell_type <- "subiculum"
SB_DOWN$exp <- "APPPS1.Il12b-/- vs APPPS1 UP"

SB_AD_UP <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/GO_csv_files/edgeR_GO/cell_type_up_0.25/subiculum_AD_Ctrl.csv_GO.csv", row.names = 1, stringsAsFactors = F)
SB_AD_UP <- SB_AD_UP[SB_AD_UP$Fisher.elim < 0.05, ]
SB_AD_UP <- SB_AD_UP[SB_AD_UP$Significant > 1, ]
SB_AD_UP$cell_type <- "subiculum"
SB_AD_UP$exp <- "APPPS1 vs WT UP"


SB_ADp40KO_UP <- read.csv("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/GO_csv_files/edgeR_GO/cell_type_up_0.25/subiculum_ADp40KO_Ctrl.csv_GO.csv", row.names = 1, stringsAsFactors = F)
SB_ADp40KO_UP <- SB_ADp40KO_UP[SB_ADp40KO_UP$Fisher.elim < 0.05, ]
SB_ADp40KO_UP <- SB_ADp40KO_UP[SB_ADp40KO_UP$Significant > 1, ]
SB_ADp40KO_UP$cell_type <- "subiculum"
SB_ADp40KO_UP$exp <- "APPPS1.Il12b-/- vs WT UP"


# Common
common <- intersect(SB_AD_UP$Term, SB_ADp40KO_UP$Term)
AD_spec <- setdiff(SB_AD_UP$Term, SB_ADp40KO_UP$Term)
ADp40KO_spec <- setdiff(SB_ADp40KO_UP$Term, SB_AD_UP$Term)

AD_spec_terms <- c("cadherin binding", "PDZ domain binding", "ligand-gated anion channel activity","excitatory extracellular ligand-gated io...",
"GTPase regulator activity", "negative regulation of cytokine-mediated...", "negative regulation of axon guidance","regulation of axonogenesis"
)

ADp40KO_spec_terms <- c("integrin binding", "voltage-gated channel activity", "ephrin receptor binding", "phospholipid binding",
                        "regulation of axon guidance", "regulation of neurotransmitter transport", "response to calcium ion")

common_terms <- c("calcium ion binding", "potassium channel activity", "axon guidance", "positive regulation of lipid biosyntheti...",
                  "transmembrane receptor protein tyrosine ...","positive regulation of MAPK cascade")

SB_UP_terms <- c("calcium ion regulated exocytosis", "release of sequestered calcium ion into ...",
                 "potassium ion transport", "voltage-gated cation channel activity", "transmitter-gated ion channel activity",
                 "regulation of synaptic transmission, GAB...",
                 "G-protein alpha-subunit binding", "GTPase regulator activity",
                 "synapse assembly", "regulation of postsynapse organization",
                 "positive regulation of Wnt signaling pat...")

SB_DOWN_terms <- c("SH2 domain binding", 
                   "potassium channel regulator activity","cellular response to calcium ion", 
                   "ion channel inhibitor activity", "calcium-dependent phospholipid binding",
                   "dendrite development", "positive regulation of dendrite extensio...",
                   "negative regulation of canonical Wnt sig...", 
                   "axon ensheathment",
                   "glutamate receptor signaling pathway",
                   "positive regulation of MAPK cascade", "regulation of MAP kinase activity",
                   "phospholipid metabolic process",  "cholesterol biosynthetic process")

SB_Sub <- rbind(SB_UP, SB_DOWN, SB_AD_UP, SB_ADp40KO_UP)

unique_terms <- c(AD_spec_terms, common_terms, ADp40KO_spec_terms, SB_UP_terms, SB_DOWN_terms) %>% unique()

SB_Sub <- SB_Sub[SB_Sub$Term %in% unique_terms, ]

SB_Sub[SB_Sub == "regulation of synaptic transmission, GAB..."] <- "regulation of synaptic transmission, GABAergic"
SB_Sub[SB_Sub == "positive regulation of Wnt signaling pat..."] <- "positive regulation of Wnt signaling pathway"
SB_Sub[SB_Sub == "release of sequestered calcium ion into ..."] <- "release of sequestered calcium ion into cytosol"
SB_Sub[SB_Sub == "negative regulation of canonical Wnt sig..."] <- "negative regulation of canonical Wnt signaling"
SB_Sub[SB_Sub == "positive regulation of dendrite extensio..."] <- "positive regulation of dendrite extension"
SB_Sub[SB_Sub == "excitatory extracellular ligand-gated io..."] <- "excitatory extracellular ligand-gated  ion channel activity"
SB_Sub[SB_Sub == "negative regulation of cytokine-mediated..."] <- "negative regulation of cytokine-mediated signaling pathway"
SB_Sub[SB_Sub == "positive regulation of lipid biosyntheti..."] <- "positive regulation of lipid biosynthetic process"
SB_Sub[SB_Sub == "transmembrane receptor protein tyrosine ..."] <- "transmembrane receptor protein tyrosine  kinase signalling pathway"


AD_spec_terms[AD_spec_terms == "excitatory extracellular ligand-gated io..."] <- "excitatory extracellular ligand-gated  ion channel activity"
AD_spec_terms[AD_spec_terms == "negative regulation of cytokine-mediated..."] <- "negative regulation of cytokine-mediated signaling pathway"

common_terms[common_terms == "positive regulation of lipid biosyntheti..."] <- "positive regulation of lipid biosynthetic process"
common_terms[common_terms == "transmembrane receptor protein tyrosine ..."] <- "transmembrane receptor protein tyrosine  kinase signalling pathway"

SB_UP_terms[SB_UP_terms == "release of sequestered calcium ion into ..."] <- "release of sequestered calcium ion into cytosol"
SB_UP_terms[SB_UP_terms == "regulation of synaptic transmission, GAB..."] <- "regulation of synaptic transmission, GABAergic"
SB_UP_terms[SB_UP_terms == "positive regulation of Wnt signaling pat..."] <- "positive regulation of Wnt signaling pathway"

SB_DOWN_terms[SB_DOWN_terms == "negative regulation of canonical Wnt sig..."] <- "negative regulation of canonical Wnt signaling"
SB_DOWN_terms[SB_DOWN_terms == "positive regulation of dendrite extensio..."] <- "positive regulation of dendrite extension"



# order
SB_Sub$Term <- factor(SB_Sub$Term, levels = c(AD_spec_terms, common_terms, ADp40KO_spec_terms, setdiff(SB_UP_terms, c(AD_spec_terms, common_terms, ADp40KO_spec_terms)), setdiff(SB_DOWN_terms, c(AD_spec_terms, common_terms, ADp40KO_spec_terms))))

SB_Sub$exp <- factor(SB_Sub$exp, levels = c("APPPS1 vs WT UP", "APPPS1.Il12b-/- vs WT UP", "APPPS1 vs APPPS1.Il12b-/- UP", "APPPS1.Il12b-/- vs APPPS1 UP"))

g1 <- ggplot(SB_Sub, 
             aes(x=exp, y=Term,
                 colour=Fisher.elim,
                 size=gene_ratio)) +
  geom_point() +
  expand_limits(x=0) + 
  labs(x="group", y="GO term", colour="p-value", size="gene_ratio") +
  theme_cowplot() + 
  theme(axis.title.y = element_blank(),
        axis.title.x =  element_blank(),
        axis.text.x = element_text(size = 12, angle = 90, vjust=0.5, family = "helvetica", color = "black"),
        axis.text.y = element_text(size = 10, family = "helvetica", color = "black"),
        strip.text.y = element_text(size = 10, angle = 0, hjust = 0, vjust= 0, family = "helvetica", color = "black"),
        strip.background = element_rect(color = "white", fill="white"),
        panel.spacing = unit(c(-.2), "lines")) + 
  scale_color_gradient(low="#003366", high="#FF9900", trans = 'reverse') + 
  guides(color = guide_colourbar(order = 1), 
         size = guide_legend(order = 2)) + guides(color = guide_colorbar(reverse=TRUE))
g1

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig5_subiculum_GO_dotplot_2nd.pdf",
       plot = g1,
       scale = 1, width = 6, height = 8.5, units = "in", device = cairo_pdf,
       dpi = 300)
