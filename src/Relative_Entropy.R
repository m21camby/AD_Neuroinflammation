#! /bin/env RScript
# written by SJK at 4. Oct. 2020
# This is for entropy analysis 
# entropy analysis by genotype

# ---------------------------- #
# function
# ---------------------------- #

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

# ---------------------------- #
# control entropy by replicate
# ---------------------------- #

entropy_ctrl.dt <- get_entropy(data9set_sub.SO, groups = c('gemgroup'), sobj_meta.data = data9set_ctrl.meta)
entropy2_ctrl.dt <- cbind(data9set_ctrl.meta, entropy_ctrl.dt)
whole2_ctrl.dt <- data.frame(entropy = entropy2_ctrl.dt[,10], test = "orig.ident", cell_type = entropy2_ctrl.dt[,9])
whole2_ctrl.dt$exp <- "WT"

# plot
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



