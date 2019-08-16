# Written by SJK 08.08 2019

# This R script is for DroNc-seq AD project
# This file aims for GO Analysis of scRNA-seq data
# Prerequisites are foreground genes which are interesting genes and background genes
# For background genes, I used FindMarkers function from Seurat with parameter logfc.threshold = 0.01 

# Ref1: https://github.com/karthikshekhar/CellTypeMIMB/blob/master/utilities.R
# Ref2: https://rpubs.com/kshekhar/349874
# Ref3: https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf

library(topGO)
library(org.Mm.eg.db)
library(ggplot2)

scRNA_topGO <- function(fg.genes = NULL,
                       bg.genes = NULL,
                       organism = "Mouse", 
                       stats.use = "fisher",
                       algorithm.use = "elim",
                       topnodes.print=20,
                       num.char=100){

if (is.null(fg.genes) | is.null(bg.genes)){
    stop("Error : Both gene lists are empty")}

# 1. Predefined list of interesting genes

n <- length(bg.genes)
geneList <- integer(n)
names(geneList) <- bg.genes
geneList[intersect(names(geneList), fg.genes)] <- 1
geneList <- factor(geneList)
# geneList object is a named factor that indicates which genes are interesting and which not.


# 2. Performing test 
# ThetopGOdataobject contains the list of genes, the list of interesting genes and the gene scores
# gene scores (optional) can be p-value for differential expression or correlation with a phenotype

onts = c( "MF", "BP", "CC" )
tab <- as.list(onts)
names(tab) <- onts

for(i in 1:3){

GOdata <- new("topGOdata",
                description = "GOanalysis",
                ontology = onts[i],
                allGenes = geneList,
                annot = annFUN.org,
                mapping = "org.Mm.eg.db",
                ID = "SYMBOL",
                nodeSize = 20) # nodeSize: prune the GO hierarchy less than nodeSize

res.result <- runTest(GOdata, statistic = stats.use, algorithm = algorithm.use)

tab[[i]] <- GenTable(GOdata, pval = res.result, topNodes = 20, numChar = 100)

}

topGOResults <- plyr::rbind.fill(tab)
topGOResults.df <- as.data.frame(topGOResults)

# Visualization
# showSigOfNodes(GOdata, score(res.result), firstSigNodes = 5, useInfo ='all')

return(topGOResults.df)

}

scRNA_topGO_plot <- function(topGOResults.df, color1 = "chocolate2"){
# "dodgerblue3" for downregulated 

topGOResults.df$pval <- as.numeric(topGOResults.df$pval)
topGOResults.df$PValue <- -log10(topGOResults.df$pval)
topGOResults.df$Term <- factor(topGOResults.df$Term, levels = rev(unique(topGOResults.df$Term)))

ggplot(topGOResults.df, aes(x = Term, y = PValue)) + geom_bar(stat="identity", color = color1, fill = color1) + coord_flip() +
  ylab("-log10(p-value)") +
    theme(axis.title.y = element_blank(), axis.text = element_text(color = "black", size = 12), 
        axis.title.x = element_text(size = 12, color = "black"), title = element_text(size = 12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        axis.ticks.y = element_blank())

}



