#! /bin/env RScript
# written by SJK at 2. May. 2020

# blue brain projet pie chart
# Neuron: 68 % (red or orange)
# Microglia 12 %
# Astrocytes 7.5%
# Oligodendrocytes 12 %
# rest 0.5 %

#.libPaths(c("/home/skim/R/usr_lib", .libPaths()))


.libPaths(c(
            "/data/rajewsky/home/skim/R/usr_lib_Seurat", .libPaths()))


.libPaths(c("/data/rajewsky/home/skim/R/usr_lib_v4", 
            "/data/rajewsky/home/skim/R/usr_lib",
            "/data/rajewsky/home/skim/R/usr_lib_Seurat", .libPaths()))


library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")



plot.data <- data.frame(cell_type = c("Neurons", "Astrocytes", "Microglia", "Oligodendrocytes", "rest"),
                        percent = c(68, 7.5, 11.99, 12.01, 0.4))

g1 <- ggplot(plot.data, aes(x = 2, y= percent, fill = factor(percent, levels = c(68, 7.5, 11.99, 12.01, 0.4)))) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y", start=0) +
  theme_void() + theme(legend.text = element_text(size = 12), legend.title = element_blank()) +
  scale_fill_manual(values=c("#CC0000", "#FFCC66", "#99CC00", "#99CCCC","#000033"),
                    labels = c("Neurons: 68%", "Astrocytes: 7.5%", "Microglia: 12%",
                               "Oligodendrocytes: 12%", "rest: 0.5%"))



ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_blue_brain_Cell_Percentage_piechart.png", 
       plot = g1, 
  scale = 1, width = 8, height = 6, units = "in", device = "png",
  dpi = 300)


##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_blue_brain_Cell_Percentage_piechart_session_info.txt")


