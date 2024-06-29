#! /bin/env RScript
# written by SJK at 3. Mar. 2020
# modified at 04. Nov. 2020
# don't use libPaths

# .libPaths(c("/data/murphy/shared_libs/R", "/usr/local/lib/R/site-librar","/usr/local/lib/R/library","/usr/lib64/R/library","/usr/share/R/library", .libPaths() ))

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

load(file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200221_figure1_for_meeting_report_Seurat_object.Robj")

data9set_cleaned.SO$cell_type <- data9set_cleaned.SO@active.ident

plot.data2 <- data9set_cleaned.SO@meta.data %>%
    dplyr::select(sample, cell_type, cluster = seurat_clusters) %>%
    mutate(cluster = cluster) %>%
    group_by(cell_type, sample) %>% 
    summarise(count = n()) %>% mutate(type_total = sum(count)) 

plot.data2 <- plot.data2 %>%
  mutate(total = sum(plot.data2$type_total)/3) %>% 
  mutate(sample_prop = count / total) %>% 
  mutate(type_prop = type_total /total) %>% 
  mutate(sample_type_prop = count /type_total)

plot.data2 <- as.data.frame(plot.data2)
plot.data2$sample <- factor(plot.data2$sample, levels = c("Ctrl", "AD", "ADp40KO"))
plot.data2$type_prop2 <- round(100*plot.data2$type_prop,2)


g1 <- ggplot(plot.data2, aes(x = 2, y= type_prop2, fill = factor(type_prop2, levels = c(52.98, 7.00, 0.58, 0.69, 4.45, 8.33, 0.30, 20.85, 3.11, 0.89, 0.83)))) + 
  geom_bar(width = 1, stat = "identity") + 
  coord_polar(theta = "y", start=0) + 
  theme_void() + theme(legend.text = element_text(size = 12, color = "black", family = "helvetica"), 
                       legend.title = element_blank()) + 
  scale_fill_manual(values=c("#CC0000", "#FF6600", "#FF9900","#CC9900", "#FFCC66", "#99CC00","#003300", "#99CCCC", "#009999", "#3399FF","#000033"), 
                    labels = c("Excitatory Neurons: 52.98%", "Inhibitory Interneurons: 7.00%", "Cajal Retzius: 0.58%",
                               "Choroid Plexus: 0.69%", "Astrocytes: 4.45%", "Microglia: 8.33%",
                               "Macrophages: 0.30%", "Oligodendrocytes: 20.85%", "OPC: 3.11%",
                               "Fibroblast: 0.89%", "Vascular cells: 0.83%"))
# +  geom_text(data =  plot.data2, aes(x = 2.1, y = 220, label = "Ex Neurons"), size = 4) + 
#  geom_text(data =  plot.data2, aes(x = 2.1, y = 130, label = "In Neurons"), size = 4) + 
#  geom_text(data =  plot.data2, aes(x = 2.2, y = 110, label = "AS"), size = 4) +
#  geom_text(data =  plot.data2, aes(x = 2.2, y = 90, label = "MG"), size = 4) + 
#  geom_text(data =  plot.data2, aes(x = 2.2, y = 50, label = "Oligo"), size = 4) + 
#  geom_text(data =  plot.data2, aes(x = 2.2, y = 10, label = "OPC"), size = 4)  

g1

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_whole_Cell_Percentage_piechart.png", 
       plot = g1, 
  scale = 1, width = 8, height = 6, units = "in", device = "png",
  dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_whole_Cell_Percentage_piechart.pdf", 
       plot = g1, 
       scale = 1, width = 8, height = 6, units = "in", device = cairo_pdf,
       dpi = 300)

ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_whole_Cell_Percentage_piechart.svg", 
       plot = g1, 
       scale = 1, width = 8, height = 6, units = "in", device = "svg",
       dpi = 300)

##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_whole_Cell_Percentage_piechart_session_info.txt")


