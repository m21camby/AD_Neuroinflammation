#! /bin/env RScript
# written by SJK at 2. Mar. 2020

.libPaths(c("/data/murphy/shared_libs/R", "/usr/local/lib/R/site-librar","/usr/local/lib/R/library","/usr/lib64/R/library","/usr/share/R/library", .libPaths() ))

library(Seurat, lib.loc = "/data/rajewsky/home/skim/R/")
library(ggplot2)
library(dplyr)
library(gridExtra)
source("/data/rajewsky/home/skim/Microglia_Heppner/Codes/DroNc_AD/src/report.sessionInfo.R")

load(file = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200221_figure1_for_meeting_report_Seurat_object.Robj")


data9set_cleaned.SO$cell_type <- data9set_cleaned.SO@active.ident

plot.data <- data9set_cleaned.SO@meta.data %>%
    dplyr::select(sample, gemgroup = gemgroup, cell_type = cell_type) %>%
    mutate(cell_type = cell_type) %>%
    group_by(gemgroup, cell_type) %>%
  summarise(count = n()) %>%
    mutate(sample_total = sum(count)) %>%
    mutate(cell_prop = count / sample_total) 

#%>%group_by(sample) %>%
#    mutate(dataset_total = sum(count)) %>%
#    ungroup() %>%
#    mutate(dataset_prop = count / dataset_total)
#sum(plot.data$dataset_prop)

#plot.data$gemgroup <- factor(plot.data$gemgroup, levels = c("1", "4", "7", "3", "6", "9", "2", "5", "8"))

plot.data$gemgroup <- factor(plot.data$gemgroup, levels = c("1", "4", "7", "2", "5", "8", "3", "6", "9"))
plot.data <- as.data.frame(plot.data)

p1 <- data.frame(expgroup = c(rep("Ctrl-1", 11), rep("ADp40KO-1", 10), rep("AD-1", 10),
                              rep("Ctrl-2", 11), rep("ADp40KO-2", 11), rep("AD-2", 11),
                              rep("Ctrl-3", 11), rep("ADp40KO-3", 11), rep("AD-3", 11)), stringsAsFactors = FALSE)


plot.data <- cbind(plot.data, p1)

plot.data$expgroup <- factor(plot.data$expgroup, levels = c("Ctrl-1", "Ctrl-2", "Ctrl-3", "AD-1", "AD-2", "AD-3", "ADp40KO-1", "ADp40KO-2", "ADp40KO-3"))

g1 <- ggplot(plot.data, aes(x = expgroup, y = cell_prop, fill = cell_type)) +
    geom_col() + theme_classic() + scale_y_continuous(expand = c(0,0)) + ylab("percentage") + 
  theme(axis.text.x = element_text(size =14, angle = 90, color = "black"), axis.title.x = element_blank(), axis.text.y = element_text(size =14, color = "black"),
        axis.title.y = element_text(size =14), legend.title = element_blank(), legend.text = element_text(size = 12)) +  
  scale_fill_manual(values=c("#CC0000", "#FF6600", "#FF9900","#CC9900", "#FFCC66", "#99CC00","#003300", "#99CCCC", "#009999", "#3399FF","#000033"))




ggsave(filename = "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Cell_Percentage.png", 
       plot = g1, 
  scale = 1, width = 8, height = 7, units = "in", device = "png",
  dpi = 300)


##################
# sessionInfo save
##################
doc = report.sessionInfo()
write.table(doc, "/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/Paper_Fig/Fig1_Cell_Percentage_session_info.txt")


