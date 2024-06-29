
load(file="/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/R_Scripts/20200210_figures_for_meeting_report_cleaned_Seurat_object_res0.8.Robj")


data9set.meta.df <- data9set_cleaned.SO@meta.data %>% as.data.frame
data9set.meta.df$CB <- rownames(data9set.meta.df)
data9set.meta.df$cell_type <- data9set_cleaned.SO@active.ident

data9set.meta.df[data9set.meta.df$gemgroup %in% 4, ]


EX.meta.df <- data9set.meta.df[data9set.meta.df$cell_type %in% c("CA1 Neurons", "CA2/CA3 Neuron", "Dentate Gyrus","Neurons", "Subiculum"), ]
IN.meta.df <- data9set.meta.df[data9set.meta.df$cell_type %in% c("Inhibitory Interneurons"), ]
MG.meta.df <- data9set.meta.df[data9set.meta.df$cell_type %in% c("Microglia"), ]
OL.meta.df <- data9set.meta.df[data9set.meta.df$cell_type %in% c("Oligo"), ]
OPC.meta.df <- data9set.meta.df[data9set.meta.df$cell_type %in% c("OPC"), ]
AS.meta.df <- data9set.meta.df[data9set.meta.df$cell_type %in% c("Astrocytes"), ]



Cell_Barcode_Extract <- function(cell_type = "EX", meta.df){
  print(cell_type) 
  meta_1.df <- meta.df[meta.df$gemgroup %in% 1, ]
  write.table(meta_1.df$CB, file = paste0("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/BigWig_Files/split_BigWig_file_by_Cell_type/Cell_barcodes_files/", cell_type, "_1.txt"), quote = FALSE)
  
  for(i in c(2:9)){
    meta_1.df <- meta.df[meta.df$gemgroup %in% i, ]
    meta_1.df$CB <- paste0(substr(meta_1.df$CB, 1, 17), 1)
    write.table(meta_1.df$CB, file = paste0("/data/rajewsky/home/skim/Microglia_Heppner/201912_10X_9sets/BigWig_Files/split_BigWig_file_by_Cell_type/Cell_barcodes_files/", cell_type, "_",i,".txt"), quote = FALSE)
    
  }
  
}


Cell_Barcode_Extract(cell_type = "EX", EX.meta.df)
Cell_Barcode_Extract(cell_type = "IN", IN.meta.df)
Cell_Barcode_Extract(cell_type = "MG", MG.meta.df)
Cell_Barcode_Extract(cell_type = "OL", OL.meta.df)
Cell_Barcode_Extract(cell_type = "OPC", OPC.meta.df)
Cell_Barcode_Extract(cell_type = "AS", AS.meta.df)


