
Remove_duplicated <- function(mouse_genes){
  
  

######################
# 1. unique genes
######################

# remove ENSEMBL duplicated rows
mouse_genes_unique.df <- mouse_genes[!mouse_genes$ensembl_gene_id_version %in% mouse_genes$ensembl_gene_id_version[duplicated(mouse_genes$ensembl_gene_id_version)], ]
# remove mgi_symbol duplicated rows
mouse_genes_unique.df <- mouse_genes_unique.df[!mouse_genes_unique.df$mgi_symbol %in% mouse_genes_unique.df$mgi_symbol[duplicated(mouse_genes_unique.df$mgi_symbol)], ]


########################
# 2. ENSEMBL duplicated
########################
# In here, I only leave the first row ENSEMBL genes and remove rest of duplicated ENSEMBL genes

ENSEMBL_duplicated_genes <- mouse_genes$ensembl_gene_id_version[duplicated(mouse_genes$ensembl_gene_id_version)]

ENSEMBL_duplicated.df <- mouse_genes[mouse_genes$ensembl_gene_id_version %in% ENSEMBL_duplicated_genes, ]

colnames(ENSEMBL_duplicated.df)

ENSEMBL_dup_removed.df <- data.frame(ensembl_gene_id_version = as.character(), mgi_symbol = as.character(), stringsAsFactors = FALSE) 

for(i in length(ENSEMBL_duplicated_genes)){
  ENSEMBL_temp.df <- ENSEMBL_duplicated.df[ENSEMBL_duplicated.df$ensembl_gene_id_version %in% ENSEMBL_duplicated_genes[i], ][1, ]
  ENSEMBL_dup_removed.df <- rbind(ENSEMBL_dup_removed.df, ENSEMBL_temp.df)
}


##########################
# 3. mgi symbol duplicated
##########################

mgi_symbol_duplicated_genes <- mouse_genes$mgi_symbol[duplicated(mouse_genes$mgi_symbol)]

mgi_symbol_duplicated.df <- mouse_genes[mouse_genes$mgi_symbol %in% mgi_symbol_duplicated_genes, ]

# add ENSEMBL to mgi_symbol as suffix
mgi_symbol_duplicated.df$mgi_symbol <- ifelse(mgi_symbol_duplicated.df$mgi_symbol == "", mgi_symbol_duplicated.df$mgi_symbol, paste0(mgi_symbol_duplicated.df$mgi_symbol, "_",mgi_symbol_duplicated.df$ensembl_gene_id_version))  
                                              
# change to NA if NA 
# this step is unnecessary
#mgi_symbol_duplicated.df$mgi_symbol <- ifelse(is.na(mgi_symbol_duplicated.df), mgi_symbol_duplicated.df$ensembl_gene_id_version, mgi_symbol_duplicated.df$mgi_symbol)

# change to ENSEMBL if ""
mgi_symbol_duplicated.df$mgi_symbol <- ifelse(mgi_symbol_duplicated.df$mgi_symbol == "", mgi_symbol_duplicated.df$ensembl_gene_id_version, mgi_symbol_duplicated.df$mgi_symbol)


mouse_genes_final.df <- rbind(mouse_genes_unique.df,ENSEMBL_dup_removed.df, mgi_symbol_duplicated.df)

return(mouse_genes_final.df)

}


Remove_duplicated_no_id_version <- function(mouse_genes){
  
  
  
  ######################
  # 1. unique genes
  ######################
  
  # remove ENSEMBL duplicated rows
  mouse_genes_unique.df <- mouse_genes[!mouse_genes$ensembl_gene_id %in% mouse_genes$ensembl_gene_id[duplicated(mouse_genes$ensembl_gene_id)], ]
  # remove mgi_symbol duplicated rows
  mouse_genes_unique.df <- mouse_genes_unique.df[!mouse_genes_unique.df$mgi_symbol %in% mouse_genes_unique.df$mgi_symbol[duplicated(mouse_genes_unique.df$mgi_symbol)], ]
  
  
  ########################
  # 2. ENSEMBL duplicated
  ########################
  # In here, I only leave the first row ENSEMBL genes and remove rest of duplicated ENSEMBL genes
  
  ENSEMBL_duplicated_genes <- mouse_genes$ensembl_gene_id[duplicated(mouse_genes$ensembl_gene_id)]
  
  ENSEMBL_duplicated.df <- mouse_genes[mouse_genes$ensembl_gene_id %in% ENSEMBL_duplicated_genes, ]
  
  colnames(ENSEMBL_duplicated.df)
  
  ENSEMBL_dup_removed.df <- data.frame(ensembl_gene_id = as.character(), mgi_symbol = as.character(), stringsAsFactors = FALSE) 
  
  for(i in length(ENSEMBL_duplicated_genes)){
    ENSEMBL_temp.df <- ENSEMBL_duplicated.df[ENSEMBL_duplicated.df$ensembl_gene_id %in% ENSEMBL_duplicated_genes[i], ][1, ]
    ENSEMBL_dup_removed.df <- rbind(ENSEMBL_dup_removed.df, ENSEMBL_temp.df)
  }
  
  
  ##########################
  # 3. mgi symbol duplicated
  ##########################
  
  mgi_symbol_duplicated_genes <- mouse_genes$mgi_symbol[duplicated(mouse_genes$mgi_symbol)]
  
  mgi_symbol_duplicated.df <- mouse_genes[mouse_genes$mgi_symbol %in% mgi_symbol_duplicated_genes, ]
  
  # add ENSEMBL to mgi_symbol as suffix
  mgi_symbol_duplicated.df$mgi_symbol <- ifelse(mgi_symbol_duplicated.df$mgi_symbol == "", mgi_symbol_duplicated.df$mgi_symbol, paste0(mgi_symbol_duplicated.df$mgi_symbol, "_",mgi_symbol_duplicated.df$ensembl_gene_id))  
  
  # change to NA if NA 
  # this step is unnecessary
  #mgi_symbol_duplicated.df$mgi_symbol <- ifelse(is.na(mgi_symbol_duplicated.df), mgi_symbol_duplicated.df$ensembl_gene_id_version, mgi_symbol_duplicated.df$mgi_symbol)
  
  # change to ENSEMBL if ""
  mgi_symbol_duplicated.df$mgi_symbol <- ifelse(mgi_symbol_duplicated.df$mgi_symbol == "", mgi_symbol_duplicated.df$ensembl_gene, mgi_symbol_duplicated.df$mgi_symbol)
  
  
  mouse_genes_final.df <- rbind(mouse_genes_unique.df,ENSEMBL_dup_removed.df, mgi_symbol_duplicated.df)
  
  return(mouse_genes_final.df)
  
}



