install.packages('BiocManager', lib = "/data/rajewsky/home/skim/R/usr_lib")
BiocManager::install('RANN', lib = "/data/rajewsky/home/skim/R/usr_lib")
install.packages('RColorBrewer', lib = "/data/rajewsky/home/skim/R/usr_lib")
install.packages('ROCR', lib = "/data/rajewsky/home/skim/R/usr_lib")
install.packages("remotes", lib = "/data/rajewsky/home/skim/R/usr_lib")

remotes::install_version("Seurat", version = "4.0.4", lib = "/data/rajewsky/home/skim/R/usr_lib")

remove.packages("Seurat",  lib = "/data/rajewsky/home/skim/R/usr_lib")
  
BiocManager::install('rsvd', version = "1.0.2", lib = "/data/rajewsky/home/skim/R/usr_lib")
# FNN, 
install.packages("rsvd", version = "1.0.0", lib = "/data/rajewsky/home/skim/R/usr_lib")

