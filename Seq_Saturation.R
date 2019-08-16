

# Written by SJK 15. 08. 2019

# In here we used non-linear regression or Loess to predict number of genes and UMIs per depth and predict if we do deep-sequencing
# pcr is predict based on linear modeling
# In every graphs, black dots are real data and the red line is prediction based on analysis.

# In here requires two path
# path1 is for subsampling directory and path2 is for output of summary file of DGE matrix by Drop-seq pipeline

library(ggplot2)

############################################
# 1. median number of genes per subsampling
############################################
DF_Genes <- function(exp, path1, path2){
  output <- matrix(ncol=2, nrow=10) # empty matrix
  for(i in c(1:9)){ # for loop for reshaping matrix from (n x 2) to (4 x n_samples)  
  DF <- read.table(paste0(path1, "out_gene_exon_tagged.dge_", exp, "_", i, ".summary.txt"), header = TRUE)
  output[i, 1] <- as.integer(sum(DF$NUM_GENIC_READS))
  output[i, 2] <- as.integer(median(DF$NUM_GENES))
  }
  output_10 <- read.table(path2, header = TRUE)
  output[10, 1] <- as.integer(sum(output_10$NUM_GENIC_READS))
  output[10, 2] <- as.integer(median(output_10$NUM_GENES))
  output.df <- as.data.frame(output)
  colnames(output.df) <- c("Read", "Genes")
  return(output.df)
}

############################################
# 2. median number of UMIs per subsampling
############################################
DF_UMIs <- function(exp, path1, path2){
  output <- matrix(ncol=2, nrow=10) # empty matrix
  for(i in c(1:9)){ # for loop for reshaping matrix from (n x 2) to (4 x n_samples)  
  DF <- read.table(paste0(path1, "out_gene_exon_tagged.dge_", exp, "_", i, ".summary.txt"), header = TRUE)
  output[i, 1] <- as.integer(sum(DF$NUM_GENIC_READS))
  output[i, 2] <- as.integer(median(DF$NUM_TRANSCRIPTS))
  }
  output_10 <- read.table(path2, header = TRUE)
  output[10, 1] <- as.integer(sum(output_10$NUM_GENIC_READS))
  output[10, 2] <- as.integer(median(output_10$NUM_TRANSCRIPTS))
  output.df <- as.data.frame(output)
  colnames(output.df) <- c("Read", "UMIs")
  return(output.df)
}

###########################################
# 3. median number of pcr per subsampling
###########################################
pcr_cal <- function(exp, path1, path2){
  output <- matrix(ncol=2, nrow=10) # empty matrix
  for(i in c(1:9)){ # for loop for reshaping matrix from (n x 2) to (4 x n_samples)  
  DF <- read.table(paste0(path1, "out_gene_exon_tagged.dge_", exp, "_", i, ".summary.txt"), header = TRUE)
  DF$pcr <- DF$NUM_GENIC_READS / DF$NUM_TRANSCRIPTS
  output[i, 1] <- as.integer(sum(DF$NUM_GENIC_READS))
  output[i, 2] <- median(DF$pcr, na.rm = TRUE)
  }
  output_10 <- read.table(path2, header = TRUE)
  output_10$pcr <- output_10$NUM_GENIC_READS / output_10$NUM_TRANSCRIPTS
  output[10, 1] <- as.integer(sum(output_10$NUM_GENIC_READS))
  output[10, 2] <- median(output_10$pcr, na.rm = TRUE)
  output.df <- as.data.frame(output)
  colnames(output.df) <- c("Read", "pcr")
  return(output.df)
}

############################################
# 4. predict based on non-linear regression
############################################
# a and b is parameters of non-linear regression calculated by nls
# c is starting prediction reads and d is ends prediction reads
Predict_nlr <- function(a,b,c,d){
  output_2 <- matrix(ncol=2, nrow=d) # empty matrix
  for(i in c(1:d)){
    k <- i+c
  mm=function(a,b,k) {a*k/(b+k)}
  output_2[i, 1] <- as.integer(mm(a, b, i+c))
  output_2[i, 2] <- as.integer(i+c)
  }
  output_2.df <- as.data.frame(output_2)
  return(output_2.df)
}

###############################################
# 5. predict based on LOESS
###############################################
# smoothes is after loess prediction object
# a is ends prediction reads - starting prediction reads and b is starting prediction reads
Predict_loess <- function(smoothed, a,b){
  output_2 <- matrix(ncol=2, nrow=a) # empty matrix
  for (i in c(1:a)){
  output_2[i, 1] <- as.integer(smoothed$fit[i])
  output_2[i, 2] <- as.integer(i+b)
}
  output_2 <- as.data.frame(output_2)
  return(output_2)
}

####################
# 6. prediction pcr
#################### 
predict_pcr <- function(predict, a, b){
  output_2 <- matrix(ncol=2, nrow=b) # empty matrix
  for (i in c(1:b)){
  #print(predict[i])
  output_2[i, 1] <- predict[i]
  output_2[i, 2] <- i+a
}
  output_2 <- as.data.frame(output_2)
  return(output_2)
}

###############
# 7. plotting
###############

predict_plot_Genes <- function(input_DF, input_subDF, lower_limit = 100, upper_limit = 120){


ggplot(input_DF, aes(x = Reads, y = Genes)) + geom_point() + 
geom_point(data = input_subDF, color = "red", size = 1) + 
theme(axis.title = element_text(face = "bold", size = 12), axis.text = element_text(face = "bold", size = 12, color = "black"), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
xlab("Number of reads (million)") + 
geom_hline(yintercept = lower_limit, linetype="dashed", color = "red") + 
geom_hline(yintercept = upper_limit, linetype="dashed", color = "red")
}


predict_plot_UMIs <- function(input_DF, input_subDF, lower_limit = 100, upper_limit = 120){

ggplot(input_DF, aes(x = Reads, y = UMIs)) + geom_point() +
geom_point(data = input_subDF, color = "red", size = 1) +
theme(axis.title = element_text(face = "bold", size = 12), axis.text = element_text(face = "bold", size = 12, color = "black"), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
xlab("Number of reads (million)") +
geom_hline(yintercept = lower_limit, linetype="dashed", color = "red") +
geom_hline(yintercept = upper_limit, linetype="dashed", color = "red")
}

predict_plot_pcr <- function(input_DF, input_subDF, lower_limit = 10, upper_limit = 20){

ggplot(input_DF, aes(x = Reads, y = pcr)) + geom_point() +
geom_point(data = input_subDF, color = "red", size = 1) +
theme(axis.title = element_text(face = "bold", size = 12), axis.text = element_text(face = "bold", size = 12, color = "black"), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
xlab("Number of reads (million)") +
geom_hline(yintercept = lower_limit, linetype="dashed", color = "red") +
geom_hline(yintercept = upper_limit, linetype="dashed", color = "red")
}









