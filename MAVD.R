# A Simulation Study of Seven Data Preprocessing Methods

#Set Working Directory
setwd("/public/linzhf")

# Load necessary libraries
library(compositions)  
library(edgeR)
library(DESeq2)
library(GUniFrac)
library(Wrench)
library(GUniFrac)
library(metagenomeSeq)

# Load custom MAVD functions
source("mavd_source.R")

############ Data Processing
#### TSS ####
res_MAVD_TSS <- function(sim=sim, res.cutoff=0.95) {
  group <- sim$group  # Extract the group information from the simulation object
  data <- sim$data    # Extract the data matrix from the simulation object
  normal_micro <- data[, group == "0"]  # Select normal microbial data
  disease_micro <- data[, group == "1"]  # Select disease microbial data
  
  # Write normal and disease microbial data to PCL files
  write.table(normal_micro, file = "/public/linzhf/TSS/normal_micro.pcl", append=FALSE, quote=FALSE, sep="\t", eol="\n", na="", row.names=TRUE, col.names=TRUE)
  write.table(disease_micro, file = "/public/linzhf/TSS/disease_micro.pcl", append=FALSE, quote=FALSE, sep="\t", eol="\n", na="", row.names=TRUE, col.names=TRUE)
  
  # Read the PCL files back into R
  TUMOR.pcl <- read.table("/public/linzhf/TSS/disease_micro.pcl", header = TRUE, sep = "\t", quote = "", fill = TRUE, comment.char = "")
  NORMAL.pcl <- read.table("/public/linzhf/TSS/normal_micro.pcl", header = TRUE, sep = "\t", quote = "", fill = TRUE, comment.char = "")
  
  # Data processing (TSS normalization)
  NORMAL.pcl <- apply(NORMAL.pcl, 2, function(x) x / sum(x))  # Normalize normal data
  TUMOR.pcl <- apply(TUMOR.pcl, 2, function(x) x / sum(x))    # Normalize tumor data
  
  # Apply MAVD Part 1 function to the normalized data
  mavd_part1 <- mavd_part1(normal.pcl = NORMAL.pcl, tumor.pcl = TUMOR.pcl)
  
  # Plot Wold's T² statistics for the normal data
  plot.wold.mavd(mavd_part1$wold.Nmat, main.extra = "Normal")
  mavd_part1$wold.Nmat  # Display Wold's T² statistics
  
  # Identify the peak index in the Wold's plot
  peak_index_l <- which(diff(mavd_part1$wold.Nmat) > 0)[1] 
  peak_index_2 <- which(diff(mavd_part1$wold.Nmat[-c(1:peak_index_l)]) < 0)[1] 
  
  K <- peak_index_l + peak_index_2  # User chooses dimension K where Wold's plot peaks (has a local maximum)
  abline(v = K, lty = "dotted", col = "red")  # Draw a vertical line at dimension K
  legend(x = K, y = 0.2 * mavd_part1$wold.Nmat[[1]], legend = paste("dim =", K), bty = "n")  # Add legend to the plot
  
  #### DATA TRANSFORMATION - 
  #### DATA IS DECOMPOSED INTO DISEASE & NORMAL COMPONENTS 
  mavd_part2 <- mavd_part2(x = mavd_part1, k = K)  # Apply MAVD Part 2 function
  D <- mavd_part2$Dc.Dmat  # Extract disease component
  N <- mavd_part2$L1.mat   # Extract normal component
  
  # Filter the components based on mean absolute difference
  diff <- data.frame(D)[apply(abs(D), 1, mean) > apply(abs(N), 1, quantile, probs=res.cutoff), ]
  gold_names_new <- rownames(diff)  # Get the names of significant features
  gold_names_new  # Return the names of significant features
  
  # return(list(Dc.Dmat=D, L1.mat=N, gold_names=gold_names_new))  # Optionally return additional data
  return(list(gold_names=gold_names_new))  # Return only the significant feature names
}

#### TMM ####
res_MAVD_TMM <- function(sim=sim, res.cutoff=0.95) {
  group <- sim$group  # Extract the group information from the simulation object
  data <- sim$data    # Extract the data matrix from the simulation object
  normal_micro <- data[, group == "0"]  # Select normal microbiome data
  disease_micro <- data[, group == "1"]  # Select disease microbiome data
  
  # Write normal and disease microbiome data to PCL files
  write.table(normal_micro, file = "/public/linzhf/TMM/1normal_micro.pcl", append=FALSE, quote=FALSE, sep="\t", eol="\n", na="", row.names=TRUE, col.names=TRUE)
  write.table(disease_micro, file = "/public/linzhf/TMM/1disease_micro.pcl", append=FALSE, quote=FALSE, sep="\t", eol="\n", na="", row.names=TRUE, col.names=TRUE)
  
  # Read the PCL files back into R
  TUMOR.pcl <- read.table("/public/linzhf/TMM/1disease_micro.pcl", header = TRUE, sep = "\t", quote = "", fill = TRUE, comment.char = "")
  NORMAL.pcl <- read.table("/public/linzhf/TMM/1normal_micro.pcl", header = TRUE, sep = "\t", quote = "", fill = TRUE, comment.char = "")
  
  # TMM normalization
  # Create DGEList objects for normal and tumor data
  NORMAL_DGEList <- DGEList(counts = NORMAL.pcl)
  TUMOR_DGEList <- DGEList(counts = TUMOR.pcl)
  
  # Perform TMM normalization
  NORMAL_TMM <- edgeR::calcNormFactors(NORMAL_DGEList, method = "TMM")
  TUMOR_TMM <- edgeR::calcNormFactors(TUMOR_DGEList, method = "TMM")
  
  # Get the TMM normalized data (counts per million)
  NORMAL.pcl <- cpm(NORMAL_TMM)
  TUMOR.pcl <- cpm(TUMOR_TMM)
  
  # Apply MAVD Part 1 function to the normalized data
  mavd_part1 <- mavd_part1(normal.pcl = NORMAL.pcl, tumor.pcl = TUMOR.pcl)
  
  # Plot Wold's T² statistics for the normal data
  plot.wold.mavd(mavd_part1$wold.Nmat, main.extra = "Normal")
  mavd_part1$wold.Nmat  # Display Wold's T² statistics
  
  # Identify the peak index in the Wold's plot
  peak_index_l <- which(diff(mavd_part1$wold.Nmat) > 0)[1] 
  peak_index_2 <- which(diff(mavd_part1$wold.Nmat[-c(1:peak_index_l)]) < 0)[1] 
  
  K <- peak_index_l + peak_index_2  # User chooses dimension K where Wold's plot peaks (has a local maximum)
  abline(v = K, lty = "dotted", col = "red")  # Draw a vertical line at dimension K
  legend(x = K, y = 0.2 * mavd_part1$wold.Nmat[[1]], legend = paste("dim =", K), bty = "n")  # Add legend to the plot
  
  #### DATA TRANSFORMATION - 
  #### DATA IS DECOMPOSED INTO DISEASE & NORMAL COMPONENTS 
  mavd_part2 <- mavd_part2(x = mavd_part1, k = K)  # Apply MAVD Part 2 function
  
  D <- mavd_part2$Dc.Dmat  # Extract disease component
  N <- mavd_part2$L1.mat   # Extract normal component
  
  # Filter the components based on mean absolute difference
  diff <- data.frame(D)[apply(abs(D), 1, mean) > apply(abs(N), 1, quantile, probs=res.cutoff), ]
  gold_names_new <- rownames(diff)  # Get the names of significant features
  gold_names_new  # Return the names of significant features
  
  # return(list(Dc.Dmat=D, L1.mat=N, gold_names=gold_names_new))  # Optionally return additional data
  return(list(gold_names=gold_names_new))  # Return only the significant feature names
}

#### RLE ####
res_MAVD_RLE <- function(sim=sim, res.cutoff=0.95) {
  group <- sim$group  # Extract the group information from the simulation object
  data <- sim$data    # Extract the data matrix from the simulation object
  normal_micro <- data[, group == "0"]  # Select normal microbiome data
  disease_micro <- data[, group == "1"]  # Select disease microbiome data
  
  # Write normal and disease microbiome data to PCL files
  write.table(normal_micro, file = "/public/linzhf/RLE/1normal_micro.pcl", append=FALSE, quote=FALSE, sep="\t", eol="\n", na="", row.names=TRUE, col.names=TRUE)
  write.table(disease_micro, file = "/public/linzhf/RLE/1disease_micro.pcl", append=FALSE, quote=FALSE, sep="\t", eol="\n", na="", row.names=TRUE, col.names=TRUE)
  
  # Read the PCL files back into R
  TUMOR.pcl <- read.table("/public/linzhf/RLE/1disease_micro.pcl", header = TRUE, sep = "\t", quote = "", fill = TRUE, comment.char = "")
  NORMAL.pcl <- read.table("/public/linzhf/RLE/1normal_micro.pcl", header = TRUE, sep = "\t", quote = "", fill = TRUE, comment.char = "")
  
  # Create DGEList objects for normal and tumor data
  NORMAL_DGEList <- DGEList(counts = NORMAL.pcl)
  TUMOR_DGEList <- DGEList(counts = TUMOR.pcl)
  
  # Perform RLE normalization
  NORMAL_RLE <- edgeR::calcNormFactors(NORMAL_DGEList, method = "RLE")
  TUMOR_RLE <- edgeR::calcNormFactors(TUMOR_DGEList, method = "RLE")
  
  # Get the RLE normalized data (counts per million)
  NORMAL.pcl <- cpm(NORMAL_RLE)
  TUMOR.pcl <- cpm(TUMOR_RLE)
  
  # Apply MAVD Part 1 function to the normalized data
  mavd_part1 <- mavd_part1(normal.pcl = NORMAL.pcl, tumor.pcl = TUMOR.pcl)
  
  # Plot Wold's T² statistics for the normal data
  plot.wold.mavd(mavd_part1$wold.Nmat, main.extra = "Normal")
  mavd_part1$wold.Nmat  # Display Wold's T² statistics
  
  # Identify the peak index in the Wold's plot
  peak_index_l <- which(diff(mavd_part1$wold.Nmat) > 0)[1] 
  peak_index_2 <- which(diff(mavd_part1$wold.Nmat[-c(1:peak_index_l)]) < 0)[1] 
  
  K <- peak_index_l + peak_index_2  # User chooses dimension K where Wold's plot peaks (has a local maximum)
  abline(v = K, lty = "dotted", col = "red")  # Draw a vertical line at dimension K
  legend(x = K, y = 0.2 * mavd_part1$wold.Nmat[[1]], legend = paste("dim =", K), bty = "n")  # Add legend to the plot
  
  #### DATA TRANSFORMATION - 
  #### DATA IS DECOMPOSED INTO DISEASE & NORMAL COMPONENTS 
  mavd_part2 <- mavd_part2(x = mavd_part1, k = K)  # Apply MAVD Part 2 function
  
  D <- mavd_part2$Dc.Dmat  # Extract disease component
  N <- mavd_part2$L1.mat   # Extract normal component
  
  # Filter the components based on mean absolute difference
  diff <- data.frame(D)[apply(abs(D), 1, mean) > apply(abs(N), 1, quantile, probs=res.cutoff), ]
  gold_names_new <- rownames(diff)  # Get the names of significant features
  gold_names_new  # Return the names of significant features
  
  # return(list(Dc.Dmat=D, L1.mat=N, gold_names=gold_names_new))  # Optionally return additional data
  return(list(gold_names=gold_names_new))  # Return only the significant feature names
}

#### GMPR ####
res_MAVD_GMPR <- function(sim=sim, res.cutoff=0.95) {
  group <- sim$group  # Extract the group information from the simulation object
  data <- sim$data    # Extract the data matrix from the simulation object
  normal_micro <- data[, group == "0"]  # Select normal microbiome data
  disease_micro <- data[, group == "1"]  # Select disease microbiome data
  
  # Write normal and disease microbiome data to PCL files
  write.table(normal_micro, file = "/public/linzhf/GMPR/1normal_micro.pcl", append=FALSE, quote=FALSE, sep="\t", eol="\n", na="", row.names=TRUE, col.names=TRUE)
  write.table(disease_micro, file = "/public/linzhf/GMPR/1disease_micro.pcl", append=FALSE, quote=FALSE, sep="\t", eol="\n", na="", row.names=TRUE, col.names=TRUE)
  
  # Read the PCL files back into R
  TUMOR.pcl <- read.table("/public/linzhf/GMPR/1disease_micro.pcl", header = TRUE, sep = "\t", quote = "", fill = TRUE, comment.char = "")
  NORMAL.pcl <- read.table("/public/linzhf/GMPR/1normal_micro.pcl", header = TRUE, sep = "\t", quote = "", fill = TRUE, comment.char = "")
  
  # Normalize normal microbiome data using GMPR
  size.factor1 <- GMPR(NORMAL.pcl)  # Calculate GMPR size factors for normal data
  normalized_GMPR1 <- apply(t(NORMAL.pcl), 2, function(x) x / size.factor1)  # Normalize data by size factors
  NORMAL.pcl <- t(normalized_GMPR1)  # Transpose back to original orientation
  
  # Normalize tumor microbiome data using GMPR
  size.factor2 <- GMPR(TUMOR.pcl)  # Calculate GMPR size factors for tumor data
  normalized_GMPR2 <- apply(t(TUMOR.pcl), 2, function(x) x / size.factor2)  # Normalize data by size factors
  TUMOR.pcl <- t(normalized_GMPR2)  # Transpose back to original orientation
  
  # Apply MAVD Part 1 function to the normalized data
  mavd_part1 <- mavd_part1(normal.pcl = NORMAL.pcl, tumor.pcl = TUMOR.pcl)
  
  # Plot Wold's T² statistics for the normal data
  plot.wold.mavd(mavd_part1$wold.Nmat, main.extra = "Normal")
  mavd_part1$wold.Nmat  # Display Wold's T² statistics
  
  # Identify the peak index in the Wold's plot
  peak_index_l <- which(diff(mavd_part1$wold.Nmat) > 0)[1] 
  peak_index_2 <- which(diff(mavd_part1$wold.Nmat[-c(1:peak_index_l)]) < 0)[1] 
  
  K <- peak_index_l + peak_index_2  # User chooses dimension K where Wold's plot peaks (has a local maximum)
  abline(v = K, lty = "dotted", col = "red")  # Draw a vertical line at dimension K
  legend(x = K, y = 0.2 * mavd_part1$wold.Nmat[[1]], legend = paste("dim =", K), bty = "n")  # Add legend to the plot
  
  #### DATA TRANSFORMATION - 
  #### DATA IS DECOMPOSED INTO DISEASE & NORMAL COMPONENTS 
  mavd_part2 <- mavd_part2(x = mavd_part1, k = K)  # Apply MAVD Part 2 function
  
  D <- mavd_part2$Dc.Dmat  # Extract disease component
  N <- mavd_part2$L1.mat   # Extract normal component
  
  # Filter the components based on mean absolute difference
  diff <- data.frame(D)[apply(abs(D), 1, mean) > apply(abs(N), 1, quantile, probs=res.cutoff), ]
  gold_names_new <- rownames(diff)  # Get the names of significant features
  gold_names_new  # Return the names of significant features
  
  # return(list(Dc.Dmat=D, L1.mat=N, gold_names=gold_names_new))  # Optionally return additional data
  return(list(gold_names=gold_names_new))  # Return only the significant feature names
}

#### Wrench ####
res_MAVD_wrench <- function(sim=sim, res.cutoff=0.95) {
  group <- sim$group  # Extract the group information from the simulation object
  data <- sim$data    # Extract the data matrix from the simulation object
  
  # Convert the data to a matrix format for Wrench analysis
  counts <- as.matrix(data)  
  group1 <- as.vector(sim$group)  # Convert group information to vector format
  
  # Apply Wrench normalization
  W <- wrench(counts, condition=group1)  # Perform Wrench normalization
  normalizationFactors <- W$nf  # Extract normalization factors
  normalized_wrench <- apply(t(data), 2, function(x) x / normalizationFactors)  # Normalize data by factors
  data <- t(normalized_wrench)  # Transpose back to original orientation
  
  # Split the normalized data into normal and disease microbiome subsets
  normal_micro <- data[, group == "0"]  
  disease_micro <- data[, group == "1"]  
  
  # Write the normal and disease microbiome data to PCL files
  write.table(normal_micro, file = "/public/linzhf/wrench/1normal_micro.pcl", append=FALSE, quote=FALSE, sep="\t", eol="\n", na="", row.names=TRUE, col.names=TRUE)
  write.table(disease_micro, file = "/public/linzhf/wrench/1disease_micro.pcl", append=FALSE, quote=FALSE, sep="\t", eol="\n", na="", row.names=TRUE, col.names=TRUE)
  
  # Read the PCL files back into R
  TUMOR.pcl <- read.table("/public/linzhf/wrench/1disease_micro.pcl", header = TRUE, sep = "\t", quote = "", fill = TRUE, comment.char = "")
  NORMAL.pcl <- read.table("/public/linzhf/wrench/1normal_micro.pcl", header = TRUE, sep = "\t", quote = "", fill = TRUE, comment.char = "")
  
  # Apply MAVD Part 1 function to the normalized data
  mavd_part1 <- mavd_part1(normal.pcl = NORMAL.pcl, tumor.pcl = TUMOR.pcl)
  
  # Plot Wold's T² statistics for the normal data
  plot.wold.mavd(mavd_part1$wold.Nmat, main.extra = "Normal")
  mavd_part1$wold.Nmat  # Display Wold's T² statistics
  
  # Identify the peak index in the Wold's plot
  peak_index_l <- which(diff(mavd_part1$wold.Nmat) > 0)[1] 
  peak_index_2 <- which(diff(mavd_part1$wold.Nmat[-c(1:peak_index_l)]) < 0)[1] 
  
  K <- peak_index_l + peak_index_2  # User chooses dimension K where Wold's plot peaks (has a local maximum)
  abline(v = K, lty = "dotted", col = "red")  # Draw a vertical line at dimension K
  legend(x = K, y = 0.2 * mavd_part1$wold.Nmat[[1]], legend = paste("dim =", K), bty = "n")  # Add legend to the plot
  
  #### DATA TRANSFORMATION - 
  #### DATA IS DECOMPOSED INTO DISEASE & NORMAL COMPONENTS 
  mavd_part2 <- mavd_part2(x = mavd_part1, k = K)  # Apply MAVD Part 2 function
  
  D <- mavd_part2$Dc.Dmat  # Extract disease component
  N <- mavd_part2$L1.mat   # Extract normal component
  
  # Filter the components based on mean absolute difference
  diff <- data.frame(D)[apply(abs(D), 1, mean) > apply(abs(N), 1, quantile, probs=res.cutoff), ]
  gold_names_new <- rownames(diff)  # Get the names of significant features
  gold_names_new  # Return the names of significant features
  
  # return(list(Dc.Dmat=D, L1.mat=N, gold_names=gold_names_new))  # Optionally return additional data
  return(list(gold_names=gold_names_new))  # Return only the significant feature names
}

#### ALR ####
res_MAVD_ALR <- function(sim=sim, res.cutoff=0.95) {
  group <- sim$group  # Extract group information from the simulation object
  data <- sim$data    # Extract the data matrix from the simulation object
  
  # Split the data into normal and disease microbiome subsets
  normal_micro <- data[, group == "0"]
  disease_micro <- data[, group == "1"]
  
  # Write the normal and disease microbiome data to PCL files
  write.table(normal_micro, file = "/public/linzhf/ALR/1normal_micro.pcl", append=FALSE, quote=FALSE, sep="\t", eol="\n", na="", row.names=TRUE, col.names=TRUE)
  write.table(disease_micro, file = "/public/linzhf/ALR/1disease_micro.pcl", append=FALSE, quote=FALSE, sep="\t", eol="\n", na="", row.names=TRUE, col.names=TRUE)
  
  # Read the PCL files back into R
  TUMOR.pcl <- read.table("/public/linzhf/ALR/1disease_micro.pcl", header = TRUE, sep = "\t", quote = "", fill = TRUE, comment.char = "")
  NORMAL.pcl <- read.table("/public/linzhf/ALR/1normal_micro.pcl", header = TRUE, sep = "\t", quote = "", fill = TRUE, comment.char = "")
  
  # ALR Normalization
  # Set pseudo-count value
  pseudo_count <- 1
  NORMAL.pcl <- NORMAL.pcl + pseudo_count  # Add pseudo-count to normal data
  TUMOR.pcl <- TUMOR.pcl + pseudo_count    # Add pseudo-count to tumor data
  
  # Transform to compositional data
  NORMAL.pcl <- apply(NORMAL.pcl, 2, function(x) x / sum(x))  # Normalize by column sums
  TUMOR.pcl <- apply(TUMOR.pcl, 2, function(x) x / sum(x))    # Normalize by column sums
  
  # Apply ALR transformation
  NORMAL.pcl <- alr(NORMAL.pcl)  # Perform ALR on normal data
  TUMOR.pcl <- alr(TUMOR.pcl)    # Perform ALR on tumor data
  
  # Apply MAVD Part 1 function to the normalized data
  mavd_part1 <- mavd_part1(normal.pcl = NORMAL.pcl, tumor.pcl = TUMOR.pcl)
  
  # Plot Wold's T² statistics for the normal data
  plot.wold.mavd(mavd_part1$wold.Nmat, main.extra = "Normal")
  mavd_part1$wold.Nmat  # Display Wold's T² statistics
  
  # Identify the peak index in the Wold's plot
  peak_index_l <- which(diff(mavd_part1$wold.Nmat) > 0)[1] 
  peak_index_2 <- which(diff(mavd_part1$wold.Nmat[-c(1:peak_index_l)]) < 0)[1] 
  
  K <- peak_index_l + peak_index_2  # User chooses dimension K where Wold's plot peaks (local maximum)
  abline(v = K, lty = "dotted", col = "red")  # Draw a vertical line at dimension K
  legend(x = K, y = 0.2 * mavd_part1$wold.Nmat[[1]], legend = paste("dim =", K), bty = "n")  # Add legend to the plot
  
  #### DATA TRANSFORMATION - 
  #### DATA IS DECOMPOSED INTO DISEASE & NORMAL COMPONENTS 
  mavd_part2 <- mavd_part2(x = mavd_part1, k = K)  # Apply MAVD Part 2 function
  
  D <- mavd_part2$Dc.Dmat  # Extract disease component
  N <- mavd_part2$L1.mat   # Extract normal component
  
  # Filter the components based on mean absolute difference
  diff <- data.frame(D)[apply(abs(D), 1, mean) > apply(abs(N), 1, quantile, probs=res.cutoff), ]
  gold_names_new <- rownames(diff)  # Get the names of significant features
  gold_names_new  # Return the names of significant features
  
  # return(list(Dc.Dmat=D, L1.mat=N, gold_names=gold_names_new))  # Optionally return additional data
  return(list(gold_names=gold_names_new))  # Return only the significant feature names
}

#### CLR ####
res_MAVD_CLR <- function(sim=sim, res.cutoff=0.95) {
  group <- sim$group  # Extract group information from the simulation object
  data <- sim$data    # Extract the data matrix from the simulation object
  
  # Split the data into normal and disease microbiome subsets
  normal_micro <- data[, group == "0"]
  disease_micro <- data[, group == "1"]
  
  # Write the normal and disease microbiome data to PCL files
  write.table(normal_micro, file = "/public/linzhf/CLR/normal_micro.pcl", append=FALSE, quote=FALSE, sep="\t", eol="\n", na="", row.names=TRUE, col.names=TRUE)
  write.table(disease_micro, file = "/public/linzhf/CLR/disease_micro.pcl", append=FALSE, quote=FALSE, sep="\t", eol="\n", na="", row.names=TRUE, col.names=TRUE)
  
  # Read the PCL files back into R
  TUMOR.pcl <- read.table("/public/linzhf/CLR/disease_micro.pcl", header = TRUE, sep = "\t", quote = "", fill = TRUE, comment.char = "")
  NORMAL.pcl <- read.table("/public/linzhf/CLR/normal_micro.pcl", header = TRUE, sep = "\t", quote = "", fill = TRUE, comment.char = "")
  
  # CLR Normalization
  # Set pseudo-count value
  pseudo_count <- 1
  NORMAL.pcl <- NORMAL.pcl + pseudo_count  # Add pseudo-count to normal data
  TUMOR.pcl <- TUMOR.pcl + pseudo_count    # Add pseudo-count to tumor data
  
  # Transform to compositional data
  NORMAL.pcl <- apply(NORMAL.pcl, 2, function(x) x / sum(x))  # Normalize by column sums
  TUMOR.pcl <- apply(TUMOR.pcl, 2, function(x) x / sum(x))    # Normalize by column sums
  
  # Apply CLR transformation
  NORMAL.pcl <- clr(NORMAL.pcl)  # Perform CLR on normal data
  TUMOR.pcl <- clr(TUMOR.pcl)    # Perform CLR on tumor data
  
  # Apply MAVD Part 1 function to the normalized data
  mavd_part1 <- mavd_part1(normal.pcl = NORMAL.pcl, tumor.pcl = TUMOR.pcl)
  
  # Plot Wold's T² statistics for the normal data
  plot.wold.mavd(mavd_part1$wold.Nmat, main.extra = "Normal")
  mavd_part1$wold.Nmat  # Display Wold's T² statistics
  
  # Identify the peak index in the Wold's plot
  peak_index_l <- which(diff(mavd_part1$wold.Nmat) > 0)[1] 
  peak_index_2 <- which(diff(mavd_part1$wold.Nmat[-c(1:peak_index_l)]) < 0)[1] 
  
  K <- peak_index_l + peak_index_2  # User chooses dimension K where Wold's plot peaks (local maximum)
  abline(v = K, lty = "dotted", col = "red")  # Draw a vertical line at dimension K
  legend(x = K, y = 0.2 * mavd_part1$wold.Nmat[[1]], legend = paste("dim =", K), bty = "n")  # Add legend to the plot
  
  #### DATA TRANSFORMATION - 
  #### DATA IS DECOMPOSED INTO DISEASE & NORMAL COMPONENTS 
  mavd_part2 <- mavd_part2(x = mavd_part1, k = K)  # Apply MAVD Part 2 function
  
  D <- mavd_part2$Dc.Dmat  # Extract disease component
  N <- mavd_part2$L1.mat   # Extract normal component
  
  # Filter the components based on mean absolute difference
  diff <- data.frame(D)[apply(abs(D), 1, mean) > apply(abs(N), 1, quantile, probs=res.cutoff), ]
  gold_names_new <- rownames(diff)  # Get the names of significant features
  gold_names_new  # Return the names of significant features
  
  # return(list(Dc.Dmat=D, L1.mat=N, gold_names=gold_names_new))  # Optionally return additional data
  return(list(gold_names=gold_names_new))  # Return only the significant feature names
}

