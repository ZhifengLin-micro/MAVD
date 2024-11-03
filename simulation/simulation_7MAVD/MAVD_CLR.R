#Set Working Directory
setwd("/public/linzhf")

# Load necessary libraries
library(compositions)  
library(edgeR)
library(DESeq2)
library(Wrench)
library(GUniFrac)
library(metagenomeSeq)

# Load custom MAVD functions
source("mavd_source.R")

############ Data Processing
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

#### MAVD_CLR simulation ####
folder_path <- "/public/linzhf/simulation_data/"
output_folder <- "/public/linzhf/simulation_output/CLR"
# Get the list of RData files in the folder
file_list <- list.files(path = folder_path, pattern = ".RData", full.names = TRUE)
# Iterate over each RData file
for (file_path in file_list) {
  tryCatch({
    load(file_path) # Load the data from the current file
    sim <- res_MAVD_CLR(sim)  # Apply the new_cut function
    file_name <- tools::file_path_sans_ext(basename(file_path)) # Extract the file name without extension
    save(sim, file = file.path(output_folder, paste0(file_name, ".RData")))
  },  error = function(e){cat(paste0("error in ",file_path))} )
}