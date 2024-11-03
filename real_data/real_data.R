setwd("/public/linzhf/real_data")

source("mavd_source.R")
source("all_function.R")
#####################################################
####MAVD_TSS####
data <- crc_zeller_genus_table
group <- crc_zeller_metadata  %>% setNames("group")
data <- as.matrix(data)
sim <- list(data=data,group=group) 

group <- sim$group
data <- sim$data
normal_micro <- data[,group == "H"]
disease_micro <- data[,group == "CRC"]

write.table(normal_micro, file = "/public/linzhf/real_data/crc_zeller/cdi_schubert/normal_micro.pcl",append=FALSE,quote=FALSE,sep="\t",eol="\n",na="",row.names=TRUE,col.names=TRUE)
write.table(disease_micro, file = "/public/linzhf/real_data/crc_zeller/cdi_schubert/disease_micro.pcl",append=FALSE,quote=FALSE,sep="\t",eol="\n",na="",row.names=TRUE,col.names=TRUE)

TUMOR.pcl <- read.table("/public/linzhf/real_data/crc_zeller/cdi_schubert/disease_micro.pcl", header = TRUE, sep = "\t", quote = "", fill = TRUE, comment.char = "")
NORMAL.pcl <- read.table("/public/linzhf/real_data/crc_zeller/cdi_schubert/normal_micro.pcl", header = TRUE, sep = "\t", quote = "", fill = TRUE, comment.char = "")

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
write.table(gold_names_new,file = "intersect/MAVD_TSS.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

##########################################
####13 DAA####
data <- crc_zeller_genus_table
group <- crc_zeller_metadata  %>% setNames("group")
data <- as.matrix(data)
sim <- list(data=data,group=group) 

physeq <- phyloseq(otu_table(sim$data, taxa_are_rows = TRUE), 
                   sample_data(data.frame(sim$group)))

tax_tab <- matrix(1,nrow = nrow(sim$data), ncol = 7)
rownames(tax_tab) = rownames(sim$data)
colnames(tax_tab) = c("Kingdom", "Phylum", "Class", "Order",
                      "Family", "Genus", "Species")
tax_tab <- as.data.frame(tax_tab)
tax_tab$Kingdom <- rownames(sim$data)

p.adj = "fdr"
p.cutoff = 0.05

aldex2 <- res_aldex2(sim, type="wi")
corncob <- res_corncob(sim)
dacomp <- res_dacomp(sim)
deseq2 <- res_deseq2(sim)
edger <- res_edger(sim)
fastancom <- res_fastancom(sim)
lefse <- res_lefse(sim)
limma <- res_limma_voom(sim)
linda <- res_linda(sim)
maaslin2 <- res_maaslin2(sim)
mseq <- res_metagenomeseq(sim)
wilcox <- res_wilcox(sim)
zicoseq <- res_zicoseq(sim)

write.table(aldex2$gold_names,file = "intersect/ALDEx2.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(corncob$gold_names,file = "intersect/Corncob.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(dacomp$gold_names,file = "intersect/DACOMP.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(deseq2$gold_names,file = "intersect/DESeq2.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(edger$gold_names,file = "intersect/edgeR.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(fastancom$gold_names,file = "intersect/fastANCOM.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(lefse$gold_names,file = "intersect/LEfSe.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(limma$gold_names,file = "intersect/limma_voom.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(linda$gold_names,file = "intersect/LinDA.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(maaslin2$gold_names,file = "intersect/maaslin2.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(mseq$gold_names,file = "intersect/metagenomeSeq.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(wilcox$gold_names,file = "intersect/Wilcox.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(zicoseq$gold_names,file = "intersect/ZicoSeq.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

# 设置文件夹路径和输出文件名
folder_path <- "/public/linzhf/real_data/crc_zeller/crc_zeller/intersect"  # 替换为你的文件夹路径
output_file <- "crc_zeller.txt"   # 合并后的文件名

# 获取文件夹中所有txt文件的路径
file_paths <- list.files(path = folder_path, pattern = "\\.txt$", full.names = TRUE)

# 打开输出文件
output_connection <- file(output_file, "w")

# 遍历所有txt文件
for (file_path in file_paths) {
  # 打开并读取txt文件内容
  file_content <- readLines(file_path)
  # 写入内容到输出文件
  writeLines(file_content, output_connection)
}

# 关闭输出文件连接
close(output_connection)

