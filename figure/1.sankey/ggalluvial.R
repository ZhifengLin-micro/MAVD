library(ggalluvial)
library(ggplot2)
library(dplyr)

inputFile="input6指标.txt"         
outFile="ggalluvial6.pdf"    
setwd("/public/linzhf/results")         
rt=read.table(inputFile, header = T, sep="\t", check.names=F)     
corLodes=to_lodes_form(rt, axes = 1:ncol(rt), id = "Cohort")

#得到输出文件
pdf(file=outFile,width=8,height=6)
mycol <- rep(c("#9B3A4D","#394A92","#E2AE79","#68AC57","#D0DCAA","#497EB2","#F0EEBB","#8E549E","#8CBDA7","#D2352C","#566CA5","#91612D","#70A0AC","#BAAFD1","#431A3D","#FFD121","#088247"),17)
ggplot(corLodes, aes(x = x, stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
  	 scale_x_discrete(expand = c(0, 0)) +  
  	 #用aes.flow控制线调颜色，forward说明颜色和前面一致，backward说明与后面一致
  	 geom_flow(width = 2/10,aes.flow = "forward") + 
	 geom_stratum(alpha = .9,width = 2/10) +
	 scale_fill_manual(values = mycol) +
	 #size = 2代表基因名字大小
	 geom_text(stat = "stratum", size = 3,color="black") +
	 xlab("") + ylab("") + theme_bw() + 
	 theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank()) + 
	 theme(panel.grid =element_blank()) + 
	 theme(panel.border = element_blank()) + 
	 ggtitle("") + guides(fill = FALSE)                            
dev.off()

