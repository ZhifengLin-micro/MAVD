
####绘制条形图OTU=100####
# 加载数据
library(dplyr)
library(ggplot2)
library(patchwork)
input <- read.delim("/public/linzhf/results/input.txt")

split_data1 <- split(input,input$effect)
effect_1 <- split_data1$"1"
effect_2 <- split_data1$"2"
effect_5 <- split_data1$"5"
effect_10 <- split_data1$"10"
# 定义自定义排序顺序
custom_order <- c("50_200", "50_500", "50_1000","200_200", "200_500", "200_1000","500_200", "500_500", "500_1000")

# 将sampleotu转换为有序因子变量
effect_1$sampleotu <- factor(effect_1$sampleotu, levels = custom_order, ordered = TRUE)
effect_2$sampleotu <- factor(effect_2$sampleotu, levels = custom_order, ordered = TRUE)
effect_5$sampleotu <- factor(effect_5$sampleotu, levels = custom_order, ordered = TRUE)
effect_10$sampleotu <- factor(effect_10$sampleotu, levels = custom_order, ordered = TRUE)

####OTU=Precision####
# 使用ggplot2绘制分面的分组条形图，旋转横坐标标签，并修改颜色
mycol <- rep(c("#9B3A4D","#394A92","#E2AE79","#68AC57","#D0DCAA","#497EB2","#F0EEBB","#8E549E","#BAAFD1","#8CBDA7","#566CA5","#91612D","#70A0AC","#D2352C","#431A3D","#FFD121","#088247"), 17)

plot_1 <- ggplot(effect_1, aes(x = name, y = mean_Precision, fill = name)) +
   geom_bar(stat = "identity", position = "dodge", color = "black") +
   geom_errorbar(aes(ymin =mean_Precision -0, ymax = mean_Precision + sd_Precision),
      position = position_dodge(width = 0.9), width = 0.25, color = "black") +
   labs(title = "Effect_size=1", y = "Precision") + 
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "bold"),
      panel.grid.major = element_blank(),  # Remove internal grid lines
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add outer border
      axis.ticks.y = element_line(color = "black", size = 0.5),
      strip.background = element_rect(fill = NA),  # Transparent strip background
      strip.text = element_text(color = "black", face = "bold")) +  # Modify strip text color if needed
   scale_fill_manual(values = mycol) +
   facet_wrap(~sampleotu, scales = "free_y", ncol = 9)+
   scale_y_continuous(limits = c(0, 1.1))
plot_1

plot_2 <- ggplot(effect_2, aes(x = name, y = mean_Precision, fill = name)) +
   geom_bar(stat = "identity", position = "dodge", color = "black") +
   geom_errorbar(aes(ymin =mean_Precision-0, ymax = mean_Precision + sd_Precision),
      position = position_dodge(width = 0.9), width = 0.25, color = "black") +
   labs(title = "Effect_size=2", y = "Precision") + 
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "bold"),
      panel.grid.major = element_blank(),  # Remove internal grid lines
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add outer border
      axis.ticks.y = element_line(color = "black", size = 0.5),
      strip.background = element_rect(fill = NA),  # Transparent strip background
      strip.text = element_text(color = "black", face = "bold")) +  # Modify strip text color if needed
   scale_fill_manual(values = mycol) +
   facet_wrap(~sampleotu, scales = "free_y", ncol = 9)+
   scale_y_continuous(limits = c(0, 1.1))

plot_2

plot_5 <- ggplot(effect_5, aes(x = name, y = mean_Precision, fill = name)) +
   geom_bar(stat = "identity", position = "dodge", color = "black") +
   geom_errorbar(aes(ymin =mean_Precision-0, ymax = mean_Precision + sd_Precision),
      position = position_dodge(width = 0.9), width = 0.25, color = "black") +
   labs(title = "Effect_size=5", y = "Precision") + 
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "bold"),
      panel.grid.major = element_blank(),  # Remove internal grid lines
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add outer border
      axis.ticks.y = element_line(color = "black", size = 0.5),
      strip.background = element_rect(fill = NA),  # Transparent strip background
      strip.text = element_text(color = "black", face = "bold")) +  # Modify strip text color if needed
   scale_fill_manual(values = mycol) +
   facet_wrap(~sampleotu, scales = "free_y", ncol = 9)+
   scale_y_continuous(limits = c(0, 1.1))
plot_5

plot_10 <- ggplot(effect_10, aes(x = name, y = mean_Precision, fill = name)) +
   geom_bar(stat = "identity", position = "dodge", color = "black") +
   geom_errorbar(aes(ymin =mean_Precision-0, ymax = mean_Precision + sd_Precision),
      position = position_dodge(width = 0.9), width = 0.25, color = "black") +
   labs(title = "Effect_size=10", y = "Precision") + 
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "bold"),
      panel.grid.major = element_blank(),  # Remove internal grid lines
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add outer border
      axis.ticks.y = element_line(color = "black", size = 0.5),
      strip.background = element_rect(fill = NA),  # Transparent strip background
      strip.text = element_text(color = "black", face = "bold")) +  # Modify strip text color if needed
   scale_fill_manual(values = mycol) +
   facet_wrap(~sampleotu, scales = "free_y", ncol = 9)+
   scale_y_continuous(limits = c(0, 1.1))
plot_10

combined_plot <- plot_1 + plot_2 + plot_5 + plot_10

# Customize the layout if needed
combined_plot <- combined_plot + plot_layout(ncol = 1)

# 保存为PDF文件
ggsave("Precision.pdf", combined_plot, width = 15, height = 12)


####TPR####
# 使用ggplot2绘制分面的分组条形图，旋转横坐标标签，并修改颜色
mycol <- rep(c("#9B3A4D","#394A92","#E2AE79","#68AC57","#D0DCAA","#497EB2","#F0EEBB","#8E549E","#BAAFD1","#8CBDA7","#566CA5","#91612D","#70A0AC","#D2352C","#431A3D","#FFD121","#088247"), 17)

plot_1 <- ggplot(effect_1, aes(x = name, y = mean_TPR, fill = name)) +
   geom_bar(stat = "identity", position = "dodge", color = "black") +
   geom_errorbar(aes(ymin =mean_TPR -0, ymax = mean_TPR + sd_TPR),
      position = position_dodge(width = 0.9), width = 0.25, color = "black") +
   labs(title = "Effect_size=1", y = "TPR") + 
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "bold"),
      panel.grid.major = element_blank(),  # Remove internal grid lines
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add outer border
      axis.ticks.y = element_line(color = "black", size = 0.5),
      strip.background = element_rect(fill = NA),  # Transparent strip background
      strip.text = element_text(color = "black", face = "bold")) +  # Modify strip text color if needed
   scale_fill_manual(values = mycol) +
   facet_wrap(~sampleotu, scales = "free_y", ncol = 9)+
   scale_y_continuous(limits = c(0, 1.1))
plot_1

plot_2 <- ggplot(effect_2, aes(x = name, y = mean_TPR, fill = name)) +
   geom_bar(stat = "identity", position = "dodge", color = "black") +
   geom_errorbar(aes(ymin =mean_TPR-0, ymax = mean_TPR + sd_TPR),
      position = position_dodge(width = 0.9), width = 0.25, color = "black") +
   labs(title = "Effect_size=2", y = "TPR") + 
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "bold"),
      panel.grid.major = element_blank(),  # Remove internal grid lines
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add outer border
      axis.ticks.y = element_line(color = "black", size = 0.5),
      strip.background = element_rect(fill = NA),  # Transparent strip background
      strip.text = element_text(color = "black", face = "bold")) +  # Modify strip text color if needed
   scale_fill_manual(values = mycol) +
   facet_wrap(~sampleotu, scales = "free_y", ncol = 9)+
   scale_y_continuous(limits = c(0, 1.1))

plot_2

plot_5 <- ggplot(effect_5, aes(x = name, y = mean_TPR, fill = name)) +
   geom_bar(stat = "identity", position = "dodge", color = "black") +
   geom_errorbar(aes(ymin =mean_TPR-0, ymax = mean_TPR + sd_TPR),
      position = position_dodge(width = 0.9), width = 0.25, color = "black") +
   labs(title = "Effect_size=5", y = "TPR") + 
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "bold"),
      panel.grid.major = element_blank(),  # Remove internal grid lines
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add outer border
      axis.ticks.y = element_line(color = "black", size = 0.5),
      strip.background = element_rect(fill = NA),  # Transparent strip background
      strip.text = element_text(color = "black", face = "bold")) +  # Modify strip text color if needed
   scale_fill_manual(values = mycol) +
   facet_wrap(~sampleotu, scales = "free_y", ncol = 9)+
   scale_y_continuous(limits = c(0, 1.1))
plot_5

plot_10 <- ggplot(effect_10, aes(x = name, y = mean_TPR, fill = name)) +
   geom_bar(stat = "identity", position = "dodge", color = "black") +
   geom_errorbar(aes(ymin =mean_TPR-0, ymax = mean_TPR + sd_TPR),
      position = position_dodge(width = 0.9), width = 0.25, color = "black") +
   labs(title = "Effect_size=10", y = "TPR") + 
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "bold"),
      panel.grid.major = element_blank(),  # Remove internal grid lines
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add outer border
      axis.ticks.y = element_line(color = "black", size = 0.5),
      strip.background = element_rect(fill = NA),  # Transparent strip background
      strip.text = element_text(color = "black", face = "bold")) +  # Modify strip text color if needed
   scale_fill_manual(values = mycol) +
   facet_wrap(~sampleotu, scales = "free_y", ncol = 9)+
   scale_y_continuous(limits = c(0, 1.1))
plot_10

combined_plot <- plot_1 + plot_2 + plot_5 + plot_10

# Customize the layout if needed
combined_plot <- combined_plot + plot_layout(ncol = 1)

# 保存为PDF文件
ggsave("TPR.pdf", combined_plot, width = 15, height = 12)

####FDR####
# 使用ggplot2绘制分面的分组条形图，旋转横坐标标签，并修改颜色
mycol <- rep(c("#9B3A4D","#394A92","#E2AE79","#68AC57","#D0DCAA","#497EB2","#F0EEBB","#8E549E","#BAAFD1","#8CBDA7","#566CA5","#91612D","#70A0AC","#D2352C","#431A3D","#FFD121","#088247"), 17)

plot_1 <- ggplot(effect_1, aes(x = name, y = mean_FDR, fill = name)) +
   geom_bar(stat = "identity", position = "dodge", color = "black") +
   geom_errorbar(aes(ymin =mean_FDR -0, ymax = mean_FDR + sd_FDR),
      position = position_dodge(width = 0.9), width = 0.25, color = "black") +
   labs(title = "Effect_size=1", y = "FDR") + 
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "bold"),
      panel.grid.major = element_blank(),  # Remove internal grid lines
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add outer border
      axis.ticks.y = element_line(color = "black", size = 0.5),
      strip.background = element_rect(fill = NA),  # Transparent strip background
      strip.text = element_text(color = "black", face = "bold")) +  # Modify strip text color if needed
   scale_fill_manual(values = mycol) +
   facet_wrap(~sampleotu, scales = "free_y", ncol = 9)+
   scale_y_continuous(limits = c(0, 1.1))
plot_1

plot_2 <- ggplot(effect_2, aes(x = name, y = mean_FDR, fill = name)) +
   geom_bar(stat = "identity", position = "dodge", color = "black") +
   geom_errorbar(aes(ymin =mean_FDR-0, ymax = mean_FDR + sd_FDR),
      position = position_dodge(width = 0.9), width = 0.25, color = "black") +
   labs(title = "Effect_size=2", y = "Precision") + 
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "bold"),
      panel.grid.major = element_blank(),  # Remove internal grid lines
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add outer border
      axis.ticks.y = element_line(color = "black", size = 0.5),
      strip.background = element_rect(fill = NA),  # Transparent strip background
      strip.text = element_text(color = "black", face = "bold")) +  # Modify strip text color if needed
   scale_fill_manual(values = mycol) +
   facet_wrap(~sampleotu, scales = "free_y", ncol = 9)+
   scale_y_continuous(limits = c(0, 1.1))

plot_2

plot_5 <- ggplot(effect_5, aes(x = name, y = mean_FDR, fill = name)) +
   geom_bar(stat = "identity", position = "dodge", color = "black") +
   geom_errorbar(aes(ymin =mean_FDR-0, ymax = mean_FDR + sd_FDR),
      position = position_dodge(width = 0.9), width = 0.25, color = "black") +
   labs(title = "Effect_size=5", y = "FDR") + 
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "bold"),
      panel.grid.major = element_blank(),  # Remove internal grid lines
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add outer border
      axis.ticks.y = element_line(color = "black", size = 0.5),
      strip.background = element_rect(fill = NA),  # Transparent strip background
      strip.text = element_text(color = "black", face = "bold")) +  # Modify strip text color if needed
   scale_fill_manual(values = mycol) +
   facet_wrap(~sampleotu, scales = "free_y", ncol = 9)+
   scale_y_continuous(limits = c(0, 1.1))
plot_5

plot_10 <- ggplot(effect_10, aes(x = name, y = mean_FDR, fill = name)) +
   geom_bar(stat = "identity", position = "dodge", color = "black") +
   geom_errorbar(aes(ymin =mean_FDR-0, ymax = mean_FDR + sd_FDR),
      position = position_dodge(width = 0.9), width = 0.25, color = "black") +
   labs(title = "Effect_size=10", y = "FDR") + 
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "bold"),
      panel.grid.major = element_blank(),  # Remove internal grid lines
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add outer border
      axis.ticks.y = element_line(color = "black", size = 0.5),
      strip.background = element_rect(fill = NA),  # Transparent strip background
      strip.text = element_text(color = "black", face = "bold")) +  # Modify strip text color if needed
   scale_fill_manual(values = mycol) +
   facet_wrap(~sampleotu, scales = "free_y", ncol = 9)+
   scale_y_continuous(limits = c(0, 1.1))
plot_10

combined_plot <- plot_1 + plot_2 + plot_5 + plot_10

# Customize the layout if needed
combined_plot <- combined_plot + plot_layout(ncol = 1)

# 保存为PDF文件
ggsave("FDR.pdf", combined_plot, width = 15, height = 12)

####F1_score####
# 使用ggplot2绘制分面的分组条形图，旋转横坐标标签，并修改颜色
mycol <- rep(c("#9B3A4D","#394A92","#E2AE79","#68AC57","#D0DCAA","#497EB2","#F0EEBB","#8E549E","#BAAFD1","#8CBDA7","#566CA5","#91612D","#70A0AC","#D2352C","#431A3D","#FFD121","#088247"), 17)

plot_1 <- ggplot(effect_1, aes(x = name, y = mean_F1_score, fill = name)) +
   geom_bar(stat = "identity", position = "dodge", color = "black") +
   geom_errorbar(aes(ymin =mean_F1_score -0, ymax = mean_F1_score + sd_F1_score),
      position = position_dodge(width = 0.9), width = 0.25, color = "black") +
   labs(title = "Effect_size=1", y = "F1_score") + 
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "bold"),
      panel.grid.major = element_blank(),  # Remove internal grid lines
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add outer border
      axis.ticks.y = element_line(color = "black", size = 0.5),
      strip.background = element_rect(fill = NA),  # Transparent strip background
      strip.text = element_text(color = "black", face = "bold")) +  # Modify strip text color if needed
   scale_fill_manual(values = mycol) +
   facet_wrap(~sampleotu, scales = "free_y", ncol = 9)+
   scale_y_continuous(limits = c(0, 1.1))
plot_1

plot_2 <- ggplot(effect_2, aes(x = name, y = mean_F1_score, fill = name)) +
   geom_bar(stat = "identity", position = "dodge", color = "black") +
   geom_errorbar(aes(ymin =mean_F1_score-0, ymax = mean_F1_score + sd_F1_score),
      position = position_dodge(width = 0.9), width = 0.25, color = "black") +
   labs(title = "Effect_size=2", y = "F1_score") + 
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "bold"),
      panel.grid.major = element_blank(),  # Remove internal grid lines
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add outer border
      axis.ticks.y = element_line(color = "black", size = 0.5),
      strip.background = element_rect(fill = NA),  # Transparent strip background
      strip.text = element_text(color = "black", face = "bold")) +  # Modify strip text color if needed
   scale_fill_manual(values = mycol) +
   facet_wrap(~sampleotu, scales = "free_y", ncol = 9)+
   scale_y_continuous(limits = c(0, 1.1))

plot_2

plot_5 <- ggplot(effect_5, aes(x = name, y = mean_F1_score, fill = name)) +
   geom_bar(stat = "identity", position = "dodge", color = "black") +
   geom_errorbar(aes(ymin =mean_F1_score-0, ymax = mean_F1_score + sd_F1_score),
      position = position_dodge(width = 0.9), width = 0.25, color = "black") +
   labs(title = "Effect_size=5", y = "F1_score") + 
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "bold"),
      panel.grid.major = element_blank(),  # Remove internal grid lines
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add outer border
      axis.ticks.y = element_line(color = "black", size = 0.5),
      strip.background = element_rect(fill = NA),  # Transparent strip background
      strip.text = element_text(color = "black", face = "bold")) +  # Modify strip text color if needed
   scale_fill_manual(values = mycol) +
   facet_wrap(~sampleotu, scales = "free_y", ncol = 9)+
   scale_y_continuous(limits = c(0, 1.1))
plot_5

plot_10 <- ggplot(effect_10, aes(x = name, y = mean_F1_score, fill = name)) +
   geom_bar(stat = "identity", position = "dodge", color = "black") +
   geom_errorbar(aes(ymin =mean_F1_score-0, ymax = mean_F1_score + sd_F1_score),
      position = position_dodge(width = 0.9), width = 0.25, color = "black") +
   labs(title = "Effect_size=10", y = "F1_score") + 
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "bold"),
      panel.grid.major = element_blank(),  # Remove internal grid lines
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add outer border
      axis.ticks.y = element_line(color = "black", size = 0.5),
      strip.background = element_rect(fill = NA),  # Transparent strip background
      strip.text = element_text(color = "black", face = "bold")) +  # Modify strip text color if needed
   scale_fill_manual(values = mycol) +
   facet_wrap(~sampleotu, scales = "free_y", ncol = 9)+
   scale_y_continuous(limits = c(0, 1.1))
plot_10

combined_plot <- plot_1 + plot_2 + plot_5 + plot_10

# Customize the layout if needed
combined_plot <- combined_plot + plot_layout(ncol = 1)

# 保存为PDF文件
ggsave("F1_score.pdf", combined_plot, width = 15, height = 12)

####FPR####
# 使用ggplot2绘制分面的分组条形图，旋转横坐标标签，并修改颜色
mycol <- rep(c("#9B3A4D","#394A92","#E2AE79","#68AC57","#D0DCAA","#497EB2","#F0EEBB","#8E549E","#BAAFD1","#8CBDA7","#566CA5","#91612D","#70A0AC","#D2352C","#431A3D","#FFD121","#088247"), 17)

plot_1 <- ggplot(effect_1, aes(x = name, y = mean_FPR, fill = name)) +
   geom_bar(stat = "identity", position = "dodge", color = "black") +
   geom_errorbar(aes(ymin =mean_FPR -0, ymax = mean_FPR + sd_FPR),
      position = position_dodge(width = 0.9), width = 0.25, color = "black") +
   labs(title = "Effect_size=1", y = "FPR") + 
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "bold"),
      panel.grid.major = element_blank(),  # Remove internal grid lines
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add outer border
      axis.ticks.y = element_line(color = "black", size = 0.5),
      strip.background = element_rect(fill = NA),  # Transparent strip background
      strip.text = element_text(color = "black", face = "bold")) +  # Modify strip text color if needed
   scale_fill_manual(values = mycol) +
   facet_wrap(~sampleotu, scales = "free_y", ncol = 9)+
   scale_y_continuous(limits = c(0, 0.5))
plot_1

plot_2 <- ggplot(effect_2, aes(x = name, y = mean_FPR, fill = name)) +
   geom_bar(stat = "identity", position = "dodge", color = "black") +
   geom_errorbar(aes(ymin =mean_FPR-0, ymax = mean_FPR + sd_FPR),
      position = position_dodge(width = 0.9), width = 0.25, color = "black") +
   labs(title = "Effect_size=2", y = "FPR") + 
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "bold"),
      panel.grid.major = element_blank(),  # Remove internal grid lines
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add outer border
      axis.ticks.y = element_line(color = "black", size = 0.5),
      strip.background = element_rect(fill = NA),  # Transparent strip background
      strip.text = element_text(color = "black", face = "bold")) +  # Modify strip text color if needed
   scale_fill_manual(values = mycol) +
   facet_wrap(~sampleotu, scales = "free_y", ncol = 9)+
   scale_y_continuous(limits = c(0, 0.5))

plot_2

plot_5 <- ggplot(effect_5, aes(x = name, y = mean_FPR, fill = name)) +
   geom_bar(stat = "identity", position = "dodge", color = "black") +
   geom_errorbar(aes(ymin =mean_FPR-0, ymax = mean_FPR + sd_FPR),
      position = position_dodge(width = 0.9), width = 0.25, color = "black") +
   labs(title = "Effect_size=5", y = "FPR") + 
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "bold"),
      panel.grid.major = element_blank(),  # Remove internal grid lines
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add outer border
      axis.ticks.y = element_line(color = "black", size = 0.5),
      strip.background = element_rect(fill = NA),  # Transparent strip background
      strip.text = element_text(color = "black", face = "bold")) +  # Modify strip text color if needed
   scale_fill_manual(values = mycol) +
   facet_wrap(~sampleotu, scales = "free_y", ncol = 9)+
   scale_y_continuous(limits = c(0, 0.5))
plot_5

plot_10 <- ggplot(effect_10, aes(x = name, y = mean_FPR, fill = name)) +
   geom_bar(stat = "identity", position = "dodge", color = "black") +
   geom_errorbar(aes(ymin =mean_FPR-0, ymax = mean_FPR + sd_FPR),
      position = position_dodge(width = 0.9), width = 0.25, color = "black") +
   labs(title = "Effect_size=10", y = "FPR") + 
   theme_minimal() +
   theme(axis.text.x = element_text(angle = 60, hjust = 1, face = "bold"),
      panel.grid.major = element_blank(),  # Remove internal grid lines
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add outer border
      axis.ticks.y = element_line(color = "black", size = 0.5),
      strip.background = element_rect(fill = NA),  # Transparent strip background
      strip.text = element_text(color = "black", face = "bold")) +  # Modify strip text color if needed
   scale_fill_manual(values = mycol) +
   facet_wrap(~sampleotu, scales = "free_y", ncol = 9)+
   scale_y_continuous(limits = c(0, 0.5))
plot_10

combined_plot <- plot_1 + plot_2 + plot_5 + plot_10

# Customize the layout if needed
combined_plot <- combined_plot + plot_layout(ncol = 1)

# 保存为PDF文件
ggsave("FPR.pdf", combined_plot, width = 15, height = 12)


