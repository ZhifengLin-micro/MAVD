# 导入 ggplot2 包
library(dplyr)
library(ggplot2)
TPRFPR <- read.delim("/public/linzhf/3.TPRFPR/TPRFPR1.txt")
data_effect_1 <- subset(TPRFPR, effect == 1)
data_effect_2 <- subset(TPRFPR, effect == 2)
data_effect_5 <- subset(TPRFPR, effect == 5)
data_effect_10 <- subset(TPRFPR, effect == 10)

sample_colors <- c("#D52527", "#614099", "#088247","#CD5A1F", "#0075C5", "#8B5549","#FF7D0D", "#F0A298", "#7D7D7D","#98C5E1", "#8481ba", "#6BB952", "#FAEE85","#234444" )

data10<- data_effect_10
data5<- data_effect_5
data2<- data_effect_2
data1<- data_effect_1
# 创建散点图并按照样本大小和OTU数量进行分面
scatter_plot10 <- ggplot(data10, aes(x = mean_FPR, y = mean_TPR, shape = name, color = name)) +
   geom_point(size = 3) +  # 绘制散点
   geom_abline(intercept = 0, slope = 1, color = "gray", linetype = "dashed") +  # 添加对角线
   labs(x = "False Positive Rate", y = "True Positive Rate", title = "Effect size = 10") +  # 添加坐标轴标签和标题
   theme_minimal() +  # 使用简洁主题
   coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +  # 设置横轴和纵轴范围
   theme(panel.border = element_rect(color = "black", fill = NA, size = 0.8),  # 设置外边框
      panel.spacing = unit(1, "lines"),  # 设置面板间距
      panel.background = element_rect(color = "gray50", fill = NA),  # 设置面板背景
      strip.background = element_blank(),  # 移除分面标题背景
      strip.text.x = element_text(size = 8, angle = 0, hjust = 0.5)) +  # 设置分面标题文本
   # facet_wrap(otu~sample)+theme_bw()+
     facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(Sample = label_both,OTU = label_both )) + theme_bw()+ # 按样本大小和OTU数量进行分面并命名
   guides(shape = guide_legend(label.position = "right", title = "model"),
      color = guide_legend(label.position = "right", title = "model")) + # 在图例中显示标签值
   scale_shape_manual(values = c(17,18,16,15,4,3,2,1,0,5,6,7,10,11)) + # 设置每种形状对应的值
   scale_color_manual(values = sample_colors) +  # 设置颜色
   theme(legend.text = element_text(colour = "black"),  # 设置图例文本颜色为黑色
      legend.title = element_text(colour = "black"),  # 设置图例标题颜色为黑色
      legend.position = "right")   # 设置图例位置为右侧

# 创建散点图并按照样本大小和OTU数量进行分面
scatter_plot5 <- ggplot(data5, aes(x = mean_FPR, y = mean_TPR, shape = name, color = name)) +
   geom_point(size = 3) +  # 绘制散点
   geom_abline(intercept = 0, slope = 1, color = "gray", linetype = "dashed") +  # 添加对角线
   labs(x = "False Positive Rate", y = "True Positive Rate", title = "Effect size = 5") +  # 添加坐标轴标签和标题
   theme_minimal() +  # 使用简洁主题
   coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +  # 设置横轴和纵轴范围
   theme(panel.border = element_rect(color = "black", fill = NA, size = 0.8),  # 设置外边框
      panel.spacing = unit(1, "lines"),  # 设置面板间距
      panel.background = element_rect(color = "gray50", fill = NA),  # 设置面板背景
      strip.background = element_blank(),  # 移除分面标题背景
      strip.text.x = element_text(size = 8, angle = 0, hjust = 0.5)) +  # 设置分面标题文本
   # facet_wrap(otu~sample)+theme_bw()+
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(Sample = label_both,OTU = label_both )) + theme_bw()+ # 按样本大小和OTU数量进行分面并命名
   guides(shape = guide_legend(label.position = "right", title = "model"),
      color = guide_legend(label.position = "right", title = "model")) + # 在图例中显示标签值
   scale_shape_manual(values = c(17,18,16,15,4,3,2,1,0,5,6,7,10,11)) + # 设置每种形状对应的值
   scale_color_manual(values = sample_colors) +  # 设置颜色
   theme(legend.text = element_text(colour = "black"),  # 设置图例文本颜色为黑色
      legend.title = element_text(colour = "black"),  # 设置图例标题颜色为黑色
      legend.position = "right")   # 设置图例位置为右侧

# 创建散点图并按照样本大小和OTU数量进行分面
scatter_plot2 <- ggplot(data2, aes(x = mean_FPR, y = mean_TPR, shape = name, color = name)) +
   geom_point(size = 3) +  # 绘制散点
   geom_abline(intercept = 0, slope = 1, color = "gray", linetype = "dashed") +  # 添加对角线
   labs(x = "False Positive Rate", y = "True Positive Rate", title = "Effect size = 2") +  # 添加坐标轴标签和标题
   theme_minimal() +  # 使用简洁主题
   coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +  # 设置横轴和纵轴范围
   theme(panel.border = element_rect(color = "black", fill = NA, size = 0.8),  # 设置外边框
      panel.spacing = unit(1, "lines"),  # 设置面板间距
      panel.background = element_rect(color = "gray50", fill = NA),  # 设置面板背景
      strip.background = element_blank(),  # 移除分面标题背景
      strip.text.x = element_text(size = 8, angle = 0, hjust = 0.5)) +  # 设置分面标题文本
   # facet_wrap(otu~sample)+theme_bw()+
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(Sample = label_both,OTU = label_both )) + theme_bw()+ # 按样本大小和OTU数量进行分面并命名
   guides(shape = guide_legend(label.position = "right", title = "model"),
      color = guide_legend(label.position = "right", title = "model")) + # 在图例中显示标签值
   scale_shape_manual(values = c(17,18,16,15,4,3,2,1,0,5,6,7,10,11)) + # 设置每种形状对应的值
   scale_color_manual(values = sample_colors) +  # 设置颜色
   theme(legend.text = element_text(colour = "black"),  # 设置图例文本颜色为黑色
      legend.title = element_text(colour = "black"),  # 设置图例标题颜色为黑色
      legend.position = "right")   # 设置图例位置为右侧

# 创建散点图并按照样本大小和OTU数量进行分面
scatter_plot1 <- ggplot(data1, aes(x = mean_FPR, y = mean_TPR, shape = name, color = name)) +
   geom_point(size = 3) +  # 绘制散点
   geom_abline(intercept = 0, slope = 1, color = "gray", linetype = "dashed") +  # 添加对角线
   labs(x = "False Positive Rate", y = "True Positive Rate", title = "Effect size = 1") +  # 添加坐标轴标签和标题
   theme_minimal() +  # 使用简洁主题
   coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +  # 设置横轴和纵轴范围
   theme(panel.border = element_rect(color = "black", fill = NA, size = 0.8),  # 设置外边框
      panel.spacing = unit(1, "lines"),  # 设置面板间距
      panel.background = element_rect(color = "gray50", fill = NA),  # 设置面板背景
      strip.background = element_blank(),  # 移除分面标题背景
      strip.text.x = element_text(size = 8, angle = 0, hjust = 0.5)) +  # 设置分面标题文本
   # facet_wrap(otu~sample)+theme_bw()+
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(Sample = label_both,OTU = label_both )) + theme_bw()+ # 按样本大小和OTU数量进行分面并命名
   guides(shape = guide_legend(label.position = "right", title = "model"),
      color = guide_legend(label.position = "right", title = "model")) + # 在图例中显示标签值
   scale_shape_manual(values = c(17,18,16,15,4,3,2,1,0,5,6,7,10,11)) + # 设置每种形状对应的值
   scale_color_manual(values = sample_colors) +  # 设置颜色
   theme(legend.text = element_text(colour = "black"),  # 设置图例文本颜色为黑色
      legend.title = element_text(colour = "black"),  # 设置图例标题颜色为黑色
      legend.position = "right")   # 设置图例位置为右侧

ggsave("FPR_TPR_effect10.pdf", plot = scatter_plot10, width =18, height = 8)
ggsave("FPR_TPR_effect5.pdf", plot = scatter_plot5, width =18, height = 8)
ggsave("FPR_TPR_effect2.pdf", plot = scatter_plot2, width =18, height = 8)
ggsave("FPR_TPR_effect1.pdf", plot = scatter_plot1, width =18, height = 8)

