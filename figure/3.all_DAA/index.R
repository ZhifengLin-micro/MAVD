###############################################
####input_data####
library(ggplot2)
library(gridExtra)
setwd("/public/linzhf/results")

input20 <- read.delim("/public/linzhf/results/input_Sample20.txt")
input20$effect <- as.numeric(ifelse(input20$effect %in% c("2", "5"), input20$effect + 2, input20$effect))
# 使用 factor 对 Samples 进行排序
input20$OTU <- factor(input20$OTU, levels = c("100", "200", "300", "400","500","800","1000", "2000"))

input50 <- read.delim("/public/linzhf/results/input_Sample50.txt")
input50$effect <- as.numeric(ifelse(input50$effect %in% c("2", "5"), input50$effect + 2, input50$effect))
# 使用 factor 对 Samples 进行排序
input50$OTU <- factor(input50$OTU, levels = c("100", "200", "300", "400","500","800","1000", "2000"))

input100 <- read.delim("/public/linzhf/results/input_Sample100.txt")
input100$effect <- as.numeric(ifelse(input100$effect %in% c("2", "5"), input100$effect + 2, input100$effect))
# 使用 factor 对 Samples 进行排序
input100$OTU <- factor(input100$OTU, levels = c("100", "200", "300", "400","500","800","1000", "2000"))

input200 <- read.delim("/public/linzhf/results/input_Sample200.txt")
input200$effect <- as.numeric(ifelse(input200$effect %in% c("2", "5"), input200$effect + 2, input200$effect))
# 使用 factor 对 Samples 进行排序
input200$OTU <- factor(input200$OTU, levels = c("100", "200", "300", "400","500","800","1000", "2000"))

input300 <- read.delim("/public/linzhf/results/input_Sample300.txt")
input300$effect <- as.numeric(ifelse(input300$effect %in% c("2", "5"), input300$effect + 2, input300$effect))
# 使用 factor 对 Samples 进行排序
input300$OTU <- factor(input300$OTU, levels = c("100", "200", "300", "400","500","800","1000", "2000"))

input500 <- read.delim("/public/linzhf/results/input_Sample500.txt")
input500$effect <- as.numeric(ifelse(input500$effect %in% c("2", "5"), input500$effect + 2, input500$effect))
# 使用 factor 对 Samples 进行排序
input500$OTU <- factor(input500$OTU, levels = c("100", "200", "300", "400","500","800","1000", "2000"))

Sample_colors <- c("#D52527", "#614099", "#088247","#CD5A1F", "#0075C5", "#8B5549","#FF7D0D", "#F0A298", "#7D7D7D","#98C5E1", "#8481ba", "#6BB952", "#FAEE85","#234444" )

####Accuracy####
# 绘制4个图
p1_20 <- ggplot(input20) +
   geom_line(aes(x = effect, y = mean_Accuracy, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_Accuracy, color = name), size = 1.6) +
   labs(y = "Accuracy") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank()  # 不显示横轴刻度
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p1_50 <- ggplot(input50) +
   geom_line(aes(x = effect, y = mean_Accuracy, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_Accuracy, color = name), size = 1.6) +
   labs(y = "Accuracy") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank(),  # 不显示横轴刻度
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p1_100 <- ggplot(input100) +
   geom_line(aes(x = effect, y = mean_Accuracy, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_Accuracy, color = name), size = 1.6) +
   labs(y = "Accuracy") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank(),  # 不显示横轴刻度
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p1_200 <- ggplot(input200) +
   geom_line(aes(x = effect, y = mean_Accuracy, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_Accuracy, color = name), size = 1.6) +
   labs(y = "Accuracy") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank(),  # 不显示横轴刻度
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p1_300 <- ggplot(input300) +
   geom_line(aes(x = effect, y = mean_Accuracy, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_Accuracy, color = name), size = 1.6) +
   labs(y = "Accuracy") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank(),  # 不显示横轴刻度
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p1_500 <- ggplot(input500) +
   geom_line(aes(x = effect, y = mean_Accuracy, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_Accuracy, color = name), size = 1.6) +
   labs(y = "Accuracy") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )


####Precision####
p2_20 <- ggplot(input20) +
   geom_line(aes(x = effect, y = mean_Precision, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_Precision, color = name), size = 1.6) +
   labs(y = "Precision") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank()  # 不显示横轴刻度
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p2_50 <- ggplot(input50) +
   geom_line(aes(x = effect, y = mean_Precision, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_Precision, color = name), size = 1.6) +
   labs(y = "Precision") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank(),  # 不显示横轴刻度
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p2_100 <- ggplot(input100) +
   geom_line(aes(x = effect, y = mean_Precision, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_Precision, color = name), size = 1.6) +
   labs(y = "Precision") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank(),  # 不显示横轴刻度
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p2_200 <- ggplot(input200) +
   geom_line(aes(x = effect, y = mean_Precision, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_Precision, color = name), size = 1.6) +
   labs(y = "Precision") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank(),  # 不显示横轴刻度
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p2_300 <- ggplot(input300) +
   geom_line(aes(x = effect, y = mean_Precision, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_Precision, color = name), size = 1.6) +
   labs(y = "Precision") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank(),  # 不显示横轴刻度
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p2_500 <- ggplot(input500) +
   geom_line(aes(x = effect, y = mean_Precision, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_Precision, color = name), size = 1.6) +
   labs(y = "Precision") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )

####TPR####
p3_20 <- ggplot(input20) +
   geom_line(aes(x = effect, y = mean_TPR , color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_TPR , color = name), size = 1.6) +
   labs(y = "TPR ") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank()  # 不显示横轴刻度
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p3_50 <- ggplot(input50) +
   geom_line(aes(x = effect, y = mean_TPR , color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_TPR , color = name), size = 1.6) +
   labs(y = "TPR ") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank(),  # 不显示横轴刻度
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p3_100 <- ggplot(input100) +
   geom_line(aes(x = effect, y = mean_TPR , color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_TPR , color = name), size = 1.6) +
   labs(y = "TPR ") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank(),  # 不显示横轴刻度
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p3_200 <- ggplot(input200) +
   geom_line(aes(x = effect, y = mean_TPR , color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_TPR , color = name), size = 1.6) +
   labs(y = "TPR ") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank(),  # 不显示横轴刻度
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p3_300 <- ggplot(input300) +
   geom_line(aes(x = effect, y = mean_TPR , color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_TPR , color = name), size = 1.6) +
   labs(y = "TPR ") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank(),  # 不显示横轴刻度
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )+
   labs(x = NULL)  # 去除 x 轴的名称


p3_500 <- ggplot(input500) +
   geom_line(aes(x = effect, y = mean_TPR , color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_TPR , color = name), size = 1.6) +
   labs(y = "TPR ") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )

####Specificity####
p4_20 <- ggplot(input20) +
   geom_line(aes(x = effect, y = mean_Specificity, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_Specificity, color = name), size = 1.6) +
   labs(y = "Specificity") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank()  # 不显示横轴刻度
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p4_50 <- ggplot(input50) +
   geom_line(aes(x = effect, y = mean_Specificity, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_Specificity, color = name), size = 1.6) +
   labs(y = "Specificity") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank(),  # 不显示横轴刻度
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p4_100 <- ggplot(input100) +
   geom_line(aes(x = effect, y = mean_Specificity, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_Specificity, color = name), size = 1.6) +
   labs(y = "Specificity") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank(),  # 不显示横轴刻度
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p4_200 <- ggplot(input200) +
   geom_line(aes(x = effect, y = mean_Specificity, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_Specificity, color = name), size = 1.6) +
   labs(y = "Specificity") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank(),  # 不显示横轴刻度
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p4_300 <- ggplot(input300) +
   geom_line(aes(x = effect, y = mean_Specificity, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_Specificity, color = name), size = 1.6) +
   labs(y = "Specificity") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank(),  # 不显示横轴刻度
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p4_500 <- ggplot(input500) +
   geom_line(aes(x = effect, y = mean_Specificity, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_Specificity, color = name), size = 1.6) +
   labs(y = "Specificity") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )


####FPR####
p5_20 <- ggplot(input20) +
   geom_line(aes(x = effect, y = mean_FPR, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_FPR, color = name), size = 1.6) +
   labs(y = "FPR") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank()  # 不显示横轴刻度
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p5_50 <- ggplot(input50) +
   geom_line(aes(x = effect, y = mean_FPR, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_FPR, color = name), size = 1.6) +
   labs(y = "FPR") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank(),  # 不显示横轴刻度
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p5_100 <- ggplot(input100) +
   geom_line(aes(x = effect, y = mean_FPR, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_FPR, color = name), size = 1.6) +
   labs(y = "FPR") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank(),  # 不显示横轴刻度
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p5_200 <- ggplot(input200) +
   geom_line(aes(x = effect, y = mean_FPR, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_FPR, color = name), size = 1.6) +
   labs(y = "FPR") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank(),  # 不显示横轴刻度
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p5_300 <- ggplot(input300) +
   geom_line(aes(x = effect, y = mean_FPR, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_FPR, color = name), size = 1.6) +
   labs(y = "FPR") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank(),  # 不显示横轴刻度
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p5_500 <- ggplot(input500) +
   geom_line(aes(x = effect, y = mean_FPR, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_FPR, color = name), size = 1.6) +
   labs(y = "FPR") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )

####FDR####
p6_20 <- ggplot(input20) +
   geom_line(aes(x = effect, y = mean_FDR, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_FDR, color = name), size = 1.6) +
   labs(y = "FDR") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank()  # 不显示横轴刻度
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p6_50 <- ggplot(input50) +
   geom_line(aes(x = effect, y = mean_FDR, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_FDR, color = name), size = 1.6) +
   labs(y = "FDR") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank(),  # 不显示横轴刻度
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p6_100 <- ggplot(input100) +
   geom_line(aes(x = effect, y = mean_FDR, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_FDR, color = name), size = 1.6) +
   labs(y = "FDR") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank(),  # 不显示横轴刻度
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p6_200 <- ggplot(input200) +
   geom_line(aes(x = effect, y = mean_FDR, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_FDR, color = name), size = 1.6) +
   labs(y = "FDR") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank(),  # 不显示横轴刻度
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )+
   labs(x = NULL)  # 去除 x 轴的名称


p6_300 <- ggplot(input300) +
   geom_line(aes(x = effect, y = mean_FDR, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_FDR, color = name), size = 1.6) +
   labs(y = "FDR") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank(),  # 不显示横轴刻度
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p6_500 <- ggplot(input500) +
   geom_line(aes(x = effect, y = mean_FDR, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_FDR, color = name), size = 1.6) +
   labs(y = "FDR") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )


####F1_score####
p7_20 <- ggplot(input20) +
   geom_line(aes(x = effect, y = mean_F1_score, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_F1_score, color = name), size = 1.6) +
   labs(y = "F1 score") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank()  # 不显示横轴刻度
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p7_50 <- ggplot(input50) +
   geom_line(aes(x = effect, y = mean_F1_score, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_F1_score, color = name), size = 1.6) +
   labs(y = "F1 score") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank(),  # 不显示横轴刻度
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p7_100 <- ggplot(input100) +
   geom_line(aes(x = effect, y = mean_F1_score, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_F1_score, color = name), size = 1.6) +
   labs(y = "F1 score") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank(),  # 不显示横轴刻度
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p7_200 <- ggplot(input200) +
   geom_line(aes(x = effect, y = mean_F1_score, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_F1_score, color = name), size = 1.6) +
   labs(y = "F1 score") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank(),  # 不显示横轴刻度
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p7_300 <- ggplot(input300) +
   geom_line(aes(x = effect, y = mean_F1_score, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_F1_score, color = name), size = 1.6) +
   labs(y = "F1 score") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      axis.text.x = element_blank(),  # 不显示横轴标签值
      axis.ticks.x = element_blank(),  # 不显示横轴刻度
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )+
   labs(x = NULL)  # 去除 x 轴的名称

p7_500 <- ggplot(input500) +
   geom_line(aes(x = effect, y = mean_F1_score, color = name), size = 1.2) +
   geom_point(aes(x = effect, y = mean_F1_score, color = name), size = 1.6) +
   labs(y = "F1 score") +
   scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c(1, 2, 5, 10)) +
   scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1)) +
   scale_color_manual(values = Sample_colors) +  # 设置颜色
   facet_grid(Sample ~ OTU, scales = "free", space = "free", labeller = labeller(OTU = label_both,Sample = label_both )) +  # Group by "name"
   theme_bw() +
   theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),  # 设置外边框
      strip.text.x = element_blank()  # 不显示分面标签（名称）
   )


#################################
####recall####
combined_plot1 <- grid.arrange(p1_20, p1_50, p1_100, p1_200, p1_300, p1_500, ncol = 1)  # You can adjust the number of columns as needed
ggsave("Accuracy_combined.pdf", plot = combined_plot1, width = 16, height = 12)
####precision####
combined_plot2 <- grid.arrange(p2_20, p2_50, p2_100, p2_200, p2_300, p2_500, ncol = 1)  # You can adjust the number of columns as needed
ggsave("Precision_combined.pdf", plot = combined_plot2, width = 16, height = 12)
####fdr####
combined_plot3 <- grid.arrange(p3_20, p3_50, p3_100, p3_200, p3_300, p3_500, ncol = 1)  # You can adjust the number of columns as needed
ggsave("TPR_combined.pdf", plot = combined_plot3, width = 16, height = 12)
####f1####
combined_plot4 <- grid.arrange(p4_20, p4_50, p4_100, p4_200, p4_300, p4_500, ncol = 1)  # You can adjust the number of columns as needed
ggsave("Specificity_combined.pdf", plot = combined_plot4, width = 16, height = 12)
####precision####
combined_plot5 <- grid.arrange(p5_20, p5_50, p5_100, p5_200, p5_300, p5_500, ncol = 1)  # You can adjust the number of columns as needed
ggsave("FPR_combined.pdf", plot = combined_plot5, width = 16, height = 12)
####fdr####
combined_plot6 <- grid.arrange(p6_20, p6_50, p6_100, p6_200, p6_300, p6_500, ncol = 1)  # You can adjust the number of columns as needed
ggsave("FDR_combined.pdf", plot = combined_plot6, width = 16, height = 12)
####f1####
combined_plot7 <- grid.arrange(p7_20, p7_50, p7_100, p7_200, p7_300, p7_500, ncol = 1)  # You can adjust the number of columns as needed
ggsave("F1_score_combined.pdf", plot = combined_plot7, width = 16, height = 12)

