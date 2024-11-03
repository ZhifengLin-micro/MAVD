library(ggplot2)
input500 <- read.delim("E:/2023biyeketi/图片/4.time/input100.txt")
otu <- input100
time <- input100$time
name <- input100$name
# 设置颜色向量
sample_colors <- c( "#614099", "#088247","#CD5A1F", "#0075C5", "#8B5549","#FF7D0D", "#F0A298", "#7D7D7D","#98C5E1", "#8481ba", "#6BB952", "#FAEE85","#234444","#D52527" )

# 自定义x轴标签
custom_labels <- c("100", "200", "500", "1000", "2000")

# 添加水平线
hlines <- data.frame(yintercept = c(-1, 1))

# 设置刻度和网格线的样式以及背景为透明
custom_theme <- theme(
   axis.ticks = element_line(size = 0.5),  # 刻度线的样式
   panel.grid.major = element_line(color = "grey"),  # 设置主要网格线的颜色
   panel.grid.minor = element_line(color = "grey", size = 0.5),  # 设置次要网格线的样式和颜色
   panel.grid.major.y = element_line(color = "grey", size = 0.5),  # 设置y轴1和-1位置的内网格线
   panel.background = element_blank(),  # 设置背景为透明
   panel.border = element_rect(color = "black", fill = NA, size = 0.5)  # 添加外边框
)

# 使用ggplot2创建折线图，并手动设置颜色、折线粗细和形状
p <- ggplot(input100, aes(x = otu, y = time, group = name, color = name, shape = name)) +
   geom_line(size = 1) +
   geom_point(size = 3) +  # 添加点层
   geom_hline(data = hlines, aes(yintercept = yintercept), color = "grey") +  # 添加水平线
   scale_color_manual(values = sample_colors) +
   scale_shape_manual(values = 1:length(unique(input100$name))) +  # 设置形状
   labs(title = "Sample: 500  Effect size: 10", x = "OTU", y = "Time") +
   custom_theme +
   scale_x_continuous(labels = custom_labels) +
   labs(x = "OTU")  # 将x轴标签修改为自定义的标签
ggsave("input100.pdf", plot = p, device = "pdf", width = 10, height =8)

####500###########

input500 <- read.delim("E:/2023biyeketi/图片/4.time/input500.txt")
otu <- input500
time <- input500$time
name <- input500$name
# 设置颜色向量
sample_colors <- c( "#614099", "#088247","#CD5A1F", "#0075C5", "#8B5549","#FF7D0D", "#F0A298", "#7D7D7D","#98C5E1", "#8481ba", "#6BB952", "#FAEE85","#234444","#D52527" )

# 自定义x轴标签
custom_labels <- c("100", "200", "500", "1000", "2000")

# 添加水平线
hlines <- data.frame(yintercept = c(-1, 1))

# 设置刻度和网格线的样式以及背景为透明
custom_theme <- theme(
   axis.ticks = element_line(size = 0.5),  # 刻度线的样式
   panel.grid.major = element_line(color = "grey"),  # 设置主要网格线的颜色
   panel.grid.minor = element_line(color = "grey", size = 0.5),  # 设置次要网格线的样式和颜色
   panel.grid.major.y = element_line(color = "grey", size = 0.5),  # 设置y轴1和-1位置的内网格线
   panel.background = element_blank(),  # 设置背景为透明
   panel.border = element_rect(color = "black", fill = NA, size = 0.5)  # 添加外边框
)

# 使用ggplot2创建折线图，并手动设置颜色、折线粗细和形状
p <- ggplot(input500, aes(x = otu, y = time, group = name, color = name, shape = name)) +
   geom_line(size = 1) +
   geom_point(size = 3) +  # 添加点层
   geom_hline(data = hlines, aes(yintercept = yintercept), color = "grey") +  # 添加水平线
   scale_color_manual(values = sample_colors) +
   scale_shape_manual(values = 1:length(unique(input100$name))) +  # 设置形状
   labs(title = "Sample: 500  Effect size: 10", x = "OTU", y = "Time") +
   custom_theme +
   scale_x_continuous(labels = custom_labels) +
   labs(x = "OTU")  # 将x轴标签修改为自定义的标签
ggsave("input500.pdf", plot = p, device = "pdf", width = 10, height =8)


