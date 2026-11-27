# 加载必要的 R 包
# 如果尚未安装，请先运行: install.packages(c("ggplot2", "dplyr", "ggpubr"))
library(ggplot2)
library(dplyr)
library(ggpubr)

# ==========================================
# 1. 数据读取与预处理
# ==========================================
# 读取丰度数据 (请确保 Carbon.txt 和 group.txt 在当前工作目录)
df_abund <- read.table("Carbon.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
df_group <- read.table("group.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 统一大小写，防止样本名不匹配 (如 p1_EV 与 P1_EV)
colnames(df_abund) <- toupper(colnames(df_abund))
df_group$sample <- toupper(df_group$sample)

# 统一处理组名称：将 D1 替换为 P，以便于与 CK 对比
df_group$treat[df_group$treat == "D1"] <- "P"

# ==========================================
# 2. 计算 Richness 指数
# ==========================================
# 统计每个样本中丰度大于 0 的基因数量
richness_values <- colSums(df_abund > 0)

df_alpha <- data.frame(
  sample = names(richness_values),
  Richness = as.numeric(richness_values)
)

# 合并多样性数据与分组数据
df_merged <- merge(df_alpha, df_group, by = "sample")

# 设定因子的水平（决定画图中从左到右的顺序）
df_merged$group2 <- factor(df_merged$group2, levels = c("M", "EV"))
df_merged$treat <- factor(df_merged$treat, levels = c("CK", "P"))

# ==========================================
# 3. 绘制 EST / WR 风格图表
# ==========================================
# 设定自然系配色：CK为青蓝色，P为陶土红
my_colors <- c("CK" = "#5C8E99", "P" = "#D47E5B")

# 构建绘图对象 p
p <- ggboxplot(
  data = df_merged, 
  x = "group2", 
  y = "Richness",
  color = "black",            # 箱线图边框颜色
  fill = "treat",             # 按处理组填充颜色
  width = 0.6,                # 箱体宽度
  outlier.shape = NA,         # 隐藏默认的异常值点（用下面的jitter代替）
  add = "jitter",             # 添加抖动散点
  add.params = list(color = "black", size = 2, alpha = 0.6, jitter = 0.2)
) +
  # 设置颜色
  scale_fill_manual(values = my_colors) +
  
  # 执行独立样本 T 检验，自动添加显著性标记
  # label = "p.signif" 会显示 ns, *, **, ***
  stat_compare_means(
    aes(group = treat), 
    method = "t.test", 
    label = "p.signif", 
    hide.ns = FALSE, 
    vjust = 0.5
  ) +
  
  # 修改坐标轴标签
  labs(x = "", y = "Richness", fill = "Treatment Group") +
  
  # EST/WR 风格主题设置
  theme_classic() +
  theme(
    text = element_text(family = "sans", size = 14),           # 字体设置
    axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
    axis.text = element_text(color = "black", size = 12),
    
    # 加粗坐标轴及边框，设置全封闭矩形框
    axis.line = element_blank(),                               # 移除经典主题默认的线
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    
    # 刻度线加粗且朝内 (向图表内部指)
    axis.ticks = element_line(linewidth = 1.2, color = "black"),
    axis.ticks.length = unit(-0.2, "cm"),
    axis.text.x = element_text(margin = margin(t = 10)),
    axis.text.y = element_text(margin = margin(r = 10)),
    
    # 图例设置在底部
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 13),
    legend.text = element_text(size = 12)
  )

# 如果需要在本地运行后保存，请取消下面两行的注释：
# ggsave("Richness_EST_style.png", plot = p, width = 6, height = 5, dpi = 300)
# ggsave("Richness_EST_style.pdf", plot = p, width = 6, height = 5)

# 您可以直接在 R 控制台中输入 `p` 来查看生成的图像。