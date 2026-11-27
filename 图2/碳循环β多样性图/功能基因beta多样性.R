# 1. 导入所需的依赖包
library(vegan)
library(ggplot2)
library(dplyr)

# 2. 读取数据 (请确保工作目录为 Carbon.txt 和 group.txt 所在的文件夹)
abund_data <- read.table("Carbon.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
group_data <- read.table("group.txt", header = TRUE, sep = "\t")

# 转置丰度矩阵 (让行是样本，列是基因)
abund_t <- t(abund_data)

# 确保分组表的样本顺序和丰度表一致
group_data <- group_data[match(rownames(abund_t), group_data$sample), ]

# 3. 计算 Bray-Curtis 距离矩阵
bc_dist <- vegdist(abund_t, method = "bray")

# 4. 主坐标分析 (PCoA)
pcoa_res <- cmdscale(bc_dist, k = 2, eig = TRUE)
pcoa_points <- as.data.frame(pcoa_res$points)
colnames(pcoa_points) <- c("PCoA1", "PCoA2")

# 计算各轴解释的方差百分比
eig_var <- pcoa_res$eig[pcoa_res$eig > 0]
explained_var <- eig_var / sum(eig_var)
pcoa1_exp <- round(explained_var[1] * 100, 1)
pcoa2_exp <- round(explained_var[2] * 100, 1)

# 将 PCoA 坐标和 metadata(分组信息) 合并
pcoa_df <- cbind(pcoa_points, group_data)
# 设定因子水平，以固定图例顺序
pcoa_df$group1 <- factor(pcoa_df$group1, levels = c("M-CK", "M-P", "EV-CK", "EV-P"))

# 5. 分别进行 PERMANOVA (Adonis) 统计检验
# 5.1 检验 M 组 (M-CK vs M-P)
idx_M <- pcoa_df$group2 == "M"
dist_M <- as.dist(as.matrix(bc_dist)[idx_M, idx_M])
adonis_M <- adonis2(dist_M ~ group1, data = pcoa_df[idx_M, ], permutations = 999)
R2_M <- adonis_M$R2[1]
P_M <- adonis_M$`Pr(>F)`[1]

# 5.2 检验 EV 组 (EV-CK vs EV-P)
idx_EV <- pcoa_df$group2 == "EV"
dist_EV <- as.dist(as.matrix(bc_dist)[idx_EV, idx_EV])
adonis_EV <- adonis2(dist_EV ~ group1, data = pcoa_df[idx_EV, ], permutations = 999)
R2_EV <- adonis_EV$R2[1]
P_EV <- adonis_EV$`Pr(>F)`[1]

# 准备添加到图中的统计结果文本
stats_label <- sprintf(
  "PERMANOVA (M-CK vs M-P):\nR² = %.3f, P = %.3f\n\nPERMANOVA (EV-CK vs EV-P):\nR² = %.3f, P = %.3f", 
  R2_M, P_M, R2_EV, P_EV
)

# 6. ES&T 期刊风格绘图 (ggplot2)
p <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, fill = group1, shape = group1)) +
  # 绘制置信椭圆 (0.866约对应Python里常用的1.5倍标准差)
  stat_ellipse(aes(color = group1), type = "norm", linetype = "dashed", linewidth = 1, level = 0.866, show.legend = FALSE) +
  # 绘制散点，设置带有黑边的实心点(shape 21和24支持填充色与边框色分离)
  geom_point(size = 4, color = "black", stroke = 1) +
  # 配置颜色和形状
  scale_fill_manual(values = c("M-CK" = "#4c97cc", "M-P" = "#f57c2e", "EV-CK" = "#4575b4", "EV-P" = "#d73027")) +
  scale_color_manual(values = c("M-CK" = "#4c97cc", "M-P" = "#f57c2e", "EV-CK" = "#4575b4", "EV-P" = "#d73027")) +
  scale_shape_manual(values = c("M-CK" = 22, "M-P" = 22, "EV-CK" = 21, "EV-P" = 21)) +
  # 坐标轴标题
  labs(
    x = paste0("PCoA 1 (", pcoa1_exp, "%)"),
    y = paste0("PCoA 2 (", pcoa2_exp, "%)"),
    fill = "Groups", shape = "Groups"
  ) +
  # 纯净的白板主题
  theme_bw() +
  theme(
    # 去除背景网格线
    panel.grid = element_blank(),
    # 字体设定：无衬线字体Arial/Helvetica风
    text = element_text(family = "sans"),
    # 加粗的坐标轴标题
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    # 坐标轴边框加粗
    panel.border = element_rect(color = "black", linewidth = 1.2),
    # 刻度线朝内(负值length)，顶端和右侧也显示刻度，符合顶级期刊要求
    axis.ticks = element_line(color = "black", linewidth = 1),
    axis.ticks.length = unit(-0.2, "cm"),
    axis.text.x = element_text(margin = margin(t = 8)),
    axis.text.y = element_text(margin = margin(r = 8)),
    # 图例设置在图内部右下侧
    legend.position = c(0.15, 0.2),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.title = element_text(face = "bold", size = 12)
  ) +
  # 嵌入统计框 (左上角位置可以通过x和y调节)
  annotate("label", x = min(pcoa_df$PCoA1), y = max(pcoa_df$PCoA2), 
           label = stats_label, hjust = 0, vjust = 1, size = 4, 
           fontface = "plain", fill = "white", color = "black", 
           label.padding = unit(0.5, "lines"), label.r = unit(0.15, "lines"))

# 7. 查看并导出高清图片
print(p)
ggsave("PCoA_EST_Style_R.png", plot = p, width = 8, height = 6, dpi = 300)

