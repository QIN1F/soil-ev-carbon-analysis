# 加载必要的包
install.packages("ggbreak")
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbreak)  # 新增：用于实现Y轴分段（截断）展示

# 1. 读取数据 (如果你跑种水平，请把 phylum.csv 换成 species.csv)
phylum <- read.csv("OTU.csv", row.names = 1, check.names = FALSE)
group <- read.table("group.txt", header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)

# 2. 数据处理：提取 Top 10 (如果是种水平可以改成 Top 15 或 Top 20)，计算相对丰度
phylum_rel <- as.data.frame(apply(phylum, 2, function(x) x/sum(x) * 100))
phylum_rel$mean_abund <- rowMeans(phylum_rel)
phylum_rel <- phylum_rel[order(phylum_rel$mean_abund, decreasing = TRUE), ]
top_10 <- rownames(phylum_rel)[1:10]

plot_data <- phylum_rel[top_10, -ncol(phylum_rel)]
plot_data["Others", ] <- 100 - colSums(plot_data)
plot_data$Phylum <- rownames(plot_data)

plot_long <- gather(plot_data, key = "Sample", value = "Abundance", -Phylum)
plot_long$Group <- group[plot_long$Sample, "group1"]

# 排序
plot_long$Phylum <- factor(plot_long$Phylum, levels = rev(c(top_10, "Others")))
plot_long$Group <- factor(plot_long$Group, levels = c("M-CK", "M-P", "EV-CK", "EV-P"))
sample_order <- group %>% arrange(factor(group1, levels = c("M-CK", "M-P", "EV-CK", "EV-P"))) %>% rownames()
plot_long$Sample <- factor(plot_long$Sample, levels = sample_order)

# 3. 自定义强主色弱辅色色卡
light_colors <- c("#cde3f2", "#fce1cd", "#e1f0c4", "#f2dcf5", "#f7f0c1", "#d6ebd1", "#f5d3d3", "#e1e6f2")
custom_palette <- c("#ebebeb", rev(c("#4c97cc", "#f57c2e", light_colors)))

# 4. 绘制基础图形
p <- ggplot(plot_long, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", width = 0.9, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = custom_palette) +
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  labs(x = "Samples", y = "Relative Abundance (%)", title = "Community Composition (Broken Y-axis)") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, color = "black"),
        strip.text = element_text(face = "bold", size = 12),
        strip.background = element_rect(fill = "#f0f0f0", color = "black"),
        panel.grid = element_blank(),
        legend.position = "right")

# 5. 添加 Y 轴截断 (分段式展示核心代码)
# 这里的 c(30, 85) 表示把Y轴 30% 到 85% 的区间折叠隐藏
# scales = 0.3 表示让上半截（被压缩后的灰色 Others 区域）的高度，只有下半截（核心差异区域）的 30%
p_broken <- p + scale_y_break(c(30, 85), scales = 0.3)

# 6. 显式出图与保存
print(p_broken) # 确保图片在 RStudio 的 Plots 窗口中弹出来
ggsave("Composition_BrokenY_LightSupport.pdf", plot = p_broken, width = 12, height = 7)
      