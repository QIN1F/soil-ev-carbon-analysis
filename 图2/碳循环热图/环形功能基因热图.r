# ==========================================
# 11个特定功能基因的环形热图绘制脚本
# ==========================================

# 0. 加载必要的包 (如未安装请取消注释下一行)
# install.packages(c("ggplot2", "dplyr", "tidyr", "patchwork"))
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# 1. 读取数据 (使用您提供的文件名)
df_raw <- read.table("Carbon.txt", header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
df_group <- read.table("group.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 清除基因名中可能存在的双引号，确保后续精准匹配
df_raw$gene_name <- gsub("\"", "", df_raw$gene_name)

# 定义您所需要的 11 个特定功能基因及其顺序
target_genes <- c(
  "fdoH",
  "aminotransferase class I and II",
  "bsdC",
  "mannan endo-1,4-beta-mannosidase",
  "pullulanase",
  "glucoamylase",
  "cellulase",
  "ack",
  "adh",
  "acdA",
  "fdhB"
)

# 在原始数据中仅过滤出这 11 个基因，并严格按照指定顺序重新排序
df_raw <- df_raw %>% filter(gene_name %in% target_genes)
df_raw <- df_raw[match(target_genes, df_raw$gene_name), ]

# 根据 group_2.txt 获取每个组对应的样本列名
cols_MCK <- df_group$sample[df_group$group1 == "M-CK"]
cols_MP  <- df_group$sample[df_group$group1 == "M-P"]
cols_EVCK <- df_group$sample[df_group$group1 == "EV-CK"]
cols_EVP  <- df_group$sample[df_group$group1 == "EV-P"]

# 2. 计算均值和显著性检验 (Welch's t-test)
mean_df <- data.frame(row.names = target_genes)
sig_MP <- character(length(target_genes))
sig_EVP <- character(length(target_genes))

for (i in seq_along(target_genes)) {
  v_MCK <- as.numeric(df_raw[i, cols_MCK])
  v_MP  <- as.numeric(df_raw[i, cols_MP])
  v_EVCK <- as.numeric(df_raw[i, cols_EVCK])
  v_EVP  <- as.numeric(df_raw[i, cols_EVP])
  
  mean_df[target_genes[i], "M-CK"] <- mean(v_MCK, na.rm = TRUE)
  mean_df[target_genes[i], "M-P"]  <- mean(v_MP, na.rm = TRUE)
  mean_df[target_genes[i], "EV-CK"] <- mean(v_EVCK, na.rm = TRUE)
  mean_df[target_genes[i], "EV-P"]  <- mean(v_EVP, na.rm = TRUE)
  
  # 计算显著性，增加判断防止因部分方差为0报错
  p_m <- if(var(v_MCK, na.rm=TRUE) + var(v_MP, na.rm=TRUE) > 0) t.test(v_MCK, v_MP, var.equal = FALSE)$p.value else 1
  p_ev <- if(var(v_EVCK, na.rm=TRUE) + var(v_EVP, na.rm=TRUE) > 0) t.test(v_EVCK, v_EVP, var.equal = FALSE)$p.value else 1
  
  sig_MP[i] <- ifelse(p_m < 0.05, "*", "")
  sig_EVP[i] <- ifelse(p_ev < 0.05, "*", "")
}

# 3. 数据标准化 (按行计算 Z-score)
z_df <- as.data.frame(t(apply(mean_df, 1, scale)))
colnames(z_df) <- c("M-CK", "M-P", "EV-CK", "EV-P")
z_df$Feature <- rownames(z_df)

# [排版核心]：添加4个空的占位符，用来在右侧强制生成一个圆环缺口 (C字形)
dummy_genes <- paste0("dummy_", 1:4)
all_genes <- c(target_genes, dummy_genes)
z_df$Feature <- factor(z_df$Feature, levels = all_genes) 

# 转为长格式以便 ggplot 绘图
df_long <- z_df %>%
  pivot_longer(cols = c("M-CK", "M-P", "EV-CK", "EV-P"), names_to = "Sample", values_to = "Value")
sample_order <- c("M-CK", "M-P", "EV-CK", "EV-P")
df_long$Sample <- factor(df_long$Sample, levels = sample_order)

# 合并显著性标记 (*)
sig_df <- data.frame(
  Feature = rep(target_genes, 2),
  Sample = c(rep("M-P", length(target_genes)), rep("EV-P", length(target_genes))),
  sig = c(sig_MP, sig_EVP)
)
df_long <- df_long %>% left_join(sig_df, by = c("Feature", "Sample"))
df_long$sig[is.na(df_long$sig)] <- ""

# 4. 计算外圈标签的精确旋转角度
N <- length(all_genes)
start_angle <- -42 # 设定起始偏移，使得缺口恰好落在右侧 (3点钟方向)

df_labels <- data.frame(
  Feature = target_genes,
  i = 1:length(target_genes)
) %>%
  mutate(
    # 计算每个标签在屏幕上的绝对角度 (逆时针排布)
    screen_angle = (90 + start_angle + 360 * (i - 0.5) / N) %% 360,
    # 如果标签在左半边 (90到270度)，文字需要靠右对齐并翻转180度，防止文字倒置
    hjust = case_when(screen_angle > 90 & screen_angle < 270 ~ 1, TRUE ~ 0),
    label_angle = case_when(screen_angle > 90 & screen_angle < 270 ~ screen_angle + 180, TRUE ~ screen_angle),
    Feature_Label = as.character(Feature)
  )

# 5. 内部组别标签 (用于填补到缺口处)
df_samples <- data.frame(
  Sample = sample_order,
  x = 1,
  y = sample_order
)
p_samples <- ggplot(df_samples, aes(x = x, y = y, label = Sample)) +
  geom_text(hjust = 0, size = 3.5, color = "black") +
  scale_y_discrete(limits = sample_order) +
  theme_void() + 
  theme(plot.margin = margin(0,0,0,0))

# 6. 绘制环形热图主图
p_circ <- ggplot(df_long, aes(x = Feature, y = Sample, fill = Value)) +
  geom_tile(color = "white", linewidth = 0.5) + 
  geom_text(aes(label = sig), size = 6, vjust = 0.75, color="black", fontface="bold") +
  
  # direction = -1 逆时针排列，与标签角度计算保持严格对应
  coord_polar(theta = "x", start = start_angle * pi / 180, direction = -1) + 
  
  scale_fill_gradientn(colors = c("#78C2B3", "#FFFFFF", "#FF8566"),
                       name = "Z-Score",
                       guide = guide_colorbar(title.theme = element_text(face = "bold", size = 10),
                                              frame.colour = "black", ticks.colour = "black")) +
  
  geom_text(data = df_labels, 
            aes(x = Feature, y = 4.7, 
                label = Feature_Label, angle = label_angle, hjust = hjust),
            inherit.aes = FALSE, size = 4) +
  
  # limits 设置确保占位符存在，扩展 y 轴制造内环空心效果
  scale_x_discrete(limits = all_genes) +
  scale_y_discrete(limits = sample_order, expand = expansion(add = c(2, 3))) +
  
  theme_void() +
  theme(
    legend.position = "right",
    plot.margin = margin(40, 40, 40, 40)
  )

# 7. 使用 patchwork 将组别标签准确插入到右侧空缺处
final_plot <- p_circ + 
  inset_element(p_samples, 
                left = 0.65, bottom = 0.45, right = 0.85, top = 0.65, 
                align_to = 'plot')

# 显示图像
print(final_plot)

# 保存图像 (您可以根据需要取消注释)
# ggsave("Filtered_Circular_Heatmap.png", final_plot, width = 9, height = 9, dpi = 300, bg = "white")