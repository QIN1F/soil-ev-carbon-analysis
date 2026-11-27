# ==============================================================================
# 脚本名称: LH4_M_Heatmap.R
# 功能描述: 绘制 MCOA 1-4 轴与细菌组 (M) 生活史策略的 Spearman 相关性热图
# ==============================================================================

# 清除环境
rm(list=ls()) 

# 加载必要的包
library(readxl)  # 用于读取 excel 文件
library(Hmisc)   # 用于计算 spearman 相关性及 P 值
library(dplyr)   # 用于数据处理
library(tidyr)   # 用于数据长宽转换
library(ggplot2) # 用于画图

# ==================== 1. 读取并准备数据 ====================
# 设置你的本地工作路径 (已修改为 lifehistory-M)
data_dir <- "E:/我的日常文件/111微塑料循环汇报/4.22生活史策略分析/我的/gene_base/lifehistory-M/"
file_path <- paste0(data_dir, "mcoa_M_samples_4axes.xlsx")

if(!file.exists(file_path)) {
  stop("错误：找不到文件 mcoa_M_samples_4axes.xlsx！请检查上一步是否成功生成。")
}

# 读取 Excel 数据
data_raw <- as.data.frame(read_excel(file_path))

# 将第一列 (name) 设置为行名，并从数据框中删除该列
rownames(data_raw) <- data_raw[, 1]
data_4axes <- data_raw[, -1]

# 自动识别列名：区分 MCOA 轴和生活史特征 (Traits)
mcoa_cols <- c("MCOA1", "MCOA2", "MCOA3", "MCOA4")
# 剩下的所有列都认为是特征列
trait_cols <- setdiff(colnames(data_4axes), mcoa_cols) 

# ==================== 2. 计算 Spearman 相关性 ====================
# rcorr 需要输入 matrix 格式
sp_res <- rcorr(as.matrix(data_4axes), type = "spearman")

# 提取 Trait (行) 与 MCOA (列) 的交叉区域
r_matrix <- sp_res$r[trait_cols, mcoa_cols]
p_matrix <- sp_res$P[trait_cols, mcoa_cols]

# ==================== 3. 数据清洗与长格式转换 ====================
# 将相关系数 (r) 矩阵转为长格式
r_long <- as.data.frame(as.table(r_matrix)) %>%
  rename(Trait = Var1, MCOA = Var2, r = Freq)

# 将显著性 (P) 矩阵转为长格式
p_long <- as.data.frame(as.table(p_matrix)) %>%
  rename(Trait = Var1, MCOA = Var2, P = Freq)

# 合并 r 和 P 数据
plot_df <- left_join(r_long, p_long, by = c("Trait", "MCOA"))

# 生成显著性星号 (*: P<0.05, **: P<0.01, ***: P<0.001)
plot_df <- plot_df %>%
  mutate(
    sig = case_when(
      P < 0.001 ~ "***",
      P < 0.01  ~ "**",
      P < 0.05  ~ "*",
      TRUE      ~ ""  
    )
  )

# 导出详细的统计结果表，方便你写文章时查阅具体数值
write.csv(plot_df, paste0(data_dir, "MCOA_Traits_Correlation_Results_M.csv"), row.names = FALSE)

# ==================== 4. 绘制热图 ====================
# 恢复原始的纵坐标顺序（默认由下至上为 MCOA1 -> MCOA4）
plot_df$MCOA <- factor(plot_df$MCOA, levels = paste0("MCOA", 1:4))

p_heatmap <- ggplot(plot_df, aes(x = Trait, y = MCOA)) +
  # 绘制色块，设置白色网格线使其更清爽
  geom_tile(aes(fill = r), color = "white", width = 0.95, height = 0.95, linewidth = 0.8) +
  # 添加显著性星号，调整位置居中偏上
  geom_text(aes(label = sig), color = "black", size = 4.5, vjust = 0.75) +
  # 设置颜色映射：负相关偏青蓝，正相关偏红橙
  scale_fill_gradient2(
    low = "#81C6C5",       
    mid = "white",         
    high = "#F08D79",      
    midpoint = 0, 
    limit = c(-1, 1),      
    name = "r",            
    breaks = c(-0.5, 0, 0.5) 
  ) +
  # 【新增】添加图表大标题
  labs(title = "M组的MCOA4轴相关性热图") + 
  theme_minimal() +
  theme(
    # 【新增】设置标题居中、加粗，并与下方图表留出间距
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    panel.grid = element_blank(),           
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10, color = "black"), 
    axis.text.y = element_text(size = 12, face = "bold", color = "black"), 
    axis.title = element_blank(),           
    legend.position = "right",              
    legend.title = element_text(face = "bold")
  ) +
  # 强制色块为正方形，保证美观度
  coord_fixed(ratio = 1)

# 在 RStudio 中预览图片
print(p_heatmap)

# ==================== 5. 导出图片 ====================
# 【修改】高度略微调大到 4，给上方的标题留出足够的空间
ggsave(paste0(data_dir, "Fig2d_MCOA_Heatmap_M.pdf"), plot = p_heatmap, width = 12, height = 4)

message("细菌组热图绘制完成！文件已保存在: ", data_dir)

