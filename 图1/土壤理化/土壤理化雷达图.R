# 安装包（首次运行时取消注释）
# install.packages(c("readxl", "fmsb", "dplyr"))

library(readxl)
library(fmsb)
library(dplyr)

# 1. 读取数据
data <- read_excel("土壤理化.xlsx", sheet = "Sheet1")
colnames(data) <- c("Sample", "Group", "SOC", "TN", "TP", "TK", "PH", "AN", "AP", "AK", "SMC")

# 2. 确保数值型
num_cols <- c("SOC", "TN", "TP", "TK", "PH", "AN", "AP", "AK", "SMC")
data[num_cols] <- lapply(data[num_cols], as.numeric)
data <- as.data.frame(data)

# 3. Wilcoxon 秩和检验（得到显著性星号）
results <- data.frame()
for (col in num_cols) {
  ck <- data[data$Group == "CK", ][[col]]
  d1 <- data[data$Group == "D1", ][[col]]
  
  test <- wilcox.test(ck, d1, alternative = "two.sided", exact = TRUE)
  
  sig <- ifelse(test$p.value < 0.001, "***", 
                ifelse(test$p.value < 0.01, "**", 
                       ifelse(test$p.value < 0.05, "*", "")))   # ns不显示星号
  
  results <- rbind(results, data.frame(
    指标 = col,
    显著性 = sig
  ))
}

# 4. 生成带星号的轴标签（与示例图完全一致）
results$custom_label <- ifelse(results$显著性 == "", 
                               results$指标, 
                               paste0(results$指标, results$显著性))
custom_vlabs <- results$custom_label

# 5. 计算各组原始均值
group_means <- aggregate(data[, num_cols], by = list(Group = data$Group), FUN = mean)
rownames(group_means) <- group_means$Group
group_means <- group_means[, -1]

# 6. 准备雷达图数据（原始均值 + 实际min/max）
min_vals <- apply(data[, num_cols], 2, min)
max_vals <- apply(data[, num_cols], 2, max)
radar_data <- rbind(max_vals, min_vals, group_means)

# 把列名换成带星号的标签
colnames(radar_data) <- custom_vlabs

# 7. 绘制雷达图 —— 关键修改：删除0%~100%（或任何刻度数字）
colors_line <- c(M0 = "#4c97cc", M1 = "#f57c2e")
colors_fill <- adjustcolor(colors_line, alpha.f = 0.25)

radarchart(radar_data,
           axistype = 1,                  # 保留网格线
           axislabcol = "transparent",    # ←←← 这里！彻底隐藏0%~100%刻度数字
           pcol = colors_line,
           pfcol = colors_fill,
           plty = 1, plwd = 3,
           cglty = 1, cglwd = 1,          # 保留蛛网网格
           cglcol = "grey",
           vlcex = 1.1,                   # 变量标签（带星号）
           title = " ",
           centerzero = TRUE)

# 8. 添加图例
legend("topright", 
       legend = c("M-CK", "M-P"),
       col = colors_line,
       lty = 1, lwd = 3, bty = "n", cex = 1.2)

