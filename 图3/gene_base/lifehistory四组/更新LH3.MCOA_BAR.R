# install.packages(c("Hmisc", "ggpubr", "tidyverse", "readr")) 
library(ggplot2)
library(dplyr)
library(ggpubr)
library(readr)
library(tidyverse)
library(grid)
library(Hmisc) 

# ==================== 1. 读取数据 ====================
mcoa_file <- "lifehistory/mcoa_240_1.csv"
if(!file.exists(mcoa_file)) stop("错误：找不到 mcoa_240_1.csv")

mcoa_240_1 <- read.csv(mcoa_file, row.names=1, check.names = F)

# ==================== 2. 相关性计算与筛选 ====================
mcoa_240_1.sp <- rcorr(as.matrix(mcoa_240_1), type="spearman")
mcoa_r <- mcoa_240_1.sp$r[-c(1,2), c(1,2)]
mcoa_p <- mcoa_240_1.sp$P[-c(1,2), c(1,2)]

# 保存结果
mcoa_r_05 <- mcoa_r
mcoa_r_05[mcoa_p > 0.05] <- 0
write.csv(mcoa_r_05, "lifehistory/life0.05.csv")

mcoa_r_01 <- mcoa_r
mcoa_r_01[mcoa_p > 0.01] <- 0
write.csv(mcoa_r_01, "lifehistory/life0.01.csv")

# ==================== 3. 绘图函数 (恢复自动排序 + 强制对齐) ====================
draw_mcoa_bar_v6 <- function(data_file, axis_name, plot_title, target_count = NULL) {
  data_plot <- read.csv(data_file)
  colnames(data_plot)[1] <- "life"
  
  # 1. 筛选并按照相关性 r 值排序 (恢复之前的排序逻辑)
  real_df <- data_plot %>%
    filter(!!sym(axis_name) != 0) %>%
    arrange(!!sym(axis_name)) %>% # 按照数值从小到大排序，coord_flip后会变成从大到小
    mutate(Correlation = ifelse(!!sym(axis_name) > 0, "Positive (r > 0)", "Negative (r < 0)"))
  
  if(nrow(real_df) == 0) return(NULL)
  
  # 2. 补齐空白行逻辑 (确保 MCOA2 的柱子和 MCOA1 一样细)
  if(!is.null(target_count) && nrow(real_df) < target_count) {
    needed <- target_count - nrow(real_df)
    for(i in 1:needed) {
      blank_label <- paste(rep(" ", i), collapse="") # 生成不可见的空白标签
      # 将空白行插入到列表的最前面，这样它们会出现在图的底部占位
      real_df <- rbind(data.frame(life=blank_label, 
                                  MCOA1=0, MCOA2=0, 
                                  Correlation="Positive (r > 0)"), real_df)
    }
  }
  
  # 锁定 Factor 顺序，确保 ggplot 按照现在的顺序画图
  real_df$life <- factor(real_df$life, levels = real_df$life)
  
  color_map <- c("Positive (r > 0)" = "#B9E44E", "Negative (r < 0)" = "#43B5D4")
  
  # 设置统一的横轴范围
  common_x_limit <- c(-1, 1) 
  
  p <- ggplot(real_df) +
    geom_bar(aes(x = life, y = !!sym(axis_name), fill = Correlation), 
             stat = "identity", color = "black", width = 0.7) +
    scale_fill_manual(values = color_map) +
    # 统一刻度范围
    scale_y_continuous(limits = common_x_limit, breaks = seq(-1, 1, by = 0.5)) + 
    labs(x = "", y = "Correlation (r)", title = plot_title) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.text.x = element_text(size = 10, color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      legend.position = "bottom"
    ) +
    coord_flip() 
  
  return(p)
}

# ==================== 4. 绘图与拼图 (0.05) ====================
data_05 <- read.csv("lifehistory/life0.05.csv")
# 以 MCOA1 的显著指标数量作为对齐基准
count_m1 <- sum(data_05$MCOA1 != 0)

p1_05 <- draw_mcoa_bar_v6("lifehistory/life0.05.csv", "MCOA1", "MCOA1 Drivers", target_count = count_m1)
p2_05 <- draw_mcoa_bar_v6("lifehistory/life0.05.csv", "MCOA2", "MCOA2 Drivers", target_count = count_m1)

if(!is.null(p1_05) & !is.null(p2_05)){
  p_combined <- ggarrange(p1_05, p2_05, ncol = 2, common.legend = TRUE, legend = "bottom")
  ggsave("lifehistory/combined_plot0.05_restored.pdf", plot = p_combined, width = 16, height = 10)
}

# ==================== 5. 处理 0.01 结果 ====================
data_01 <- read.csv("lifehistory/life0.01.csv")
count_m1_01 <- sum(data_01$MCOA1 != 0)

p1_01 <- draw_mcoa_bar_v6("lifehistory/life0.01.csv", "MCOA1", "MCOA1 Drivers (P < 0.01)", target_count = count_m1_01)

if(!is.null(p1_01)) {
  ggsave("lifehistory/MCOA1_only_0.01_restored.pdf", plot = p1_01, width = 9, height = 10)
}

message("运行成功！纵坐标已恢复自动排序，且 MCOA2 柱子厚度已对齐。")

