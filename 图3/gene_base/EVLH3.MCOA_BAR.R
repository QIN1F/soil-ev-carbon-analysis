# install.packages(c("Hmisc", "ggpubr", "tidyverse", "readr")) 
library(ggplot2)
library(dplyr)
library(ggpubr)
library(readr)
library(tidyverse)
library(grid)
library(Hmisc) 

# ==================== 1. 读取数据 (对接 LH2_EV) ====================
mcoa_file <- "lifehistory/mcoa_EV_samples.csv"
if(!file.exists(mcoa_file)) stop("错误：找不到 mcoa_EV_samples.csv")

mcoa_EV_data <- read.csv(mcoa_file, row.names=1, check.names = F)
colnames(mcoa_EV_data)[1:2] <- c("MCOA1", "MCOA2")

# ==================== 2. 相关性计算与筛选 ====================
mcoa_sp <- rcorr(as.matrix(mcoa_EV_data), type="spearman")
mcoa_r <- mcoa_sp$r[-c(1,2), c(1,2)]
mcoa_p <- mcoa_sp$P[-c(1,2), c(1,2)]

mcoa_r_05 <- mcoa_r
mcoa_r_05[mcoa_p > 0.05] <- 0
write.csv(mcoa_r_05, "lifehistory/life0.05_EV.csv")

mcoa_r_01 <- mcoa_r
mcoa_r_01[mcoa_p > 0.01] <- 0
write.csv(mcoa_r_01, "lifehistory/life0.01_EV.csv")

# ==================== 3. 绘图函数 (柱子调细版) ====================
draw_mcoa_bar_EV_slim <- function(data_file, axis_name, plot_title) {
  data_plot <- read.csv(data_file)
  
  if(ncol(data_plot) >= 3) {
    colnames(data_plot)[2:3] <- c("MCOA1", "MCOA2")
  }
  colnames(data_plot)[1] <- "life"
  
  real_df <- data_plot %>%
    filter(!!sym(axis_name) != 0) %>%
    arrange(!!sym(axis_name)) %>% 
    mutate(Correlation = ifelse(!!sym(axis_name) > 0, "Positive (r > 0)", "Negative (r < 0)"))
  
  if(nrow(real_df) == 0) return(NULL)
  
  real_df$life <- factor(real_df$life, levels = real_df$life)
  color_map <- c("Positive (r > 0)" = "#B9E44E", "Negative (r < 0)" = "#43B5D4")
  
  p <- ggplot(real_df) +
    # 【修改点】：width 从 0.5 调为 0.3，柱子会变得很清秀
    geom_bar(aes(x = life, y = !!sym(axis_name), fill = Correlation), 
             stat = "identity", color = "black", width = 0.3) + 
    scale_fill_manual(values = color_map) +
    scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.5)) + 
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

# ==================== 4. 绘图与拼图 (P < 0.05_EV) ====================
data_05_EV <- read.csv("lifehistory/life0.05_EV.csv")
colnames(data_05_EV)[2:3] <- c("MCOA1", "MCOA2")

p1_05_EV <- draw_mcoa_bar_EV_slim("lifehistory/life0.05_EV.csv", "MCOA1", "MCOA1 Drivers (P < 0.05)")
p2_05_EV <- draw_mcoa_bar_EV_slim("lifehistory/life0.05_EV.csv", "MCOA2", "MCOA2 Drivers (P < 0.05)")

if(!is.null(p1_05_EV) | !is.null(p2_05_EV)){
  plot_list <- list(p1_05_EV, p2_05_EV)
  plot_list <- plot_list[!sapply(plot_list, is.null)]
  p_combined <- ggarrange(plotlist = plot_list, ncol = length(plot_list), common.legend = TRUE, legend = "bottom")
  
  num_bars <- max(sum(data_05_EV$MCOA1 != 0), sum(data_05_EV$MCOA2 != 0))
  # 稍微拉高基础高度到 5，让细柱子有更多呼吸空间
  dynamic_height <- max(5, num_bars * 0.6 + 2)
  
  ggsave("lifehistory/combined_plot0.05_EV_slim.pdf", plot = p_combined, 
         width = 6 * length(plot_list), height = dynamic_height)
}

# ==================== 5. 极显著 (P < 0.01_EV) ====================
if(file.exists("lifehistory/life0.01_EV.csv")){
  data_01_EV <- read.csv("lifehistory/life0.01_EV.csv")
  colnames(data_01_EV)[2:3] <- c("MCOA1", "MCOA2")
  
  p1_01_EV <- draw_mcoa_bar_EV_slim("lifehistory/life0.01_EV.csv", "MCOA1", "MCOA1 Drivers (P < 0.01)")
  if(!is.null(p1_01_EV)) {
    num_bars_01 <- sum(data_01_EV$MCOA1 != 0)
    # 针对指标极少的情况，强制给一个舒展的高度
    dynamic_height_01 <- max(4, num_bars_01 * 1.0 + 1)
    ggsave("lifehistory/MCOA1_only_0.01_EV_slim.pdf", plot = p1_01_EV, width = 8, height = dynamic_height_01)
  }
}

message("调细版柱状图已完成！请查看带有 _slim 后缀的文件。")

