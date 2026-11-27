# install.packages(c("Hmisc", "ggpubr", "tidyverse", "readr")) 
library(ggplot2)
library(dplyr)
library(ggpubr)
library(readr)
library(tidyverse)
library(grid)
library(Hmisc) 

# ==================== 1. 读取数据 (指向 lifehistory-M 文件夹) ====================
path_in <- "lifehistory-M/"
mcoa_file <- paste0(path_in, "mcoa_M_samples_final.csv")
if(!file.exists(mcoa_file)) stop("错误：找不到 mcoa_M_samples_final.csv，请先运行 MCoA 代码")

# 读取数据，保留第一列 name 作为 index
mcoa_M_data_raw <- read.csv(mcoa_file, check.names = F)
rownames(mcoa_M_data_raw) <- mcoa_M_data_raw$name

# 核心修正：移除非数值列 (name, group) 以便进行相关性计算
# 假设前四列是 name, Axis1, Axis2, group
mcoa_M_numeric <- mcoa_M_data_raw %>% select(-name, -group)

# 强制将前两列坐标重命名为 MCOA1 和 MCOA2
colnames(mcoa_M_numeric)[1:2] <- c("MCOA1", "MCOA2")

# ==================== 2. 相关性计算与筛选 ====================
mcoa_sp <- rcorr(as.matrix(mcoa_M_numeric), type="spearman")

# 提取性状与前两个坐标轴的相关系数和 P 值
mcoa_r <- mcoa_sp$r[-c(1,2), c(1,2)]
mcoa_p <- mcoa_sp$P[-c(1,2), c(1,2)]

# 保存 P < 0.05 结果
mcoa_r_05 <- mcoa_r
mcoa_r_05[mcoa_p > 0.05] <- 0
write.csv(mcoa_r_05, paste0(path_in, "life0.05_M.csv"))

# 保存 P < 0.01 结果
mcoa_r_01 <- mcoa_r
mcoa_r_01[mcoa_p > 0.01] <- 0
write.csv(mcoa_r_01, paste0(path_in, "life0.01_M.csv"))

# ==================== 3. 绘图函数 (比例优化版) ====================
draw_mcoa_bar_M_v8 <- function(data_path, axis_name, plot_title) {
  if(!file.exists(data_path)) return(NULL)
  data_plot <- read.csv(data_path)
  
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
    geom_bar(aes(x = life, y = !!sym(axis_name), fill = Correlation), 
             stat = "identity", color = "black", width = 0.5) + 
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

# ==================== 4. 绘 plot 与拼图 (P < 0.05) ====================
file_05 <- paste0(path_in, "life0.05_M.csv")
if(file.exists(file_05)){
  p1_05 <- draw_mcoa_bar_M_v8(file_05, "MCOA1", "MCOA1 Drivers (P < 0.05)")
  p2_05 <- draw_mcoa_bar_M_v8(file_05, "MCOA2", "MCOA2 Drivers (P < 0.05)")
  
  plot_list <- list(p1_05, p2_05)
  plot_list <- plot_list[!sapply(plot_list, is.null)]
  
  if(length(plot_list) > 0){
    p_combined <- ggarrange(plotlist = plot_list, ncol = length(plot_list), common.legend = TRUE, legend = "bottom")
    
    data_05_check <- read.csv(file_05)
    num_bars <- max(sum(data_05_check[,2] != 0, na.rm=T), sum(data_05_check[,3] != 0, na.rm=T))
    dynamic_height <- max(6, num_bars * 0.5 + 2)
    
    ggsave(paste0(path_in, "combined_plot0.05_M_final.pdf"), plot = p_combined, 
           width = 6 * length(plot_list), height = dynamic_height)
  }
}

# ==================== 5. 绘 plot 与拼图 (P < 0.01) ====================
file_01 <- paste0(path_in, "life0.01_M.csv")
if(file.exists(file_01)){
  p1_01 <- draw_mcoa_bar_M_v8(file_01, "MCOA1", "MCOA1 Drivers (P < 0.01)")
  
  if(!is.null(p1_01)) {
    data_01_check <- read.csv(file_01)
    num_bars_01 <- sum(data_01_check[,2] != 0, na.rm=T)
    dynamic_height_01 <- max(5, num_bars_01 * 0.8 + 1)
    
    ggsave(paste0(path_in, "MCOA1_only_0.01_M_final.pdf"), plot = p1_01, 
           width = 8, height = dynamic_height_01)
  }
}

message("✅ 比例优化版处理完成！驱动因子分析结果已存入 lifehistory-M/")

