rm(list=ls()) # 清除环境

library(omicade4)
library(ggplot2)
library(readxl)
library(openxlsx) 

# ==================== 0. 你的专属本地路径设置 ====================
# 数据文件夹路径 (已更改为 lifehistory-M)
data_dir <- "E:/我的日常文件/111微塑料循环汇报/4.22生活史策略分析/我的/gene_base/lifehistory-M/"
# 分组文件路径
group_file <- "E:/我的日常文件/111微塑料循环汇报/4.22生活史策略分析/我的/group_lifehistory.txt"

# ==================== 1. 样品与视觉设置 (细菌组专版) ====================
# 【注意】请确保 group_lifehistory.txt 中细菌组的名字是这两者，如果不同请自行修改！
group_levels <- c("M0", "M1") 
# 换了一组区别于 EV 组的蓝红色，方便你把两张图放在一起对比时不混淆
mycolors <- c("#3498DB", "#E74C3C") 
myshapes <- c(rep(21, 12)) 

# ==================== 2. 数据载入 ====================
data1 <- read.delim(paste0(data_dir, "1.Life_history_trait_M.txt"), sep="\t", header=T, check.names=F, row.names = 1)
data2 <- read.delim(paste0(data_dir, "2.eggNOG_M.txt"), sep="\t", header=T, check.names=F, row.names = 1)
data3 <- read.delim(paste0(data_dir, "3.CAZy_M.txt"), sep="\t", header=T, check.names=F, row.names = 1)

data4_path <- paste0(data_dir, "4.KEGG_M.xlsx")
if(file.exists(data4_path)){
  temp_data4 <- as.data.frame(read_excel(data4_path))
  rownames(temp_data4) <- temp_data4[,1]
  data4 <- temp_data4[,-1]
}

Label <- read.delim(group_file, sep="\t", header=T, check.names=F, row.names = 1)

# ==================== 2.5 数据清洗 ====================
data1 <- data1[rowSums(data1) > 0, ]
data2 <- data2[rowSums(data2) > 0, ]
data3 <- data3[rowSums(data3) > 0, ]
data4 <- data4[rowSums(data4) > 0, ]

# ==================== 3. MCoA 计算 ====================
NCI60_4arrays <- list(data1, data2, data3, data4)
names(NCI60_4arrays) <- c('Life_history_trait','eggNOG','CAZy','KEGG')

# 提取前 4 个轴
mcoin <- mcia(NCI60_4arrays, cia.nf = 4, cia.scan = FALSE, nsc = TRUE)

# ==================== 4. 提取与保存统计结果 ====================
total_eig <- sum(mcoin$mcoa$pseudoeig, na.rm = TRUE)
var1 <- mcoin$mcoa$pseudoeig[1] / total_eig * 100
var2 <- mcoin$mcoa$pseudoeig[2] / total_eig * 100

write.table(as.data.frame(mcoin$mcoa$pseudoeig), paste0(data_dir, "pseudoeig_M.txt"), sep="\t", col.names = NA)
write.table(mcoin$mcoa$RV, paste0(data_dir, "RV_M.txt"), sep="\t", col.names = NA)

# 提取用于画散点图的前两个轴
coordinate2 <- mcoin$mcoa$SynVar[, 1:2]
colnames(coordinate2) <- c("Axis1","Axis2")
coordinate2$name <- row.names(coordinate2)
Label$name <- row.names(Label)

coordinate2 <- merge(coordinate2, Label, by="name")
if(nrow(coordinate2) == 0) stop("错误：样本名匹配失败，请检查分组文件和 M 样本名是否一致！")
coordinate2$group <- factor(coordinate2$group, levels = group_levels)

# ==================== 5. 绘图 (细菌组特写 MCOA1 vs MCOA2) ====================
p_m <- ggplot(data=coordinate2, aes(x=Axis1, y=Axis2, fill=group, shape=group)) + 
  geom_point(size = 5, stroke = 0.6, colour = "gray20") +
  theme_bw() +
  geom_hline(yintercept=0, linetype=2, color='gray70') + 
  geom_vline(xintercept=0, linetype=2, color='gray70') +
  xlab(paste0("MCoA1 (", sprintf("%.2f", var1), "%)")) +
  ylab(paste0("MCoA2 (", sprintf("%.2f", var2), "%)")) +
  scale_fill_manual(values = mycolors) +
  scale_shape_manual(values = myshapes) +
  theme(aspect.ratio = 1, panel.grid = element_blank())

print(p_m)
ggsave(paste0(data_dir, "MCOA_MOnly_Analysis.pdf"), width = 6, height = 5)

# ==================== 6. 导出 4 轴数据用于图 2d 热图 ====================
coordinate_4axes <- mcoin$mcoa$SynVar[, 1:4]
colnames(coordinate_4axes) <- c("MCOA1", "MCOA2", "MCOA3", "MCOA4")
coordinate_4axes$name <- row.names(coordinate_4axes)

data11 <- as.data.frame(t(data1))
data11$name <- rownames(data11)

all_combined_M_4axes <- merge(coordinate_4axes, data11, by="name")

# 导出带 4 个轴的专属文件
write.xlsx(all_combined_M_4axes, paste0(data_dir, "mcoa_M_samples_4axes.xlsx"), rowNames = FALSE)

# 保留旧版的 2 轴数据
write.xlsx(all_combined_M_4axes[, !names(all_combined_M_4axes) %in% c("MCOA3", "MCOA4")], paste0(data_dir, "mcoa_M_samples.xlsx"), rowNames = FALSE)

message("运行成功！细菌组4轴数据已保存在: ", paste0(data_dir, "mcoa_M_samples_4axes.xlsx"))

