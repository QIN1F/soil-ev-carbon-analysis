rm(list=ls())
library(omicade4)
library(ggplot2)
library(readxl)
library(dplyr)

# ==================== 1. 设置与路径 ====================
# 指向第一个代码生成的文件夹
path_in <- "lifehistory-M/"
group_file <- "../group_lifehistory.txt" 

# 视觉设置 (针对 M0 和 M1 两个分组)
group_levels <- c("M0", "M1")
mycolors <- c("#3C8DBC", "#D33724")
myshapes <- c(21, 21) # 只需要 2 个值对应两个 level

# ==================== 2. 载入第一个代码生成的数据 ====================
# 读取 1.txt, 2.txt, 3.txt (注意 row.names = 1)
data1 <- read.delim(paste0(path_in, "1.Life_history_trait_M.txt"), sep="\t", header=T, check.names=F, row.names = 1)
data2 <- read.delim(paste0(path_in, "2.eggNOG_M.txt"), sep="\t", header=T, check.names=F, row.names = 1)
data3 <- read.delim(paste0(path_in, "3.CAZy_M.txt"), sep="\t", header=T, check.names=F, row.names = 1)

# 读取 4.xlsx
data4_path <- paste0(path_in, "4.KEGG_M.xlsx")
temp_data4 <- as.data.frame(read_excel(data4_path))
rownames(temp_data4) <- temp_data4[,1]
data4 <- temp_data4[,-1]

# 读取位于上一级的分组文件
Label <- read.delim(group_file, sep="\t", header=T, check.names=F, row.names = 1)

# ==================== 3. 数据清洗与对齐 ====================
# 过滤全 0 行
data1 <- data1[rowSums(data1) > 0, ]
data2 <- data2[rowSums(data2) > 0, ]
data3 <- data3[rowSums(data3) > 0, ]
data4 <- data4[rowSums(data4) > 0, ]

# 关键：确保四个表的列名（样本）顺序完全一致
common_samples <- colnames(data1)
data2 <- data2[, common_samples]
data3 <- data3[, common_samples]
data4 <- data4[, common_samples]

# ==================== 4. MCoA 计算 ====================
NCI60_4arrays <- list(data1, data2, data3, data4)
names(NCI60_4arrays) <- c('Life_history_trait','eggNOG','CAZy','KEGG')

# nsc = TRUE 进行标准化处理
mcoin <- mcia(NCI60_4arrays, cia.nf = 2, cia.scan = FALSE, nsc = TRUE)

# ==================== 5. 结果提取与保存 ====================
# 计算解释率
total_eig <- sum(mcoin$mcoa$pseudoeig, na.rm = TRUE)
var1 <- mcoin$mcoa$pseudoeig[1] / total_eig * 100
var2 <- mcoin$mcoa$pseudoeig[2] / total_eig * 100

# 保存 RV 系数和特征值 (存入 lifehistory-M)
write.table(as.data.frame(mcoin$mcoa$pseudoeig), paste0(path_in, "pseudoeig_M.txt"), sep="\t", col.names = NA)
write.table(mcoin$mcoa$RV, paste0(path_in, "RV_M.txt"), sep="\t", col.names = NA)

# 整理坐标数据
coordinate2 <- mcoin$mcoa$SynVar
colnames(coordinate2) <- c("Axis1","Axis2")
coordinate2$name <- row.names(coordinate2)
Label$name <- row.names(Label)

# 合并分组信息
coordinate2 <- merge(coordinate2, Label, by="name")
coordinate2$group <- factor(coordinate2$group, levels = group_levels)

# ==================== 6. 绘图与导出 ====================
p2 <- ggplot(data=coordinate2, aes(x=Axis1, y=Axis2, fill=group, shape=group)) + 
  geom_point(size = 5, stroke = 0.6, colour = "gray20") +
  theme_bw() +
  geom_hline(yintercept=0, linetype=2, color='gray70') + 
  geom_vline(xintercept=0, linetype=2, color='gray70') +
  xlab(paste0("MCoA1 (", sprintf("%.2f", var1), "%)")) +
  ylab(paste0("MCoA2 (", sprintf("%.2f", var2), "%)")) +
  scale_fill_manual(values = mycolors) +
  scale_shape_manual(values = myshapes) +
  theme(aspect.ratio = 1, 
        panel.grid = element_blank(),
        legend.title = element_blank())

print(p2)

# 保存 PDF 和供下个代码使用的 CSV
ggsave(paste0(path_in, "MCOA_BacOnly_Analysis.pdf"), width = 6, height = 5)

# 整合 MCoA 坐标与原始生活史性状丰度
data1_t <- as.data.frame(t(data1))
data1_t$name <- rownames(data1_t)
all_combined_M <- merge(coordinate2[,c("name","Axis1","Axis2", "group")], data1_t, by="name")

# 导出 CSV
write.csv(all_combined_M, paste0(path_in, "mcoa_M_samples_final.csv"), row.names = F)

message("✅ 细菌组 MCoA 专题分析完成！文件已存入 lifehistory-M/")

