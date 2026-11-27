
rm(list=ls()) # 清除环境


library(omicade4)
library(ggplot2)
library(readxl)

# ==================== 1. 样品与视觉设置 ====================
group_levels <- c("M0", "M1", "EV0", "EV1")
mycolors <- c("#3C8DBC", "#D33724", "#00A65A", "#F39C12") 
myshapes <- c(rep(21, 6), rep(22, 6), rep(24, 6), rep(25, 6))

# ==================== 2. 数据载入 ====================
data1 <- read.delim("lifehistory/1.Life_history_trait.txt", sep="\t", header=T, check.names=F, row.names = 1)
data2 <- read.delim("lifehistory/2.eggNOG.txt", sep="\t", header=T, check.names=F, row.names = 1)
data3 <- read.delim("lifehistory/3.CAZy.txt", sep="\t", header=T, check.names=F, row.names = 1)

data4_path <- "lifehistory/4.KEGG.xlsx"
if(file.exists(data4_path)){
  temp_data4 <- as.data.frame(read_excel(data4_path))
  rownames(temp_data4) <- temp_data4[,1]
  data4 <- temp_data4[,-1]
}

Label <- read.delim("../group_lifehistory.txt", sep="\t", header=T, check.names=F, row.names = 1)

# ==================== 2.5 核心修正：数据清洗 (解决报错) ====================
# 过滤掉那些在所有样本中丰度之和为 0 的行
data1 <- data1[rowSums(data1) > 0, ]
data2 <- data2[rowSums(data2) > 0, ]
data3 <- data3[rowSums(data3) > 0, ]
data4 <- data4[rowSums(data4) > 0, ]

# 检查是否有数据集被删空了
if(nrow(data4) == 0) stop("错误：KEGG 数据集在过滤全0行后没有剩余数据，请检查原始文件！")

# ==================== 3. MCoA 计算 ====================
NCI60_4arrays <- list(data1, data2, data3, data4)
names(NCI60_4arrays) <- c('Life_history_trait','eggNOG','CAZy','KEGG')

# 现在运行 mcia 应该就不会报错了
mcoin <- mcia(NCI60_4arrays, cia.nf = 4, cia.scan = FALSE, nsc = TRUE)

# ... (后面绘图和导出的代码保持不变)
# ==================== 4. 提取与保存结果 ====================
total_eig <- sum(mcoin$mcoa$pseudoeig, na.rm = TRUE)
var1 <- mcoin$mcoa$pseudoeig[1] / total_eig * 100
var2 <- mcoin$mcoa$pseudoeig[2] / total_eig * 100

write.table(as.data.frame(mcoin$mcoa$pseudoeig), "lifehistory/pseudoeig.txt", sep="\t", col.names = NA)
write.table(mcoin$mcoa$RV, "lifehistory/RV.txt", sep="\t", col.names = NA)

coordinate2 <- mcoin$mcoa$SynVar
colnames(coordinate2) <- c("Axis1","Axis2","Axis3","Axis4")

# ==================== 5. 绘图 (针对 sample 列名的修正) ====================

# 1. 处理坐标表的样本名
coordinate2$name <- row.names(coordinate2)

# 2. 处理分组表的样本名
# 因为读取时设置了 row.names = 1，sample 变成了行名，现在把它取出来存入 name 列
Label$name <- row.names(Label)

# 3. 合并两个表
# 这一步会自动通过 name 列对齐 24 个样本
coordinate2 <- merge(coordinate2, Label, by="name")

# 检查合并是否成功，如果不成功通常是样本名大小写或特殊符号不一致
if(nrow(coordinate2) == 0) stop("错误：样本名匹配失败，请检查数据表和分组表的样本名是否一致！")

# 4. 设置分组因子顺序
coordinate2$group <- factor(coordinate2$group, levels = group_levels)

# 5. 绘图
p2 <- ggplot(data=coordinate2, aes(x=Axis1, y=Axis2, fill=group, shape=group)) + 
  geom_point(size = 5, stroke = 0.6, colour = "gray20") +
  theme_bw() +
  geom_hline(yintercept=0, linetype=2, color='gray70') + 
  geom_vline(xintercept=0, linetype=2, color='gray70') +
  xlab(paste0("MCoA1 (", sprintf("%.2f", var1), "%)")) +
  ylab(paste0("MCoA2 (", sprintf("%.2f", var2), "%)")) +
  scale_fill_manual(values = mycolors) +
  scale_shape_manual(values = myshapes) +
  theme(aspect.ratio = 1, panel.grid = element_blank())

print(p2)
ggsave("lifehistory/MCOA_24samples_Analysis.pdf", width = 7, height = 6)


# ==================== 6. 导出 downstream 文件 (用于 LH3) ====================

# 提取关键列
coordinate22 <- coordinate2[,c("name","Axis1","Axis2")]
colnames(coordinate22) <- c("name","MCOA1","MCOA2")

# 转换 data1 (Life history trait) 的格式
data11 <- as.data.frame(t(data1))
data11$name <- rownames(data11)

# 合并坐标与功能指标
all_combined <- merge(coordinate22, data11, by="name")

# 最终导出到 lifehistory 文件夹
write.csv(all_combined, "lifehistory/mcoa_240_1.csv", row.names = F)

message("LH2 全部运行成功！请查看 lifehistory 文件夹下的 PDF 和 CSV 文件。")

