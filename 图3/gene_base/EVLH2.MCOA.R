rm(list=ls()) # 清除环境

library(omicade4)
library(ggplot2)
library(readxl)

# ==================== 1. 样品与视觉设置 (囊泡组专版) ====================
# 注意：确保你的 group_lifehistory.txt 分组文件中对应的组名是 EV0 和 EV1
group_levels <- c("EV0", "EV1")
mycolors <- c("#00A65A", "#F39C12") # EV0绿色，EV1橙色（保持与总图一致）
# 12个样本，统一使用 21 号形状（带填充的圆圈），看起来更整齐
myshapes <- c(rep(21, 12)) 

# ==================== 2. 数据载入 (对接 LH1_EV 生成的文件) ====================
# 读取带 _EV 后缀的文件
data1 <- read.delim("lifehistory/1.Life_history_trait_EV.txt", sep="\t", header=T, check.names=F, row.names = 1)
data2 <- read.delim("lifehistory/2.eggNOG_EV.txt", sep="\t", header=T, check.names=F, row.names = 1)
data3 <- read.delim("lifehistory/3.CAZy_EV.txt", sep="\t", header=T, check.names=F, row.names = 1)

data4_path <- "lifehistory/4.KEGG_EV.xlsx"
if(file.exists(data4_path)){
  temp_data4 <- as.data.frame(read_excel(data4_path))
  rownames(temp_data4) <- temp_data4[,1]
  data4 <- temp_data4[,-1]
}

# 载入分组信息
Label <- read.delim("../group_lifehistory.txt", sep="\t", header=T, check.names=F, row.names = 1)

# ==================== 2.5 数据清洗 ====================
data1 <- data1[rowSums(data1) > 0, ]
data2 <- data2[rowSums(data2) > 0, ]
data3 <- data3[rowSums(data3) > 0, ]
data4 <- data4[rowSums(data4) > 0, ]

# ==================== 3. MCoA 计算 ====================
NCI60_4arrays <- list(data1, data2, data3, data4)
names(NCI60_4arrays) <- c('Life_history_trait','eggNOG','CAZy','KEGG')

# 计算 MCoA
mcoin <- mcia(NCI60_4arrays, cia.nf = 2, cia.scan = FALSE, nsc = TRUE)

# ==================== 4. 提取与保存统计结果 ====================
total_eig <- sum(mcoin$mcoa$pseudoeig, na.rm = TRUE)
var1 <- mcoin$mcoa$pseudoeig[1] / total_eig * 100
var2 <- mcoin$mcoa$pseudoeig[2] / total_eig * 100

# 保存统计文件，带 EV 标识
write.table(as.data.frame(mcoin$mcoa$pseudoeig), "lifehistory/pseudoeig_EV.txt", sep="\t", col.names = NA)
write.table(mcoin$mcoa$RV, "lifehistory/RV_EV.txt", sep="\t", col.names = NA)

coordinate2 <- mcoin$mcoa$SynVar
colnames(coordinate2) <- c("Axis1","Axis2")
coordinate2$name <- row.names(coordinate2)
Label$name <- row.names(Label)

# 合并分组信息：这会自动保留 EV0 和 EV1 的 12 个样本
coordinate2 <- merge(coordinate2, Label, by="name")

if(nrow(coordinate2) == 0) stop("错误：样本名匹配失败，请检查分组文件和 EV 样本名是否一致！")

coordinate2$group <- factor(coordinate2$group, levels = group_levels)

# ==================== 5. 绘图 (囊泡组特写) ====================
p_ev <- ggplot(data=coordinate2, aes(x=Axis1, y=Axis2, fill=group, shape=group)) + 
  geom_point(size = 5, stroke = 0.6, colour = "gray20") +
  theme_bw() +
  geom_hline(yintercept=0, linetype=2, color='gray70') + 
  geom_vline(xintercept=0, linetype=2, color='gray70') +
  xlab(paste0("MCoA1 (", sprintf("%.2f", var1), "%)")) +
  ylab(paste0("MCoA2 (", sprintf("%.2f", var2), "%)")) +
  scale_fill_manual(values = mycolors) +
  scale_shape_manual(values = myshapes) +
  theme(aspect.ratio = 1, panel.grid = element_blank())

print(p_ev)
ggsave("lifehistory/MCOA_EVOnly_Analysis.pdf", width = 6, height = 5)

# ==================== 6. 导出用于 LH3_EV 的数据 ====================
coordinate22 <- coordinate2[,c("name","Axis1","Axis2")]
colnames(coordinate22) <- c("name","MCOA1","MCOA2")

data11 <- as.data.frame(t(data1))
data11$name <- rownames(data11)

# 合并坐标与功能指标
all_combined_EV <- merge(coordinate22, data11, by="name")

# 最终导出到 lifehistory 文件夹，文件名带 EV
write.csv(all_combined_EV, "lifehistory/mcoa_EV_samples.csv", row.names = F)

message("LH2_EV 运行成功！已生成 MCOA_EVOnly_Analysis.pdf 和 mcoa_EV_samples.csv")

