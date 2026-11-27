library(readxl)
library(dplyr)
library(openxlsx)

# ==================== 1. 路径设置 (基于您在 gene-base 文件夹内运行) ====================
# 结果导出文件夹：自动在 gene-base 下创建
path_out <- "lifehistory-M/"
if(!dir.exists(path_out)) dir.create(path_out)

# 外部文件夹文件 (使用 ../ 退回上一级目录)
index      <- "../annotation/database_lifehistory_G2.txt" 
group_file <- "../group_lifehistory.txt" 

# 内部数据文件 (直接在当前 gene-base 目录下读取)
KO_file         <- "KO.profile.entry.xlsx"
eggNOG_file      <- "eggNOG.Category.tpm.xlsx"
eggNOG_ID_file   <- "eggNOG.profile.entry.xlsx"
CAZy_file        <- "CAZy.tpm.xlsx"
LEVELC_file      <- "LEVELC.tpm.xlsx"

# 目标样本清单
bac_samples <- c("CK1", "CK2", "CK3", "CK4", "CK5", "CK6", "P1", "P2", "P3", "P4", "P5", "P6")

# 定义数据过滤函数：保留第一列并提取 12 个目标样本
filter_bac_data <- function(df) {
  existing_bac_cols <- bac_samples[bac_samples %in% colnames(df)]
  return(df[, c(1, which(colnames(df) %in% existing_bac_cols))])
}

# ==================== 2. 数据读取与预检查 ====================
# 检查路径是否正确
if(!file.exists(index)) stop("❌ 找不到数据库文件！请确认您在 gene-base 文件夹下运行 R。")
if(!file.exists(group_file)) stop("❌ 找不到分组文件！")

# 读取基础数据
Lifestyle <- read.delim(index) %>% unique()
Group_info <- read.delim(group_file) # 读取分组文件供后续分析参考 

# 1. 处理 KO [cite: 192, 313]
Data_3 <- read_excel(KO_file) %>% filter_bac_data()

# 2. 处理 eggNOG_Category (整理分类标签并合并求和) [cite: 194, 327]
Data_4 <- read_excel(eggNOG_file)
colnames(Data_4)[1] <- "eggNOG_Category"
Data_4$eggNOG_Category <- gsub("^[^:]*:", "", Data_4$eggNOG_Category) 
Data_4_M <- Data_4 %>%
  filter_bac_data() %>%
  group_by(eggNOG_Category) %>%
  summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))

# 3. 处理 eggNOG_ID [cite: 313, 316]
Data_5 <- read_excel(eggNOG_ID_file)
Data_5 <- Data_5[, -ncol(Data_5)] %>% filter_bac_data()
colnames(Data_5)[1] <- "eggNOG_ID"

# 4. 处理 CAZy (整合描述信息) [cite: 194, 313]
Data_6 <- read_excel(CAZy_file) %>% filter_bac_data()
Lifestyle6_M <- Lifestyle %>%
  select(Function, CAZy) %>%
  merge(Data_6, by = "CAZy") %>%
  rename(Description = CAZy)

# ==================== 3. 整合综合生活史表 [cite: 329, 347] ====================
Lifestyle3 <- Lifestyle %>% select(Function, KO) %>% merge(Data_3, by = "KO")
Lifestyle4 <- Lifestyle %>% select(Function, eggNOG_Category) %>% 
  mutate(eggNOG_Category = gsub("^[^:]*:", "", eggNOG_Category)) %>%
  merge(Data_4_M, by = "eggNOG_Category")

# 合并所有功能来源
results <- bind_rows(Lifestyle3, Lifestyle4, Lifestyle6_M %>% mutate(Origin="CAZy"))

# 生成最终的生活史策略矩阵 (Function 为行名) [cite: 304, 305]
df2_M <- results %>%
  group_by(Function) %>%
  summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))

# ==================== 4. 结果导出 (至 lifehistory-M/) ====================
# 1. 生活史主表 (用于后续 MCOA/SEM 分析) [cite: 327, 329]
write.table(df2_M, paste0(path_out, "1.Life_history_trait_M.txt"), sep="\t", row.names = FALSE, quote = FALSE)

# 2. 详细功能表
write.table(Data_4_M, paste0(path_out, "2.eggNOG_M.txt"), sep="\t", row.names = FALSE, quote = FALSE)
write.table(CAZy_Export_M <- Lifestyle6_M %>% select(-Function), 
            paste0(path_out, "3.CAZy_M.txt"), sep="\t", row.names = FALSE, quote = FALSE)

# 3. KEGG 数据 [cite: 193]
Data_C <- read_excel(LEVELC_file) %>% filter_bac_data()
write.xlsx(Data_C, paste0(path_out, "4.KEGG_M.xlsx"))

print(paste0("✅ LH1_M 处理完成！结果已存入：", getwd(), "/", path_out))

