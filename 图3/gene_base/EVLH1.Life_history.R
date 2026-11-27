# install.packages("openxlsx") 
library(readxl)
library(dplyr)
library(openxlsx)

# ==================== 环境设置 ====================
path_out = "lifehistory/"
if(!dir.exists(path_out)) dir.create(path_out)

index <- "../annotation/database_lifehistory_G2.txt" 

# 美吉结果文件名 (xlsx)
KO_file <- "KO.profile.entry.xlsx"
eggNOG_file <- "eggNOG.Category.tpm.xlsx"
eggNOG_ID_file <- "eggNOG.profile.entry.xlsx"
CAZy_file <- "CAZy.tpm.xlsx"
LEVELC_file <- "LEVELC.tpm.xlsx"

# 定义囊泡组样本列名 (EV0: CK1_EV-6_EV, EV1: p1_EV-6_EV)
# 注意：请务必核对你的 xlsx 文件中这些样本的大小写是否完全一致
ev_samples <- c("CK1_EV", "CK2_EV", "CK3_EV", "CK4_EV", "CK5_EV", "CK6_EV", 
                "p1_EV", "p2_EV", "p3_EV", "p4_EV", "p5_EV", "p6_EV")

# 辅助函数：只保留第一列 ID 和 12 个囊泡样本列
filter_ev_data <- function(df) {
  existing_ev_cols <- ev_samples[ev_samples %in% colnames(df)]
  return(df[, c(1, which(colnames(df) %in% existing_ev_cols))])
}

# 读取生活史数据库
Lifestyle <- read.delim(index) %>% unique()

# ==================== 1. KO 处理 ====================
Data_3 <- read_excel(KO_file) %>% filter_ev_data()

Lifestyle3 <- Lifestyle %>%
  select(Function, KO) %>%
  merge(Data_3, by = "KO") %>%
  rename(Description = KO)

# ==================== 2. eggNOG_Category 处理 (合并重复项) ====================
Data_4 <- read_excel(eggNOG_file)
colnames(Data_4)[1] <- "eggNOG_Category"
Data_4$eggNOG_Category <- gsub("^[^:]*:", "", Data_4$eggNOG_Category) 

Data_4_EV <- Data_4 %>%
  filter_ev_data() %>%
  group_by(eggNOG_Category) %>%
  summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))

Lifestyle4 <- Lifestyle %>%
  select(Function, eggNOG_Category) %>%
  mutate(eggNOG_Category = gsub("^[^:]*:", "", eggNOG_Category)) %>%
  merge(Data_4_EV, by = "eggNOG_Category") %>%
  rename(Description = eggNOG_Category)

# ==================== 3. eggNOG_ID 处理 ====================
Data_5 <- read_excel(eggNOG_ID_file)
Data_5 <- Data_5[, -ncol(Data_5)] %>% filter_ev_data()
colnames(Data_5)[1] <- "eggNOG_ID"

# ==================== 4. CAZy 处理 (确保 ID 唯一) ====================
Data_6 <- read_excel(CAZy_file) %>% filter_ev_data()

Lifestyle6 <- Lifestyle %>%
  select(Function, CAZy) %>%
  merge(Data_6, by = "CAZy") %>%
  rename(Description = CAZy)

# ==================== 数据整合与求和 (1.Life_history_trait_EV.txt) ====================
results <- bind_rows(
  Lifestyle3 %>% mutate(Origin="KO"), 
  Lifestyle4 %>% mutate(Origin="eggNOG_Category"), 
  Lifestyle6 %>% mutate(Origin="CAZy")
)

df2_EV <- results %>%
  group_by(Function) %>%
  summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE))) %>%
  replace(is.na(.), 0)

# ==================== 导出结果 (带 _EV 后缀) ====================

# 1. 导出生活史主表
write.table(df2_EV, paste0(path_out, "1.Life_history_trait_EV.txt"), 
            sep="\t", row.names = FALSE, quote = FALSE)

# 2. 导出 eggNOG 分类表
write.table(Data_4_EV, paste0(path_out, "2.eggNOG_EV.txt"), 
            sep="\t", row.names = FALSE, quote = FALSE)

# 3. 导出 CAZy 表 (用 Description 确保唯一行名)
CAZy_Export_EV <- Lifestyle6 %>% select(-Function)
write.table(CAZy_Export_EV, paste0(path_out, "3.CAZy_EV.txt"), 
            sep="\t", row.names = FALSE, quote = FALSE)

# 4. 处理 LEVELC (KEGG) 文件
Data_C <- read_excel(LEVELC_file) %>% filter_ev_data()
write.xlsx(Data_C, paste0(path_out, "4.KEGG_EV.xlsx"))

print("LH1_EV (囊泡组专题) 分析完成！生成文件带有 _EV 后缀。")

