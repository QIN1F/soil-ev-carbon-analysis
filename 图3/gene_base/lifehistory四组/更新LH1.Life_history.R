install.packages("openxlsx")
# Load necessary libraries
library(readxl)
library(dplyr)
library(openxlsx)

# 环境设置
# 确保你的 .xlsx 文件都在当前 gene_base 文件夹下
path <- "../annotation/lifehistory/"
if(!dir.exists("lifehistory")) dir.create("lifehistory")
path_out="lifehistory/"
index <- "../annotation/database_lifehistory_G2.txt"  # 数据库路径（保持 .txt）

# 美吉结果文件名定义 (已改为 .xlsx)
KO_file <- "KO.profile.entry.xlsx"
eggNOG_file <- "eggNOG.Category.tpm.xlsx"
eggNOG_ID_file <- "eggNOG.profile.entry.xlsx"
CAZy_file <- "CAZy.tpm.xlsx"
LEVELC_file <- "LEVELC.tpm.xlsx"

out_dirs <- path_out

# 读取生活史数据库 (这是个文本文件，保留 read.delim)
Lifestyle <- read.delim(index)
Lifestyle <- unique(Lifestyle)

################# 1. KO #####################
# 使用 read_excel 读取 xlsx 文件
Data_3 <- read_excel(KO_file)

Lifestyle3 <- Lifestyle %>%
  select(Function, KO) %>%
  merge(Data_3, by = "KO") %>%
  rename(Description = KO) %>%
  mutate(Origin = "KO")

################# 2. eggNOG_Category #####################
Data_4 <- read_excel(eggNOG_file)
colnames(Data_4)[1] <- "eggNOG_Category"
# 去除可能存在的前缀
Data_4$eggNOG_Category <- gsub("^[^:]*:", "", Data_4$eggNOG_Category) 

Lifestyle4 <- Lifestyle %>%
  select(Function, eggNOG_Category) %>%
  merge(Data_4, by = "eggNOG_Category") %>%
  rename(Description = eggNOG_Category) %>%
  mutate(Origin = "eggNOG_Category")

################# 3. eggNOG_ID #####################
Data_5 <- read_excel(eggNOG_ID_file)
colnames(Data_5)[1] <- "eggNOG_ID"
# 去除最后一列描述列
Data_5 <- Data_5[,-ncol(Data_5)]

Lifestyle5 <- Lifestyle %>%
  select(Function, eggNOG_ID) %>%
  merge(Data_5, by = "eggNOG_ID") %>%
  rename(Description = eggNOG_ID) %>%
  mutate(Origin = "eggNOG_ID")

################# 4. CAZy #####################
Data_6 <- read_excel(CAZy_file)

Lifestyle6 <- Lifestyle %>%
  select(Function, CAZy) %>%
  merge(Data_6, by = "CAZy") %>%
  rename(Description = CAZy) %>%
  mutate(Origin = "CAZy")

################# Combine Data #####################
results <- bind_rows(Lifestyle3, Lifestyle4, Lifestyle5, Lifestyle6)

################# Summing by Function #####################
# 自动计算从第3列到倒数第2列的丰度总和
df2 <- results %>%
  select(-Origin) %>%
  group_by(Function) %>%
  # 使用 \(x) 定义匿名函数，这是新版 R 的推荐写法
  summarise(across(3:(ncol(results)-2), \(x) sum(x, na.rm = TRUE))) %>%
  replace(is.na(.), 0)

Lifestyle66 <- Lifestyle6[,-c(2,ncol(Lifestyle6))]

################# Exporting Data #####################
# 输出结果依然为 .txt (制表符分隔)，方便后续代码 LH2/LH3 直接读取
write.table(Data_4, paste0(out_dirs, "2.eggNOG.txt"), sep="\t", row.names = FALSE, quote = FALSE)
write.table(Lifestyle66, paste0(out_dirs, "3.CAZy.txt"), sep="\t", row.names = FALSE, quote = FALSE)
write.table(df2, paste0(out_dirs, "1.Life_history_trait.txt"), sep="\t", row.names = FALSE, quote = FALSE)

# 复制 LEVELC 文件
file.copy(
  from = LEVELC_file,
  to   = paste0(out_dirs, "4.KEGG.xlsx") # 这里输出也对应改为 xlsx
)

print("LH1 分析完成！请检查 lifehistory 文件夹下的输出文件。")

