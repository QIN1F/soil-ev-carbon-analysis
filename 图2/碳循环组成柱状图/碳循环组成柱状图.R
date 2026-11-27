# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)


# 🎨 自定义颜色设置区 (参考提取配色)
# ==========================================
custom_colors <- c(
  "Step 6: Central Carbon & Others"    = "#F2D59D",  # 灰色 (背景属性强，放最上层)
  "Step 5: Carbon Fixation"            = "#E1E5A5",  # 淡土黄色
  "Step 4: Fermentation & Acetate"     = "#c4e2b3",  # 橘黄色
  "Step 3: C1 & Methane Metabolism"    = "#7AB69E",  # 砖红色 (突出重点)
  "Step 2: Aromatics Degradation"      = "#588FA8",  # 灰绿色
  "Step 1: Organic Carbon Degradation" = "#9990B6"   # 灰紫色 (深色压底，放最下层)
)
# ==========================================
# ==========================================

# 1. Read data files
carbon_data <- read.table("Carbon.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
group_data <- read.table("group.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 2. Categorize genes into 6 classic scientific steps
categorize_gene_by_step <- function(gene) {
  gene <- tolower(gene)
  if (grepl("amylase|arabinosidase|galactosidase|glucosidase|glucuronidase|mannosidase|xylosidase|cellobiosidase|cellulase|chitiniase|pullulanase|fae|hexosaminidase", gene)) {
    return("Step 1: Organic Carbon Degradation")
  } else if (grepl("cata|ubix|bsdc|bcr", gene)) {
    return("Step 2: Aromatics Degradation")
  } else if (grepl("pmo|mxa|mau|fdh|fdo|fgh|frm", gene)) {
    return("Step 3: C1 & Methane Metabolism")
  } else if (grepl("acetaldehyde|acetate|ack|acs|pta|adh|ldh|pfld|acda|por", gene)) {
    return("Step 4: Fermentation & Acetate")
  } else if (grepl("form ii|acla|aclb|cdh|coos", gene)) {
    return("Step 5: Carbon Fixation")
  } else {
    return("Step 6: Central Carbon & Others") 
  }
}

# 3. Apply categorization
carbon_data$Step_Category <- sapply(carbon_data$gene_name, categorize_gene_by_step)

# 4. Calculate relative abundance (%)
sample_cols <- setdiff(names(carbon_data), c("gene_name", "Step_Category"))
relative_data <- carbon_data
relative_data[sample_cols] <- sweep(relative_data[sample_cols], 2, colSums(relative_data[sample_cols]), "/") * 100

# 5. Aggregate & Merge
step_data <- relative_data %>%
  group_by(Step_Category) %>%
  summarise(across(all_of(sample_cols), sum), .groups = "drop")

long_data <- step_data %>% pivot_longer(cols = all_of(sample_cols), names_to = "sample", values_to = "Relative_Abundance")
merged_data <- long_data %>% inner_join(group_data, by = "sample")

final_summary <- merged_data %>%
  group_by(group1, Step_Category) %>%
  summarise(Mean_Abundance = mean(Relative_Abundance), .groups = "drop")

# 6. Reorder factors (Important for stacking order)
final_summary$Step_Category <- factor(final_summary$Step_Category, 
                                      levels = c("Step 6: Central Carbon & Others", 
                                                 "Step 5: Carbon Fixation", 
                                                 "Step 4: Fermentation & Acetate", 
                                                 "Step 3: C1 & Methane Metabolism", 
                                                 "Step 2: Aromatics Degradation", 
                                                 "Step 1: Organic Carbon Degradation"))

# 7. Plotting with strict EST/WR composition style and MANUAL colors
p <- ggplot(final_summary, aes(x = group1, y = Mean_Abundance, fill = Step_Category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6, color = "white", linewidth = 0.4) +
  # ✨ 这里是关键修改：调用您自定义的颜色 ✨
  scale_fill_manual(values = custom_colors) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  labs(
    x = "Experimental Group",
    y = "Mean Relative Abundance (%)",
    fill = "Carbon Cycle Pathways"
  ) +
  theme_classic() + 
  theme(
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text = element_text(color = "black", size = 12, face = "bold"),
    axis.title = element_text(color = "black", size = 14, face = "bold"),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    legend.position = "right",
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p)

# Save with parameters ideal for journal submission
ggsave("Custom_Colors_Carbon_Cycle.pdf", plot = p, width = 7.5, height = 5.5, device = "pdf")

