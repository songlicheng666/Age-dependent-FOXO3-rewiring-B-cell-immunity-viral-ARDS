# =============================================================================
# Script: 01_cell_proportion_calculation.R
# Description: Cell-type proportion calculation across four groups (major lineages and
# B-cell subtypes); Wilcoxon tests with Bonferroni correction; boxplot visualization
# Author: Licheng Song
# Institution: Chinese PLA General Hospital (Eighth Medical Center)
# Study: Age-dependent B-cell immune dysregulation in virus-associated ARDS
# Year: 2025
# Input: meta_data (Seurat metadata data.frame with orig.ident, group, celltype_final, celltype_fine)
# Output: 1_Control_Major_Raw.csv, 2_ARDS_T1_Major_Raw.csv, *_Boxplot.pdf
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)

# 1. 基础设置
setwd("D:/博士毕业/毕业文章/F2单细胞分群及细胞比例/细胞比例计算") 
# meta_data <- read.csv("计算细胞比例的metadata.csv")
sample_col <- "orig.ident" 

# ==========================================
# 任务 1: Control 组的大群 (Child vs Adult)
# ==========================================
df_1 <- meta_data %>% 
  filter(group %in% c("childcontrol", "adultcontrol"))

prop_1 <- df_1 %>%
  group_by(!!sym(sample_col), group, celltype_final) %>%
  summarise(count = n(), .groups = 'drop') %>%
  complete(!!sym(sample_col), celltype_final, fill = list(count = 0)) %>%
  group_by(!!sym(sample_col)) %>%
  mutate(group = unique(group[!is.na(group)])) %>% # 补全因complete丢失的group
  mutate(proportion = count / sum(count) * 100) %>%
  ungroup()

write_csv(prop_1, "1_Control_Major_Raw.csv")
print("已导出: 1_Control_Major_Raw.csv")

# ==========================================
# 任务 2: ARDS T1 组的大群 (Child vs Adult)
# ==========================================
df_2 <- meta_data %>% 
  filter(group %in% c("childARDS)T1", "adultARDS)T1"))

prop_2 <- df_2 %>%
  group_by(!!sym(sample_col), group, celltype_final) %>%
  summarise(count = n(), .groups = 'drop') %>%
  complete(!!sym(sample_col), celltype_final, fill = list(count = 0)) %>%
  group_by(!!sym(sample_col)) %>%
  mutate(group = unique(group[!is.na(group)])) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  ungroup()

write_csv(prop_2, "2_ARDS_T1_Major_Raw.csv")
print("已导出: 2_ARDS_T1_Major_Raw.csv")

# ==========================================
# 准备 B 细胞数据 (用于任务 3 和 4)
# ==========================================
# 筛选 B 细胞谱系
b_cell_meta <- meta_data %>% 
  filter(celltype_final %in% c("B_cell", "Plasma", "Plasmablast", "Memory_B", "Naive_B"))

# ==========================================
# 任务 3: Control 组的 B 细胞亚群
# ==========================================
df_3 <- b_cell_meta %>% 
  filter(group %in% c("childcontrol", "adultcontrol"))

prop_3 <- df_3 %>%
  group_by(!!sym(sample_col), group, celltype_fine) %>% # 注意这里是 celltype_fine
  summarise(count = n(), .groups = 'drop') %>%
  complete(!!sym(sample_col), celltype_fine, fill = list(count = 0)) %>%
  group_by(!!sym(sample_col)) %>%
  mutate(group = unique(group[!is.na(group)])) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  ungroup()

write_csv(prop_3, "3_Control_B_Subtypes_Raw.csv")
print("已导出: 3_Control_B_Subtypes_Raw.csv")

# ==========================================
# 任务 4: ARDS T1 组的 B 细胞亚群
# ==========================================
df_4 <- b_cell_meta %>% 
  filter(group %in% c("childARDS)T1", "adultARDS)T1"))

prop_4 <- df_4 %>%
  group_by(!!sym(sample_col), group, celltype_fine) %>% # 注意这里是 celltype_fine
  summarise(count = n(), .groups = 'drop') %>%
  complete(!!sym(sample_col), celltype_fine, fill = list(count = 0)) %>%
  group_by(!!sym(sample_col)) %>%
  mutate(group = unique(group[!is.na(group)])) %>%
  mutate(proportion = count / sum(count) * 100) %>%
  ungroup()

write_csv(prop_4, "4_ARDS_T1_B_Subtypes_Raw.csv")
print("已导出: 4_ARDS_T1_B_Subtypes_Raw.csv")













library(ggpubr)
library(readr)
library(dplyr)
library(tidyr)
library(rstatix)
library(ggplot2)

# 定义颜色
group_colors <- c(
  "childcontrol" = "#8DD3C7", 
  "adultcontrol" = "#BEBADA", 
  "childARDS)T1" = "#FB8072", 
  "adultARDS)T1" = "#80B1D3"
)

# 定义一个绘图函数 (仅用于画图和统计，不涉及数据计算)
run_stat_and_plot <- function(csv_path, group1, group2, cell_col, output_prefix, y_lab) {
 
  prop_df <- read_csv(csv_path, show_col_types = FALSE)
  
  # 2. 关键：重新设置 Group 因子顺序 (否则画图顺序会乱)
  prop_df$group <- factor(prop_df$group, levels = c(group1, group2))
  
  # 3. 统计检验
  stat_test <- prop_df %>%
    group_by(!!sym(cell_col)) %>%
    wilcox_test(proportion ~ group) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj") %>%
    add_xy_position(x = cell_col, dodge = 0.6)
  
  # 4. 描述性统计 (Mean/SD)
  desc_stats <- prop_df %>%
    group_by(!!sym(cell_col), group) %>%
    summarise(Mean = mean(proportion), SD = sd(proportion), .groups = 'drop') %>%
    pivot_wider(names_from = group, values_from = c(Mean, SD), names_glue = "{group}_{.value}")
  
  # 5. 导出统计结果表
  final_stats <- left_join(desc_stats, stat_test, by = cell_col) %>%
    select(!!sym(cell_col), contains("Mean"), contains("SD"), p, p.adj, p.adj.signif)
  
  write_csv(final_stats, paste0(output_prefix, "_Stats_Result.csv"))
  
  # 6. 绘图
  p <- ggplot(prop_df, aes(x = !!sym(cell_col), y = proportion)) + 
    geom_boxplot(aes(fill = group), width = 0.6, outlier.shape = NA, alpha = 0.8) +
    geom_point(aes(fill = group), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), 
               size = 1.2, alpha = 0.6, color = "black") +
    scale_fill_manual(values = group_colors) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      panel.grid.major = element_line(colour = "grey92"),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    ) +
    labs(x = "", y = y_lab, fill = "Group") +
    stat_pvalue_manual(stat_test, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE)
  
  pdf(paste0(output_prefix, "_Boxplot.pdf"), width = 6, height = 4)
  print(p)
  dev.off()
  
  message(paste0("完成: ", output_prefix))
}

# 1. Control 大群
run_stat_and_plot(
  csv_path = "1_Control_Major_Raw.csv",
  group1 = "childcontrol", group2 = "adultcontrol",
  cell_col = "celltype_final",
  output_prefix = "1_Control_Major_Final",
  y_lab = "Proportion (%)"
)

# 2. ARDS T1 大群
run_stat_and_plot(
  csv_path = "2_ARDS_T1_Major_Raw.csv",
  group1 = "childARDS)T1", group2 = "adultARDS)T1",
  cell_col = "celltype_final",
  output_prefix = "2_ARDS_T1_Major_Final",
  y_lab = "Proportion (%)"
)

# 3. Control B细胞亚群
run_stat_and_plot(
  csv_path = "3_Control_B_Subtypes_Raw.csv",
  group1 = "childcontrol", group2 = "adultcontrol",
  cell_col = "celltype_fine",
  output_prefix = "3_Control_B_Subtypes_Final",
  y_lab = "Proportion in B lineage (%)"
)

# 4. ARDS T1 B细胞亚群
run_stat_and_plot(
  csv_path = "4_ARDS_T1_B_Subtypes_Raw.csv",
  group1 = "childARDS)T1", group2 = "adultARDS)T1",
  cell_col = "celltype_fine",
  output_prefix = "4_ARDS_T1_B_Subtypes_Final",
  y_lab = "Proportion in B lineage (%)"
)
