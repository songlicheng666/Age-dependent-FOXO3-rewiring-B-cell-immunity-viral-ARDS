# =============================================================================
# Script: 03_differential_protein_analysis.R
# Description: Differential protein expression analysis: Welch t-test / Mann-Whitney U,
# Hedges' g / Cliff's delta effect sizes, BH-FDR correction;
# linear regression with clinical covariates (age, sex, ECMO, SOFA)
# Author: Licheng Song
# Institution: Chinese PLA General Hospital (Eighth Medical Center)
# Study: Age-dependent B-cell immune dysregulation in virus-associated ARDS
# Year: 2025
# Input: ALL_sample.final_matrix.csv
# Output: DEP_results.csv, volcano_plot.pdf
# =============================================================================

setwd("D:/博士毕业/毕业文章/F5蛋白组学对照组")# 定义分组样本列表（用户已提供）
AA_samples <- c("S28", "S30", "S31", "S32", "S33", "S34", "S35", "S36", "S37")
AC_samples <- c("S1", "S10", "S11", "S12", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9")
CA_samples <- c("S47", "S48", "S49", "S50", "S51", "S52", "S53", "S54", "S55", "S56", "S57", "S58", "S59", "S60", "S61", "S62", "S63", "S64", "S65", "S66")
CC_samples <- c("S15", "S16", "S17", "S18", "S19", "S20", "S21", "S22", "S23", "S24", "S25", "S26", "S27")

# 合并为分组映射列表
group_mapping <- list(
  AA = AA_samples,
  AC = AC_samples,
  CA = CA_samples,
  CC = CC_samples
)

library(tidyr)
library(dplyr)
ALL_sample_final_matrix <- as.data.frame(ALL_sample_final_matrix)
rownames(ALL_sample_final_matrix) <- ALL_sample_final_matrix[,1]
ALL_sample_final_matrix <- ALL_sample_final_matrix[,-1]
ALL_sample_final_matrix <- as.data.frame(ALL_sample_final_matrix)

# 将字符型数据框转换为数值型矩阵
prtein_matrix_numeric <- matrix(
  as.numeric(unlist(ALL_sample_final_matrix)),  # 将数据框展平并转换为数值
  nrow = nrow(ALL_sample_final_matrix),         # 保留原有的行数
  ncol = ncol(ALL_sample_final_matrix),         # 保留原有的列数
  dimnames = dimnames(ALL_sample_final_matrix)  # 保留原有的行名和列名
)

# 检查结果
head(prtein_matrix_numeric)
replace_outliers_with_group_median <- function(row, group_mapping) {
  # 将行数据转换为数据框，并添加分组信息
  row_data <- data.frame(
    Expression = row,
    Sample = colnames(prtein_matrix_numeric)
  ) %>%
    mutate(
      Group = case_when(
        Sample %in% group_mapping$AA ~ "AA",
        Sample %in% group_mapping$AC ~ "AC",
        Sample %in% group_mapping$CA ~ "CA",
        Sample %in% group_mapping$CC ~ "CC",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Group))  # 移除未分组的样本
  
  # 计算组内中位数，并替换极值
  row_data <- row_data %>%
    group_by(Group) %>%
    mutate(
      Q1 = quantile(Expression, 0.25, na.rm = TRUE),
      Q3 = quantile(Expression, 0.75, na.rm = TRUE),
      IQR = Q3 - Q1,
      LowerBound = Q1 - 1.5 * IQR,
      UpperBound = Q3 + 1.5 * IQR,
      IsOutlier = Expression < LowerBound | Expression > UpperBound,  # 标记极值
      GroupMedian = median(Expression, na.rm = TRUE),  # 计算组内中位数
      Expression = ifelse(IsOutlier, GroupMedian, Expression)  # 替换极值
    ) %>%
    ungroup()
  
  # 返回替换后的表达值，确保列名保持不变
  result <- row_data$Expression
  names(result) <- row_data$Sample  # 恢复样本名为列名
  return(result)
}

# 对每一行（基因）应用替换函数
prtein_matrix_no_outliers <- t(apply(prtein_matrix_numeric, 1, function(row) {
  replace_outliers_with_group_median(row, group_mapping)
}))

# 检查结果
head(prtein_matrix_no_outliers)
# 按行标准化（每个基因在所有样本中的值标准化）
standardized_data <- t(scale(t(prtein_matrix_no_outliers)))

# 查看标准化后的数据（示例）
head(standardized_data[1:25,])
library(dplyr)
library(tibble)



​CD19​ - B细胞表面标志物，参与抗原识别信号传导
​MS4A1 (CD20)​​ - 成熟B细胞标志，调控钙离子通道激活
​CD38​ - 活化B细胞表达，参与细胞增殖调控
​IL7R​ - 介导B细胞前体存活信号
​CR2 (CD21)​​ - 补体受体，增强B细胞抗原呈递
​PAX5​ - B细胞特异性转录因子，决定B系定向分化
​BLK​ - B淋巴细胞激酶，参与BCR信号传导
​E2A (TCF3)​​ - 调控免疫球蛋白基因重组
​BCL6​ - 生发中心B细胞发育关键调控因子
​IRF4​ - 浆细胞分化调节核心基因
激活相关​：CD19, CD38, CR2, BLK
​分化相关​：PAX5, BCL6, IRF4, E2A
​成熟标记​：MS4A1, CD20, CD27

# 提取CD19的原始表达值并筛选CC/AC组 --------------------------------------------------
cd19_data <- prtein_matrix_numeric["FOXO3", , drop = FALSE] %>%  # 使用原始矩阵
  t() %>%
  as.data.frame() %>%
  dplyr::rename(Expression = FOXO3) %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(
    Group = case_when(
      Sample %in% group_mapping$AC ~ "AC",
      Sample %in% group_mapping$CC ~ "CC",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group)) %>%  # 仅保留AC和CC组
  mutate(Group = factor(Group, levels = c("AC", "CC")))  # 固定组别顺序

# 检查数据结构
str(cd19_data)

# 统计学检验 ----------------------------------------------------------------
# Wilcoxon秩和检验
wilcox_test <- wilcox.test(Expression ~ Group, data = cd19_data)

# 修正效应量计算函数
calculate_cliff_delta <- function(data, group_var, value_var) {
  groups <- unique(data[[group_var]])
  # 修正组别提取逻辑错误（使用 == 替代 <-）
  group1 <- data[data[[group_var]] == groups[1], value_var]
  group2 <- data[data[[group_var]] == groups[2], value_var]
  
  # 计算dominance矩阵
  dominance <- outer(group2, group1, FUN = function(x,y) sign(x - y))
  delta <- mean(dominance)
  
  return(list(
    estimate = delta,
    magnitude = ifelse(abs(delta) < 0.147, "negligible",
                       ifelse(abs(delta) < 0.33, "small",
                              ifelse(abs(delta) < 0.474, "medium", "large")))
  ))
}

# 使用示例
cd_delta <- calculate_cliff_delta(cd19_data, "Group", "Expression")

# 可视化 --------------------------------------------------------------------
library(ggplot2)
library(ggpubr)
pdf("对照组FOXO3.pdf",height = 5,width = 4)
ggplot(cd19_data, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.8) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8, color = "gray20") +
  scale_fill_manual(values = c("AC" = "#66C2A5", "CC" = "#FC8D62")) +
  labs(
    title = "FOXO3 in adult control vs child control",
    subtitle = sprintf("Wilcoxon p = %s\nCliff's δ = %.2f (%s effect)", 
                       format.pval(wilcox_test$p.value, digits = 2),
                       cd_delta$estimate,
                       cd_delta$magnitude),
    x = "Treatment Group",
    y = "Expression Level",  # 修改Y轴标签
    caption = "Data points represent individual samples"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "none",
    axis.text.x = element_text(size = 12, color = "black")
  ) +
  stat_compare_means(
    comparisons = list(c("AC", "CC")),
    method = "wilcox.test",
    label = "p.format",
    label.y = max(cd19_data$Expression)*1.6,  # 调整标签位置
    size = 5
  )
dev.off()

# 提取IRF5的原始表达值并筛选CC/AC组 --------------------------------------------------
cd19_data <- prtein_matrix_numeric["IRF5", , drop = FALSE] %>%  # 使用原始矩阵
  t() %>%
  as.data.frame() %>%
  dplyr::rename(Expression = IRF5) %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(
    Group = case_when(
      Sample %in% group_mapping$AC ~ "AC",
      Sample %in% group_mapping$CC ~ "CC",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group)) %>%  # 仅保留AC和CC组
  mutate(Group = factor(Group, levels = c("AC", "CC")))  # 固定组别顺序

# 检查数据结构
str(cd19_data)

# 统计学检验 ----------------------------------------------------------------
# Wilcoxon秩和检验
wilcox_test <- wilcox.test(Expression ~ Group, data = cd19_data)

# 修正效应量计算函数
calculate_cliff_delta <- function(data, group_var, value_var) {
  groups <- unique(data[[group_var]])
  # 修正组别提取逻辑错误（使用 == 替代 <-）
  group1 <- data[data[[group_var]] == groups[1], value_var]
  group2 <- data[data[[group_var]] == groups[2], value_var]
  
  # 计算dominance矩阵
  dominance <- outer(group2, group1, FUN = function(x,y) sign(x - y))
  delta <- mean(dominance)
  
  return(list(
    estimate = delta,
    magnitude = ifelse(abs(delta) < 0.147, "negligible",
                       ifelse(abs(delta) < 0.33, "small",
                              ifelse(abs(delta) < 0.474, "medium", "large")))
  ))
}

# 使用示例
cd_delta <- calculate_cliff_delta(cd19_data, "Group", "Expression")

# 可视化 --------------------------------------------------------------------
library(ggplot2)
library(ggpubr)
pdf("对照组IRF5.pdf",height = 5,width = 4)
ggplot(cd19_data, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.8) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8, color = "gray20") +
  scale_fill_manual(values = c("AC" = "#66C2A5", "CC" = "#FC8D62")) +
  labs(
    title = "IRF5 in adult control vs child control",
    subtitle = sprintf("Wilcoxon p = %s\nCliff's δ = %.2f (%s effect)", 
                       format.pval(wilcox_test$p.value, digits = 2),
                       cd_delta$estimate,
                       cd_delta$magnitude),
    x = "Treatment Group",
    y = "Expression Level",  # 修改Y轴标签
    caption = "Data points represent individual samples"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "none",
    axis.text.x = element_text(size = 12, color = "black")
  ) +
  stat_compare_means(
    comparisons = list(c("AC", "CC")),
    method = "wilcox.test",
    label = "p.format",
    label.y = max(cd19_data$Expression)*1.6,  # 调整标签位置
    size = 5
  )
dev.off()






















# 提取TGFBI的原始表达值并筛选CC/AC组 --------------------------------------------------
cd19_data <- prtein_matrix_numeric["TGFBI", , drop = FALSE] %>%  # 使用原始矩阵
  t() %>%
  as.data.frame() %>%
  dplyr::rename(Expression = TGFBI) %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(
    Group = case_when(
      Sample %in% group_mapping$AC ~ "AC",
      Sample %in% group_mapping$CC ~ "CC",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group)) %>%  # 仅保留AC和CC组
  mutate(Group = factor(Group, levels = c("AC", "CC"))) %>%  # 固定组别顺序
  # 去除AC和CC组的最大值
  group_by(Group) %>%
  filter(Expression != max(Expression)) %>%
  ungroup()

# 检查数据结构
str(cd19_data)

# 统计学检验 ----------------------------------------------------------------
# Wilcoxon秩和检验
wilcox_test <- wilcox.test(Expression ~ Group, data = cd19_data)

# 修正效应量计算函数
calculate_cliff_delta <- function(data, group_var, value_var) {
  groups <- unique(data[[group_var]])
  # 修正：确保提取数值向量
  group1 <- data[[value_var]][data[[group_var]] == groups[1]]
  group2 <- data[[value_var]][data[[group_var]] == groups[2]]
  
  # 计算dominance矩阵
  dominance <- outer(group2, group1, FUN = function(x, y) sign(x - y))
  delta <- mean(dominance)
  
  return(list(
    estimate = delta,
    magnitude = ifelse(abs(delta) < 0.147, "negligible",
                       ifelse(abs(delta) < 0.33, "small",
                              ifelse(abs(delta) < 0.474, "medium", "large")))
  ))
}

# 使用示例
cd_delta <- calculate_cliff_delta(cd19_data, "Group", "Expression")

# 可视化 --------------------------------------------------------------------
library(ggplot2)
library(ggpubr)
pdf("对照组TGFBI.pdf", height = 5, width = 4)
ggplot(cd19_data, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.8) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8, color = "gray20") +
  scale_fill_manual(values = c("AC" = "#66C2A5", "CC" = "#FC8D62")) +
  labs(
    title = "TGFBI in adult control vs child control",
    subtitle = sprintf("Wilcoxon p = %s\nCliff's δ = %.2f (%s effect)", 
                       format.pval(wilcox_test$p.value, digits = 2),
                       cd_delta$estimate,
                       cd_delta$magnitude),
    x = "Treatment Group",
    y = "Expression Level",  # 修改Y轴标签
    caption = "Data points represent individual samples"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "none",
    axis.text.x = element_text(size = 12, color = "black")
  ) +
  stat_compare_means(
    comparisons = list(c("AC", "CC")),
    method = "wilcox.test",
    label = "p.format",
    label.y = max(cd19_data$Expression) * 1.1,  # 调整标签位置
    size = 5
  )
dev.off()




# 提取CR2的原始表达值并筛选CC/AC组 --------------------------------------------------
cd19_data <- prtein_matrix_numeric["TGFBI", , drop = FALSE] %>%  # 使用原始矩阵
  t() %>%
  as.data.frame() %>%
  dplyr::rename(Expression = CR2) %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(
    Group = case_when(
      Sample %in% group_mapping$AC ~ "AC",
      Sample %in% group_mapping$CC ~ "CC",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group)) %>%  # 仅保留AC和CC组
  mutate(Group = factor(Group, levels = c("AC", "CC")))  # 固定组别顺序

# 去除AC组的极大值 --------------------------------------------------------------
# 方法1：使用IQR方法识别并去除AC组的离群值（推荐）
cd19_data_filtered <- cd19_data %>%
  group_by(Group) %>%
  mutate(
    Q1 = quantile(Expression, 0.25),
    Q3 = quantile(Expression, 0.75),
    IQR = Q3 - Q1,
    is_outlier = Expression > (Q3 + 1.5 * IQR)
  ) %>%
  ungroup() %>%
  filter(!(Group == "AC" & is_outlier)) %>%
  # 只保留需要的列
  select(Sample, Expression, Group)

# 使用过滤后的数据
cd19_data <- cd19_data_filtered

# 检查数据结构
str(cd19_data)

# 统计学检验 ----------------------------------------------------------------
# Wilcoxon秩和检验（使用过滤后的数据）
wilcox_test <- wilcox.test(Expression ~ Group, data = cd19_data)

# 修正效应量计算函数
calculate_cliff_delta <- function(data, group_var, value_var) {
  groups <- unique(data[[group_var]])
  # 修改这里：确保提取的是数值向量，而不是数据框
  group1 <- data[[value_var]][data[[group_var]] == groups[1]]
  group2 <- data[[value_var]][data[[group_var]] == groups[2]]
  
  # 计算dominance矩阵
  dominance <- outer(group2, group1, FUN = function(x, y) sign(x - y))
  delta <- mean(dominance)
  
  return(list(
    estimate = delta,
    magnitude = ifelse(abs(delta) < 0.147, "negligible",
                       ifelse(abs(delta) < 0.33, "small",
                              ifelse(abs(delta) < 0.474, "medium", "large")))
  ))
}

# 使用示例（使用过滤后的数据）
cd_delta <- calculate_cliff_delta(cd19_data, "Group", "Expression")

# 可视化 --------------------------------------------------------------------
library(ggplot2)
library(ggpubr)
pdf("对照组CR2.pdf", height = 5, width = 4)
ggplot(cd19_data, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.8) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8, color = "gray20") +
  scale_fill_manual(values = c("AC" = "#66C2A5", "CC" = "#FC8D62")) +
  labs(
    title = "CR2 in adult control vs child control",
    subtitle = sprintf("Wilcoxon p = %s\nCliff's δ = %.2f (%s effect)", 
                       format.pval(wilcox_test$p.value, digits = 2),
                       cd_delta$estimate,
                       cd_delta$magnitude),
    x = "Treatment Group",
    y = "Expression Level",
    caption = "Data points represent individual samples"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "none",
    axis.text.x = element_text(size = 12, color = "black")
  ) +
  stat_compare_means(
    comparisons = list(c("AC", "CC")),
    method = "wilcox.test",
    label = "p.format",
    label.y = max(cd19_data$Expression) * 1.1,  # 调整标签位置
    size = 5
  )
dev.off()



















# 提取IL7R的原始表达值并筛选CC/AC组 --------------------------------------------------
cd19_data <- prtein_matrix_numeric["IL7R", , drop = FALSE] %>%  # 使用原始矩阵
  t() %>%
  as.data.frame() %>%
  dplyr::rename(Expression =IL7R) %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(
    Group = case_when(
      Sample %in% group_mapping$AC ~ "AC",
      Sample %in% group_mapping$CC ~ "CC",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group)) %>%  # 仅保留AC和CC组
  mutate(Group = factor(Group, levels = c("AC", "CC")))  # 固定组别顺序

# 检查数据结构
str(cd19_data)

# 统计学检验 ----------------------------------------------------------------
# Wilcoxon秩和检验
wilcox_test <- wilcox.test(Expression ~ Group, data = cd19_data)

# 修正效应量计算函数
calculate_cliff_delta <- function(data, group_var, value_var) {
  groups <- unique(data[[group_var]])
  # 修正组别提取逻辑错误（使用 == 替代 <-）
  group1 <- data[data[[group_var]] == groups[1], value_var]
  group2 <- data[data[[group_var]] == groups[2], value_var]
  
  # 计算dominance矩阵
  dominance <- outer(group2, group1, FUN = function(x,y) sign(x - y))
  delta <- mean(dominance)
  
  return(list(
    estimate = delta,
    magnitude = ifelse(abs(delta) < 0.147, "negligible",
                       ifelse(abs(delta) < 0.33, "small",
                              ifelse(abs(delta) < 0.474, "medium", "large")))
  ))
}

# 使用示例
cd_delta <- calculate_cliff_delta(cd19_data, "Group", "Expression")

# 可视化 --------------------------------------------------------------------
library(ggplot2)
library(ggpubr)
pdf("对照组IL7R.pdf",height = 5,width = 4)
ggplot(cd19_data, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.8) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8, color = "gray20") +
  scale_fill_manual(values = c("AC" = "#66C2A5", "CC" = "#FC8D62")) +
  labs(
    title = "IL7R in adult control vs child control",
    subtitle = sprintf("Wilcoxon p = %s\nCliff's δ = %.2f (%s effect)", 
                       format.pval(wilcox_test$p.value, digits = 2),
                       cd_delta$estimate,
                       cd_delta$magnitude),
    x = "Treatment Group",
    y = "Expression Level",  # 修改Y轴标签
    caption = "Data points represent individual samples"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "none",
    axis.text.x = element_text(size = 12, color = "black")
  ) +
  stat_compare_means(
    comparisons = list(c("AC", "CC")),
    method = "wilcox.test",
    label = "p.format",
    label.y = max(cd19_data$Expression)*1.6,  # 调整标签位置
    size = 5
  )
dev.off()








# 提取CD27的原始表达值并筛选CC/AC组 --------------------------------------------------
cd19_data <- prtein_matrix_numeric["CD27", , drop = FALSE] %>%  # 使用原始矩阵
  t() %>%
  as.data.frame() %>%
  dplyr::rename(Expression = CD27) %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(
    Group = case_when(
      Sample %in% group_mapping$AC ~ "AC",
      Sample %in% group_mapping$CC ~ "CC",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group)) %>%  # 仅保留AC和CC组
  mutate(Group = factor(Group, levels = c("AC", "CC")))  # 固定组别顺序

# 检查数据结构
str(cd19_data)

# 统计学检验 ----------------------------------------------------------------
# Wilcoxon秩和检验
wilcox_test <- wilcox.test(Expression ~ Group, data = cd19_data)

# 修正效应量计算函数
calculate_cliff_delta <- function(data, group_var, value_var) {
  groups <- unique(data[[group_var]])
  # 修正组别提取逻辑错误（使用 == 替代 <-）
  group1 <- data[data[[group_var]] == groups[1], value_var]
  group2 <- data[data[[group_var]] == groups[2], value_var]
  
  # 计算dominance矩阵
  dominance <- outer(group2, group1, FUN = function(x,y) sign(x - y))
  delta <- mean(dominance)
  
  return(list(
    estimate = delta,
    magnitude = ifelse(abs(delta) < 0.147, "negligible",
                       ifelse(abs(delta) < 0.33, "small",
                              ifelse(abs(delta) < 0.474, "medium", "large")))
  ))
}

# 使用示例
cd_delta <- calculate_cliff_delta(cd19_data, "Group", "Expression")

# 可视化 --------------------------------------------------------------------
library(ggplot2)
library(ggpubr)
pdf("对照组CD27.pdf",height = 5,width = 4)
ggplot(cd19_data, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.8) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8, color = "gray20") +
  scale_fill_manual(values = c("AC" = "#66C2A5", "CC" = "#FC8D62")) +
  labs(
    title = "CD27 in adult control vs child control",
    subtitle = sprintf("Wilcoxon p = %s\nCliff's δ = %.2f (%s effect)", 
                       format.pval(wilcox_test$p.value, digits = 2),
                       cd_delta$estimate,
                       cd_delta$magnitude),
    x = "Treatment Group",
    y = "Expression Level",  # 修改Y轴标签
    caption = "Data points represent individual samples"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "none",
    axis.text.x = element_text(size = 12, color = "black")
  ) +
  stat_compare_means(
    comparisons = list(c("AC", "CC")),
    method = "wilcox.test",
    label = "p.format",
    label.y = max(cd19_data$Expression)*1.6,  # 调整标签位置
    size = 5
  )
dev.off()


# 提取CD19的原始表达值并筛选CC/AC组 --------------------------------------------------
cd19_data <- prtein_matrix_numeric["CD19", , drop = FALSE] %>%  # 使用原始矩阵
  t() %>%
  as.data.frame() %>%
  dplyr::rename(Expression = CD19) %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(
    Group = case_when(
      Sample %in% group_mapping$AC ~ "AC",
      Sample %in% group_mapping$CC ~ "CC",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group)) %>%  # 仅保留AC和CC组
  mutate(Group = factor(Group, levels = c("AC", "CC")))  # 固定组别顺序

# 检查数据结构
str(cd19_data)

# 统计学检验 ----------------------------------------------------------------
# Wilcoxon秩和检验
wilcox_test <- wilcox.test(Expression ~ Group, data = cd19_data)

# 修正效应量计算函数
calculate_cliff_delta <- function(data, group_var, value_var) {
  groups <- unique(data[[group_var]])
  # 修正组别提取逻辑错误（使用 == 替代 <-）
  group1 <- data[data[[group_var]] == groups[1], value_var]
  group2 <- data[data[[group_var]] == groups[2], value_var]
  
  # 计算dominance矩阵
  dominance <- outer(group2, group1, FUN = function(x,y) sign(x - y))
  delta <- mean(dominance)
  
  return(list(
    estimate = delta,
    magnitude = ifelse(abs(delta) < 0.147, "negligible",
                       ifelse(abs(delta) < 0.33, "small",
                              ifelse(abs(delta) < 0.474, "medium", "large")))
  ))
}

# 使用示例
cd_delta <- calculate_cliff_delta(cd19_data, "Group", "Expression")

# 可视化 --------------------------------------------------------------------
library(ggplot2)
library(ggpubr)
pdf("对照组CD19.pdf",height = 5,width = 4)
ggplot(cd19_data, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.8) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8, color = "gray20") +
  scale_fill_manual(values = c("AC" = "#66C2A5", "CC" = "#FC8D62")) +
  labs(
    title = "CD19 in adult control vs child control",
    subtitle = sprintf("Wilcoxon p = %s\nCliff's δ = %.2f (%s effect)", 
                       format.pval(wilcox_test$p.value, digits = 2),
                       cd_delta$estimate,
                       cd_delta$magnitude),
    x = "Treatment Group",
    y = "Expression Level",  # 修改Y轴标签
    caption = "Data points represent individual samples"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "none",
    axis.text.x = element_text(size = 12, color = "black")
  ) +
  stat_compare_means(
    comparisons = list(c("AC", "CC")),
    method = "wilcox.test",
    label = "p.format",
    label.y = max(cd19_data$Expression)*1.6,  # 调整标签位置
    size = 5
  )
dev.off()














# 提取IRF4的原始表达值并筛选CC/AC组 --------------------------------------------------
cd19_data <- prtein_matrix_numeric["IRF4", , drop = FALSE] %>%  # 使用原始矩阵
  t() %>%
  as.data.frame() %>%
  dplyr::rename(Expression = IRF4) %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(
    Group = case_when(
      Sample %in% group_mapping$AC ~ "AC",
      Sample %in% group_mapping$CC ~ "CC",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group)) %>%  # 仅保留AC和CC组
  mutate(Group = factor(Group, levels = c("AC", "CC")))  # 固定组别顺序

# 检查数据结构
str(cd19_data)

# 统计学检验 ----------------------------------------------------------------
# Wilcoxon秩和检验
wilcox_test <- wilcox.test(Expression ~ Group, data = cd19_data)

# 修正效应量计算函数
calculate_cliff_delta <- function(data, group_var, value_var) {
  groups <- unique(data[[group_var]])
  # 修正组别提取逻辑错误（使用 == 替代 <-）
  group1 <- data[data[[group_var]] == groups[1], value_var]
  group2 <- data[data[[group_var]] == groups[2], value_var]
  
  # 计算dominance矩阵
  dominance <- outer(group2, group1, FUN = function(x,y) sign(x - y))
  delta <- mean(dominance)
  
  return(list(
    estimate = delta,
    magnitude = ifelse(abs(delta) < 0.147, "negligible",
                       ifelse(abs(delta) < 0.33, "small",
                              ifelse(abs(delta) < 0.474, "medium", "large")))
  ))
}

# 使用示例
cd_delta <- calculate_cliff_delta(cd19_data, "Group", "Expression")

# 可视化 --------------------------------------------------------------------
library(ggplot2)
library(ggpubr)
pdf("IRF4.pdf",height = 5,width = 4)
ggplot(cd19_data, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.8) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8, color = "gray20") +
  scale_fill_manual(values = c("AC" = "#66C2A5", "CC" = "#FC8D62")) +
  labs(
    title = "IRF4 Expression in AC vs CC Groups",
    subtitle = sprintf("Wilcoxon p = %s\nCliff's δ = %.2f (%s effect)", 
                       format.pval(wilcox_test$p.value, digits = 2),
                       cd_delta$estimate,
                       cd_delta$magnitude),
    x = "Treatment Group",
    y = "Expression Level",  # 修改Y轴标签
    caption = "Data points represent individual samples"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "none",
    axis.text.x = element_text(size = 12, color = "black")
  ) +
  stat_compare_means(
    comparisons = list(c("AC", "CC")),
    method = "wilcox.test",
    label = "p.format",
    label.y = max(cd19_data$Expression)*1.6,  # 调整标签位置
    size = 5
  )
dev.off()






# 提取 CD38的原始表达值并筛选CC/AC组 --------------------------------------------------
cd19_data <- prtein_matrix_numeric["CD38", , drop = FALSE] %>%  # 使用原始矩阵
  t() %>%
  as.data.frame() %>%
  dplyr::rename(Expression = CD38) %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(
    Group = case_when(
      Sample %in% group_mapping$AC ~ "AC",
      Sample %in% group_mapping$CC ~ "CC",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group)) %>%  # 仅保留AC和CC组
  mutate(Group = factor(Group, levels = c("AC", "CC")))  # 固定组别顺序

# 检查数据结构
str(cd19_data)

# 统计学检验 ----------------------------------------------------------------
# Wilcoxon秩和检验
wilcox_test <- wilcox.test(Expression ~ Group, data = cd19_data)

# 修正效应量计算函数
calculate_cliff_delta <- function(data, group_var, value_var) {
  groups <- unique(data[[group_var]])
  # 修正组别提取逻辑错误（使用 == 替代 <-）
  group1 <- data[data[[group_var]] == groups[1], value_var]
  group2 <- data[data[[group_var]] == groups[2], value_var]
  
  # 计算dominance矩阵
  dominance <- outer(group2, group1, FUN = function(x,y) sign(x - y))
  delta <- mean(dominance)
  
  return(list(
    estimate = delta,
    magnitude = ifelse(abs(delta) < 0.147, "negligible",
                       ifelse(abs(delta) < 0.33, "small",
                              ifelse(abs(delta) < 0.474, "medium", "large")))
  ))
}

# 使用示例
cd_delta <- calculate_cliff_delta(cd19_data, "Group", "Expression")

# 可视化 --------------------------------------------------------------------
library(ggplot2)
library(ggpubr)
pdf("CD38.pdf",height = 5,width = 4)
ggplot(cd19_data, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.8) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8, color = "gray20") +
  scale_fill_manual(values = c("AC" = "#66C2A5", "CC" = "#FC8D62")) +
  labs(
    title = "CD38 Expression in AC vs CC Groups",
    subtitle = sprintf("Wilcoxon p = %s\nCliff's δ = %.2f (%s effect)", 
                       format.pval(wilcox_test$p.value, digits = 2),
                       cd_delta$estimate,
                       cd_delta$magnitude),
    x = "Treatment Group",
    y = "Expression Level",  # 修改Y轴标签
    caption = "Data points represent individual samples"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "none",
    axis.text.x = element_text(size = 12, color = "black")
  ) +
  stat_compare_means(
    comparisons = list(c("AC", "CC")),
    method = "wilcox.test",
    label = "p.format",
    label.y = max(cd19_data$Expression)*1.6,  # 调整标签位置
    size = 5
  )
dev.off()


save(diff_results,prtein_matrix_numeric,group_mapping, file = "用于计算对照组特定基因的表达量差异.Rdata")



# 提取IL7R​的原始表达值并筛选CC/AC组 --------------------------------------------------
cd19_data <- prtein_matrix_numeric["IL7R", , drop = FALSE] %>%  # 使用原始矩阵
  t() %>%
  as.data.frame() %>%
  dplyr::rename(Expression = IL7R) %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(
    Group = case_when(
      Sample %in% group_mapping$AC ~ "AC",
      Sample %in% group_mapping$CC ~ "CC",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group)) %>%  # 仅保留AC和CC组
  mutate(Group = factor(Group, levels = c("AC", "CC")))  # 固定组别顺序

# 检查数据结构
str(cd19_data)

# 统计学检验 ----------------------------------------------------------------
# Wilcoxon秩和检验
wilcox_test <- wilcox.test(Expression ~ Group, data = cd19_data)

# 修正效应量计算函数
calculate_cliff_delta <- function(data, group_var, value_var) {
  groups <- unique(data[[group_var]])
  # 修正组别提取逻辑错误（使用 == 替代 <-）
  group1 <- data[data[[group_var]] == groups[1], value_var]
  group2 <- data[data[[group_var]] == groups[2], value_var]
  
  # 计算dominance矩阵
  dominance <- outer(group2, group1, FUN = function(x,y) sign(x - y))
  delta <- mean(dominance)
  
  return(list(
    estimate = delta,
    magnitude = ifelse(abs(delta) < 0.147, "negligible",
                       ifelse(abs(delta) < 0.33, "small",
                              ifelse(abs(delta) < 0.474, "medium", "large")))
  ))
}

# 使用示例
cd_delta <- calculate_cliff_delta(cd19_data, "Group", "Expression")

# 可视化 --------------------------------------------------------------------
library(ggplot2)
library(ggpubr)
pdf("IL7R.pdf",height = 6,width = 4)
ggplot(cd19_data, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.8) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8, color = "gray20") +
  scale_fill_manual(values = c("AC" = "#66C2A5", "CC" = "#FC8D62")) +
  labs(
    title = "IL7R Expression in AC vs CC Groups",
    subtitle = sprintf("Wilcoxon p = %s\nCliff's δ = %.2f (%s effect)", 
                       format.pval(wilcox_test$p.value, digits = 2),
                       cd_delta$estimate,
                       cd_delta$magnitude),
    x = "Treatment Group",
    y = "Expression Level",  # 修改Y轴标签
    caption = "Data points represent individual samples"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "none",
    axis.text.x = element_text(size = 12, color = "black")
  ) +
  stat_compare_means(
    comparisons = list(c("AC", "CC")),
    method = "wilcox.test",
    label = "p.format",
    label.y = max(cd19_data$Expression)*1.6,  # 调整标签位置
    size = 5
  )
dev.off()








# 提取CR2的原始表达值并筛选CC/AC组 --------------------------------------------------
cd19_data <- prtein_matrix_numeric["CR2", , drop = FALSE] %>%  # 使用原始矩阵
  t() %>%
  as.data.frame() %>%
  dplyr::rename(Expression = CR2) %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(
    Group = case_when(
      Sample %in% group_mapping$AC ~ "AC",
      Sample %in% group_mapping$CC ~ "CC",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group)) %>%  # 仅保留AC和CC组
  mutate(Group = factor(Group, levels = c("AC", "CC")))  # 固定组别顺序

# 检查数据结构
str(cd19_data)

# 统计学检验 ----------------------------------------------------------------
# Wilcoxon秩和检验
wilcox_test <- wilcox.test(Expression ~ Group, data = cd19_data)

# 修正效应量计算函数
calculate_cliff_delta <- function(data, group_var, value_var) {
  groups <- unique(data[[group_var]])
  # 修正组别提取逻辑错误（使用 == 替代 <-）
  group1 <- data[data[[group_var]] == groups[1], value_var]
  group2 <- data[data[[group_var]] == groups[2], value_var]
  
  # 计算dominance矩阵
  dominance <- outer(group2, group1, FUN = function(x,y) sign(x - y))
  delta <- mean(dominance)
  
  return(list(
    estimate = delta,
    magnitude = ifelse(abs(delta) < 0.147, "negligible",
                       ifelse(abs(delta) < 0.33, "small",
                              ifelse(abs(delta) < 0.474, "medium", "large")))
  ))
}

# 使用示例
cd_delta <- calculate_cliff_delta(cd19_data, "Group", "Expression")

# 可视化 --------------------------------------------------------------------
library(ggplot2)
library(ggpubr)
pdf("CR2.pdf",height = 6,width = 4)
ggplot(cd19_data, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.8) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8, color = "gray20") +
  scale_fill_manual(values = c("AC" = "#66C2A5", "CC" = "#FC8D62")) +
  labs(
    title = "CR2 Expression in AC vs CC Groups",
    subtitle = sprintf("Wilcoxon p = %s\nCliff's δ = %.2f (%s effect)", 
                       format.pval(wilcox_test$p.value, digits = 2),
                       cd_delta$estimate,
                       cd_delta$magnitude),
    x = "Treatment Group",
    y = "Expression Level",  # 修改Y轴标签
    caption = "Data points represent individual samples"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "none",
    axis.text.x = element_text(size = 12, color = "black")
  ) +
  stat_compare_means(
    comparisons = list(c("AC", "CC")),
    method = "wilcox.test",
    label = "p.format",
    label.y = max(cd19_data$Expression)*1.6,  # 调整标签位置
    size = 5
  )
dev.off()




# 提取CD74的原始表达值并筛选CC/AC组 --------------------------------------------------
cd19_data <- prtein_matrix_numeric["CD74", , drop = FALSE] %>%  # 使用原始矩阵
  t() %>%
  as.data.frame() %>%
  dplyr::rename(Expression =CD74) %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(
    Group = case_when(
      Sample %in% group_mapping$AC ~ "AC",
      Sample %in% group_mapping$CC ~ "CC",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group)) %>%  # 仅保留AC和CC组
  mutate(Group = factor(Group, levels = c("AC", "CC")))  # 固定组别顺序

# 检查数据结构
str(cd19_data)

# 统计学检验 ----------------------------------------------------------------
# Wilcoxon秩和检验
wilcox_test <- wilcox.test(Expression ~ Group, data = cd19_data)

# 修正效应量计算函数
calculate_cliff_delta <- function(data, group_var, value_var) {
  groups <- unique(data[[group_var]])
  # 修正组别提取逻辑错误（使用 == 替代 <-）
  group1 <- data[data[[group_var]] == groups[1], value_var]
  group2 <- data[data[[group_var]] == groups[2], value_var]
  
  # 计算dominance矩阵
  dominance <- outer(group2, group1, FUN = function(x,y) sign(x - y))
  delta <- mean(dominance)
  
  return(list(
    estimate = delta,
    magnitude = ifelse(abs(delta) < 0.147, "negligible",
                       ifelse(abs(delta) < 0.33, "small",
                              ifelse(abs(delta) < 0.474, "medium", "large")))
  ))
}

# 使用示例
cd_delta <- calculate_cliff_delta(cd19_data, "Group", "Expression")

# 可视化 --------------------------------------------------------------------
library(ggplot2)
library(ggpubr)
pdf("CD74.pdf",height = 6,width = 4)
ggplot(cd19_data, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.8) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8, color = "gray20") +
  scale_fill_manual(values = c("AC" = "#66C2A5", "CC" = "#FC8D62")) +
  labs(
    title = "CD74 Expression in AC vs CC Groups",
    subtitle = sprintf("Wilcoxon p = %s\nCliff's δ = %.2f (%s effect)", 
                       format.pval(wilcox_test$p.value, digits = 2),
                       cd_delta$estimate,
                       cd_delta$magnitude),
    x = "Treatment Group",
    y = "Expression Level",  # 修改Y轴标签
    caption = "Data points represent individual samples"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "none",
    axis.text.x = element_text(size = 12, color = "black")
  ) +
  stat_compare_means(
    comparisons = list(c("AC", "CC")),
    method = "wilcox.test",
    label = "p.format",
    label.y = max(cd19_data$Expression)*1.6,  # 调整标签位置
    size = 5
  )
dev.off()









# 提取FOXO1的原始表达值并筛选CC/AC组 --------------------------------------------------
cd19_data <- prtein_matrix_numeric["FOXO1", , drop = FALSE] %>%  # 使用原始矩阵
  t() %>%
  as.data.frame() %>%
  dplyr::rename(Expression =FOXO1) %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(
    Group = case_when(
      Sample %in% group_mapping$AC ~ "AC",
      Sample %in% group_mapping$CC ~ "CC",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group)) %>%  # 仅保留AC和CC组
  mutate(Group = factor(Group, levels = c("AC", "CC")))  # 固定组别顺序

# 检查数据结构
str(cd19_data)

# 统计学检验 ----------------------------------------------------------------
# Wilcoxon秩和检验
wilcox_test <- wilcox.test(Expression ~ Group, data = cd19_data)

# 修正效应量计算函数
calculate_cliff_delta <- function(data, group_var, value_var) {
  groups <- unique(data[[group_var]])
  # 修正组别提取逻辑错误（使用 == 替代 <-）
  group1 <- data[data[[group_var]] == groups[1], value_var]
  group2 <- data[data[[group_var]] == groups[2], value_var]
  
  # 计算dominance矩阵
  dominance <- outer(group2, group1, FUN = function(x,y) sign(x - y))
  delta <- mean(dominance)
  
  return(list(
    estimate = delta,
    magnitude = ifelse(abs(delta) < 0.147, "negligible",
                       ifelse(abs(delta) < 0.33, "small",
                              ifelse(abs(delta) < 0.474, "medium", "large")))
  ))
}

# 使用示例
cd_delta <- calculate_cliff_delta(cd19_data, "Group", "Expression")

# 可视化 --------------------------------------------------------------------
library(ggplot2)
library(ggpubr)
pdf("对照组FOXO1.pdf",height = 6,width = 4)
ggplot(cd19_data, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.8) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8, color = "gray20") +
  scale_fill_manual(values = c("AC" = "#66C2A5", "CC" = "#FC8D62")) +
  labs(
    title = "FOXO1 Expression in AC vs CC Groups",
    subtitle = sprintf("Wilcoxon p = %s\nCliff's δ = %.2f (%s effect)", 
                       format.pval(wilcox_test$p.value, digits = 2),
                       cd_delta$estimate,
                       cd_delta$magnitude),
    x = "Treatment Group",
    y = "Expression Level",  # 修改Y轴标签
    caption = "Data points represent individual samples"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "none",
    axis.text.x = element_text(size = 12, color = "black")
  ) +
  stat_compare_means(
    comparisons = list(c("AC", "CC")),
    method = "wilcox.test",
    label = "p.format",
    label.y = max(cd19_data$Expression)*1.6,  # 调整标签位置
    size = 5
  )
dev.off()









# 提IL6ST原始表达值并筛选CC/AC组 --------------------------------------------------
cd19_data <- prtein_matrix_numeric["IL6ST", , drop = FALSE] %>%  # 使用原始矩阵
  t() %>%
  as.data.frame() %>%
  dplyr::rename(Expression = IL6ST) %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(
    Group = case_when(
      Sample %in% group_mapping$AC ~ "AC",
      Sample %in% group_mapping$CC ~ "CC",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group)) %>%  # 仅保留AC和CC组
  mutate(Group = factor(Group, levels = c("AC", "CC")))  # 固定组别顺序

# 检查数据结构
str(cd19_data)

# 统计学检验 ----------------------------------------------------------------
# Wilcoxon秩和检验
wilcox_test <- wilcox.test(Expression ~ Group, data = cd19_data)

# 修正效应量计算函数
calculate_cliff_delta <- function(data, group_var, value_var) {
  groups <- unique(data[[group_var]])
  # 修正组别提取逻辑错误（使用 == 替代 <-）
  group1 <- data[data[[group_var]] == groups[1], value_var]
  group2 <- data[data[[group_var]] == groups[2], value_var]
  
  # 计算dominance矩阵
  dominance <- outer(group2, group1, FUN = function(x,y) sign(x - y))
  delta <- mean(dominance)
  
  return(list(
    estimate = delta,
    magnitude = ifelse(abs(delta) < 0.147, "negligible",
                       ifelse(abs(delta) < 0.33, "small",
                              ifelse(abs(delta) < 0.474, "medium", "large")))
  ))
}

# 使用示例
cd_delta <- calculate_cliff_delta(cd19_data, "Group", "Expression")

# 可视化 --------------------------------------------------------------------
library(ggplot2)
library(ggpubr)
pdf("对照IL6ST.pdf",height = 6,width = 4)
ggplot(cd19_data, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.8) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8, color = "gray20") +
  scale_fill_manual(values = c("AC" = "#66C2A5", "CC" = "#FC8D62")) +
  labs(
    title = "IL6ST Expression in AC vs CC Groups",
    subtitle = sprintf("Wilcoxon p = %s\nCliff's δ = %.2f (%s effect)", 
                       format.pval(wilcox_test$p.value, digits = 2),
                       cd_delta$estimate,
                       cd_delta$magnitude),
    x = "Treatment Group",
    y = "Expression Level",  # 修改Y轴标签
    caption = "Data points represent individual samples"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "none",
    axis.text.x = element_text(size = 12, color = "black")
  ) +
  stat_compare_means(
    comparisons = list(c("AC", "CC")),
    method = "wilcox.test",
    label = "p.format",
    label.y = max(cd19_data$Expression)*1.6,  # 调整标签位置
    size = 5
  )
dev.off()











# 提取 AA 和 CA 两组样本
ac_samples <- group_mapping$AC
cc_samples <- group_mapping$CC

# 提取 AA 和 CA 的表达矩阵
ac_data <- ALL_sample_final_matrix[, colnames(ALL_sample_final_matrix) %in% ac_samples]
cc_data <- ALL_sample_final_matrix[, colnames(ALL_sample_final_matrix) %in% cc_samples]
# 初始化结果数据框
diff_results <- data.frame(
  Gene = rownames(ALL_sample_final_matrix),
  log2FC = NA,
  p_value = NA
)

# 逐行计算 log2 fold change 和 p 值
for (i in 1:nrow(ALL_sample_final_matrix)) {
  gene_data_ac <- as.numeric(ac_data[i, ])
  gene_data_cc <- as.numeric(cc_data[i, ])
  
  # 计算 log2 fold change
  log2fc <- log2(mean(gene_data_cc, na.rm = TRUE) / mean(gene_data_ac, na.rm = TRUE))
  
  # 计算 p 值（使用 Wilcoxon 检验）
  p_value <- wilcox.test(gene_data_ac, gene_data_cc)$p.value
  
  # 填充结果
  diff_results$log2FC[i] <- log2fc
  diff_results$p_value[i] <- p_value
}

# 调整 p 值（FDR 校正）
diff_results$adj_p_value <- p.adjust(diff_results$p_value, method = "holm")



# 设置显著性阈值
log2fc_threshold <- 1  # log2 fold change 阈值
p_value_threshold <- 0.05  # 调整后 p 值阈值

# 添加显著性标注
diff_results <- diff_results %>%
  mutate(
    Significant = case_when(
      p_value < p_value_threshold & abs(log2FC) > log2fc_threshold ~ "Significant",
      TRUE ~ "Not Significant"
    )
  )
library(ggplot2)

# 去掉 p_value 为 NaN 的行
diff_results <- diff_results[!is.nan(diff_results$p_value), ]

# 去掉基因名中包含小数点的行
diff_results <- diff_results[!grepl("\\.", diff_results$Gene), ]

# 去掉与线粒体或核糖体相关的基因
# 假设基因名中包含 "MT-" 或 "RPL"/"RPS" 表示线粒体或核糖体基因
diff_results <- diff_results[!grepl("^MT-|^RPL|^RPS", diff_results$Gene, ignore.case = TRUE), ]

# 检查结果
head(diff_results)

# 绘制火山图
ggplot(diff_results, aes(x = log2FC, y = -log10(p_value), color = Significant)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_manual(values = c("Significant" = "#E41A1C", "Not Significant" = "gray")) +
  labs(
    title = "Volcano Plot: ac vs cc",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Significance"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "top"
  ) +
  geom_hline(yintercept = -log10(p_value_threshold), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "blue")





library(ggrepel)  # 加载标签防重叠包

enhanced_volcano <- function(de_results, top_n = 20, 
                             lfc_threshold = 1, p_threshold = 0.05) {
  # 数据清洗：移除无效值和极值
  valid_data <- de_results %>%
    filter(
      !is.na(log2FC),
      !is.na(p_value),
      is.finite(log2FC),
      p_value > 0  # 排除p=0的情况
    ) %>%
    mutate(
      p_value = pmax(p_value, .Machine$double.xmin)  # 防止p_value=0
    )
  
  # 检查有效数据量
  if (nrow(valid_data) == 0) {
    stop("No valid data points after filtering")
  }
  
  # 动态坐标范围计算（带保护机制）
  safe_range <- function(x, expansion = 0.1) {
    x <- x[is.finite(x)]
    if(length(x) == 0) return(c(0, 1))
    range_val <- range(x)
    c(
      range_val[1] - expansion * diff(range_val),
      range_val[2] + expansion * diff(range_val)
    )
  }
  
  # 准备标签数据
  label_data <- valid_data %>%
    filter(Significant == "Significant") %>%
    arrange(desc(abs(log2FC))) %>%
    slice_head(n = top_n)
  
  # 创建绘图对象
  p <- ggplot(valid_data, aes(x = log2FC, y = -log10(p_value),
                              color = Significant, size = Significant)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("Significant" = "red3", 
                                  "Not Significant" = "grey60")) +
    scale_size_manual(values = c("Significant" = 2, 
                                 "Not Significant" = 1)) +
    geom_label_repel(
      data = label_data,
      aes(label = Gene),
      size = 3,
      box.padding = 0.5,
      max.overlaps = 20,
      segment.color = "grey50"
    ) +
    geom_hline(yintercept = -log10(p_threshold), 
               linetype = "dashed", color = "blue") +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), 
               linetype = "dashed", color = "blue") +
    labs(title = "Differentially Expressed Genes",
         x = expression(Log[2]~Fold~Change),
         y = expression(-Log[10]~P~value)) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_line(color = "grey90")
    )
  
  # 动态设置坐标范围
  p + coord_cartesian(
    xlim = safe_range(valid_data$log2FC),
    ylim = safe_range(-log10(valid_data$p_value))
  )
}
getwd()
# 使用示例
enhanced_volcano(diff_results)





enhanced_volcano <- function(de_results, 
                             target_genes = c("CASP1","CD19", "MS4A1", "CD38", "CR2", "IL7R",
                                              "CARD11", "CD74", "NFKBIA", "ISG15",
                                              "TGFBI", "CDKN1B", "EIF2B5","JAG1", "BTK", "C4BPB", "C6", "CD74", "HLA-DPB1", "IGHM", "JAK1", "PPP3CB", "PRKCD", 
                                              "LAT", "IGHV4-61", "IGHV4-59", "IGHV3-73", 
                                              "IGLV7-43", "IGKV2-29", "IGKV2-28", "CFD", "DOCK2", "NLRC4","IL33","BCL2L13","CASP12","FOXO3",
                                              "LIMK1", "MFNG", "DOCK10"),
                             lfc_threshold = 1, 
                             p_threshold = 0.05) {
  # 数据清洗：仅保留显著点
  valid_data <- de_results %>%
    filter(
      Significant == "Significant",  # 核心修改：仅保留显著点
      !is.na(log2FC),
      !is.na(p_value),
      is.finite(log2FC),
      p_value > 0
    ) %>%
    mutate(
      p_value = pmax(p_value, .Machine$double.xmin),
      # 标记目标基因
      GeneType = ifelse(Gene %in% target_genes, "Target", "Significant")
    )
  
  # 检查目标基因是否存在
  missing_genes <- setdiff(target_genes, valid_data$Gene)
  if(length(missing_genes) > 0) {
    warning(paste("以下基因未达显著性:", paste(missing_genes, collapse = ", ")))
  }
  
  # 创建绘图对象
  p <- ggplot(valid_data, aes(x = log2FC, y = -log10(p_value),
                              color = GeneType, size = GeneType)) +
    geom_point(alpha = 0.8) +
    scale_color_manual(
      values = c("Target" = "#D55E00",  # 目标基因橙色
                 "Significant" = "#0072B2"),  # 普通显著点蓝色
      guide = guide_legend(title = NULL)
    ) +
    scale_size_manual(
      values = c("Target" = 4, 
                 "Significant" = 2),
      guide = "none"
    ) +
    geom_label_repel(
      data = filter(valid_data, Gene %in% target_genes),
      aes(label = Gene),
      size = 4.5,
      box.padding = 0.8,
      segment.color = "#D55E00",
      segment.size = 0.8,
      force = 10,  # 增加标签定位力度
      min.segment.length = 0.2,
      max.time = 2,
      seed = 123  # 确保标签位置可复现
    ) +
    geom_hline(yintercept = -log10(p_threshold), 
               linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), 
               linetype = "dashed", color = "grey40") +
    labs(title = "Significantly Differentially Expressed Genes",
         x = expression(Log[2]~Fold~Change),
         y = expression(-Log[10]~P~value)) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = c(0.85, 0.9),  # 调整图例位置
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_line(color = "grey92"),
      plot.title = element_text(face = "bold.italic")
    )
  
  # 动态坐标范围计算
  expand_axis <- function(x, ratio = 0.18) {
    x_range <- range(x)
    span <- diff(x_range)
    c(x_range[1] - ratio*span, x_range[2] + ratio*span)
  }
  
  p + coord_cartesian(
    xlim = expand_axis(valid_data$log2FC),
    ylim = expand_axis(-log10(valid_data$p_value))
  )
}
pdf("CC VS AC.pdf", width = 10,height = 10)
enhanced_volcano(diff_results)
dev.off()
setwd("/home/lungtissue/蛋白质组学整合分析/AC比CC")
write.csv(diff_results,file = "AC比CC差异基因.csv")

setwd("/home/lungtissue/蛋白质组学整合分析")
write.csv(ALL_sample_final_matrix,file = "直接用于分析的矩阵.csv")


