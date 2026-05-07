# =============================================================================
# Script: 07_ARDS_protein_expression_boxplot.R
# Description: Protein-level expression visualization for candidate prognostic markers
# (FOXO3, TLR7, IL6ST, etc.) across all five groups with statistical annotations
# Author: Licheng Song
# Institution: Chinese PLA General Hospital (Eighth Medical Center)
# Study: Age-dependent B-cell immune dysregulation in virus-associated ARDS
# Year: 2025
# Input: ALL_sample.final_matrix.csv
# Output: candidate_protein_boxplots.pdf
# =============================================================================

library(dplyr)
library(tidyr)
library(tibble)
setwd("D:/博士毕业/毕业文章/F6蛋白组学ARDS")# 定义分组样本列表（用户已提供）
# ==============================================================================
# 1. 读取数据 (保持原始列名)
# ==============================================================================
# 假设第一列是蛋白/基因ID，将其设为行名 (row.names = 1)
# check.names = FALSE 非常重要，防止 R 自动修改样本名 (如将 S-1 改为 S.1)
raw_data <- read.csv("ALL_sample.final_matrix.csv", 
                     check.names = FALSE, 
                     row.names = 1,
                     stringsAsFactors = FALSE)

# 获取数据中实际存在的所有样本名
available_samples <- colnames(raw_data)

cat("数据加载完成。\n")
cat("数据矩阵包含行数 (蛋白数):", nrow(raw_data), "\n")
cat("数据矩阵包含列数 (样本数):", ncol(raw_data), "\n\n")

# ==============================================================================
# 2. 定义原始分组 (用户提供的完整列表)
# ==============================================================================

# AC组 (S1 - S14)
AC_samples_def <- c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12", "S13", "S14")

# CC组 (S15 - S27)
CC_samples_def <- c("S15", "S16", "S17", "S18", "S19", "S20", "S21", "S22", "S23", "S24", "S25", "S26", "S27")

# AA组 (S28-S47, S67-S76, S87-S108, BF1-BF5)
AA_samples_def <- c(
  "S28", "S29", "S30", "S31", "S32", "S33", "S34", "S35", "S36", "S37", 
  "S38", "S39", "S40", "S41", "S42", "S43", "S44", "S45", "S46", "S47",
  "S67", "S68", "S69", "S70", "S71", "S72", "S73", "S74", "S75", "S76",
  "S87", "S88", "S89", "S90", "S91", "S92", "S93", "S94", "S95", "S96", 
  "S97", "S98", "S99", "S100", "S101", "S102", "S103", "S104", "S105", "S106", "S107", "S108",
  "BF1", "BF2", "BF3", "BF4", "BF5"
)

# CA1组 (S48-S57, S77-S86)
CA1_samples_def <- c(
  "S48", "S49", "S50", "S51", "S52", "S53", "S54", "S55", "S56", "S57",
  "S77", "S78", "S79", "S80", "S81", "S82", "S83", "S84", "S85", "S86"
)

# CA2组 (S58 - S66)
CA2_samples_def <- c("S58", "S59", "S60", "S61", "S62", "S63", "S64", "S65", "S66")

# ==============================================================================
# 3. 自动校对：只保留数据中存在的样本 (取交集)
# ==============================================================================

# 定义一个辅助函数来过滤和汇报
filter_samples <- function(defined_list, available_list, group_name) {
  valid_samples <- intersect(defined_list, available_list)
  cat(sprintf("[%s] 定义数量: %d -> 实际匹配数量: %d\n", 
              group_name, length(defined_list), length(valid_samples)))
  
  # 如果有丢失的样本，打印出来（可选）
  missing <- setdiff(defined_list, available_list)
  if(length(missing) > 0) {
    cat(sprintf("   警告: %d 个样本未在数据中找到: %s ...\n", 
                length(missing), paste(head(missing, 3), collapse=", ")))
  }
  return(valid_samples)
}

cat("正在校对样本...\n")
AC_samples <- filter_samples(AC_samples_def, available_samples, "AC")
CC_samples <- filter_samples(CC_samples_def, available_samples, "CC")
AA_samples <- filter_samples(AA_samples_def, available_samples, "AA")
CA1_samples <- filter_samples(CA1_samples_def, available_samples, "CA1")
CA2_samples <- filter_samples(CA2_samples_def, available_samples, "CA2")

# 合并 CA 总组 (如果需要)
CA_samples <- c(CA1_samples, CA2_samples)

# ==============================================================================
# 4. 构建最终的分组映射列表 (供后续分析使用)
# ==============================================================================
group_mapping <- list(
  AC = AC_samples,
  CC = CC_samples,
  AA = AA_samples,
  CA1 = CA1_samples,
  CA2 = CA2_samples,
  CA_ALL = CA_samples # 包含所有 CA
)

# 转换为数值矩阵 (用于后续计算)
# 注意：这里我们只保留所有分组中涉及到的样本，去掉可能存在的无关列
all_target_samples <- c(AC_samples, CC_samples, AA_samples, CA1_samples, CA2_samples)
prtein_matrix_numeric <- data.matrix(raw_data[, all_target_samples])

library(tidyr)
library(dplyr)
calculate_cliff_delta <- function(data, group_var, value_var) {
  groups <- unique(data[[group_var]])
  # 确保有且只有两个组
  if (length(groups) != 2) {
    stop("Cliff's delta 目前只支持两组比较")
  }
  
  # 提取两个组的数值向量
  group1 <- data[[value_var]][data[[group_var]] == groups[1]]
  group2 <- data[[value_var]][data[[group_var]] == groups[2]]
  
  # 计算 dominance 矩阵
  dominance <- outer(group2, group1, FUN = function(x, y) sign(x - y))
  delta <- mean(dominance, na.rm = TRUE)
  
  magnitude <- ifelse(
    abs(delta) < 0.147, "negligible",
    ifelse(
      abs(delta) < 0.33, "small",
      ifelse(
        abs(delta) < 0.474, "medium", "large"
      )
    )
  )
  
  list(
    estimate  = delta,
    magnitude = magnitude
  )
}

get_gene_ac_cc_from_original <- function(gene_symbol) {
  # gene_symbol：例如 "CD19", "BTK" 等
  # 使用你原来的做法：从 prtein_matrix_numeric 抽一行，然后标注 AC/CC
  
  if (!gene_symbol %in% rownames(prtein_matrix_numeric)) {
    stop(paste("基因不在矩阵中:", gene_symbol))
  }
  
  df <- prtein_matrix_numeric[gene_symbol, , drop = FALSE] %>%
    t() %>%
    as.data.frame() %>%
    # 把这一列重命名为 Expression（不再用基因名做列名，避免后面写很多 rename）
    setNames("Expression") %>%
    tibble::rownames_to_column("Sample") %>%
    dplyr::mutate(
      Group = case_when(
        Sample %in% group_mapping$CA1 ~ "CA1",
        Sample %in% group_mapping$AA ~ "AA",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(Group)) %>%
    # 关键：这里把顺序改成 CC 在左，AC 在右
    dplyr::mutate(Group = factor(Group, levels = c("CA1", "AA"))) 
  
  df$Gene <- gene_symbol
  return(df)
}
ALL_sample_final_matrix <- as.data.frame(raw_data)
rownames(ALL_sample_final_matrix) <- ALL_sample_final_matrix[,1]
ALL_sample_final_matrix <- ALL_sample_final_matrix[,-1]
#ALL_sample_final_matrix <- as.data.frame(ALL_sample_final_matrix)

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
        Sample %in% group_mapping$CA1 ~ "CA1",
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

# 2) 按功能模块指定基因顺序（你可以按自己想法改）
gene_groups <- list(
  "B-cell activation / markers" = c("CD19", "CD27", "CD38", "CR2"),
  "Development / signaling"     = c("BTK", "JAK1", "IL7R", "IL6ST", "CD74", "IRF4", "FOXO1", "FOXO3"),
  "Inflammation / cell death"   = c("CASP1", "TNFRSF10B", "NLRC4", "IRF5", "OAS1"),
  "Matrix / tissue"             = c("TGFBI")
)

gene_order <- unique(unlist(gene_groups))
# 防止写错名字：只保留矩阵里真实存在的基因
gene_order <- gene_order[gene_order %in% rownames(prtein_matrix_numeric)]

# 2.1 拼接所有基因的长表数据
violin_data_list <- lapply(gene_order, function(g) {
  get_gene_ac_cc_from_original(g)
})
violin_data <- bind_rows(violin_data_list)
# 按基因缩放表达量：使每个基因在自己的 panel 中数值在 0~1 左右
violin_data <- violin_data %>%
  dplyr::group_by(Gene) %>%
  dplyr::mutate(
    Expression_scaled = Expression / (max(Expression, na.rm = TRUE) + 1e-6)
  ) %>%
  dplyr::ungroup()

# 2.2 计算每个基因的 Wilcoxon p 值
get_gene_stats <- function(gene_symbol) {
  df <- get_gene_ac_cc_from_original(gene_symbol)
  
  # Wilcoxon
  w <- wilcox.test(Expression ~ Group, data = df)
  
  # Cliff's delta
  cd <- calculate_cliff_delta(df, group_var = "Group", value_var = "Expression")
  
  tibble(
    Gene         = gene_symbol,
    p_value      = w$p.value,
    p_label      = format.pval(w$p.value, digits = 2, eps = 1e-4),
    cliff_delta  = cd$estimate,
    magnitude    = cd$magnitude
  )
}

stat_df <- bind_rows(lapply(gene_order, get_gene_stats)) %>%
  # 构造你要显示的完整文字
  mutate(
    label = sprintf(
      "Wilcoxon p = %s\nCliff's \u03b4 = %.2f (%s effect)",
      p_label, cliff_delta, magnitude
    )
  )

library(ggplot2)

# 3.1 为 p 值标注准备一个数据表：每个基因一个点，y 在该基因表达的顶部稍微往上
# 计算每个基因的 y 坐标（这里假设你已经改成 Expression_scaled，如果还在用 Expression，就把字段名换回去）
label_pos <- violin_data %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(y_pos = max(Expression_scaled, na.rm = TRUE) * 1.05, .groups = "drop")

# 合并统计结果和 y 位置
stat_df_for_plot <- stat_df %>%
  dplyr::left_join(label_pos, by = "Gene")

# 3.2 颜色（注意 factor 顺序在 get_gene_ac_cc_from_original 里是 AC,CC；
#       但我们画 x 轴要 CC 左 AC 右，因此这里重设一次 factor）
violin_data$Group <- factor(violin_data$Group, levels = c("CA1", "AA"))
group_colors <- c("CA1" = "#FC8D62", "AA" = "#66C2A5")

# 3.3 确保 facet 顺序按 gene_order 来
violin_data$Gene <- factor(violin_data$Gene, levels = gene_order)
stat_df_for_plot$Gene <- factor(stat_df_for_plot$Gene, levels = gene_order)
# 先基于 violin_data 计算每个基因的平均表达量，用于排序
gene_mean_expr <- violin_data %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(mean_expr = mean(Expression, na.rm = TRUE), .groups = "drop")

# 按平均表达量从低到高排序：低表达排在前面
gene_order_by_expr <- gene_mean_expr %>%
  dplyr::arrange(mean_expr) %>%
  dplyr::pull(Gene)

# 用新的顺序重设 factor 水平
violin_data$Gene     <- factor(violin_data$Gene,     levels = gene_order_by_expr)
stat_df_for_plot$Gene <- factor(stat_df_for_plot$Gene, levels = gene_order_by_expr)

# 3.4 画图并导出 PDF
pdf("CA1_AA_selected_genes_violin_combined_with_p_and_cliff.pdf", width = 7, height = 5)

ggplot(violin_data, aes(x = Group, y = Expression_scaled, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.8) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, size = 0.25) +
  geom_jitter(width = 0.08, size = 0.9, alpha = 0.8, color = "gray20") +
  scale_fill_manual(values = group_colors) +
  facet_wrap(~ Gene, nrow = 3) +
  # 在每个 panel 顶部写完整的 p 和 Cliff's δ
  geom_text(
    data = stat_df_for_plot,
    aes(x = 1.5, y = y_pos, label = label),
    inherit.aes = FALSE,
    size = 1,      # 字稍微小一点
    lineheight = 0.9 # 行距紧一点
  ) +
  labs(
    x = NULL,
    y = "Relative protein expression"
  ) +
  theme_classic(base_size = 9.5) +
  theme(
    strip.text = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 7, color = "black"),
    legend.position = "none",
    plot.margin = margin(3, 3, 3, 3)
  )

dev.off()
save(ALL_sample_final_matrix,gene_groups,group_mapping,prtein_matrix_no_outliers,prtein_matrix_numeric,raw_data,standardized_data,AA_samples,AC_samples,CA_samples,CA1_samples,CA2_samples,CC,CC_samples,file = "最全的用于计算蛋白组学的数据.Rdata")
