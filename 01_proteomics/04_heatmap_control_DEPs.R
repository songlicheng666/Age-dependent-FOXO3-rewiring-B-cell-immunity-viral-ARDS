# =============================================================================
# Script: 04_heatmap_control_DEPs.R
# Description: Heatmap of significant differentially expressed proteins (DEPs)
# between pediatric controls (CC) and adult controls (AC)
# Author: Licheng Song
# Institution: Chinese PLA General Hospital (Eighth Medical Center)
# Study: Age-dependent B-cell immune dysregulation in virus-associated ARDS
# Year: 2025
# Input: prtein_matrix_numeric (from 03_differential_protein_analysis.R), diff_results
# Output: 差异基因热图.pdf
# =============================================================================


head(diff_results)
# 筛选 p_value ≤ 0.05 的基因
diff_genes <- diff_results %>%
  dplyr::filter(p_value <= 0.05) %>%
  dplyr::pull(Gene)  # 提取基因名列表
# 获取CC和AC组的样本名
CC_samples <- group_mapping$CC
AC_samples <- group_mapping$AC

# 筛选表达矩阵中的对应样本（并转换为数值矩阵）
prtein_matrix_CC_AC <- prtein_matrix_numeric[diff_genes, c(CC_samples, AC_samples), drop = FALSE]
prtein_matrix_CC_AC <- as.matrix(prtein_matrix_CC_AC)  # 确保为矩阵格式
# 标准化数据（行方向：每个基因在样本间的Z-score）
heatmap_data <- t(scale(t(prtein_matrix_CC_AC)))
# 创建分组注释数据框
annotation_col <- data.frame(
  Group = factor(c(rep("CC", length(CC_samples)), rep("AC", length(AC_samples))))
)
rownames(annotation_col) <- c(CC_samples, AC_samples)  # 样本名为行名

# 自定义分组颜色
group_colors <- list(Group = c(CC = "#E69F00", AC = "#56B4E9"))  # 橙色-青色
library(pheatmap)
pdf("差异基因热图.pdf",width = 7,height = 6)
pheatmap(
  mat = heatmap_data,
  annotation_col = annotation_col,
  annotation_colors = group_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),  # 颜色梯度
  show_colnames = FALSE,     # 不显示样本名
  show_rownames = F,      # 显示基因名（若过多可设为FALSE）
  cluster_rows = TRUE,       # 对行（基因）聚类
  cluster_cols = TRUE,       # 对列（样本）聚类
  clustering_method = "ward.D2",  # 聚类算法
  main = "Differential Genes (CC vs AC, p ≤ 0.05)"
)
dev.off()
