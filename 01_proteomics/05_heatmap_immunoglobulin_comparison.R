# =============================================================================
# Script: 05_heatmap_immunoglobulin_comparison.R
# Description: Heatmap of immunoglobulin family proteins with directional annotation
# (AC_UP / CC_UP) between adult and pediatric controls
# Author: Licheng Song
# Institution: Chinese PLA General Hospital (Eighth Medical Center)
# Study: Age-dependent B-cell immune dysregulation in virus-associated ARDS
# Year: 2025
# Input: prtein_matrix_no_outliers (from 03_differential_protein_analysis.R), diff_results
# Output: 免疫球蛋白AC和CC两组比较热图.pdf
# =============================================================================

# 定义免疫球蛋白基因标识符（根据实际基因名调整正则表达式）
ig_genes <- grep("^IG[HKL]|^IGJ|^IGV", rownames(prtein_matrix_no_outliers), value = TRUE, ignore.case = TRUE)

# 提取CA1和CC组的差异结果（假设diff_results已包含两组比较）
AC_vs_CC <- diff_results %>%
  filter(Gene %in% ig_genes)

# 标记显著上调基因
AC_vs_CC <- AC_vs_CC %>%
  mutate(
    Regulation = case_when(
      log2FC > 0 & p_value <= 0.05 ~ "CC_UP",
      log2FC < 0 & p_value <= 0.05 ~ "AC_UP",
      TRUE ~ "Non_sig"
    )
  )

# 生成热图数据 -----------------------------------------------------------
heatmap_data <- prtein_matrix_no_outliers[AC_vs_CC$Gene, c(group_mapping$AC, group_mapping$CC)]
heatmap_zscore <- t(scale(t(heatmap_data)))  # 行标准化

# 构建行注释（基因调控方向）
row_anno <- data.frame(
  Regulation = factor(AC_vs_CC$Regulation, 
                      levels = c("AC_UP", "CC_UP", "Non_sig")),
  row.names = AC_vs_CC$Gene
)

# 自定义颜色方案
anno_colors <- list(
  Regulation = c("AC_UP" = "#E41A1C",  # 红色
                 "CC_UP" = "#4DAF4A",    # 绿色
                 "Non_sig" = "grey90")
)


# 检查标准化后的数据是否存在无效值
invalid_values <- sum(is.na(heatmap_zscore) | is.infinite(heatmap_zscore))
if (invalid_values > 0) {
  warning(paste("Found", invalid_values, "NA/NaN/Inf values. Removing affected genes."))
  
  # 删除包含无效值的行（基因）
  valid_rows <- apply(heatmap_zscore, 1, function(row) {
    !any(is.na(row) | is.infinite(row))
  })
  heatmap_zscore <- heatmap_zscore[valid_rows, ]
  
  # 同步更新注释信息
  row_anno <- row_anno[rownames(row_anno) %in% rownames(heatmap_zscore), , drop = FALSE]
}
# 检查并移除标准差为零的基因
zero_var_genes <- apply(heatmap_zscore, 1, function(row) sd(row, na.rm = TRUE) == 0)
if (sum(zero_var_genes) > 0) {
  warning(paste("Removing", sum(zero_var_genes), "genes with zero variance."))
  heatmap_zscore <- heatmap_zscore[!zero_var_genes, ]
  row_anno <- row_anno[rownames(row_anno) %in% rownames(heatmap_zscore), , drop = FALSE]
}
# 绘制热图 --------------------------------------------------------------
library(pheatmap)
setwd("C:/Users/songl/OneDrive/桌面/多组学分析/对照组比较")
# 确保数据有效后绘图
pdf("免疫球蛋白AC和CC两组比较热图.pdf",width = 4.5,height = 13)
  pheatmap(
    mat = heatmap_zscore,
    annotation_col = data.frame(
      Group = rep(c("AC", "CC"), 
                  c(length(group_mapping$AC), length(group_mapping$CC))),
      row.names = colnames(heatmap_zscore)
    ),
    annotation_row = row_anno,
    annotation_colors = anno_colors,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    main = "Immunoglobulin Genes: AC vs CC Expression",
    fontsize_row = 8,
    gaps_col = length(group_mapping$AC),
    clustering_distance_rows = "euclidean",
    clustering_method = "ward.D2",
    silent = FALSE  # 显示计算过程
  )
dev.off()
save(prtein_matrix_no_outliers,diff_results,file = "用于计算免疫球蛋白AC和CC两组比较热图.Rdata")
