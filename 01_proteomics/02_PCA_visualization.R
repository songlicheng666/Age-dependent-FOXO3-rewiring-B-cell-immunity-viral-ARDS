# =============================================================================
# Script: 02_PCA_visualization.R
# Description: Principal component analysis (PCA) on standardized protein matrix
# with ellipse-annotated scatter plot for five sample groups
# Author: Licheng Song
# Institution: Chinese PLA General Hospital (Eighth Medical Center)
# Study: Age-dependent B-cell immune dysregulation in virus-associated ARDS
# Year: 2025
# Input: ALL_sample.final_matrix.csv (FOT-normalized protein abundance matrix)
# Output: PCA_Plot.pdf
# =============================================================================

setwd("C:/Users/songl/OneDrive/桌面/多组学分析")# 定义分组样本列表（用户已提供）
# 加载必要包
library(ggplot2)
library(ggrepel)

# 合并为分组映射列表
group_mapping <- list(
  AA = AA_samples,
  AC = AC_samples,
  CA1 = CA1_samples,
  CA2 = CA2_samples,
  CC = CC_samples
)
# 准备分组信息 ---------------------------------------------------------------
sample_groups <- data.frame(
  Sample = colnames(prtein_matrix_no_outliers),
  Group = case_when(
    colnames(prtein_matrix_no_outliers) %in% group_mapping$AA ~ "AA",
    colnames(prtein_matrix_no_outliers) %in% group_mapping$AC ~ "AC",
    colnames(prtein_matrix_no_outliers) %in% group_mapping$CA1 ~ "CA1",
    colnames(prtein_matrix_no_outliers) %in% group_mapping$CA2 ~ "CA2",
    colnames(prtein_matrix_no_outliers) %in% group_mapping$CC ~ "CC",
    TRUE ~ "Unknown"
  )
)
head(standardized_data)
# 检查并处理无效特征 ---------------------------------------------------------
# 检测包含NA/NaN的基因
na_genes <- apply(standardized_data, 1, function(x) any(is.na(x)))
# 检测零方差基因
zero_var_genes <- apply(standardized_data, 1, var, na.rm = TRUE) == 0

# 生成过滤后的数据矩阵
filtered_data <- standardized_data[!na_genes & !zero_var_genes, ]

# 再次检查数据维度
cat("原始基因数:", nrow(standardized_data), "\n过滤后基因数:", nrow(filtered_data))

# 重新执行PCA分析 ----------------------------------------------------------
pca_result <- prcomp(t(filtered_data), scale. = FALSE)  # 注意转置矩阵

# 提取主成分数据
pca_data <- data.frame(
  Sample = rownames(pca_result$x),
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2]
) %>% 
  left_join(sample_groups, by = "Sample") %>%
  mutate(Group = factor(Group, levels = c("AA", "AC", "CA1", "CA2", "CC")))
table(pca_data$Group)
# 生成方差解释比例（处理可能出现的负值）
variance <- pca_result$sdev^2
variance_percent <- round(variance[1:2]/sum(variance)*100, 1)

# 可视化优化 --------------------------------------------------------------
library(ggforce)
ggplot(pca_data, aes(PC1, PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.7) +  # 调整点透明度增强重叠区可视化
  geom_mark_ellipse(
    aes(fill = Group, label = Group),
    expand = unit(2, "mm"),  # 扩大椭圆区域
    label.buffer = unit(10, "mm"),  # 增加标签间距
    con.type = "none"  # 移除连接线
  ) +
  scale_color_manual(
    name = "Experimental Group",
    values = c(AA="#F8766D", AC="#00BA38", 
               CA1="#619CFF", CA2="#6A3D9A",  # 修改CA2颜色增强区分度
               CC="#FF69B4")) +
  scale_fill_manual(
    name = "Experimental Group",  # 添加填充图例
    values = c(AA="#F8766D", AC="#00BA38", 
               CA1="#619CFF", CA2="#6A3D9A", 
               CC="#FF69B4")) +
  labs(
    x = paste0("PC1 (", variance_percent[1], "%)"),
    y = paste0("PC2 (", variance_percent[2], "%)"),
    caption = "Ellipses represent 95% confidence regions"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.grid.major = element_line(color = "grey90"),
    plot.caption = element_text(color = "gray40", size = 10)
  ) +
  guides(
    color = guide_legend(nrow = 1),  # 水平排列图例
    fill = "none"  # 隐藏填充图例（已通过颜色图例统一）
  )

# 保存高分辨率图像
ggsave("PCA_Plot.pdf", width = 6.5, height = 7, dpi = 300)
save(diff_results,filtered_data,group_mapping,pca_data,pca_result,prtein_matrix_no_outliers,prtein_matrix_numeric,sample_groups,standardized_data,file = "PCA绘图数据.Rdata")
