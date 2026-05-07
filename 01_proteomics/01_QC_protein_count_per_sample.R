# =============================================================================
# Script: 01_QC_protein_count_per_sample.R
# Description: Per-sample protein detection count, outlier flagging (3-sigma rule),
# and group-level summary statistics for DIA-MS plasma proteomics
# Author: Licheng Song
# Institution: Chinese PLA General Hospital (Eighth Medical Center)
# Study: Age-dependent B-cell immune dysregulation in virus-associated ARDS
# Year: 2025
# Input: ALL_sample.raw_matrix.csv (raw protein intensity matrix)
# Output: protein_count_summary.csv, protein_count_boxplot.pdf
# =============================================================================

# 加载必要包
library(dplyr)
library(tidyr)
setwd("C:/Users/songl/OneDrive/桌面/多组学分析")
# 数据读取与预处理 ---------------------------------------------------------------
raw_data <- read.csv("C:/Users/songl/OneDrive/桌面/多组学分析/ALL_sample.raw_matrix.csv",
                     stringsAsFactors = FALSE, 
                     check.names = FALSE)  # 保持原始列名格式

# 提取样本列（使用正则表达式精确匹配S/L开头+纯数字的列）
sample_cols <- grep("^(S|L)\\d+$", colnames(raw_data), value = TRUE)
prtein_matrix <- raw_data[, sample_cols] 

# 转换为数值矩阵（带异常值检测）
prtein_matrix <- apply(prtein_matrix, 2, function(x) {
  x_numeric <- as.numeric(x)
  if(any(is.na(x_numeric))) {
    invalid_count <- sum(is.na(x_numeric))
    warning(paste("发现", invalid_count, "个非数值数据点，已转换为NA"))
  }
  x_numeric
})
rownames(prtein_matrix) <- raw_data$Gene_symbol

# 分组定义（精确到具体样本名）----------------------------------------------------
group_definition <- list(
  AA = c("S28", "S30", "S31", "S32", "S33", "S34", "S35", "S36", "S37", "S38", "S39",
         "S40", "S41", "S42", "S43", "S44", "S45", "S46", "S47"),
  AC = c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9",
         "S10", "S11", "S12", "S13", "S14"),
  CA1 = c("S47", "S48", "S49", "S50", "S51", "S52", "S53", "S54", "S55", "S56", "S57"),
  CA2 = c("S58", "S59", "S60", "S61", "S62", "S63", "S64", "S65", "S66"),
  CC = c("S15", "S16", "S17", "S18", "S19", "S20", "S21", "S22", 
         "S23", "S24", "S25", "S26", "S27")
)

# 创建分组因子并过滤未分组样本 ----------------------------------------------------
sample_group <- factor(
  case_when(
    colnames(prtein_matrix) %in% group_definition$AA ~ "AA",
    colnames(prtein_matrix) %in% group_definition$AC ~ "AC",
    colnames(prtein_matrix) %in% group_definition$CA1 ~ "CA1",
    colnames(prtein_matrix) %in% group_definition$CA2 ~ "CA2",
    colnames(prtein_matrix) %in% group_definition$CC ~ "CC",
    TRUE ~ NA_character_  # 明确标记未分组
  ),
  levels = c("AA", "AC", "CA1", "CA2", "CC")
)

valid_samples <- !is.na(sample_group)
if(sum(!valid_samples) > 0){
  removed_samples <- colnames(prtein_matrix)[!valid_samples]
  message("已移除未分组样本：", paste(removed_samples, collapse = ", "))
  prtein_matrix <- prtein_matrix[, valid_samples]
  sample_group <- sample_group[valid_samples]
}

# 计算检测到的蛋白数量 -----------------------------------------------------------
protein_counts <- apply(prtein_matrix, 2, function(x) {
  sum(x > 0, na.rm = TRUE)  # 排除NA和0值
})

# 创建分析数据集 -----------------------------------------------------------------
df <- data.frame(
  Sample = colnames(prtein_matrix),
  Group = sample_group,
  ProteinCount = protein_counts,
  stringsAsFactors = FALSE
) %>%
  # 过滤可能的异常值（按3σ原则）
  mutate(
    z_score = scale(ProteinCount),
    ProteinCount = ifelse(abs(z_score) > 3, NA, ProteinCount)
  ) %>%
  drop_na(ProteinCount)

# 分组统计报告 -------------------------------------------------------------------
group_stats <- df %>%
  group_by(Group) %>%
  summarise(
    N = n(),
    Mean = round(mean(ProteinCount), 1),
    SD = round(sd(ProteinCount), 1),
    Median = median(ProteinCount),
    Min = min(ProteinCount),
    Max = max(ProteinCount),
    .groups = "drop"
  ) %>%
  mutate(
    CI = paste0(round(Mean - 1.96*SD/sqrt(N), 1), " - ", 
                round(Mean + 1.96*SD/sqrt(N), 1))
  )

print(group_stats)

# 统计检验与可视化 --------------------------------------------------------------
# Kruskal-Wallis检验
kruskal_result <- kruskal.test(ProteinCount ~ Group, data = df)
cat("\nKruskal-Wallis检验结果:\n")
print(kruskal_result)

if (!require("FSA")) install.packages("FSA", repos = "https://cloud.r-project.org/")# 事后检验（使用Holm校正）
library(FSA)
dunn_test <- dunnTest(ProteinCount ~ Group, data = df, method = "holm")
cat("\nDunn事后检验结果:\n")
print(dunn_test)

# 可视化
library(ggplot2)
ggplot(df, aes(x = Group, y = ProteinCount, fill = Group)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1.5) +
  labs(title = "各组蛋白质检测数量分布",
       subtitle = paste("Kruskal-Wallis p =", format.pval(kruskal_result$p.value, digits = 3)),
       y = "检测到的蛋白质数量",
       x = "实验分组") +
  theme_bw(base_size = 12) +
  scale_fill_brewer(palette = "Set2")









library(ggpubr)  # 用于添加统计标注
# 提取AC和CC组的数据
cc_ac_data <- df %>% 
  filter(Group %in% c("CC", "AC")) %>%
  mutate(Group = factor(Group, levels = c("AC", "CC")))  # 调整因子顺序

# Mann-Whitney U检验（Wilcoxon秩和检验）
wilcox_test <- wilcox.test(ProteinCount ~ Group, data = cc_ac_data, exact = FALSE)
cat("Wilcoxon秩和检验结果:\n")
print(wilcox_test)
pdf("对照组比较蛋白数量.pdf",height = 4,width = 4)
# 创建带统计标注的小提琴图
ggplot(cc_ac_data, aes(x = Group, y = ProteinCount, fill = Group)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1.5) +
  stat_compare_means(
    comparisons = list(c("AC", "CC")),
    method = "wilcox.test",
    label = "p.signif",
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                       symbols = c("​**​*", "​**​", "*", "ns")),
    tip.length = 0.01,
    vjust = 0.5
  ) +
  scale_fill_manual(values = c("#1f78b4", "#33a02c")) +
  labs(title = "AC vs CC",
       y = "Number of proteins",
       x = "group") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")
dev.off()









# 提取AC和CC组的数据
ca_aa_data <- df %>% 
  filter(Group %in% c("AA", "CA1")) %>%
  mutate(Group = factor(Group, levels = c("AA", "CA1")))  # 调整因子顺序

# Mann-Whitney U检验（Wilcoxon秩和检验）
wilcox_test <- wilcox.test(ProteinCount ~ Group, data = ca_aa_data, exact = FALSE)
cat("Wilcoxon秩和检验结果:\n")
print(wilcox_test)
pdf("ardst1比较蛋白数量.pdf",height = 4,width = 4)
# 创建带统计标注的小提琴图
ggplot(ca_aa_data, aes(x = Group, y = ProteinCount, fill = Group)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1.5) +
  stat_compare_means(
    comparisons = list(c("AA", "CA1")),
    method = "wilcox.test",
    label = "p.signif",
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                       symbols = c("​**​*", "​**​", "*", "ns")),
    tip.length = 0.01,
    vjust = 0.5
  ) +
  scale_fill_manual(values = c("#1f78b4", "#33a02c")) +
  labs(title = "AA vs CA1",
       y = "Number of proteins",
       x = "group") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")
dev.off()







# 提取AC和CC组的数据
ca_aa_data <- df %>% 
  filter(Group %in% c("AA", "CA1")) %>%
  mutate(Group = factor(Group, levels = c("AA", "CA1")))  # 调整因子顺序

# Mann-Whitney U检验（Wilcoxon秩和检验）
wilcox_test <- wilcox.test(ProteinCount ~ Group, data = ca_aa_data, exact = FALSE)
cat("Wilcoxon秩和检验结果:\n")
print(wilcox_test)
pdf("ardst1比较蛋白数量.pdf",height = 4,width = 4)
# 创建带统计标注的小提琴图
ggplot(ca_aa_data, aes(x = Group, y = ProteinCount, fill = Group)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1.5) +
  stat_compare_means(
    comparisons = list(c("AA", "CA1")),
    method = "wilcox.test",
    label = "p.signif",
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                       symbols = c("​**​*", "​**​", "*", "ns")),
    tip.length = 0.01,
    vjust = 0.5
  ) +
  scale_fill_manual(values = c("#1f78b4", "#33a02c")) +
  labs(title = "AA vs CA1",
       y = "Number of proteins",
       x = "group") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")
dev.off()











# 调整分组顺序 ----------------------------------------------------------------
df <- df %>%
  mutate(Group = factor(Group, 
                        levels = c("AC", "CC", "AA", "CA1", "CA2"),
                        labels = c("AC", "CC", "AA", "CA1", "CA2")))

# 定义需要比较的组对（根据您的要求）-------------------------------------------
comparison_list <- list(
  c("AC", "CC"),
  c("AA", "CA1"),
  c("AC", "AA"),
  c("CC", "CA1"),
  c("CA1", "CA2")
)

# 执行两两比较（带Holm校正）-----------------------------------------------------
stat_results <- ggpubr::compare_means(
  ProteinCount ~ Group, 
  data = df,
  method = "wilcox.test",
  p.adjust.method = "holm",
  comparisons = comparison_list
)
head(df)
print(unique(stat_results$group1))
print(unique(stat_results$group2))
print(levels(df$Group))
# 生成可视化图表 ---------------------------------------------------------------
# 预先计算显著性标记的垂直位置
y_positions <- max(df$ProteinCount) * (1 + seq(0.05, by = 0.15, length.out = length(comparison_list)))
pdf("所有组的蛋白数量比较.pdf",height = 3,width = 6)
# 生成可视化图表
ggplot(df, aes(x = Group, y = ProteinCount)) +
  geom_violin(aes(fill = Group), alpha = 0.7, trim = FALSE) +
  geom_boxplot(aes(fill = Group), width = 0.2, outlier.shape = NA, show.legend = FALSE) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 2, color = "gray30") +
    scale_fill_manual(
    values = c("#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", "#59A14F")  # 手动指定颜色
  ) +
  labs(
    title = "Number of protein",
    subtitle = "Wilcoxon rank-sum test (Holm correction)",
    y = "Protein detected",
    x = ""
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "gray40")
  ) +
  coord_cartesian(ylim = c(0, max(df$ProteinCount)*1.3))  # 扩展Y轴显示范围
dev.off()



protein_ranks <- data.frame(
  Group = rep(c("AA", "AC", "CA1", "CA2", "CC"), each=2500),
  Rank = rep(1:2500, 5),
  Intensity = c(rnorm(2500, mean=3), rnorm(2500, mean=2.8), 
                rnorm(2500, mean=2.5), rnorm(2500, mean=2.3), 
                rnorm(2500, mean=2)) 
)

top_proteins <- list(
  AA = c("APOA1","HPX","GC","C3","SERP1O3","VTN","CP","AGT","A1BG","APOC1"),
  AC = c("ProteinA","ProteinB","ProteinC","ProteinD","ProteinE","ProteinF","ProteinG","ProteinH","ProteinI","ProteinJ"),
  CA1 = c("Biomarker1","Biomarker2","Biomarker3","Biomarker4","Biomarker5","Biomarker6","Biomarker7","Biomarker8","Biomarker9","Biomarker10"),
  CA2 = c("MarkerA","MarkerB","MarkerC","MarkerD","MarkerE","MarkerF","MarkerG","MarkerH","MarkerI","MarkerJ"),
  CC = c("GeneX","GeneY","GeneZ","GeneA1","GeneB1","GeneC1","GeneD1","GeneE1","GeneF1","GeneG1")
)

# 可视化代码
library(ggplot2)
library(ggrepel)

ggplot(protein_ranks, aes(x=Rank, y=Intensity, color=Group)) +
  geom_smooth(se=FALSE, method = "loess", span=0.3) +
  geom_label_repel(
    data = do.call(rbind, lapply(names(top_proteins), function(g){
      data.frame(Group=g, 
                 Protein=top_proteins[[g]],
                 x=seq(50, by=250, length.out=10),
                 y=max(protein_ranks$Intensity)*0.95 - 0.2*as.numeric(factor(g)))
    })),
    aes(x=x, y=y, label=Protein, fill=Group),
    color="white",
    box.padding = 0.5,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")) +
  labs(x = "Rank of quantified proteins", 
       y = "Quantitative protein intensity (log10)",
       title = "Protein Abundance Distribution Across Five Groups") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid.major = element_line(color="gray90"),
        axis.text = element_text(size=10))















# 生成蛋白质丰度分布数据 --------------------------------------------------
# 计算每个蛋白质在各组的中位强度（log10转换）
protein_ranks <- raw_data %>%
  select(Gene_symbol, all_of(colnames(prtein_matrix))) %>%  # 选择所有样本列
  pivot_longer(cols = -Gene_symbol, names_to = "Sample", values_to = "Intensity") %>%
  mutate(
    Intensity = as.numeric(Intensity),
    log_Intensity = log10(Intensity + 1e-3)  # 添加极小值避免log(0)
  ) %>%
  left_join(data.frame(Sample = colnames(prtein_matrix), Group = sample_group), 
            by = "Sample") %>%
  group_by(Group, Gene_symbol) %>%
  summarise(
    Median_Intensity = median(log_Intensity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Group) %>%
  mutate(
    Rank = dense_rank(desc(Median_Intensity))  # 按中位强度降序排名
  ) %>%
  ungroup()

# 获取每个组前10高丰度蛋白
top_proteins <- protein_ranks %>%
  filter(Rank <= 10) %>%
  arrange(Group, Rank) %>%
  split(.$Group) %>%
  lapply(function(x) head(x$Gene_symbol, 10))

# 可视化代码 --------------------------------------------------------------------
library(ggplot2)
library(ggrepel)
head(protein_ranks)
# 检查Median_Intensity分布
summary(protein_ranks$Median_Intensity)
# 应看到最小值 >= log10(0 + 1e-3) = -3

# 检查各组数据量
table(protein_ranks$Group)
# 确保每组至少有几十个观测值

# 检查前10名蛋白的强度值
protein_ranks %>% 
  filter(Rank <= 10) %>%
  group_by(Group) %>%
  summarise(Min = min(Median_Intensity), 
            Max = max(Median_Intensity))
library(ggplot2)
library(ggrepel)
# 数据预处理优化 --------------------------------------------------------------
# 过滤低质量数据（保留至少在一个组中强度>10的蛋白）
valid_proteins <- protein_ranks %>%
  group_by(Gene_symbol) %>%
  filter(max(Median_Intensity) > 1) %>%  # log10(10) = 1
  ungroup()

# 计算分位数截断点（去除极端低值）
intensity_q <- quantile(valid_proteins$Median_Intensity, probs = 0.01)
valid_proteins <- valid_proteins %>% 
  filter(Median_Intensity > intensity_q)

# 可视化代码优化 --------------------------------------------------------------
library(ggplot2)
library(ggrepel)
library(mgcv)  # 需要加载GAM支持

# 强制设置分组顺序
df$Group <- factor(df$Group, levels = c("AC", "CC", "AA", "CA1", "CA2"))

# 计算分面布局参数
n_groups <- length(levels(df$Group))
ncol <- 5  # 每行显示3组
nrow <- ceiling(n_groups/ncol)
pdf("protein_distribution.pdf", 
    width = 15, 
    height = 5)
# 优化后的可视化代码
ggplot(protein_ranks, aes(x = Rank, y = Median_Intensity, color = Group)) +
  geom_smooth(
    aes(group = Group),
    method = "gam",
    formula = y ~ s(x, bs = "cs", k = 40),
    se = FALSE,
    linewidth = 1.2
  ) +
  geom_text_repel(
    data = protein_ranks %>% 
      filter(Rank <= 10) %>%
      group_by(Group) %>%
      mutate(
        y_pos = Median_Intensity + seq(-0.2, 0.2, length.out = n())
      ),
    aes(label = Gene_symbol, y = y_pos),
    size = 3.5,
    box.padding = 0.15,
    max.overlaps = 20,
    min.segment.length = 0.3,
    direction = "both"
  ) +
  scale_color_manual(
    values = c(
      AC = "#377EB8",   # 蓝色
      CC = "#4DAF4A",   # 绿色
      AA = "#E41A1C",   # 红色
      CA1 = "#984EA3",  # 紫色
      CA2 = "#FF7F00"   # 橙色
    )
  ) +
  labs(
    title = "Protein Abundance Distribution",
    x = "Rank of quantified proteins",
    y = expression(log[10]("Quantitative protein intensity")),
    caption = "Top 10 abundant proteins labeled"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(color = "grey92"),
    strip.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  facet_wrap(
    ~ Group, 
    nrow = nrow,
    labeller = labeller(Group = function(x) paste0(x, " (n = 107)"))
  ) +
  coord_cartesian(ylim = c(-1, 4)) +
  scale_x_continuous(
    breaks = seq(0, 2500, 500),
    labels = function(x) format(x, big.mark = ",", scientific = FALSE)
  )
dev.off()
