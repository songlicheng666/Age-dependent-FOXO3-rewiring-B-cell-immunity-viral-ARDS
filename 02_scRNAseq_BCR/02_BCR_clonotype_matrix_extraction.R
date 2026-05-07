# =============================================================================
# Script: 02_BCR_clonotype_matrix_extraction.R
# Description: BCR clonotype abundance matrix construction, isotype proportion per sample,
# clonal homeostasis summary; integration with Seurat cell-state annotations
# Author: Licheng Song
# Institution: Chinese PLA General Hospital (Eighth Medical Center)
# Study: Age-dependent B-cell immune dysregulation in virus-associated ARDS
# Year: 2025
# Input: scBCR_RNA_PB (scRepertoire object), CountsSeuratB_umap_J (Seurat object)
# Output: clonotype_proportion_matrix.csv, isotype_proportion_matrix.csv
# =============================================================================

head(breast.TCGA$data.train)
table(scBCR_RNA_PB$cloneType)
table(scBCR_RNA_PB$Isotype)
table(scBCR_RNA_PB$group1)
getwd()
setwd("/home/lungtissue/蛋白质组学整合分析")
table(CountsSeuratB_umap_J$group)
table(CountsSeuratB_umap_J$group1)
# 合并CA和CAT2假设scBCR_RNA_PB是你的Seurat对象
# 首先复制group1到group2
CountsSeuratB_umap_J$group2 <- CountsSeuratB_umap_J$group

# 使用gsub函数将CAT2替换为CA
CountsSeuratB_umap_J$group2 <- gsub("C_ECMO_off", "C_ECMO_on", CountsSeuratB_umap_J$group2)

# 确认group2中只有AA, AC, CA, CC四组
unique_groups <- unique(CountsSeuratB_umap_J$group2)
print(unique_groups)
table(CountsSeuratB_umap_J$sample)


# 提取样本和克隆类型信息
samples <- scBCR_RNA_PB$sample
clone_types <- scBCR_RNA_PB$cloneType

# 创建一个空的数据框，用于存储结果
result_df <- data.frame(matrix(ncol = length(unique(clone_types)), nrow = length(unique(samples))))
rownames(result_df) <- unique(samples)
colnames(result_df) <- unique(clone_types)

# 计算每个样本中各个克隆类型的比例
for (sample in unique(samples)) {
  # 获取当前样本的克隆类型
  current_clone_types <- clone_types[samples == sample]
  
  # 计算每个克隆类型的数量
  clone_type_counts <- table(current_clone_types)
  
  # 计算比例
  clone_type_proportions <- clone_type_counts / sum(clone_type_counts)
  
  # 将结果存储到数据框中
  result_df[sample, names(clone_type_proportions)] <- clone_type_proportions
}

# 用0填充NA值（如果某个样本没有某个克隆类型）
result_df[is.na(result_df)] <- 0

# 打印结果数据框
print(result_df)
result_cloneType_proportion <- result_df




# 提取样本和克隆类型信息
samples <- CountsSeuratB_umap_J$sample
Isotype <- scBCR_RNA_PB$Isotype

# 创建一个空的数据框，用于存储结果
result_Isotype <- data.frame(matrix(ncol = length(unique(Isotype)), nrow = length(unique(samples))))
rownames(result_Isotype) <- unique(samples)
colnames(result_Isotype) <- unique(Isotype)

# 计算每个样本中各个克隆类型的比例
for (sample in unique(samples)) {
  # 获取当前样本的克隆类型
  current_clone_types <- Isotype[samples == sample]
  
  # 计算每个克隆类型的数量
  clone_type_counts <- table(current_clone_types)
  
  # 计算比例
  clone_type_proportions <- clone_type_counts / sum(clone_type_counts)
  
  # 将结果存储到数据框中
  result_Isotype[sample, names(clone_type_proportions)] <- clone_type_proportions
}

# 用0填充NA值（如果某个样本没有某个克隆类型）
result_Isotype[is.na(result_Isotype)] <- 0

# 打印结果数据框
print(result_Isotype)
result_Isotype_proportion <- result_Isotype



#四组综合的差异基因
Idents(CountsSeuratB_umap_J) <- CountsSeuratB_umap_J$group2

# 计算所有组之间的差异基因
all_differential_genes <- FindAllMarkers(CountsSeuratB_umap_J)

# 查看结果
head(all_differential_genes)
table(CountsSeuratB_umap_J$sample)





# 原始样本名称
original_samples <- c("A_CTRL_34", "A_CTRL_35", "A_CTRL_36", "A_CTRL_37", "A_CTRL_38", "A_CTRL_39", "A_CTRL_40", "A_CTRL_41", "A_CTRL_42", "A_CTRL_43", "A_CTRL_44", "A_CTRL_45",
                      "A_ECMO_on_46", "A_ECMO_on_47", "A_ECMO_on_48", "A_ECMO_on_49", "A_ECMO_on_50", "A_ECMO_on_51", "A_ECMO_on_52", "A_ECMO_on_53", "A_ECMO_on_54",
                      "C_CTRL_1", "C_CTRL_10", "C_CTRL_11", "C_CTRL_12", "C_CTRL_13", "C_CTRL_2", "C_CTRL_3", "C_CTRL_4", "C_CTRL_5", "C_CTRL_6", "C_CTRL_7", "C_CTRL_8", "C_CTRL_9",
                      "C_ECMO_off_24", "C_ECMO_off_25", "C_ECMO_off_26", "C_ECMO_off_27", "C_ECMO_off_28", "C_ECMO_off_29", "C_ECMO_off_30", "C_ECMO_off_31", "C_ECMO_off_32", "C_ECMO_off_33",
                      "C_ECMO_on_14", "C_ECMO_on_15", "C_ECMO_on_16", "C_ECMO_on_17", "C_ECMO_on_18", "C_ECMO_on_19", "C_ECMO_on_20", "C_ECMO_on_21", "C_ECMO_on_22", "C_ECMO_on_23")

# 新的样本名称
new_samples <- c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12",
                 "S28", "S30", "S31", "S32", "S33", "S34", "S35", "S36", "S37", 
                 "S15", "S16", "S17", "S18", "S19", "S20", "S21", "S22", "S23", "S24", "S25", "S26", "S27",
                 "S47", "S48", "S49", "S50", "S51", "S52", "S53", "S54", "S55", "S56",
                 "S57", "S58", "S59", "S60", "S61", "S62", "S63", "S64", "S65", "S66")

# 方法1：使用 plyr::mapvalues
library(plyr)
CountsSeuratB_umap_J$sample2 <- plyr::mapvalues(CountsSeuratB_umap_J$sample,
                                        from = names(sample_mapping),
                                        to = sample_mapping)

# 方法2：使用数据框方式进行映射
sample_df <- data.frame(
  old = names(sample_mapping),
  new = sample_mapping,
  stringsAsFactors = FALSE
)
CountsSeuratB_umap_J$sample2 <- sample_df$new[match(CountsSeuratB_umap_J$sample, sample_df$old)]
table(CountsSeuratB_umap_J$sample2,CountsSeuratB_umap_J$group2)
#提取基因表达矩阵
# 默认使用标准化后的数据 (normalized data)
expr_matrix <- GetAssayData(CountsSeuratB_umap_J, slot = "data")

# 获取细胞的样本信息
cell_sample_info <- CountsSeuratB_umap_J$sample2

# 按样本聚合表达量（计算每个样本中每个基因的平均表达量）
# 使用 aggregate 函数按样本进行聚合
expr_by_sample <- aggregate(t(as.matrix(expr_matrix)), 
                            list(Sample = cell_sample_info), 
                            mean)

# 转置矩阵并整理行列名
final_matrix <- t(expr_by_sample[,-1])
colnames(final_matrix) <- expr_by_sample$Sample

# 现在 final_matrix 就是你想要的矩阵：
# - rownames 是基因名
# - colnames 是 sample2 的名称
# - 值是每个样本中每个基因的平均表达量

# 可以查看矩阵的维度
dim(final_matrix)

# 查看矩阵的前几行和前几列
head(final_matrix[, 1:5])
# 获取不以 "RP" 和 "ATP" 开头的基因的行名
keep_genes <- rownames(final_matrix)[!grepl("^(RP|ATP)", rownames(final_matrix))]

# 使用这些基因名过滤矩阵
filtered_matrix <- final_matrix[keep_genes, ]

# 查看过滤前后的维度
cat("原始矩阵维度：", dim(final_matrix), "\n")
cat("过滤后矩阵维度：", dim(filtered_matrix), "\n")

# 可以查看过滤后矩阵的前几行
head(filtered_matrix[, 1:5])
# 获取不包含"."的基因的行名
keep_genes <- rownames(final_matrix)[!grepl("\\.", rownames(final_matrix))]

# 使用这些基因名过滤矩阵
filtered_matrix <- final_matrix[keep_genes, ]

# 查看过滤前后的维度
cat("原始矩阵维度：", dim(final_matrix), "\n")
cat("过滤后矩阵维度：", dim(filtered_matrix), "\n")

# 可以查看过滤后矩阵的前几行
head(filtered_matrix[, 1:5])
# 带有更多参数控制的导出版本
write.csv(filtered_matrix, 
          file = "sc_plasmacell_genematrix.csv",
          quote = FALSE,  # 不给字符串加引号
          row.names = TRUE  # 保留行名（基因名）
)
head(result_cloneType_proportion)
# 修正重复的样本名称
new_samples_fixed <- new_samples

name_mapping_fixed <- setNames(new_samples_fixed, original_samples)

# 替换行名
new_rownames <- name_mapping_fixed[current_rownames]
rownames(result_cloneType_proportion) <- new_rownames

# 使用 !is.na 筛选列名
result_cloneType_proportion <- result_cloneType_proportion[, !is.na(colnames(result_cloneType_proportion))]

# 查看结果
head(result_cloneType_proportion)
result_cloneType_proportion <- t(result_cloneType_proportion)




rownames(result_Isotype_proportion) <- new_rownames

# 查看结果
head(result_Isotype_proportion)
result_Isotype_proportion <- t(result_Isotype_proportion)


# 带有更多参数控制的导出版本
write.csv(result_Isotype_proportion, 
          file = "result_Isotype_proportion.csv",
          quote = FALSE,  # 不给字符串加引号
          row.names = TRUE  # 保留行名（基因名）
)
write.csv(result_cloneType_proportion, 
          file = "result_cloneType_proportion.csv",
          quote = FALSE,  # 不给字符串加引号
          row.names = TRUE  # 保留行名（基因名）
)

# 加载 readxl 包
library(readxl)

# 读取 Excel 文件
prtein_matrix <- read_excel("ALL_sample.final_matrix.xlsx")

# 查看数据的前几行
head(prtein_matrix)

# 查看数据的基本信息
str(prtein_matrix)
# 将第一列设置为行名并删除该列
rownames(prtein_matrix) <- prtein_matrix[[1]]
# 查看结果
head(prtein_matrix)
# 将 filtered_matrix 转换为数据框
filtered_matrix <- as.data.frame(filtered_matrix)
result_cloneType_proportion <- as.data.frame(result_cloneType_proportion)
result_Isotype_proportion <- as.data.frame(result_Isotype_proportion)

# 获取所有数据框的列名
cols_protein <- colnames(prtein_matrix)
cols_cloneType <- colnames(result_cloneType_proportion)
cols_Isotype <- colnames(result_Isotype_proportion)
cols_filtered <- colnames(filtered_matrix)

# 找出共同的列名
common_cols <- Reduce(intersect, list(cols_protein, 
                                      cols_cloneType, 
                                      cols_Isotype, 
                                      cols_filtered))

# 使用共同的列名从每个数据框中提取相应的列
protein_filtered <- prtein_matrix[, common_cols]
cloneType_filtered <- result_cloneType_proportion[, common_cols]
Isotype_filtered <- result_Isotype_proportion[, common_cols]
filtered_filtered <- filtered_matrix[, common_cols]

# 查看结果
print(paste("共同的列数：", length(common_cols)))
print("共同的列名：")
print(common_cols)









# 首先创建一个按分组的样本顺序列表
AA_samples <- c("S28", "S30", "S31", "S32", "S33", "S34", "S35", "S36", "S37")
AC_samples <- c("S1", "S10", "S11", "S12", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9")
CA_samples <- c("S47", "S48", "S49", "S50", "S51", "S52", "S53", "S54", "S55", "S56", "S57", "S58", "S59", "S60", "S61", "S62", "S63", "S64", "S65", "S66")
CC_samples <- c("S15", "S16", "S17", "S18", "S19", "S20", "S21", "S22", "S23", "S24", "S25", "S26", "S27")

# 合并所有样本，按照 AA, AC, CA, CC 的顺序
ordered_samples <- c(AA_samples, AC_samples, CA_samples, CC_samples)

# 获取四个数据框共同的列
common_cols <- Reduce(intersect, list(colnames(protein_filtered), 
                                      colnames(cloneType_filtered), 
                                      colnames(Isotype_filtered), 
                                      colnames(filtered_filtered)))

# 找出在 common_cols 中的样本列（与 ordered_samples 相交的列）
sample_cols <- intersect(ordered_samples, common_cols)

# 找出非样本列（可能是基因名或其他特征）
other_cols <- setdiff(common_cols, sample_cols)

# 创建最终的列顺序（先放其他列，再放样本列，样本列按分组顺序排列）
final_order <- c(other_cols, sample_cols)

# 重新排序每个数据框的列
protein_filtered <- protein_filtered[, final_order]
cloneType_filtered <- cloneType_filtered[, final_order]
Isotype_filtered <- Isotype_filtered[, final_order]
filtered_filtered <- filtered_filtered[, final_order]

# 验证列顺序是否一致
all.equal(colnames(protein_filtered), colnames(cloneType_filtered))
all.equal(colnames(cloneType_filtered), colnames(Isotype_filtered))
all.equal(colnames(Isotype_filtered), colnames(filtered_filtered))

# 查看新的列顺序
print(colnames(protein_filtered))






# 根据 scBCR_RNA_PB$group2 和 scBCR_RNA_PB$sample2 的对应关系创建 factor

# 首先创建一个样本到分组的映射
sample_group_mapping <- c(
  # AA samples
  setNames(rep("AA", length(c("S28", "S30", "S31", "S32", "S33", "S34", "S35", "S36"))),
           c("S28", "S30", "S31", "S32", "S33", "S34", "S35", "S36")),
  
  # AC samples
  setNames(rep("AC", length(c("S1", "S10", "S11", "S12", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9"))),
           c("S1", "S10", "S11", "S12", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9")),
  
  # CA samples
  setNames(rep("CA", length(c("S47", "S48", "S49", "S50", "S51", "S52", "S53", "S54", "S55", "S56", "S57", "S58", "S59", "S60", "S61", "S62", "S63", "S64", "S65", "S66"))),
           c("S47", "S48", "S49", "S50", "S51", "S52", "S53", "S54", "S55", "S56", "S57", "S58", "S59", "S60", "S61", "S62", "S63", "S64", "S65", "S66")),
  
  # CC samples
  setNames(rep("CC", length(c("S15", "S16", "S17", "S18", "S19", "S20", "S21", "S22", "S23", "S24", "S25", "S26", "S27"))),
           c("S15", "S16", "S17", "S18", "S19", "S20", "S21", "S22", "S23", "S24", "S25", "S26", "S27"))
)

# 创建 factor 变量
group_factor <- factor(sample_group_mapping[colnames(protein_filtered)], 
                       levels = c("AA", "AC", "CA", "CC"))

# 查看结果
print("factor 变量：")
print(table(group_factor))  # 显示各个水平及其频数
print("前几个值：")
print(head(group_factor))



# 创建新的 list
scBCR_data <- list()

# 将数据框和 factor 添加到 list 中
scBCR_data$data.train <- list(
  protein_filtered,          # [[1]] protein matrix
  cloneType_filtered,       # [[2]] cloneType proportion
  Isotype_filtered,         # [[3]] Isotype proportion
  group_factor,             # [[4]] group factor
  filtered_filtered         # [[5]] filtered matrix
)

# 为 list 中的元素添加名称
names(scBCR_data$data.train) <- c(
  "protein_matrix",
  "cloneType_proportion",
  "Isotype_proportion",
  "group_factor",
  "filtered_matrix"
)

# 验证 list 的结构
print("List 结构：")
str(scBCR_data)

# 验证 factor 变量
print("\nFactor 变量：")
print(table(scBCR_data$data.train$group_factor))

# 验证数据框维度
print("\n各数据框维度：")
print(sapply(scBCR_data$data.train[c(1,2,3,5)], dim))












# set a list of all the X dataframes
data = list(
  protein = t(as.matrix(sapply(scBCR_data$data.train$protein_matrix, as.numeric))),
  cloneType = t(as.matrix(scBCR_data$data.train$cloneType_proportion)),
  Isotype = t(as.matrix(scBCR_data$data.train$Isotype_proportion)),
  gene = t(as.matrix(scBCR_data$data.train$filtered_matrix))
)
# 检查数据维度
lapply(data, dim)

# 检查数据类型
lapply(data, class)
lapply(data, function(x) typeof(head(x[,1])))

# 设置响应变量
Y = scBCR_data$data.train$group_factor
summary(Y)


list.keepX = c(25, 25) # select arbitrary values of features to keep
list.keepY = c(25, 25)

# generate three pairwise PLS models
pls1 <- spls(data[["protein"]], data[["gene"]], 
             keepX = list.keepX, keepY = list.keepY)
list.keepX = c(25, 25) # select arbitrary values of features to keep
list.keepY = c(5, 5)
pls2 <- spls(data[["protein"]], data[["cloneType"]], 
             keepX = list.keepX, keepY = list.keepY)
pls3 <- spls(data[["protein"]], data[["Isotype"]], 
             keepX = list.keepX, keepY = list.keepY)

# plot features of first PLS
plotVar(pls1, cutoff = 0.5, title = "(a) protein vs gene", 
        legend = c("protein", "gene"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

# plot features of second PLS
plotVar(pls2, cutoff = 0.1, title = "(b) protein vs cloneType", 
        legend = c("protein", "cloneType"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

# plot features of third PLS
plotVar(pls3, cutoff = 0.1, title = "(c) protein vs Isotype", 
        legend = c("protein", "Isotype"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

# calculate correlation of miRNA and mRNA
cor(pls1$variates$X, pls1$variates$Y) 
# calculate correlation of miRNA and proteins
cor(pls2$variates$X, pls2$variates$Y) 
# calculate correlation of mRNA and proteins
cor(pls3$variates$X, pls3$variates$Y) 
#The moderate to high correlation between these features means that an appropriate value for the design matrix would be about ~0.8 - 0.9. However, as discussed in the N-Integration Methods page, values above 0.5 will cause a reduction in predictive ability of the model - and prediction is what is desired in this context. Hence, a value of 0.1 will be used to prioritise the discriminative ability of the model.
# for square matrix filled with 0.1s
# 完整的检查流程

# 详细检查数据结构
# 1. 首先检查数据的基本统计信息
print("每列的方差：")
col_vars = apply(data$protein, 2, var)
summary(col_vars)










# 2. 使用正确的方式调用 nearZeroVar
# 确保加载 caret 包
library(caret)
protein_nzv = nearZeroVar(data$protein)
print("近零方差特征的位置：")
print(protein_nzv)

# 3. 使用正索引方式选择列
keep_cols = setdiff(1:ncol(data$protein), protein_nzv)
data$protein = data$protein[, keep_cols]

gene_nzv = nearZeroVar(data$gene)
print("近零方差特征的位置：")
print(gene_nzv)

# 3. 使用正索引方式选择列
keep_cols = setdiff(1:ncol(data$gene), gene_nzv)
data$gene = data$gene[, keep_cols]










# 4. 验证结果
#cat("\n处理后的蛋白质数据维度：", dim(data$protein)[2], "列\n")
# 对所有数据块进行近零方差特征处理
#for(block in names(data)) {
#  cat("\n处理", block, "数据块：\n")
  # 获取近零方差特征
#  nzv = nearZeroVar(data[[block]])
#  if(length(nzv) > 0) {
    # 使用正索引方式移除近零方差特征
#    keep_cols = setdiff(1:ncol(data[[block]]), nzv)
#    data[[block]] = data[[block]][, keep_cols, drop = FALSE]
 #   cat("移除了", length(nzv), "个近零方差特征\n")
#    cat("剩余特征数量：", ncol(data[[block]]), "\n")
#  } else {
#    cat("未发现近零方差特征\n")
#  }
#}
# 显示处理后每个数据块的维度
#cat("\n处理后各数据块的维度：\n")
#for(block in names(data)) {
#  cat(block, ": ", dim(data[[block]])[2], " 特征\n", sep="")
#}
# 然后构建模型
# 对每个数据块进行预处理
#for(block in names(data)) {
  # 跳过小型数据块（如 cloneType 和 Isotype）
#  if(ncol(data[[block]]) <= 10) next
  # 计算相关性矩阵
#  cor_matrix = cor(data[[block]])
  # 找出高度相关的特征 (先使用较高的阈值 0.98)
 # high_cor = findCorrelation(cor_matrix, cutoff = 0.80)
  
 # if(length(high_cor) > 0) {
    # 移除高度相关的特征
 #   data[[block]] = data[[block]][, -high_cor, drop = FALSE]
 #   cat(block, "数据块移除了", length(high_cor), 
 #       "个高度相关特征，剩余", ncol(data[[block]]), "个特征\n")
#  }
#}
# 显示处理后的维度
#cat("\n处理后各数据块的维度：\n")
#for(block in names(data)) {
#  cat(block, ": ", dim(data[[block]])[2], " 特征\n", sep="")
#}

# 重新构建模型
basic.diablo.model = block.splsda(X = data, 
                                  Y = Y, 
                                  ncomp = 4, 
                                  design = design,
                                  scale = TRUE)

# 运行交叉验证
# To choose the number of components for the final DIABLO model, the function perf() is run with 10-fold cross-validation repeated 10 times.    run component number tuning with repeated CV
perf.diablo = perf(basic.diablo.model, validation = 'Mfold', 
                   folds = 10, nrepeat = 10) 






#上述方法仍然不起作用
# 使用 PCA 进行降维
library(stats)

for(block in names(data)) {
  # 只处理大型数据块
  if(ncol(data[[block]]) > 1000) {
    # 进行PCA
    pca_result = prcomp(data[[block]], scale. = TRUE)
    
    # 选择解释95%方差所需的主成分数
    var_explained = cumsum(pca_result$sdev^2/sum(pca_result$sdev^2))
    n_comp = which(var_explained >= 0.95)[1]
    
    # 使用PCA转换后的数据
    data[[block]] = pca_result$x[, 1:n_comp, drop = FALSE]
    
    cat(block, "数据块通过PCA降维到", n_comp, "个特征\n")
  }
}

# 重新构建模型
basic.diablo.model = block.splsda(X = data, 
                                  Y = Y, 
                                  ncomp = 3, 
                                  design = design,
                                  scale = TRUE)

# 运行交叉验证
perf.diablo = perf(basic.diablo.model, 
                   validation = 'Mfold', 
                   folds = 10, 
                   nrepeat = 10)
















data = list(
  protein = t(as.matrix(sapply(scBCR_data$data.train$protein_matrix, as.numeric))),
  cloneType = t(as.matrix(scBCR_data$data.train$cloneType_proportion)),
  Isotype = t(as.matrix(scBCR_data$data.train$Isotype_proportion)),
  gene = t(as.matrix(scBCR_data$data.train$filtered_matrix))
)



rm(CountsSeuratB_umap_J)


plot(perf.diablo) # plot output of tuning
#From the performance plot above (Figure 2) we observe that both overall and balanced error rate (BER) decrease from 1 to 2 components. The standard deviation indicates a potential slight gain in adding more components. The centroids.dist distance seems to give the best accuracy (see Supplemental Material in [6]). Considering this distance and the BER, the output $choice.ncomp indicates the optimal number of components for the final DIABLO model.
# set the optimal ncomp value
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"] 
# show the optimal choice for ncomp for each dist metric
perf.diablo$choice.ncomp$WeightedVote 
# set grid of values for each component to test
test.keepX = list (gene = c(5:9, seq(10, 18, 2), seq(20,30,5)), 
                   Isotype = c(2:4), 
                   cloneType = c(2:4),
                   protein = c(5:9, seq(10, 18, 2), seq(20,30,5)))

# run the feature selection tuning
tune.ARDS = tune.block.splsda(X = data, Y = Y, ncomp = ncomp, 
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 10, nrepeat = 1,
                              dist = "centroids.dist")


save(tune.ARDS, data, file = "所有单细胞数据加蛋白组学加BCR两个参数.Rdata")
list.keepX = tune.ARDS$choice.keepX # set the optimal values of features to retain
list.keepX
# set the optimised DIABLO model
final.diablo.model = block.splsda(X = data, Y = Y, ncomp = ncomp, 
                                  keepX = list.keepX, design = design)
final.diablo.model$design # design matrix for the final model
#可以使用函数 selectVar() 提取所选变量，例如在 gene 块中，如下所示。请注意，可以从perf()函数的输出中提取所选变量的稳定性。
# the features selected to form the first component
selectVar(final.diablo.model, block = 'gene', comp = 3)$gene$name 
##  [1] "ZNF552"  "KDM4B"   "CCNA2"   "LRIG1"   "PREX1"   "FUT8"    "C4orf34"
##  [8] "TTC39A"  "ASPM"    "SLC43A3" "MEX3A"   "SEMA3C"  "E2F1"    "STC2"   
## [15] "FMNL2"   "LMO4"    "MED13L"  "DTWD2"   "CSRP2"   "NTN4"    "KIF13B" 
## [22] "NCAPG2"  "SLC19A2" "EPHB3"   "FAM63A"
selectVar(final.diablo.model, block = 'protein', comp = 3)$protein$name
#plotDIABLO() is a diagnostic plot to check whether the correlation between components from each data set has been maximised as specified in the design matrix. We specify which dimension to be assessed with the ncomp argument.
plotDiablo(final.diablo.model, ncomp = 4)
#如图 3 所示，每个数据集的第一个分量彼此高度相关（由左下角的大数字表示）。与样本亚型相关的颜色和省略号表示每个成分区分不同肿瘤亚型的鉴别能力。对于第一个组件，每个子类型的质心都是不同的，但每个样本组在其置信椭圆中存在适度的重叠。
#带有该函数的plotIndiv()样本图将每个样本投影到每个模块的元件所跨越的空间中（图 4）。使用此图可以更好地评估样本的聚类。mRNA 数据的聚类质量似乎最高，而 miRNA 聚类质量最低。这表明 mRNA 在模型中可能具有更多的鉴别能力。

plotIndiv(final.diablo.model, ind.names = FALSE, legend = TRUE, 
          title = 'DIABLO Sample Plots')
#In the arrow plot below (Figure), the start of the arrow indicates the centroid between all data sets for a given sample and the tips of the arrows indicate the location of that sample in each block. Such graphics highlight the agreement between all data sets at the sample level. While somewhat difficult to interpret, even qualitatively, Figure 5 shows that the agreement within the LumA group seems to be the highest and lowest in the Her2 group.
#在下面的箭头图（图）中，箭头的开头表示给定样品的所有数据集之间的质心，箭头的尖端表示该样品在每个块中的位置。此类图形突出显示了样本级别所有数据集之间的一致性。虽然有些难以解释，甚至在定性上也很难解释，但图 5 显示LumA，组内的一致性似乎是Her2组中最高和最低的。

plotArrow(final.diablo.model, ind.names = FALSE, legend = TRUE, 
          title = 'DIABLO')

#Several graphical outputs are available to visualise and mine the associations between the selected variables.
#可以使用多个图形输出来可视化和挖掘所选变量之间的关联。

#The best starting point to evaluate the correlation structure between variables is with the correlation circle plot, depicted in Figure 6. A majority of the miRNA variables are positively correlated with the first component while the mRNA variables seem to separate along this dimension. These first two components correlate highly with the selected variables from the proteomics dataset. From this, the correlation of each selected feature from all three datasets can be evaluated based on their proximity. s
#评估变量之间相关结构的最佳起点是使用相关圆图，如图 6 所示。大多数 miRNA 变量与第一个组分呈正相关，而 mRNA 变量似乎沿此维度分离。前两个组成部分与蛋白质组学数据集中选定的变量高度相关。由此，可以根据所有三个数据集中每个选定要素的接近度来评估它们的相关性。s

plotVar(final.diablo.model, var.names = FALSE, 
        style = 'graphics', legend = TRUE,
        pch = c(16, 17, 15,18), 
        cex = c(2,2,2,2), 
        col = c('darkorchid', 'brown1', 'lightgreen',"orange"))

#The circos plot is exclusive to integrative frameworks and represents the correlations between variables of different types, represented on the side quadrants. From Figure 7, it seems that the miRNA variables are almost entirely negatively correlated with the other two dataframes. The proteomics features are the opposite, such that they display primarily positive correlations while the mRNA variables are more mixed. Note that these correlations are above a value of 0.7 (cutoff = 0.7). All the interpretations made above are only relevant for features with very strong correlations.
#circos 图是综合框架独有的，表示不同类型变量之间的相关性，在侧象限上表示。从图 7 中可以看出，miRNA 变量似乎几乎完全与其他两个数据帧呈负相关。蛋白质组学特征正好相反，因此它们主要显示正相关，而 mRNA 变量则更加混合。请注意，这些相关性高于 0.7 （cutoff = 0.7 ） 的值。上面所做的所有解释仅与具有非常强相关性的特征相关。

circosPlot(final.diablo.model, cutoff = 0.7, line = TRUE,
           color.blocks= c('darkorchid', 'brown1', 'lightgreen',"orange"),
           color.cor = c("chocolate3","grey20"), size.labels = 1.5)

#Another visualisation of the correlations between the different types of variables is the relevance network, which is also built on the similarity matrix (as is the circos plot). Each colour represents a type of variable. Figure 8 shows this network which has a lower cutoff an Figure 7 (cutoff = 0.4). Two distinct clusters can be observed, though due to the density of the plot the relationships within the cluster are hard to determine. The interactive version of this plot would be useful here.
#不同类型变量之间相关性的另一种可视化方式是相关性网络，它也建立在相似性矩阵之上（就像 circos 图一样）。每种颜色代表一种类型的变量。图 8 显示了这个网络，其截止值较低，图 7 （cutoff = 0.4 ）。可以观察到两个不同的集群，但由于图的密度，集群内的关系很难确定。此图的交互式版本在此处将很有用。

network(final.diablo.model, blocks = c(1,2,3),
        color.node = c('darkorchid', 'brown1', 'lightgreen'), cutoff = 0.4)

#The network can be saved in a .gml format to be input into the software Cytoscape, using the R package igraph. An example is shown directly below.
#可以使用R软件包igraph将网络保存为输入到 Cytoscape 软件中的.gml格式。下面直接显示了一个示例。

library(igraph)
my.network = network(final.diablo.model, blocks = c(1,2,3),
                     color.node = c('darkorchid', 'brown1', 'lightgreen'), cutoff = 0.4)
write.graph(my.network$gR, file = "myNetwork.gml", format = "gml")
#The function plotLoadings() visualises the loading weights of each selected variable on each component and each data set (Figure 9). The colour indicates the class in which the variable has the maximum level of expression (contrib = 'max') using the median (method = 'median'). Figure 9 depicts the loading values for the second dimension.
#该函数plotLoadings()可视化每个组件和每个数据集上每个选定变量的加载权重（图 9）。颜色表示变量具有最大表达式级别 （contrib = 'max' ） 的类，使用中位数 （method = 'median' ）。图 9 描述了第二个维度的载荷值。
plotLoadings(final.diablo.model, comp = 2, contrib = 'max', method = 'median')

#该cimDIABLO()函数是一个聚类图像图，专门用于表示每个样品的多组学分子特征表达。从图 10 中，可以确定一组样品在一组特征中的均一表达水平区域。例如，样品Her2是唯一一组特定蛋白质表现出极高表达水平的组（如图 10 中部底部的红色小块所示）。这表明这些特征对于此子类型具有相当的区别。

cimDiablo(final.diablo.model)


#我们使用函数 perf() .该方法在对象final.diablo.model输入的预先指定的参数上运行block.splsda()模型，但在交叉验证的样本上运行模型。然后，我们评估对遗漏样本的预测准确性。

#In addition to the usual (balanced) classification error rates, predicted dummy variables and variates, as well as the stability of the selected features, the perf() function for DIABLO outputs the performance based on Majority Vote (each data set votes for a class for a particular test sample) or a weighted vote, where the weight is defined according to the correlation between the latent component associated to a particular data set and the outcome.
#除了通常的（平衡的）分类错误率、预测的虚拟变量和变量以及所选特征的稳定性外，DIABLO perf() 的函数还根据多数投票（每个数据集投票支持特定测试样本的一个类）或加权投票来输出性能，其中权重是根据与特定数据集相关的潜在成分与结果之间的相关性定义的。

#Since the tune() function was used with the centroid.dist argument, the outputs of the perf() function for that same distance are examined. The following code may take a few minutes to run.
#由于该tune()函数与centroid.dist参数一起使用，因此将检查相同距离的perf()函数输出。以下代码可能需要几分钟才能运行。

# run repeated CV performance evaluation
perf.diablo = perf(final.diablo.model, validation = 'Mfold', 
                   M = 10, nrepeat = 10, 
                   dist = 'centroids.dist') 

perf.diablo$MajorityVote.error.rate

#From the above output, it can be seen that the error rate across the board is quite low, indicating the constructed DIABLO model does a fairly good job of classifying novel samples.
#从上面的输出中可以看出，整体错误率相当低，这表明构建的 DIABLO 模型在对新样本的分类方面做得相当好。

#An AUC plot per block can also be obtained using the function auroc(). The interpretation of this output may not be particularly insightful in relation to the performance evaluation of our methods, but can complement the statistical analysis..
#也可以使用函数 auroc() 获得每个块的 AUC 图 。对于我们方法的性能评估，这个输出的解释可能并不特别有洞察力，但可以补充统计分析。

auc.splsda = auroc(final.diablo.model, roc.block = "miRNA", 
                   roc.comp = 2, print = FALSE)

#The predict() function predicts the class of samples from a test set. In our specific case, one data set is missing in the test set but the method can still be applied. Make sure the name of the blocks correspond exactly.
#该predict()函数预测测试集中的样本类。在我们的特定情况下，测试集中缺少一个数据集，但仍然可以应用该方法。确保块的名称完全对应。

data.test.TCGA = list(mRNA = breast.TCGA$data.test$mrna,
                      miRNA = breast.TCGA$data.test$mirna)

predict.diablo = predict(final.diablo.model, newdata = data.test.TCGA)
#The confusion table compares the real subtypes with the predicted subtypes for a 2-component model, for the distance of interest. This model performs quite well as it makes only two errors:
#  混淆表将 2 分量模型的真实子类型与预测子类型进行比较，以了解感兴趣距离。这个模型表现得相当不错，因为它只犯了两个错误：

confusion.mat = get.confusion_matrix(truth = breast.TCGA$data.test$subtype,
                                     predicted = predict.diablo$WeightedVote$centroids.dist[,2])
confusion.mat
##       predicted.as.Basal predicted.as.Her2 predicted.as.LumA
## Basal                 20                 1                 0
## Her2                   0                14                 0
## LumA                   0                 1                34
#These two errors correspond to a very low balanced error rate.
#这两个误差对应于非常低的平衡误差率。

get.BER(confusion.mat)
## [1] 0.02539683