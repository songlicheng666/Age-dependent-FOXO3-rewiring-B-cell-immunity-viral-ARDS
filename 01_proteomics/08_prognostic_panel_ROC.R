# =============================================================================
# Script: 08_prognostic_panel_ROC.R
# Description: Integration of single-cell major lineage proportions with proteomics;
# logistic regression and ROC curve analysis for prognostic protein panels
# Author: Licheng Song
# Institution: Chinese PLA General Hospital (Eighth Medical Center)
# Study: Age-dependent B-cell immune dysregulation in virus-associated ARDS
# Year: 2025
# Input: ALL_sample.final_matrix.xlsx, filtered_matrix (single-cell proportions)
# Output: ROC_curves.pdf, AUC_comparison.csv
# =============================================================================



# 加载 readxl 包
library(readxl)

# 读取 Excel 文件
prtein_matrix <- read_excel("/home/lungtissue/蛋白质组学整合分析/ALL_sample.final_matrix.xlsx")

# 查看数据的前几行
head(prtein_matrix)

# 查看数据的基本信息
str(prtein_matrix)
# 将第一列设置为行名并删除该列
rownames(prtein_matrix) <- prtein_matrix[[1]]
# 查看结果
head(prtein_matrix)
prtein_matrix <- as.matrix(prtein_matrix)
prtein_matrix <- prtein_matrix[,-1]
# 将 filtered_matrix 转换为数据框
filtered_matrix <- as.data.frame(filtered_matrix)

prtein_matrix <- as.data.frame(prtein_matrix)
# 获取所有数据框的列名
cols_protein <- colnames(prtein_matrix)
cols_filtered <- colnames(filtered_matrix)

# 找出共同的列名
common_cols <- Reduce(intersect, list(cols_protein,    cols_filtered))

# 使用共同的列名从每个数据框中提取相应的列
protein_filtered <- prtein_matrix[, common_cols]
filtered_filtered <- filtered_matrix[, common_cols]

# 查看结果
print(paste("共同的列数：", length(common_cols)))
print("共同的列名：")
print(common_cols)

ordered_samples <- c(AA_samples, AC_samples, CA1_samples, CA2_samples, CC_samples)

# 获取四个数据框共同的列
common_cols <- Reduce(intersect, list(colnames(protein_filtered), 
                                      colnames(filtered_filtered)))

# 找出在 common_cols 中的样本列（与 ordered_samples 相交的列）
sample_cols <- intersect(ordered_samples, common_cols)

# 找出非样本列（可能是基因名或其他特征）
other_cols <- setdiff(common_cols, sample_cols)

# 创建最终的列顺序（先放其他列，再放样本列，样本列按分组顺序排列）
final_order <- c(other_cols, sample_cols)

# 重新排序每个数据框的列
protein_filtered <- protein_filtered[, final_order]
filtered_filtered <- filtered_filtered[, final_order]

# 验证列顺序是否一致
all.equal(colnames(protein_filtered), colnames(filtered_filtered))
# 查看新的列顺序
print(colnames(protein_filtered))






# 根据 scBCR_RNA_PB$group2 和 scBCR_RNA_PB$sample2 的对应关系创建 factor

# 创建 factor 变量
group_factor <- factor(sample_group_mapping[colnames(protein_filtered)], 
                       levels = c("AA", "AC", "CA1", "CA2", "CC"))

# 查看结果
print("factor 变量：")
print(table(group_factor))  # 显示各个水平及其频数
print("前几个值：")
print(head(group_factor))

head(filtered_filtered)
# 提取基因名称列表
gene_names <- rownames(filtered_filtered)
# 定义需要去除的基因类型规则
mito_genes <- grepl("^MT-", gene_names, ignore.case = TRUE)   # 线粒体基因（以MT-开头）
ribo_genes <- grepl("^RPL|^RPS", gene_names, ignore.case = TRUE) # 核糖体基因（RPL/RPS开头）
linc_genes <- grepl("^LINC", gene_names)
PCDH_genes <- grepl("^PCDH", gene_names)# LINC开头的基因
TRAV_genes <- grepl("^TRAV", gene_names)# LINC开头的基因
TRBV_genes <- grepl("^TRBV", gene_names)# LINC开头的基因
has_dots <- grepl("\\.", gene_names)   
has_line <- grepl("\\-", gene_names)   # 包含小数点的基因
# 综合过滤条件（使用逻辑或|组合所有需要排除的条件）
remove_genes <- mito_genes | ribo_genes | linc_genes | has_dots| PCDH_genes| TRAV_genes |  has_line| TRBV_genes
# 保留符合要求的基因（取反操作）
filtered_filtered <- filtered_filtered[!remove_genes, ]



head(protein_filtered)
# 提取基因名称列表
gene_names <- rownames(protein_filtered)
# 定义需要去除的基因类型规则
mito_genes <- grepl("^MT-", gene_names, ignore.case = TRUE)   # 线粒体基因（以MT-开头）
ribo_genes <- grepl("^RPL|^RPS", gene_names, ignore.case = TRUE) # 核糖体基因（RPL/RPS开头）
linc_genes <- grepl("^LINC", gene_names)                     # LINC开头的基因
PCDH_genes <- grepl("^PCDH", gene_names)# LINC开头的基因
TRAV_genes <- grepl("^TRAV", gene_names)# LINC开头的基因
TRBV_genes <- grepl("^TRBV", gene_names)# LINC开头的基因
has_dots <- grepl("\\.", gene_names)   
has_line <- grepl("\\-", gene_names)   # 包含小数点的基因
# 综合过滤条件（使用逻辑或|组合所有需要排除的条件）
remove_genes <- mito_genes | ribo_genes | linc_genes | has_dots| PCDH_genes| TRAV_genes |  has_line| TRBV_genes
# 保留符合要求的基因（取反操作）
protein_filtered <- protein_filtered[!remove_genes, ]






# 创建新的 list
scBCR_data <- list()

# 将数据框和 factor 添加到 list 中
scBCR_data$data.train <- list(
  protein_filtered,          # [[1]] protein matrix
  group_factor,             # [[4]] group factor
  filtered_filtered         # [[5]] filtered matrix
)

# 为 list 中的元素添加名称
names(scBCR_data$data.train) <- c(
  "protein_matrix","group_factor",
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
print(sapply(scBCR_data$data.train[c(1,2,3)], dim))












# set a list of all the X dataframes
data = list(
  protein = t(as.matrix(scBCR_data$data.train$protein_matrix)),
  gene = t(as.matrix(scBCR_data$data.train$filtered_matrix))
)
data$protein <- apply(data$protein, 2, as.numeric)
# 检查数据维度
lapply(data, dim)

# 检查数据类型
lapply(data, class)
lapply(data, function(x) typeof(head(x[,1])))


#两种数据的sPLS
X <- data$gene # use the gene expression data as the X matrix
Y <- data$protein # use the clinical data as the Y matrix



#Preliminary Analysis with PCA
pca.gene <- pca(X, ncomp = 10, center = TRUE, scale = TRUE)
pca.protein <- pca(Y, ncomp = 10, center = TRUE, scale = TRUE)

plot(pca.gene)
plot(pca.protein)


plotIndiv(pca.gene, comp = c(1, 2), 
          group = scBCR_data$data.train$group_factor, 
          ind.names = F, 
          legend = TRUE, title = 'ARDS gene, PCA comp 1 - 2')

plotIndiv(pca.protein, comp = c(1, 2), 
          group = scBCR_data$data.train$group_factor, 
          ind.names = F, 
          legend = TRUE, title = 'ARDS protein, PCA comp 1 - 2')


spls.liver <- spls(X = X, Y = Y, ncomp = 5, mode = 'regression')



# repeated CV tuning of component count
perf.spls.liver <- perf(spls.liver, validation = 'Mfold',
                        folds = 10, nrepeat = 5) 
#看官网解释
plot(perf.spls.liver, criterion = 'Q2.total')

# set range of test values for number of variables to use from X dataframe
list.keepX <- c(seq(20, 50, 5))
# set range of test values for number of variables to use from Y dataframe
list.keepY <- c(20, 50, 5) 


tune.spls.liver <- tune.spls(X, Y, ncomp = 2,
                             test.keepX = list.keepX,
                             test.keepY = list.keepY,
                             nrepeat = 1, folds = 10, # use 10 folds
                             mode = 'regression', measure = 'cor') 
plot(tune.spls.liver)         # use the correlation measure for tuning
#The optimal number of features to use for both datasets can be extracted through the below calls.
tune.spls.liver$choice.keepX
## comp1 comp2 
##    25    20
tune.spls.liver$choice.keepY
## comp1 comp2 
##     3     6
#These values will be stored to form the final model.

# extract optimal number of variables for X dataframe
optimal.keepX <- tune.spls.liver$choice.keepX 

# extract optimal number of variables for Y datafram
optimal.keepY <- tune.spls.liver$choice.keepY

optimal.ncomp <-  length(optimal.keepX) # extract optimal number of components


# use all tuned values from above
final.spls.liver <- spls(X, Y, ncomp = optimal.ncomp, 
                         keepX = optimal.keepX,
                         keepY = optimal.keepY,
                         mode = "regression") # explanitory approach being used, 
# hence use regression mode



plotIndiv(final.spls.liver, ind.names = FALSE, 
          rep.space = "X-variate", # plot in X-variate subspace
          group = liver.toxicity$treatment$Time.Group, # colour by time group
          pch = as.factor(liver.toxicity$treatment$Dose.Group), 
          col.per.group = color.mixo(1:4), 
          legend = TRUE, legend.title = 'Time', legend.title.pch = 'Dose')

plotIndiv(final.spls.liver, ind.names = FALSE,
          rep.space = "Y-variate", # plot in Y-variate subspace
          group = liver.toxicity$treatment$Time.Group, # colour by time group
          pch = as.factor(liver.toxicity$treatment$Dose.Group), 
          col.per.group = color.mixo(1:4), 
          legend = TRUE, legend.title = 'Time', legend.title.pch = 'Dose')


plotIndiv(final.spls.liver, ind.names = FALSE, 
          rep.space = "XY-variate", # plot in averaged subspace
          group = liver.toxicity$treatment$Time.Group, # colour by time group
          pch = as.factor(liver.toxicity$treatment$Dose.Group), # select symbol
          col.per.group = color.mixo(1:4),                      # by dose group
          legend = TRUE, legend.title = 'Time', legend.title.pch = 'Dose')


#The plots can also be represented in 3D using style = '3d', as seen in Figure 7.

col.tox <- color.mixo(as.numeric(as.factor(liver.toxicity$treatment[, 4]))) # create set of colours
plotIndiv(final.spls.liver, ind.names = FALSE, 
          rep.space = "XY-variate", # plot in averaged subspace
          axes.box = "both", col = col.tox, style = '3d')



#The plotArrow() option is useful in this context to visualise the level of agreement between data sets. It can be seen in Figure 8 that specific groups of samples seem to be located far apart from one data set to the other, indicating a potential discrepancy between the information extracted.

plotArrow(final.spls.liver, ind.names = FALSE,
          group = liver.toxicity$treatment$Time.Group, # colour by time group
          col.per.group = color.mixo(1:4),
          legend.title = 'Time.Group')

#The stability of a given feature is defined as the proportion of cross validation folds (across repeats) where it was selected for to be used for a given component. Stability values (for the X component) can be extracted via perf.spls.liver$features$stability.X. Figure 9(a) and (b) depict these stabilities as histograms for the first two components respectively. Both components use the same set of features fairly consistently across repeated folds, meaning the variance of the data can be attributed to a specific set of features. The first component displays this quality much more than the second.

# form new perf() object which utilises the final model
perf.spls.liver <- perf(final.spls.liver, 
                        folds = 5, nrepeat = 10, # use repeated cross-validation
                        validation = "Mfold", 
                        dist = "max.dist",  # use max.dist measure
                        progressBar = FALSE)

# plot the stability of each feature for the first two components, 
# 'h' type refers to histogram
par(mfrow=c(1,2)) 
plot(perf.spls.liver$features$stability.X[[1]], type = 'h',
     ylab = 'Stability',
     xlab = 'Features',
     main = '(a) Comp 1', las =2,
     xlim = c(0, 150))
plot(perf.spls.liver$features$stability.X$comp2, type = 'h',
     ylab = 'Stability',
     xlab = 'Features',
     main = '(b) Comp 2', las =2,
     xlim = c(0, 300))


#The relationship between the features and components can be explored using a correlation circle plot. This highlights the contributing variables that together explain the covariance between the two datasets. Specific subsets of molecules can be further investigated. Figure 10 shows the correlations between selected genes (names not shown), between selected clinical parameters and the relationship between sets of genes and certain clinical parameters.

plotVar(final.spls.liver, cex = c(3,4), var.names = c(FALSE, TRUE))


#Two substructures can be observed from this network. First, the small cluster to the top right, where clinical feature ALB.g.dL is shown to be negatively correlated with three genetic features and none else. The second substructure includes three clinical features which are (mostly) positively correlated with a large set of genetic features.

color.edge <- color.GreenRed(50)  # set the colours of the connecting lines

# X11() # To open a new window for Rstudio
network(final.spls.liver, comp = 1:2,
        cutoff = 0.7, # only show connections with a correlation above 0.7
        shape.node = c("rectangle", "circle"),
        color.node = c("cyan", "pink"),
        color.edge = color.edge,
        save = 'png', # save as a png to the current working directory
        name.save = 'sPLS Liver Toxicity Case Study Network Plot')



#Another complementary plot used for exploration of feature structure in sPLS is the Cluster Image Map (CIM). The same technique of using the X11() function or save/save.name parameters may be required here too. Figure 12 shows that the clinical variables can be separated into three clusters, each of them either positively or negatively associated with two groups of genes. This is similar to what we have observed in Figure 11. The large red cluster corresponds to the largest substructure in the network while the large blue cluster was not depicted in Figure 11 due to the use of the cutoff parameter.

cim(final.spls.liver, comp = 1:2, xlab = "clinic", ylab = "genes")
