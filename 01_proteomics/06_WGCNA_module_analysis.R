# =============================================================================
# Script: 06_WGCNA_module_analysis.R
# Description: Weighted gene co-expression network analysis (WGCNA):
# soft-threshold selection, module detection (dynamic tree cut),
# module-trait correlation (biweight midcorrelation), hub protein identification
# Author: Licheng Song
# Institution: Chinese PLA General Hospital (Eighth Medical Center)
# Study: Age-dependent B-cell immune dysregulation in virus-associated ARDS
# Year: 2025
# Input: ALL_sample.final_matrix.xlsx (log2-transformed)
# Output: WGCNA_module_trait_heatmap.pdf, hub_proteins.csv
# =============================================================================

library(WGCNA)
library(reshape2)
library(stringr)
library(DESeq2)
setwd("C:/Users/songl/OneDrive/桌面/多组学分析/WGCNA")
# 加载包
library(DESeq2)
library(readxl)
# 读取 Excel 文件
prtein_matrix <- read_excel("C:/Users/songl/OneDrive/桌面/多组学分析/ALL_sample.final_matrix.xlsx")

# 查看数据的前几行
head(prtein_matrix)

# 查看数据的基本信息
str(prtein_matrix)
# 将第一列设置为行名并删除该列
rownames(prtein_matrix) <- prtein_matrix[[1]]
# 查看结果
head(prtein_matrix)
# 数据准备
gene_names <- prtein_matrix[, 1]
data_matrix <- as.matrix(prtein_matrix[, -1])
data_matrix <- apply(data_matrix, 2, as.numeric)
rownames(data_matrix) <- rownames( prtein_matrix)
# 首先创建一个按分组的样本顺序列表
AA_samples <- c("S28", "S30", "S31", "S32", "S33", "S34", "S35", "S36", "S37")
AC_samples <- c("S1", "S10", "S11", "S12", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9")
CA1_samples <- c("S47", "S48", "S49", "S50", "S51", "S52", "S53", "S54", "S55", "S56", "S57")
CA2_samples <- c("S58", "S59", "S60", "S61", "S62", "S63", "S64", "S65", "S66")
CC_samples <- c("S15", "S16", "S17", "S18", "S19", "S20", "S21", "S22", "S23", "S24", "S25", "S26", "S27")
ordered_samples <- c(AA_samples, AC_samples, CA1_samples, CA2_samples, CC_samples)
# 创建分组信息
group <- c(rep("AA", length(AA_samples)),
           rep("AC", length(AC_samples)),
           rep("CA1", length(CA1_samples)),
           rep("CA2", length(CA2_samples)),
           rep("CC", length(CC_samples)))
colData <- data.frame(condition = factor(group))
rownames(colData) <- ordered_samples

# 确保样本顺序一致
data_matrix <- data_matrix[, ordered_samples]
data_matrix <- round(data_matrix)
# 构建 DESeq 数据对象
dds <- DESeqDataSetFromMatrix(countData = data_matrix,
                              colData = colData,
                              design = ~ condition)

# 标准化
dds <- DESeq(dds)

# 方差稳定化变换
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
vsd_data <- assay(vsd)

# Log2(x + 1) 变换
normalized_counts <- counts(dds, normalized = TRUE)
log2_data <- log2(normalized_counts + 1)

# 比较结果
par(mfrow = c(1, 2))
boxplot(vsd_data, main = "VST Transformation", las = 2)
boxplot(log2_data, main = "Log2(x+1) Transformation", las = 2)

# 保存结果
write.csv(vsd_data, "vsd_transformed_data.csv", row.names = TRUE)
write.csv(log2_data, "log2_transformed_data.csv", row.names = TRUE)

# 常规表达矩阵，log2转换后或
# Deseq2的varianceStabilizingTransformation转换的数据
# 如果有批次效应，需要事先移除，可使用removeBatchEffect
# 如果有系统偏移(可用boxplot查看基因表达分布是否一致)，
# 需要quantile normalization

exprMat <- "ALL_sample.final_matrix.csv"

# 官方推荐 "signed" 或 "signed hybrid"
# 为与原文档一致，故未修改
type = "unsigned"

# 相关性计算
# 官方推荐 biweight mid-correlation & bicor
# corType: pearson or bicor
# 为与原文档一致，故未修改
corType = "pearson"

corFnc = ifelse(corType=="pearson", cor, bicor)
# 对二元变量，如样本性状信息计算相关性时，
# 或基因表达严重依赖于疾病状态时，需设置下面参数
maxPOutliers = ifelse(corType=="pearson",1,0.05)

# 关联样品性状的二元变量时，设置
robustY = ifelse(corType=="pearson",T,F)

##导入数据##
dataExpr <- log2_data

dim(dataExpr)
head(dataExpr)[,1:8]
## [1] 3600  134
## 筛选中位绝对偏差前75%的基因，至少MAD大于0.01
## 筛选后会降低运算量，也会失去部分信息
## 也可不做筛选，使MAD大于0即可
m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad >
                                max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]

## 转换为样品在行，基因在列的矩阵
dataExpr <- as.data.frame(t(dataExprVar))

## 检测缺失值
gsg = goodSamplesGenes(dataExpr, verbose = 3)

##  Flagging genes and samples with too many missing values...
##   ..step 1

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:",
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)

dim(dataExpr)

## [1]  134 2697

head(dataExpr)[,1:8]

##       MMT00000051 MMT00000080 MMT00000102 MMT00000149 MMT00000159
## F2_2  -0.02260000 -0.04870000  0.17600000  0.07680000 -0.14800000
## F2_3   0.06170000  0.05820000 -0.18900000  0.18600000  0.17700000
## F2_14 -0.12900000 -0.04830000 -0.06500000  0.21400000 -0.13200000
## F2_15  0.08710000 -0.03710000 -0.00846000  0.12000000  0.10700000
## F2_19 -0.11500000  0.02510000 -0.00574000  0.02100000 -0.11900000
## F2_20 -0.06502607  0.08504274 -0.01807182  0.06222751 -0.05497686
##       MMT00000207 MMT00000212 MMT00000241
## F2_2   0.06870000  0.06090000 -0.01770000
## F2_3   0.10100000  0.05570000 -0.03690000
## F2_14  0.10900000  0.19100000 -0.15700000
## F2_15 -0.00858000 -0.12100000  0.06290000
## F2_19  0.10500000  0.05410000 -0.17300000
## F2_20 -0.02441415  0.06343181  0.06627665

## 查看是否有离群样品
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers,
                        networkType=type, verbose=5)

par(mfrow = c(1,2))
cex1 = 0.9
# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
# 网络越符合无标度特征 (non-scale)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# 筛选标准。R-square=0.85
abline(h=0.85,col="red")

# Soft threshold与平均连通性
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,
     cex=cex1, col="red")
power = sft$powerEstimate
power
# 无向网络在power小于15或有向网络power小于30内，没有一个power值可以使
# 无标度网络图谱结构R^2达到0.8，平均连接度较高如在100以上，可能是由于
# 部分样品与其他样品差别太大。这可能由批次效应、样品异质性或实验条件对
# 表达影响太大等造成。可以通过绘制样品聚类查看分组信息和有无异常样品。
# 如果这确实是由有意义的生物变化引起的，也可以使用下面的经验power值。
if (is.na(power)){
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))      
                 )
  )
}




# 检查数据中是否有缺失值
sum(is.na(dataExpr))  # 如果结果不是0，需要处理缺失值

# 去除全为零或接近零的行
dataExpr <- dataExpr[rowSums(dataExpr) > 0, ]
# 填补缺失值（用列的中位数替代）
dataExpr[is.na(dataExpr)] <- apply(dataExpr, 2, median, na.rm = TRUE)

##一步法网络构建：One-step network construction and module detection##
# power: 上一步计算的软阈值
# maxBlockSize: 计算机能处理的最大模块的基因数量 (默认5000)；
#  4G内存电脑可处理8000-10000个，16G内存电脑可以处理2万个，32G内存电脑可
#  以处理3万个
#  计算资源允许的情况下最好放在一个block里面。
# corType: pearson or bicor
# numericLabels: 返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
# saveTOMs：最耗费时间的计算，存储起来，供后续使用
# mergeCutHeight: 合并模块的阈值，越大模块越少
# 使用 Pearson 相关性
#corType <- "pearson"

# 或者使用 bicor，但调整 maxPOutliers
corType <- "bicor"
maxPOutliers <- 0.1  # 允许一定比例的异常值
net = blockwiseModules(dataExpr, power = 2, maxBlockSize = 10769,
                       TOMType = type, minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType,
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = paste0(exprMat, ".tom"),
                       verbose = 3)

table(net$colors)
## 灰色的为**未分类**到模块的基因。
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
# Plot the dendrogram and the module colors underneath
# 如果对结果不满意，还可以recutBlockwiseTrees，节省计算时间
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)




# module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
MEs = net$MEs

### 不需要重新计算，改下列名字就好
### 官方教程是重新计算的，起始可以不用这么麻烦
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)

# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap",
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T,
                      xLabelsAngle = 90)
## 如果有表型数据，也可以跟ME数据放一起，一起出图
#MEs_colpheno = orderMEs(cbind(MEs_col, traitData))
#plotEigengeneNetworks(MEs_colpheno, "Eigengene adjacency heatmap",
#                      marDendro = c(3,3,2,4),
#                      marHeatmap = c(3,4,2,2), plotDendrograms = T,
#                      xLabelsAngle = 90)













# 如果采用分步计算，或设置的blocksize>=总基因数，直接load计算好的TOM结果
# 否则需要再计算一遍，比较耗费时间
TOM = TOMsimilarityFromExpr(dataExpr, power=power, corType=corType, networkType=type)
load(net$TOMFiles[1], verbose=T)

## Loading objects:
##   TOM

TOM <- as.matrix(TOM)

dissTOM = 1-TOM
# Transform dissTOM with a power to make moderately strong
# connections more visible in the heatmap
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function

# 这一部分特别耗时，行列同时做层级聚类
TOMplot(plotTOM, net$dendrograms, moduleColors,
        main = "Network heatmap plot, all genes")
probes = colnames(dataExpr)
dimnames(TOM) <- list(probes, probes)

# Export the network into edge and node list files Cytoscape can read
# threshold 默认为0.5, 可以根据自己的需要调整，也可以都导出后在
# cytoscape中再调整
cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = paste(exprMat, ".edges.txt", sep=""),
                               nodeFile = paste(exprMat, ".nodes.txt", sep=""),
                               weighted = TRUE, threshold = 0,
                               nodeNames = probes, nodeAttr = moduleColors)






trait <- "TraitsClean.txt"
# 读入表型数据，不是必须的
if(trait != "") {
  traitData <- read.table(file=trait, sep='\t', header=T, row.names=1,
                          check.names=FALSE, comment='',quote="")
  sampleName = rownames(dataExpr)
  traitData = traitData[match(sampleName, rownames(traitData)), ]
}

### 模块与表型数据关联
if (corType=="pearsoon") {
  modTraitCor = cor(MEs_col, traitData, use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}

## Warning in bicor(x, y, use = use, ...): bicor: zero MAD in variable 'y'.
## Pearson correlation was used for individual columns with zero (or missing)
## MAD.

# signif表示保留几位小数
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
# 设置图形边距
par(mar = c(5, 5, 4, 2))  # 下、左、上、右的边距

# 绘制热图
pdf("GCWNA与临床信息结合.pdf",width = 10,height = 5)
labeledHeatmap(Matrix = modTraitCor, 
               xLabels = colnames(traitData),
               yLabels = colnames(MEs_col),
               cex.lab = 1,
               ySymbols = colnames(MEs_col), 
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix, 
               setStdMargins = FALSE,
               cex.text = 0.5, 
               zlim = c(-1, 1),
               main = paste("Module-trait relationships"))
dev.off()


















## 从上图可以看到MEmagenta与Insulin_ug_l相关

## 模块内基因与表型数据关联

# 性状跟模块虽然求出了相关性，可以挑选最相关的那些模块来分析，
# 但是模块本身仍然包含非常多的基因，还需进一步的寻找最重要的基因。
# 所有的模块都可以跟基因算出相关系数，所有的连续型性状也可以跟基因的表达
# 值算出相关系数。
# 如果跟性状显著相关基因也跟某个模块显著相关，那么这些基因可能就非常重要
# 。

### 计算模块与基因的相关性矩阵

if (corType=="pearsoon") {
  geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}

# 计算性状与基因的相关性矩阵

## 只有连续型性状才能进行计算，如果是离散变量，在构建样品表时就转为0-1矩阵。

if (corType=="pearsoon") {
  geneTraitCor = as.data.frame(cor(dataExpr, traitData, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(dataExpr, traitData, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}

## Warning in bicor(x, y, use = use, ...): bicor: zero MAD in variable 'y'.
## Pearson correlation was used for individual columns with zero (or missing)
## MAD.

# 最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析
# 确保模块和表型名称一致
module = "blue"  # 注意这里是模块颜色，而不是 "MEyellow"
pheno = "Age"

# 获取模块名称
modNames = substring(colnames(MEs_col), 3)

# 获取对应列索引
module_column = match(module, modNames)
pheno_column = match(pheno, colnames(traitData))

# 检查列索引是否有效
if (is.na(module_column)) stop("Module not found in MEs_col")
if (is.na(pheno_column)) stop("Phenotype not found in traitData")

# 获取模块内的基因
moduleGenes = moduleColors == module

# 检查是否有基因被选中
if (sum(moduleGenes) == 0) stop("No genes found in the specified module")

# 绘制散点图
sizeGrWindow(7, 7)
par(mfrow = c(1, 1))

# 绘制散点图，并为点添加阴影
verboseScatterplot(
  abs(geneModuleMembership[moduleGenes, module_column]),
  abs(geneTraitCor[moduleGenes, pheno_column]),
  xlab = paste("Module Membership in", module, "module"),
  ylab = paste("Gene significance for", pheno),
  main = paste("Module membership vs. gene significance\n"),
  cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = NA  # 不直接绘制点
)

# 添加阴影效果
points(
  abs(geneModuleMembership[moduleGenes, module_column]),
  abs(geneTraitCor[moduleGenes, pheno_column]),
  col = adjustcolor(module, alpha.f = 0.5),  # 调整透明度，模拟阴影
  pch = 16,  # 实心圆点
  cex = 1.5  # 放大点的大小
)

# 添加主点（覆盖阴影）
points(
  abs(geneModuleMembership[moduleGenes, module_column]),
  abs(geneTraitCor[moduleGenes, pheno_column]),
  col = module,  # 使用模块颜色
  pch = 16,  # 实心圆点
  cex = 1.2  # 正常大小
)
save(colData,cyt,data_matrix,dataExpr,dds,dissTOM,log2_data,MEs, moduleGenes,moduleColors,moduleLabels,ordered_samples,probes,geneModuleMembership,exprMat,robustY,net,MEs_col,modTraitCor,plotTOM,prtein_matrix,sft,TOM,vsd,traitData,group,power,file = "计算WGCNA的必要文件.Rdata")

























# 确保 moduleColors 的长度与基因名称一致
# 基因名称是列名，因此使用 colnames(dataExpr)
moduleGenesList <- split(colnames(dataExpr), moduleColors)

# 查看某个模块的基因，例如黄色模块
yellowModuleGenes <- moduleGenesList[["turquoise"]]
print(yellowModuleGenes)

# 查看所有模块及其包含的基因
View(moduleGenesList)
write.csv(yellowModuleGenes,file = "turquoise-FOXO1FOXO3所在模块的所有基因.csv")
# 确保 moduleColors 的长度与基因名称一致
# 基因名称是列名，因此使用 colnames(dataExpr)
moduleGenesList <- split(colnames(dataExpr), moduleColors)

# 查看某个模块的基因，例如黄色模块
yellowModuleGenes <- moduleGenesList[["yellow"]]
print(yellowModuleGenes)

# 查看所有模块及其包含的基因
View(moduleGenesList)
write.csv(yellowModuleGenes,file = "yellow-FOXO1FOXO3所在模块的所有基因.csv")



# 查看某个模块的基因，例如黄色模块
yellowModuleGenes <- moduleGenesList[["brown"]]
print(yellowModuleGenes)

# 查看所有模块及其包含的基因
View(moduleGenesList)
write.csv(yellowModuleGenes,file = "brown-FOXO1FOXO3所在模块的所有基因.csv")


# 查看某个模块的基因，例如黄色模块
yellowModuleGenes <- moduleGenesList[["red"]]
print(yellowModuleGenes)

# 查看所有模块及其包含的基因
View(moduleGenesList)
write.csv(yellowModuleGenes,file = "red-FOXO1FOXO3所在模块的所有基因.csv")

# 查看某个模块的基因，例如黄色模块
yellowModuleGenes <- moduleGenesList[["green"]]
print(yellowModuleGenes)

# 查看所有模块及其包含的基因
View(moduleGenesList)
write.csv(yellowModuleGenes,file = "green-FOXO1FOXO3所在模块的所有基因.csv")



yellowModuleGenes <- moduleGenesList[["blue"]]
print(yellowModuleGenes)

# 查看所有模块及其包含的基因
View(moduleGenesList)
write.csv(yellowModuleGenes,file = "blue-FOXO1FOXO3所在模块的所有基因.csv")


yellowModuleGenes <- moduleGenesList[["grey"]]
print(yellowModuleGenes)

# 查看所有模块及其包含的基因
View(moduleGenesList)
write.csv(yellowModuleGenes,file = "grey-FOXO1FOXO3所在模块的所有基因.csv")

# 定义感兴趣的基因列表
genes_of_interest <- c("FOXO3", "FOXO1", "TLR7", "CD19", "CD38", "CD48", "IL6ST","CR2","TRAF2")

# 确保 moduleColors 的长度与基因名称一致
# 基因名称是列名，因此使用 colnames(dataExpr)
moduleGenesList <- split(colnames(dataExpr), moduleColors)

# 创建一个空的数据框用于存储结果
gene_module_mapping <- data.frame(Gene = character(),
                                  Module = character(),
                                  stringsAsFactors = FALSE)

# 遍历感兴趣的基因，查找它们所属的模块
for (gene in genes_of_interest) {
  for (module in names(moduleGenesList)) {
    if (gene %in% moduleGenesList[[module]]) {
      gene_module_mapping <- rbind(gene_module_mapping, data.frame(Gene = gene, Module = module))
    }
  }
}

# 查看结果
print(gene_module_mapping)

# 保存结果到CSV文件
write.csv(gene_module_mapping, file = "B细胞调节基因gene_module_mapping.csv", row.names = FALSE)






# 确保 moduleColors 的名称是基因名称
names(moduleColors) <- rownames(dataExpr)

# 按模块分组基因
moduleGenesList <- split(names(moduleColors), moduleColors)

# 提取感兴趣模块的基因
modules_of_interest <- c("turquoise","green","red", "blue", "brown", "yellow")  # 替换为你的模块名称
module_genes <- unlist(moduleGenesList[modules_of_interest])

# 检查是否有基因
if (length(module_genes) == 0) {
  stop("感兴趣模块中没有基因，请检查模块名称或数据！")
}
# 提取差异基因的表达数据
module_genes_data <- dataExpr[module_genes, , drop = FALSE]

# 标准化数据
module_genes_data_scaled <- t(scale(t(module_genes_data)))
# 为每个基因添加模块名称
row_annotation <- data.frame(
  Module = moduleColors[module_genes],
  row.names = module_genes
)

# 检查注释是否正确
print(head(row_annotation))
# 检查数据是否正确
print(dim(module_genes_data_scaled))  # 确保有基因和样本
# 加载绘图包
library(pheatmap)

# 绘制热图
pheatmap(
  module_genes_data_scaled,
  cluster_rows = TRUE,  # 行聚类
  cluster_cols = TRUE,  # 列聚类
  annotation_row = row_annotation,  # 添加模块名称注释
  show_rownames = FALSE,  # 隐藏基因名称
  show_colnames = TRUE,  # 显示样本名称
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Standardized Module Differential Genes Expression"
)


genes_to_mark <- c("FOXO3", "FOXO1", "TLR7", "CD38", "CD48", "IL6ST","IRF4","IGHG1","IL7R","IGHD","IGHM","TRAF2","CR2")
row_annotation <- data.frame(
  Module = moduleColors[module_genes],
  Marked = ifelse(module_genes %in% genes_to_mark, "Marked", "Unmarked"),
  row.names = module_genes
)

# 检查注释是否正确
print(head(row_annotation))
# 加载必要的包
library(ComplexHeatmap)
library(circlize)

# 标准化数据
module_genes_data_scaled <- t(scale(t(module_genes_data)))

# 创建模块颜色注释
module_colors <- moduleColors[module_genes]  # 提取模块颜色
row_annotation <- rowAnnotation(
  Module = module_colors,
  col = list(Module = structure(unique(module_colors), names = unique(module_colors))),
  show_annotation_name = FALSE  # 隐藏注释标题
)

# 创建热图对象
heatmap <- Heatmap(
  module_genes_data_scaled,
  name = "Expression",
  show_row_names = FALSE,  # 隐藏基因名称
  show_column_names = TRUE,  # 显示样本名称
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
)

# 添加右侧标注
row_anno_mark <- rowAnnotation(
  Marked = anno_mark(
    at = which(rownames(module_genes_data_scaled) %in% genes_to_mark),
    labels = genes_to_mark,
    side = "right"
  )
)

# 绘制热图
pdf("module基因的热图并标记.pdf",height = 6,width = 10)
draw(row_annotation + heatmap + row_anno_mark, merge_legend = TRUE)
dev.off()
