# =============================================================================
# Script: 04_SCENIC_transcription_factor_analysis.R
# Description: SCENIC regulon inference pipeline for B-cell transcription factor analysis:
# GENIE3 (co-expression network) -> RcisTarget (motif enrichment) -> AUCell (activity scoring);
# FOXO3 regulon activity comparison between pediatric and adult ARDS
# Author: Licheng Song
# Institution: Chinese PLA General Hospital (Eighth Medical Center)
# Study: Age-dependent B-cell immune dysregulation in virus-associated ARDS
# Year: 2025
# Input: CountsSeuratB_umap_J (Seurat object, Naive_B subset); cisTarget_databases/ (local)
# Output: SCENIC regulon activity matrices, FOXO3_regulon_activity_comparison.pdf
# =============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GENIE3", "AUCell", "RcisTarget"))
# 安装 devtools（如果尚未安装）
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCENIC") 
packageVersion("SCENIC")
getwd()
# 加载必要的库
library(SCENIC)
library(Seurat)
library(dplyr)
table(CountsSeuratB_umap_J$group1)
setwd("/home/lungtissue/ARDSpbmc/新分析/B细胞转录因子分析")
# 从 Seurat 对象中提取表达矩阵
BC <- CountsSeuratB_umap_J[,CountsSeuratB_umap_J$celltype2%in%c("Naive_B")]
#只看CA和AA组的NAIVE B
BC <- BC[,BC$group1%in%c("CA","AA","CC","AC")]
#BC <- BCauc
exprMat <- as.matrix(GetAssayData(object = BC, slot = "data"))

# 提取分组信息（假设分组信息存储在 meta.data 中的 "group" 列）
group_info <- BC@meta.data$group1
names(group_info) <- rownames(BC@meta.data)


data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
# 初始化 SCENIC 设置
scenicOptions <- initializeScenic(org = "hgnc", dbDir = "cisTarget_databases", nCores = 44)

genesKept <- geneFiltering(exprMat, scenicOptions, minCountsPerGene = 3 * 0.01 * ncol(exprMat), minSamples = ncol(exprMat) * 0.01)
head(genesKept)
length(genesKept)
# 设置基因表达矩阵和基因集
exprMat_filtered <- exprMat[rownames(exprMat) %in% genesKept, ]
#exprMat_filtered <- exprMat[rownames(exprMat) %in% genesKept(scenicOptions), ]
# 使用 GENIE3 推断基因调控网络
library(GENIE3)
rm(exprMat)
#weightMatrix <- GENIE3(exprMat_filtered, regulators = rownames(exprMat_filtered), nCores = 20)

#这样会把每个基因都当做调控基因，如果已经知道调控基因，可以放在一个向量里,如：
# Genes that are used as candidate regulators
regulators <- c("SP1","STAT5A","AP1","KLF6","FOXO3","AKT1","AMPK")
#weightMatrix <- GENIE3(exprMat_filtered, regulators=regulators)










#抽样来计算
# 1. 启用20核并行
# 2. 对细胞抽样至2000个
# 3. 调整算法参数
# 从每个分组中抽取相同数量的细胞
set.seed(123)  # 确保可重复性

# 1. 获取每个分组的细胞
groups <- unique(group_info)
cells_by_group <- lapply(groups, function(grp) {
  names(group_info)[group_info == grp]
})
names(cells_by_group) <- groups

# 2. 找出最小分组的大小
min_group_size <- min(sapply(cells_by_group, length))

# 3. 从每个分组中抽取相同数量的细胞（最小分组大小）
sampled_cells <- unlist(lapply(cells_by_group, function(cells) {
  sample(cells, size = min_group_size, replace = FALSE)
}))

# 4. 创建抽样后的表达矩阵
exprMat_sampled <- exprMat_filtered[, sampled_cells]

# 5. 验证抽样结果
table(group_info[sampled_cells])
exprMat_opt <- exprMat_filtered[, sampled_cells]
runGenie3(exprMat_opt, scenicOptions, nParts = 10)
#weightMatrix <- GENIE3(exprMat_opt,
#                       regulators = regulators,
#                       nCores = 20,
#                       treeMethod = "ET",
#                       nTrees = 200)


exprMat_log <- log2(exprMat+1) #log标准话原始矩阵
#scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 1
scenicOptions@settings$seed <- 123
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions) #1. 获取共表达模块
save(scenicOptions,exprMat_opt,exprMat,file = "跑runSCENIC_2_createRegulons之前.Rdata")
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions,coexMethod=c("top5perTarget"))  #2. 获取regulons
#可使用参数coexMethod=c("top5perTarget")，coexMethod可供选择，作者尝试了多种策略过滤低相关性TF-Target,并建议是6种过滤标准都用,其实也并没有耗时多少。默认都计算。
#w001：以每个TF为核心保留weight>0.001的基因形成共表达模块；#w005，#top50
library(doParallel)#top5perTarget：每个基因保留weight值top5的TF得到精简的TF-Target关联表，然后把基因分配给TF构建共表达模块；#top10perTarget；#top50perTarget
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log) #3. 对细胞中的 GRN（调节子）评分


# 检查 scenicOptions 对象
str(scenicOptions)
list.files("int")
regulonAUC <- readRDS("int/3.4_regulonAUC.Rds")
# 提取 AUC 数据
#regulonAUC <- getAUC(scenicOptions)
head(regulonAUC)
str(regulonAUC)
# 加载必要的包
library(AUCell)
library(pheatmap)
library(ComplexHeatmap)
library(SummarizedExperiment)
# 提取 AUC 矩阵
aucMatrix <- as.matrix(assay(regulonAUC, "AUC"))

# 检查 AUC 矩阵的维度
dim(aucMatrix)

# 使用 pheatmap 绘制热图
pheatmap(aucMatrix, scale = "row", clustering_distance_rows = "euclidean", show_rownames = FALSE)

# 使用 ComplexHeatmap 绘制热图
Heatmap(aucMatrix, name = "AUC", show_row_names = FALSE, show_column_names = FALSE)
rownames(aucMatrix)
# 筛选部分数据进行绘图
subsetMatrix <- aucMatrix[1:30, 1:100]
#regulators <- c("SP1","STAT5A","AP1","KLF6","FOXO3","AKT1","AMPK")
pheatmap(subsetMatrix, scale = "row", clustering_distance_rows = "euclidean", show_rownames = TRUE,show_colnames =F)






saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat_log)

nPcs <- c(5,15,50)
scenicOptions@settings$seed <- 123 # same seed for all of them
cellInfo <- data.frame(BC@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=="orig.ident")] <- "sample"
colnames(cellInfo)[which(colnames(cellInfo)=="nFeature_RNA")] <- "nGene"
colnames(cellInfo)[which(colnames(cellInfo)=="nCount_RNA")] <- "nUMI"
colnames(cellInfo)[which(colnames(cellInfo)=="seurat_clusters")] <- "cluster"
colnames(cellInfo)[which(colnames(cellInfo)=="group1")] <- "group"
cellInfo <- cellInfo[,c("sample","nGene","nUMI","cluster","group")]
saveRDS(cellInfo, file="int/cellInfo.Rds")
colVars <- list(group=c("AA"="forestgreen",  "AC"="darkorange",  "CA"="magenta4", 
                           "CC"="hotpink"))
colVars$group <- colVars$group[intersect(names(colVars$group), cellInfo$group)]
saveRDS(colVars, file="int/colVars.Rds")

# Plot as pdf (individual files in int/):
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
# Run t-SNE with different settings:
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
pdf("int/AUC_tsne.pdf",width=10,height=10)
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="group", cex=.5)
dev.off()


##导入原始regulonAUC矩阵，也就是运行完runSCENIC_3_scoreCells后产生的AUC矩阵，可以查看每个GRN在每个细胞中的AUC活性打分
AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- data.frame(t(AUCmatrix@assays@data@listData$AUC), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC) #把(替换成_
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC) #把)去掉，最后如把 KLF3_extended (79g) 替换成KLF3_extended_79g
colnames(AUCmatrix) <- RegulonName_AUC
pbmcauc <- AddMetaData(BC, AUCmatrix) #把AUC矩阵添加到pbmc的metadata信息中
pbmcauc@assays$integrated <- NULL
saveRDS(pbmcauc,'BCauc.rds')

##导入二进制regulonAUC矩阵
BINmatrix <- readRDS("int/4.1_binaryRegulonActivity.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(','_',RegulonName_BIN)
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN
pbmcbin <- AddMetaData(BC, BINmatrix)
pbmcbin@assays$integrated <- NULL
saveRDS(pbmcbin, 'BCbin.rds')
library(ggplot2)
##利用Seurat可视化AUC
dir.create('scenic_seurat')
#FeaturePlot
colnames(AUCmatrix)
colnames(BINmatrix)
GRNs <-intersect(colnames(AUCmatrix),colnames(BINmatrix))
head(GRNs)
for(i in 1:length(GRNs)){
  p1 = FeaturePlot(pbmcauc, features=GRNs[i], label=T, reduction = 'umap')
  p2 = FeaturePlot(pbmcbin, features=GRNs[i], label=T, reduction = 'umap')
  p3 = DimPlot(BC, reduction = 'umap', group.by = "group1", label=T)
  plotc = p1|p2|p3
  ggsave(paste("scenic_seurat/",GRNs[i],".png",sep=""), plotc, width=14 ,height=4)
}















#RidgePlot&VlnPlot ,和单细胞展示基因表达一样，这儿可以利用小提琴图展示GRNs的活性得分
for(i in 1:length(GRNs)){
  p1 = RidgePlot(pbmcauc, features =GRNs[i], group.by="group1") + theme(legend.position='none')
  p2 = VlnPlot(pbmcauc, features =GRNs[i], pt.size = 0, group.by="group1") + theme(legend.position='none')
  plotc = p1 + p2
  ggsave(paste("scenic_seurat/","Ridge-Vln_",GRNs[i],".png",sep=""),plotc, width=10, height=8)
}
library(ComplexHeatmap)
library(pheatmap)
library(circlize)
cellInfo <- readRDS("int/cellInfo.Rds")
group = as.data.frame(subset(cellInfo,select = 'group'))
AUCmatrix <- t(AUCmatrix)
BINmatrix <- t(BINmatrix)
#挑选部分感兴趣的regulons
my.regulons <- intersect(rownames(AUCmatrix),rownames(BINmatrix))
myAUCmatrix <- t(AUCmatrix[rownames(AUCmatrix)%in%my.regulons,])
myBINmatrix <- t(BINmatrix[rownames(BINmatrix)%in%my.regulons,])
# 数据聚合：按组计算转录因子AUC的平均值
group_mean_AUC <- aggregate(myAUCmatrix, by = list(group$group), FUN = mean)
rownames(group_mean_AUC) <- group_mean_AUC[, 1]  # 设置行名为组名
group_mean_AUC <- group_mean_AUC[, -1]  # 删除聚合列

# 绘制热图
ht_opt$raster_temp_image_max_width <- 10000  # 设置临时图像宽度
ht_opt$raster_temp_image_max_height <- 10000  # 设置临时图像高度
#使用regulon原始AUC值绘制热图
pdf("scenic_seurat/AUC_heatmap.pdf",width=14,height=4)
Heatmap(
  group_mean_AUC,
  name = "AUC",
  show_row_names = TRUE,
  show_column_names = T,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  col = colorRamp2(c(min(group_mean_AUC), max(group_mean_AUC)), c("blue", "red"))
)
dev.off()
#使用regulon二进制BIN值绘制热图
group_mean_BIN <- aggregate(myBINmatrix, by = list(group$group), FUN = mean)
rownames(group_mean_BIN) <- group_mean_BIN[, 1]  # 设置行名为组名
group_mean_BIN <- group_mean_BIN[, -1]  # 删除聚合列

pdf("scenic_seurat/BIN_heatmap.pdf",width=14,height=4)
Heatmap(group_mean_BIN,
name = "BIN",
show_row_names = TRUE,
show_column_names = T,
cluster_rows = TRUE,
cluster_columns = TRUE,
col = colorRamp2(c(min(group_mean_BIN), max(group_mean_BIN)), c("blue", "red"))
)
dev.off()

#如果觉得regulons太多了，有一些是没有意义的，可以使用RSS来识别挑选各个细胞类型代表性regulons，进行分析
regulonAUC <-loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "group"])
rssPlot <-plotRSS(rss) #大小 rss评分，颜色 Z-score 
rssPlot
ggsave('int/rss.pdf', rssPlot$plot, width=4 ,height=5)

ggsave('int/rss.pdf', rssPlot$plot, width=5 ,height=6)



save(aucMatrix,BC,BINmatrix,AUCmatrix,cellInfo,colVars,exprMat,exprMat_log,group,myAUCmatrix,myBINmatrix,pbmcauc,pbmcbin,rss,regulonAUC,rssPlot,scenicOptions,GRNs,file = "全部计算完毕后.Rdata")


# 提取FOXO3及其上游调控因子的AUC数据
subset_data <- group_mean_AUC[, c("ETS1_extended_105g","KLF6_16g", "SPIB_45g","SPI1_245g", "IKZF1_extended_10g", "FOS_extended_82g", "NFKB1_extended_33g")]

# 绘制热图
pdf("FOXO3及上游调控因子热图.pdf",width=8,height=4)
Heatmap(subset_data,
        name = "AUC",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        col = colorRamp2(c(min(subset_data), max(subset_data)), c("blue", "red")),
        column_title = "FOXO3及上游调控因子",
        row_title = "样本组（CC, CA, AC, AA）")
dev.off()

#AUC 高 + RSS 高：这个 regulon 在该 group 中既强又特异，很可能是这个群的关键调控因子
#AUC 高但 RSS 不高：在很多不同群里都比较高，是“普遍活跃”的调控网络，不是某一群特异的

library(ggplot2)
rownames(rss)
colnames(rss)
# 提取RSS评分数据
rss_data <- data.frame(
  Group = colnames(rss),
  RSS = rss["KLF6 (16g)", ]
)

# 绘制RSS评分图
ggplot(rss_data, aes(x = Group, y = RSS)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  labs(title = "转录因子的RSS评分",
       x = "调控因子",
       y = "RSS评分") +
  coord_flip()

library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

# rss: 行 = regulon, 列 = group (CC, CA, AC, AA)
dim(rss)
rownames(rss)[1:5]
colnames(rss)

# 每个 group 选 Top10 regulon
topN <- 10

top_regulons_list <- lapply(colnames(rss), function(grp) {
  # 对该列从大到小排序，选前 topN 个 regulon
  ord <- order(rss[, grp], decreasing = TRUE, na.last = NA)
  regs <- rownames(rss)[ord][1:min(topN, length(ord))]
  data.frame(
    group   = grp,
    regulon = regs,
    RSS     = rss[regs, grp],
    stringsAsFactors = FALSE
  )
})

top_regulons_df <- do.call(rbind, top_regulons_list)

# 看一下结果
head(top_regulons_df)
table(top_regulons_df$group)
# 所有入选的 regulon 去重
top_regulons_unique <- unique(top_regulons_df$regulon)
length(top_regulons_unique)

# 从 rss 中取出这些 regulon 在所有 group 的 RSS
rss_top_mat <- rss[top_regulons_unique, colnames(rss), drop = FALSE]

# 可选：按行最大 RSS 排序（让热图更好看）
max_order <- order(apply(rss_top_mat, 1, max), decreasing = TRUE)
rss_top_mat <- rss_top_mat[max_order, ]
# 颜色映射：蓝-白-红（低到高）
rss_min <- min(rss_top_mat, na.rm = TRUE)
rss_max <- max(rss_top_mat, na.rm = TRUE)

col_fun <- colorRamp2(
  c(rss_min, (rss_min + rss_max) / 2, rss_max),
  c("navy", "white", "firebrick")
)

pdf("RSS_top10_regulons_heatmap.pdf", width = 8, height = 10)

Heatmap(
  rss_top_mat,
  name = "RSS",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,  # group 顺序按 CC, CA, AC, AA
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_title = "Per-group top regulons (by RSS)",
  row_title = "Regulons"
)

dev.off()


library(ggplot2)
# regulonAUC 里的 regulon 名称
rownames(rss)[grep("FOXO3", rownames(rss), ignore.case = TRUE)]
#策略调整：把分析改成“FOXO3 轴（axis）”
#可以把研究构建成三层：
#
#FOXO3 基因表达在四组中的变化
#与 FOXO3 紧密相关 / 上游调控的 TF regulons（你已有的那 94 个）
#FOXO3 相关下游凋亡/应激基因的表达 + 这些基因的调控 regulons
BC$group1<-factor(BC$group1,levels = c("CC","AC","CA","AA"))
colnames(group_mean_AUC)[grep("FOXO3", colnames(group_mean_AUC), ignore.case = TRUE)]
FeaturePlot(BC, features = "FOXO3", reduction = "umap")
FeaturePlot(BC, features = "FOXO3", reduction = "umap", split.by = "group1")
FeaturePlot(BC, features = "FOS", reduction = "umap", split.by = "group1")
#2. 4 组间小提琴图 + 统计
#复制
library(ggpubr)

Idents(BC) <- "group1"  # CC/CA/AC/AA

VlnPlot(BC, features = "FOXO3", pt.size = 0.1, group.by = "group1") +
  stat_compare_means(
    comparisons = list(
      c("CC","CA"),
      c("CC","AC"),
      c("CC","AA"),
      c("CA","AC"),
      c("CA","AA"),
      c("AC","AA")
    ),
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE
  ) +
  theme_minimal() +
  labs(title = "FOXO3 expression in Naive B across groups",
       x = "Group", y = "FOXO3 (normalized expression)")
#这部分是基因表达层面，不是 regulon。


#四、围绕 FOXO3 选择“上游 TF regulons”，重点看它们的 AUC & RSS
#既然没有 FOXO3 regulon，可以：

#用文献/KEGG/STRING 挑出 “可能调控 FOXO3 的 TF”：
#比如你前面提过的：KLF6, SP1, ETS1, NFKB1, STAT1/3, CEBPB, 等
#在你已有的 regulon 列表中，找到这些 TF 的 regulon 名
#（你 group_mean_AUC 的列名里已经有）：
#例如，基于你贴出的列名，可以先选一批：
tf_axis <- c(
  "KLF6_16g",
  "KLF6_extended_69g",
  "SP1_extended_28g",
  "ETS1_extended_105g",
  "NFKB1_extended_33g",
  "STAT1_extended_50g",
  "STAT3_extended_15g",
  "CEBPB_extended_149g",
  "SPI1_245g"
)

# 保证它们确实存在于 group_mean_AUC
tf_axis <- intersect(tf_axis, colnames(group_mean_AUC))
tf_axis
library(ComplexHeatmap)
library(circlize)

tf_auc_mat <- as.matrix(group_mean_AUC[, tf_axis, drop = FALSE])

# 行是 group，列是 TF regulon
# 为方便阅读，可以转置成 regulon 在行
tf_auc_mat_t <- t(tf_auc_mat)

col_fun_auc <- colorRamp2(
  c(min(tf_auc_mat_t), (min(tf_auc_mat_t)+max(tf_auc_mat_t))/2, max(tf_auc_mat_t)),
  c("navy", "white", "firebrick")
)

pdf("可能调控 FOXO3 的 TF_FOXO3_axis_TF_regulons_meanAUC_heatmap.pdf", width = 6, height = 3)
Heatmap(
  tf_auc_mat_t,
  name = "Mean AUC",
  col = col_fun_auc,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 7),
  column_title = "Groups (CC, CA, AC, AA)",
  row_title = "Candidate upstream TF regulons"
)
dev.off()

##2. 对这些 TF regulons 各自看 RSS（特异性），挑出“在哪一组最特点”
#复制
# rss: 行 = regulon, 列 = group
tf_rss_mat <- rss[ gsub("_", " ", gsub("_([0-9]+g)$", " (\\1)", tf_axis)), , drop = FALSE ]
# 上面这行是从 AUC 列名（KLF6_16g）转换回 RSS 的行名（KLF6 (16g)）
# 如果不想这么复杂，可以手动写一个命名对应表。

# 简化方式：先看所有 regulon 名，手动对应：
rownames(rss)[1:20]
# 手工匹配后，例如：
rss_names_for_tf <- c(
  "KLF6 (16g)"            = "KLF6_16g",
  "KLF6_extended (69g)"   = "KLF6_extended_69g",
  "SP1_extended (28g)"    = "SP1_extended_28g",
  "ETS1_extended (105g)"  = "ETS1_extended_105g",
  "NFKB1_extended (33g)"  = "NFKB1_extended_33g",
  "STAT1_extended (50g)"  = "STAT1_extended_50g",
  "STAT3_extended (15g)"  = "STAT3_extended_15g",
  "CEBPB_extended (149g)" = "CEBPB_extended_149g",
  "SPI1 (245g)"           = "SPI1_245g"
)

# 只保留确实存在的行
rss_names_for_tf <- rss_names_for_tf[names(rss_names_for_tf) %in% rownames(rss)]

tf_rss_mat <- rss[names(rss_names_for_tf), , drop = FALSE]

col_fun_rss <- colorRamp2(
  c(min(tf_rss_mat), (min(tf_rss_mat)+max(tf_rss_mat))/2, max(tf_rss_mat)),
  c("navy", "white", "firebrick")
)

pdf("rss层面看调控FOXO3_axis_TF_regulons_RSS_heatmap.pdf", width = 6, height = 3)
Heatmap(
  tf_rss_mat,
  name = "RSS",
  col = col_fun_rss,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 7),
  column_title = "Groups",
  row_title = "Candidate upstream TF regulons"
)
dev.off()

#这张图回答的问题是：

#“在四组中，哪些 TF regulons 对特定组（如 AA / AC）的特异性最高”，

#可以结合 FOXO3 表达高/低的组，推测“上游驱动”。



#五、在细胞层面看这些 TF regulon 与 FOXO3 表达的关系
#你已经有 pbmcauc（Naive B 的 Seurat 对象 + regulon AUC meta.data）。

#1. FeaturePlot：比对 FOXO3 expression 和某个 TF regulon AUC
#例如看SPI1_245g：
FeaturePlot(pbmcauc, features = c("FOXO3", "KLF6_16g"),split.by = "group1",
            reduction = "umap", ncol = 2)


Idents(pbmcauc) <- "group1"

VlnPlot(pbmcauc, features = c("KLF6_16g","FOXO3"),
        pt.size = 0.05, group.by = "group1", ncol = 2)
#虽然没有 FOXO3 regulon，但你可以从文献选定一批 FOXO3 经典靶基因/通路基因，比如：
foxo3_downstream_genes <- c("BID","DDIT3","MAP3K5","BCL2L11","GADD45A","CDKN1B","SP1", "STAT5A", "AKT1")
#1. 在 Naive B 中直接看这些基因的表达（VlnPlot / DotPlot）
pdf("FOXO3调控的下游的基因的变化.pdf",width = 10,height = 6)
VlnPlot(BC,
        features = foxo3_downstream_genes,
        group.by = "group1",
        pt.size = 0.01, ncol = 3)
dev.off()
#看这些基因是否落入某些 TF regulons 的 target 中
regulons <- loadInt(scenicOptions, "regulons")  # 名 = "TF (Ng)" ；值 = target gene 向量

# 看每个候选 TF regulon 里，是否包含 FOXO3 下游基因
lapply(names(rss_names_for_tf), function(reg_name) {
  targets <- regulons[[reg_name]]
  intersect(targets, foxo3_downstream_genes)
})









str(group_mean_AUC)
head(group_mean_AUC)
# 提取FOXO3的表达数据
foxo3_data <- data.frame(
  Group = rownames(group_mean_AUC),  # 提取组别
  FOXO3_AUC = group_mean_AUC[, "KLF6_16g"]  # 提取 KLF6_16g 列的数据
)

head(foxo3_data)
# 绘制箱线图
ggplot(foxo3_data, aes(x = Group, y = FOXO3_AUC, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "FOXO3在四组中的差异表达",
       x = "组别",
       y = "FOXO3 AUC") +
  scale_fill_manual(values = c("CC" = "red", "CA" = "blue", "AC" = "green", "AA" = "purple"))








library(igraph)
library(igraph)

nodes <- data.frame(
  name  = c("FOXO3", "BID","BCL2L11","DDIT3","GADD45A",
            "MAP3K5","FAS","TNFSF10","SOD2","CAT","GPX1","SESN1"),
  group = c("Target", "Regulator", "Regulator", "Regulator", "Regulator", "Regulator", "Regulator", "Regulator", "Regulator", "Regulator", "Regulator", "Regulator"),
  stringsAsFactors = FALSE
)

edges <- data.frame(
  from = c("FOXO3", "FOXO3", "FOXO3", "FOXO3","FOXO3", "FOXO3", "FOXO3",  "FOXO3", "FOXO3", "FOXO3", "FOXO3"),
  to   = c( "BID","BCL2L11","DDIT3","GADD45A",
            "MAP3K5","FAS","TNFSF10","SOD2","CAT","GPX1","SESN1"),
  stringsAsFactors = FALSE
)

g <- graph_from_data_frame(edges, vertices = nodes, directed = TRUE)

vertex_cols <- ifelse(nodes$group == "Regulator", "yellow", "red")

plot(g,
     vertex.size  = 30,
     vertex.color = vertex_cols,
     edge.arrow.size = 0.5,
     main = "FOXO3 调控轴示意图")


#可以，把这张“示意图”升级成一个带权重、带表达信息、带分组信息的“半定量网络图”。思路是：

#用 真实数据 给边加“调控强度”（相关系数）
#用 表达/活性 给节点加“大小/颜色”
#再用 group（AA/AC/CA/CC） 做分面或颜色编码，让图既有生物意义，又不至于假装是 SCENIC 算出来的网络。
#下面给你一套可以直接改用的方案。
#你现在有：

#Seurat 对象 BC（Naive B，含 4 个 group）
#基因表达矩阵（FOXO3、BID、DDIT3…）
#regulon AUC 也有，但 FOXO3 没 regulon，所以这里先用基因表达做相关性。
#1. 计算 FOXO3 与下游基因的相关系数（全体 Naive B）
library(Seurat)
library(dplyr)
#你现在的核心科学问题其实是：

#“为什么儿童 ARDS 凋亡低，但 FOXO3 反而高？

#FOXO3 在儿童到底更偏向 ‘保护/适应’ 还是 ‘促凋亡’？”

#那对应就需要把 FOXO3 的靶基因拆成两个功能维度，而你现在这组基因里：

#促凋亡 / 应激：BID, BCL2L11, DDIT3, GADD45A, MAP3K5
#细胞周期/调控：CDKN1B
#信号/上游/共调控：AKT1, SP1, STAT5A（其实很多不是典型“FOXO3 下游”）
#这些基因混得很杂，且保护性/抗氧化基因几乎没选到，不利于讲你现在这个“保护 vs 凋亡”的故事。


# 确保这些基因在对象里确实存在
genes_all <- c("FOXO3","BID","BCL2L11","DDIT3",
               "MAP3K5","FAS","TNFSF10","GADD45A","SOD2","CAT","GPX1","SESN1")
genes_all <- intersect(genes_all, rownames(BC))
genes_all


# 此时行 = 基因, 列 = 细胞
expr_df <- FetchData(BC, vars = genes_all)  # 行 = 细胞, 列 = 基因
colnames(expr_df)

foxo3_cor <- sapply(genes_all[genes_all != "FOXO3"], function(g) {
  if(!g %in% colnames(expr_df)) return(NA_real_)
  cor(expr_df$FOXO3, expr_df[[g]], method = "pearson", use = "complete.obs")
})

foxo3_cor

# 构建 edges 数据框（带权重）
targets <- genes_all[genes_all != "FOXO3"]

edges <- data.frame(
  from  = "FOXO3",
  to    = targets,
  weight = foxo3_cor[targets],
  stringsAsFactors = FALSE
)
edges
# 计算基因在所有 Naive B 细胞中的平均表达
gene_mean_expr <- rowMeans(expr_mat)

nodes <- data.frame(
  name  = genes_all,
  role  = c("Center", rep("Downstream", length(genes_all) - 1)),
  mean_expr = gene_mean_expr[genes_all],
  stringsAsFactors = FALSE
)

nodes
#假设你已经有 Naive B 内部的差异分析结果（比如 AA vs CC），如果没有，下面是一个简单示例（注意：真实使用时请根据你之前 DEG 分析结果替换）：
Idents(BC) <- "group1"

# 示例：AA vs CC 做差异分析
deg_CA_AA <- FindMarkers(BC, ident.1 = "CA", ident.2 = "AA", features = genes_all,
                         logfc.threshold = 0, min.pct = 0)

deg_CA_AA <- deg_CA_AA %>%
  mutate(gene = rownames(.))

# 把 log2FC 映射到节点
nodes <- nodes %>%
  left_join(deg_CA_AA[, c("gene","avg_log2FC")],
            by = c("name" = "gene"))

# 缺失的设为 0
nodes$avg_log2FC[is.na(nodes$avg_log2FC)] <- 0
library(igraph)

g <- graph_from_data_frame(edges, vertices = nodes, directed = TRUE)

# Node size: mean expression
nodes$mean_expr[is.na(nodes$mean_expr)] <- 0
vsize <- scales::rescale(nodes$mean_expr, to = c(15, 40))

# Edge width: |correlation|
edge_weight <- E(g)$weight
edge_lwd <- scales::rescale(abs(edge_weight), to = c(1, 6))

# Edge color: positive = red, negative = green
edge_cols <- ifelse(edge_weight >= 0, "red", "green")

# Node color: consistent with edges, log2FC > 0 = red, <= 0 = green
logfc <- nodes$avg_log2FC
node_cols <- ifelse(logfc > 0, "red", "green")

set.seed(123)
pdf("FOXO3调控网络.pdf",height = 9,width = 9)
plot(
  g,
  vertex.size      = vsize,
  vertex.label.cex = 0.8,
  vertex.color     = node_cols,
  edge.width       = edge_lwd,
  edge.color       = edge_cols,
  edge.arrow.size  = 0.6,
  layout           = layout_in_circle(g, order = V(g)$name),
  main             = "FOXO3 regulatory axis (correlation & expression)"
)

# Legend
legend(
  "topleft",
  legend = c("log2FC > 0 / positive", "log2FC ≤ 0 / negative"),
  col    = c("red", "green"),
  pch    = 16,
  pt.cex = 1.5,
  bty    = "n"
)
legend(
  "bottomleft",
  legend = c("Edge width = |correlation|"),
  lwd    = 3,
  col    = "grey20",
  bty    = "n"
)
dev.off()
#你可以、也很值得“从 FOXO3 角度做区别分析”，但前提是把它定位为“一个 stress‑adaptation TF 轴中的关键节点”，而不是简化成“凋亡开关”。
#一、FOXO3_apoptosis_score vs FOXO3_protective_score（基因表达层面）
#1. 定义两个 FOXO3 下游基因集
#放在你已经定义 genes_all 那一段附近即可：
# FOXO3 促凋亡 / 应激模块（根据文献 + 你当前数据）
foxo3_apoptosis_genes  <- c("BCL2L11", "BID", "FAS", "TNFSF10", "DDIT3", "MAP3K5")

# FOXO3 保护 / 应激适应模块
foxo3_protective_genes <- c("GADD45A", "SOD2", "CAT", "GPX1", "SESN1")

# 只保留在 BC 中存在的基因
foxo3_apoptosis_genes  <- intersect(foxo3_apoptosis_genes,  rownames(BC))
foxo3_protective_genes <- intersect(foxo3_protective_genes, rownames(BC))
foxo3_apoptosis_genes
foxo3_protective_genes
#2. 计算每个细胞的模块 score（用 Seurat 的 AddModuleScore）
# 使用归一化表达矩阵（RNA@data）
DefaultAssay(BC) <- "RNA"

BC <- AddModuleScore(
  object = BC,
  features = list(foxo3_apoptosis_genes),
  name = "FOXO3_apoptosis"
)

BC <- AddModuleScore(
  object = BC,
  features = list(foxo3_protective_genes),
  name = "FOXO3_protective"
)

# Seurat 会自动生成列名 like "FOXO3_apoptosis1"、"FOXO3_protective1"
colnames(BC@meta.data)[grepl("FOXO3_", colnames(BC@meta.data))]
#3. 在 CA vs AA 比较这两个 score + 画图
library(ggpubr)

Idents(BC) <- "group1"
# 只看 CA 和 AA
BC_sub <- subset(BC, idents = c("CA", "AA"))

# 促凋亡模块
p_apop <- VlnPlot(
  BC_sub,
  features = "FOXO3_apoptosis1",
  group.by = "group1",
  pt.size = 0
) +
  stat_compare_means(
    comparisons = list(c("CA", "AA")),
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() +
  labs(title = "FOXO3-apoptosis module score (CA vs AA)",
       x = "Group", y = "Module score")

# 保护模块
p_prot <- VlnPlot(
  BC_sub,
  features = "FOXO3_protective1",
  group.by = "group1",
  pt.size = 0
) +
  stat_compare_means(
    comparisons = list(c("CA", "AA")),
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal() +
  labs(title = "FOXO3-protective module score (CA vs AA)",
       x = "Group", y = "Module score")

pdf("FOXO3_apoptosis_vs_protective_scores_CA_AA.pdf", width = 8, height = 4)
p_apop + p_prot
dev.off()
#如果你想看四组（CC/AC/CA/AA）的整体趋势，也可以再来一个四组版本：
p_apop_all <- VlnPlot(
  BC,
  features = "FOXO3_apoptosis1",
  group.by = "group1",
  pt.size = 0
) + theme_minimal() +
  labs(title = "FOXO3-apoptosis module score (all groups)",
       x = "Group", y = "Module score")

p_prot_all <- VlnPlot(
  BC,
  features = "FOXO3_protective1",
  group.by = "group1",
  pt.size = 0
) + theme_minimal() +
  labs(title = "FOXO3-protective module score (all groups)",
       x = "Group", y = "Module score")

pdf("FOXO3_apoptosis_vs_protective_scores_all_groups.pdf", width = 8, height = 4)
p_apop_all + p_prot_all
dev.off()
#你已经把 AUC/BIN 加到 Seurat（pbmcauc, pbmcbin）并且 group1/cluster 信息也在 meta.data 里，我们只需要聚焦你关心的 TF 轴。

#1. 定义 FOXO3-axis TF regulons，并从 AUC/BIN 矩阵取出
#你前面已经写了一版，这里稍微整理为“轴”：
# FOXO3-axis TF regulons（根据你前面的列名）
foxo3_axis_regs <- c(
  "KLF6_16g",
  "KLF6_extended_69g",
  "ETS1_extended_105g",
  "SPIB_45g",
  "SPI1_245g",
  "NFKB1_extended_33g",
  "FOS_extended_82g",
  "STAT1_extended_50g",
  "STAT3_extended_15g",
  "CEBPB_extended_149g"
)

foxo3_axis_regs <- intersect(foxo3_axis_regs, colnames(AUCmatrix))
foxo3_axis_regs
#2. 画 UMAP：这些 regulon 的 AUC + FOXO3 表达
# 确保 pbmcauc 有 group1、seurat_clusters 等
Idents(pbmcauc) <- "group1"

# FOXO3 表达 + 一个代表性 TF regulon（比如 KLF6_16g）
pdf("FOXO3_axis_FOXO3_and_KLF6_UMAP_by_group.pdf", width = 10, height = 6)
FeaturePlot(pbmcauc,
            features = c("FOXO3", "KLF6_16g"),
            split.by = "group1",
            reduction = "umap", ncol = 2)
dev.off()

# 所有 FOXO3-axis regulons 的 AUC UMAP（按 group 分面）
pdf("FOXO3_axis_TF_regulons_AUC_UMAP_by_group.pdf", width = 12, height = 3 * ceiling(length(foxo3_axis_regs)/3))
for(f in foxo3_axis_regs){
  print(
    FeaturePlot(pbmcauc,
                features = f,
                split.by = "group1",
                reduction = "umap") +
      ggtitle(paste(f, "AUC"))
  )
}
dev.off()
#3. 看这些 regulons 在 cluster（B 细胞 cluster）中的平均活性 + CA 富集
#如果你在 BC 中已有 seurat_clusters，我们可以对 regulon AUC 按 cluster & group 平均：
# 从 pbmcauc@meta.data 中取 AUC + cluster + group
meta_auc <- pbmcauc@meta.data[, c("group1", "seurat_clusters", foxo3_axis_regs)]
colnames(meta_auc)[1:2] <- c("group", "cluster")

# 按 cluster 计算平均 AUC（所有细胞），再看在哪个 group 富集
library(dplyr)

cluster_mean_auc <- meta_auc %>%
  group_by(cluster) %>%
  summarise(across(all_of(foxo3_axis_regs), mean, na.rm = TRUE))

cluster_mean_auc
#如果你想可视化“哪个 cluster FOXO3-high / regulon-high，并且在 CA 中富集”，可以做：
# 先找 FOXO3-high 的细胞（表达层面）
DefaultAssay(BC) <- "RNA"
foxo3_expr <- FetchData(BC, vars = "FOXO3")
BC$FOXO3_high <- ifelse(foxo3_expr$FOXO3 > quantile(foxo3_expr$FOXO3, 0.75), "high", "low")

pdf("FOXO3_high_cells_UMAP_by_group.pdf", width = 6, height = 5)
DimPlot(BC, reduction = "umap", group.by = "FOXO3_high", split.by = "group1")
dev.off()

# 再结合一个轴 TF（如 KLF6_16g）：
auc_k <- FetchData(pbmcauc, vars = "KLF6_16g")
pbmcauc$KLF6_high <- ifelse(auc_k$KLF6_16g > quantile(auc_k$KLF6_16g, 0.75), "high", "low")

pdf("FOXO3_KLF6_high_cells_UMAP_by_group.pdf", width = 10, height = 5)
DimPlot(pbmcauc, reduction = "umap", group.by = "KLF6_high", split.by = "group1")
dev.off()
#三、把 FOXO3 放进更大的 TF 轴模型（AUC + RSS）
#1. 你已有：可能调控 FOXO3 的 TF_FOXO3_axis_TF_regulons_meanAUC_heatmap.pdf
#这里再补一张“同一批 TF regulons 的 RSS 热图 + FOXO3 表达的箱线图”。

#1）RSS 热图（你已经有一版，稍微加个标题说明）
# rss_names_for_tf 已经在你代码里定义过
# 这里再画一次并带中文标题
tf_rss_mat <- rss[names(rss_names_for_tf), , drop = FALSE]

col_fun_rss <- colorRamp2(
  c(min(tf_rss_mat), (min(tf_rss_mat)+max(tf_rss_mat))/2, max(tf_rss_mat)),
  c("navy", "white", "firebrick")
)

pdf("FOXO3_axis_TF_regulons_RSS_heatmap_中文标注版.pdf", width = 6, height = 3)
Heatmap(
  tf_rss_mat,
  name = "RSS",
  col = col_fun_rss,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 7),
  column_title = "Groups (CC, CA, AC, AA)",
  row_title = "Candidate FOXO3-axis TF regulons"
)
dev.off()
#2）FOXO3 表达 + 轴上其中 1–2 个 TF regulon（如 KLF6_16g）的 AUC 小提琴图
library(ggpubr)

Idents(pbmcauc) <- "group1"

p_foxo3_expr <- VlnPlot(
  BC,
  features = "FOXO3",
  group.by = "group1",
  pt.size = 0.1
) +
  theme_minimal() +
  labs(title = "FOXO3 expression across groups",
       x = "Group", y = "FOXO3 (normalized)")

p_klf6_auc <- VlnPlot(
  pbmcauc,
  features = "KLF6_16g",
  group.by = "group1",
  pt.size = 0.1
) +
  theme_minimal() +
  labs(title = "KLF6 regulon AUC across groups",
       x = "Group", y = "AUC score")

pdf("FOXO3_and_KLF6_axis_expression_AUC_vln.pdf", width = 8, height = 4)
p_foxo3_expr + p_klf6_auc
dev.off()
save(aucMatrix, BC, BINmatrix, AUCmatrix,
     cellInfo, colVars, exprMat, exprMat_log, group,
     myAUCmatrix, myBINmatrix, pbmcauc, pbmcbin,
     rss, regulonAUC, rssPlot, scenicOptions, GRNs,
     foxo3_apoptosis_genes, foxo3_protective_genes,
     file = "全部计算完毕后.Rdata")







library(ggplot2)
library(dplyr)

# 1. 取数据
df_cor <- FetchData(pbmcauc, vars = c("FOXO3", "KLF6_16g", "group1"))
df_cor$group1 <- factor(df_cor$group1, levels = c("CC","AC","CA","AA"))

# 2. 每组筛选：KLF6_16g、FOXO3 都 > 0
df_cor_filt <- df_cor %>%
  dplyr::group_by(group1) %>%
  dplyr::filter(KLF6_16g > 0 & FOXO3 > 0) %>%
  dplyr::ungroup()
write.csv(df_cor_filt,file="df_cor_filt.csv")
df_cor_filt <- read.csv(file="df_cor_filt.csv")
# 3. 分组算相关
cor_by_group <- df_cor_filt %>%
  dplyr::group_by(group1) %>%
  dplyr::summarise(
    n = length(FOXO3),
    r = cor(KLF6_16g, FOXO3, method = "spearman"),
    .groups = "drop"
  )
cor_by_group

# 4. 分面散点 + 阴影回归线 + 轴线
p <- ggplot(df_cor_filt, aes(x = KLF6_16g, y = FOXO3)) +
  geom_point(alpha = 0.5, size = 0.6, color = "grey30") +
  geom_smooth(method = "lm", se = TRUE,
              color = "steelblue", fill = "lightblue") +
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = 0, color = "black") +
  facet_wrap(~ group1, nrow = 2) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  labs(
    title = "Correlation between KLF6 regulon activity and FOXO3 expression\n(cells with both > 0)",
    x = "KLF6 regulon AUC",
    y = "FOXO3 expression"
  )

p
ggsave("FOXO3_KLF6_correlation_by_group_filtered_nonzero.pdf",
       p, width = 6, height = 5)
