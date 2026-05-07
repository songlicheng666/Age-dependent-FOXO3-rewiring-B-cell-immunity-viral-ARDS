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

AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- data.frame(t(AUCmatrix@assays@data@listData$AUC), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC) 
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC) 
colnames(AUCmatrix) <- RegulonName_AUC
pbmcauc <- AddMetaData(BC, AUCmatrix) 
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

my.regulons <- intersect(rownames(AUCmatrix),rownames(BINmatrix))
myAUCmatrix <- t(AUCmatrix[rownames(AUCmatrix)%in%my.regulons,])
myBINmatrix <- t(BINmatrix[rownames(BINmatrix)%in%my.regulons,])

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
regulonAUC <-loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "group"])
rssPlot <-plotRSS(rss) 
rssPlot
ggsave('int/rss.pdf', rssPlot$plot, width=4 ,height=5)

ggsave('int/rss.pdf', rssPlot$plot, width=5 ,height=6)

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

BC$group1<-factor(BC$group1,levels = c("CC","AC","CA","AA"))
colnames(group_mean_AUC)[grep("FOXO3", colnames(group_mean_AUC), ignore.case = TRUE)]
FeaturePlot(BC, features = "FOXO3", reduction = "umap")
FeaturePlot(BC, features = "FOXO3", reduction = "umap", split.by = "group1")
FeaturePlot(BC, features = "FOS", reduction = "umap", split.by = "group1")

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

tf_axis <- intersect(tf_axis, colnames(group_mean_AUC))
tf_axis
library(ComplexHeatmap)
library(circlize)

tf_auc_mat <- as.matrix(group_mean_AUC[, tf_axis, drop = FALSE])

tf_auc_mat_t <- t(tf_auc_mat)

col_fun_auc <- colorRamp2(
  c(min(tf_auc_mat_t), (min(tf_auc_mat_t)+max(tf_auc_mat_t))/2, max(tf_auc_mat_t)),
  c("navy", "white", "firebrick")
)

pdf("TF_FOXO3_axis_TF_regulons_meanAUC_heatmap.pdf", width = 6, height = 3)
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

rss_names_for_tf <- rss_names_for_tf[names(rss_names_for_tf) %in% rownames(rss)]

tf_rss_mat <- rss[names(rss_names_for_tf), , drop = FALSE]

col_fun_rss <- colorRamp2(
  c(min(tf_rss_mat), (min(tf_rss_mat)+max(tf_rss_mat))/2, max(tf_rss_mat)),
  c("navy", "white", "firebrick")
)

pdf("FOXO3_axis_TF_regulons_RSS_heatmap.pdf", width = 6, height = 3)
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

lapply(names(rss_names_for_tf), function(reg_name) {
  targets <- regulons[[reg_name]]
  intersect(targets, foxo3_downstream_genes)
})


str(group_mean_AUC)
head(group_mean_AUC)

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

#1. 计算 FOXO3 与下游基因的相关系数（全体 Naive B）
library(Seurat)
library(dplyr)

genes_all <- c("FOXO3","BID","BCL2L11","DDIT3",
               "MAP3K5","FAS","TNFSF10","GADD45A","SOD2","CAT","GPX1","SESN1")
genes_all <- intersect(genes_all, rownames(BC))
genes_all

expr_df <- FetchData(BC, vars = genes_all) 
colnames(expr_df)

foxo3_cor <- sapply(genes_all[genes_all != "FOXO3"], function(g) {
  if(!g %in% colnames(expr_df)) return(NA_real_)
  cor(expr_df$FOXO3, expr_df[[g]], method = "pearson", use = "complete.obs")
})

foxo3_cor

targets <- genes_all[genes_all != "FOXO3"]

edges <- data.frame(
  from  = "FOXO3",
  to    = targets,
  weight = foxo3_cor[targets],
  stringsAsFactors = FALSE
)
edges

gene_mean_expr <- rowMeans(expr_mat)

nodes <- data.frame(
  name  = genes_all,
  role  = c("Center", rep("Downstream", length(genes_all) - 1)),
  mean_expr = gene_mean_expr[genes_all],
  stringsAsFactors = FALSE
)

nodes

Idents(BC) <- "group1"

# 示例：AA vs CC 做差异分析
deg_CA_AA <- FindMarkers(BC, ident.1 = "CA", ident.2 = "AA", features = genes_all,
                         logfc.threshold = 0, min.pct = 0)

deg_CA_AA <- deg_CA_AA %>%
  mutate(gene = rownames(.))

nodes <- nodes %>%
  left_join(deg_CA_AA[, c("gene","avg_log2FC")],
            by = c("name" = "gene"))

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

foxo3_apoptosis_genes  <- c("BCL2L11", "BID", "FAS", "TNFSF10", "DDIT3", "MAP3K5")

foxo3_protective_genes <- c("GADD45A", "SOD2", "CAT", "GPX1", "SESN1")

foxo3_apoptosis_genes  <- intersect(foxo3_apoptosis_genes,  rownames(BC))
foxo3_protective_genes <- intersect(foxo3_protective_genes, rownames(BC))
foxo3_apoptosis_genes
foxo3_protective_genes

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

colnames(BC@meta.data)[grepl("FOXO3_", colnames(BC@meta.data))]

library(ggpubr)

Idents(BC) <- "group1"

BC_sub <- subset(BC, idents = c("CA", "AA"))

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


