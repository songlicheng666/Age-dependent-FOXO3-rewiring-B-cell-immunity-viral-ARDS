# =============================================================================
# Script: 02_spatial_Bcell_niche_analysis.R
# Description: B-cell spot identification (NNLS/RCTD weights or marker-based scoring),
# FOXO3 module scoring, region-level enrichment, spatial correlation with
# apoptosis/inflammation scores, LIVE vs DEATH outcome comparison
# Author: Licheng Song
# Institution: Chinese PLA General Hospital (Eighth Medical Center)
# Study: Age-dependent B-cell immune dysregulation in virus-associated ARDS
# Year: 2025
# Input: Spatial_integrated / Spatial_integrated_Bcell_region_analysis (Seurat Visium object)
# Output: B_cell_spatial_distribution.pdf, FOXO3_niche_correlation.pdf
# =============================================================================

# 前置依赖与参数
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(Matrix)
library(presto)          # 快速差异分析（可选）
library(fgsea)           # 富集分析（可选）
library(msigdbr)         # MSigDB通路（可选）
library(reshape2)

# 设定默认assay为表达层用于打分（优先RNA或SCT，其次Spatial）
if("RNA" %in% names(Spatial_integrated@assays)) {
  DefaultAssay(Spatial_integrated) <- "RNA"
} else if("SCT" %in% names(Spatial_integrated@assays)) {
  DefaultAssay(Spatial_integrated) <- "SCT"
} else {
  DefaultAssay(Spatial_integrated) <- "Spatial"
}
Spatial_integrated <- Spatial_integrated_Bcell_region_analysis
# 1 定位B细胞相关spot
# 如果已经有celltypeprops assay（NNLS或RCTD权重），挑选B细胞权重较高的spot；否则用标记基因近似筛选
DefaultAssay(Spatial_integrated) <- "celltypeprops"
available_props <- tryCatch({rownames(Spatial_integrated[["celltypeprops"]]@data)}, error=function(e) NULL)
b_prop_name <- NULL
if(!is.null(available_props)) {
  # 在celltypeprops中搜索可能的B细胞条目名
  candidates <- grep("B|Bcell|B_cell|B.cells|Bcells", available_props, value = TRUE, ignore.case = TRUE)
  if(length(candidates) > 0) {
    b_prop_name <- candidates[1]
    message(paste("检测到B细胞权重特征:", b_prop_name))
  } else {
    message("未在celltypeprops中检测到B细胞权重特征 将采用标记基因筛选法")
  }
} else {
  message("celltypeprops assay不可用 将采用标记基因筛选法")
}

# 切回表达assay用于后续打分
if("RNA" %in% names(Spatial_integrated@assays)) {
  DefaultAssay(Spatial_integrated) <- "RNA"
} else if("SCT" %in% names(Spatial_integrated@assays)) {
  DefaultAssay(Spatial_integrated) <- "SCT"
} else {
  DefaultAssay(Spatial_integrated) <- "Spatial"
}

# B细胞标记集合（可根据物种与数据调整）
b_markers_core <- c("MS4A1","CD79A","CD79B","BANK1","CD74","MZB1")         # B细胞与浆细胞核心
b_markers_naive <- c("MS4A1","CD19","IGHM","IGHD","TNFRSF13C")             # 初始/成熟B
b_markers_gc <- c("BCL6","AICDA","LMO2","NEK6","POU2AF1")                  # 生发中心
b_markers_plasma <- c("MZB1","XBP1","JCHAIN","SDC1","PRDM1")               # 浆细胞
b_markers_activation <- c("CD86","CD40","HLA-DRA","HLA-DPA1","HLA-DPB1")   # 激活与抗原呈递
b_markers_isotype <- c("IGHG1","IGHG3","IGHA1","IGHA2","IGHE")             # 类转换
b_markers_chemokine <- c("CXCR4","CXCR5","CCR7")                           # 归巢与定位
b_markers_prolif <- c("MKI67","TOP2A","PCNA")                              # 增殖
b_markers_regulatory <- c("IL10","TGFB1","PDCD1LG2")                       # 免疫调节

# 确保基因存在
gene_exists <- function(genes, obj) {
  genes[genes %in% rownames(obj)]
}
b_sets <- list(
  B_core = b_markers_core,
  B_naive = b_markers_naive,
  B_GC = b_markers_gc,
  B_plasma = b_markers_plasma,
  B_activation = b_markers_activation,
  B_isotype = b_markers_isotype,
  B_chemokine = b_markers_chemokine,
  B_proliferation = b_markers_prolif,
  B_regulatory = b_markers_regulatory
)
b_sets <- lapply(b_sets, gene_exists, obj = Spatial_integrated)

# 2 定义B细胞spot集合
if(!is.null(b_prop_name) && "celltypeprops" %in% names(Spatial_integrated@assays)) {
  DefaultAssay(Spatial_integrated) <- "celltypeprops"
  b_weight <- FetchData(Spatial_integrated, vars = b_prop_name)[,1]
  Spatial_integrated$B_weight <- b_weight
  # 阈值可调 例如top 30%或固定阈值
  thr <- quantile(b_weight, 0.7, na.rm = TRUE)
  Spatial_integrated$B_enriched <- b_weight >= thr
  message(paste("按NNLS/RCTD权重定义B_enriched 阈值:", round(thr,3)))
} else {
  # 用核心标记module score近似定义B_enriched
  DefaultAssay(Spatial_integrated) <- if("RNA" %in% names(Spatial_integrated@assays)) "RNA" else if("SCT" %in% names(Spatial_integrated@assays)) "SCT" else "Spatial"
  b_core_genes <- b_sets$B_core
  Spatial_integrated <- AddModuleScore(Spatial_integrated, features = list(b_core_genes), name = "BcoreScore", assay = DefaultAssay(Spatial_integrated))
  # AddModuleScore生成 BcoreScore1 列
  core_thr <- quantile(Spatial_integrated$BcoreScore1, 0.7, na.rm = TRUE)
  Spatial_integrated$B_enriched <- Spatial_integrated$BcoreScore1 >= core_thr
  Spatial_integrated$B_weight <- Spatial_integrated$BcoreScore1
  message(paste("按标记打分定义B_enriched 阈值:", round(core_thr,3)))
}

# 恢复表达assay
DefaultAssay(Spatial_integrated) <- if("RNA" %in% names(Spatial_integrated@assays)) "RNA" else if("SCT" %in% names(Spatial_integrated@assays)) "SCT" else "Spatial"

# 3 计算功能module scores
for(set_name in names(b_sets)) {
  genes <- b_sets[[set_name]]
  if(length(genes) >= 3) {  # 至少3个基因更稳健
    Spatial_integrated <- AddModuleScore(Spatial_integrated, features = list(genes), name = paste0(set_name,"_Score"), assay = DefaultAssay(Spatial_integrated))
    # 生成形如 B_core_Score1 的列名
    colname <- paste0(set_name,"_Score1")
    Spatial_integrated@meta.data[[colname]] <- Spatial_integrated@meta.data[[colname]]
  } else {
    message(paste("集合", set_name, "可用基因不足 跳过"))
  }
}

# 4 在LIVE与DEATH间比较B细胞功能分数
# 仅在B_enriched spot中比较
md <- Spatial_integrated@meta.data
md$sample <- md$orig.ident
score_cols <- grep("_Score1$", colnames(md), value = TRUE)
regions <- sort(unique(Spatial_integrated$regionLevel))
compare_fun <- function(df, score_cols) {
  res_list <- lapply(score_cols, function(sc){
    d_sub <- df %>% filter(B_enriched)
    # 计算每个样本的均值
    means <- d_sub %>% group_by(sample) %>% summarise(mean_score = mean(.data[[sc]], na.rm = TRUE), .groups = "drop")
    # t检验与效应量
    x <- d_sub %>% filter(sample == "LIVE") %>% pull(sc)
    y <- d_sub %>% filter(sample == "DEATH") %>% pull(sc)
    if(length(x) > 2 && length(y) > 2) {
      tt <- tryCatch(t.test(x, y), error=function(e) NULL)
      # Cohen's d
      cohen_d <- (mean(x, na.rm=TRUE) - mean(y, na.rm=TRUE)) / sqrt(((sd(x,na.rm=TRUE)^2 + sd(y,na.rm=TRUE)^2))/2)
      pval <- if(!is.null(tt)) tt$p.value else NA
    } else {
      pval <- NA; cohen_d <- NA
    }
    data.frame(module = sub("_Score1","", sc),
               LIVE_mean = means$mean_score[means$sample=="LIVE"],
               DEATH_mean = means$mean_score[means$sample=="DEATH"],
               diff = means$mean_score[means$sample=="LIVE"] - means$mean_score[means$sample=="DEATH"],
               cohen_d = cohen_d,
               pval = pval)
  })
  do.call(rbind, res_list)
}
b_comp <- compare_fun(md, score_cols)
print("B细胞功能模块在LIVE与DEATH的比较：")
print(b_comp)

# 5 可视化模块分数差异
plot_df <- md %>% filter(B_enriched) %>% select(sample, all_of(score_cols)) %>%
  melt(id.vars = "sample", variable.name = "module", value.name = "score")
plot_df$module <- sub("_Score1","", plot_df$module)

p_violin <- ggplot(plot_df, aes(x = module, y = score, fill = sample)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, position = position_dodge(width = 0.75)) +
  scale_fill_manual(values = c("LIVE" = "forestgreen", "DEATH" = "red3")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "B细胞功能模块分数在LIVE与DEATH的分布", x = "模块", y = "Module score")
print(p_violin)

# 6 空间可视化比较B细胞权重与关键功能
DefaultAssay(Spatial_integrated) <- "celltypeprops"
if(!is.null(b_prop_name)) {
  p_sp1 <- SpatialFeaturePlot(Spatial_integrated, features = b_prop_name, pt.size.factor = 6500, images = names(Spatial_integrated@images)) + ggtitle("B细胞权重空间分布")
  print(p_sp1)
}
DefaultAssay(Spatial_integrated) <- if("RNA" %in% names(Spatial_integrated@assays)) "RNA" else if("SCT" %in% names(Spatial_integrated@assays)) "SCT" else "Spatial"
key_modules_to_show <- c("B_naive","B_plasma","B_isotype","B_GC")
key_cols <- paste0(key_modules_to_show, "_Score1")
key_cols <- key_cols[key_cols %in% colnames(Spatial_integrated@meta.data)]
for(col in key_cols) {
  Spatial_integrated@meta.data[[col]] <- Spatial_integrated@meta.data[[col]]
  p_sp <- SpatialFeaturePlot(Spatial_integrated, features = col, pt.size.factor = 6500, images = names(Spatial_integrated@images)) +
    ggtitle(paste("空间分布", col))
  print(p_sp)
}

# 7 B细胞相关基因差异表达比较
# 在B_enriched spot之间比较LIVE vs DEATH的表达（以RNA或SCT为默认assay）
DefaultAssay(Spatial_integrated) <- if("RNA" %in% names(Spatial_integrated@assays)) "RNA" else if("SCT" %in% names(Spatial_integrated@assays)) "SCT" else "Spatial"
Idents(Spatial_integrated) <- Spatial_integrated$orig.ident
b_cells <- WhichCells(Spatial_integrated, expression = B_enriched)
if(length(b_cells) > 10) {
  de_markers <- unique(unlist(b_sets))
  de_markers <- de_markers[de_markers %in% rownames(Spatial_integrated)]
  # 只在这些基因上做差异
  de_res <- FindMarkers(Spatial_integrated, ident.1 = "LIVE", ident.2 = "DEATH",
                        subset.ident = NULL, features = de_markers,
                        test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0.1,
                        assay = DefaultAssay(Spatial_integrated), slot = "data",
                        subset.cells = b_cells)
  de_res$gene <- rownames(de_res)
  de_res <- de_res[order(de_res$p_val_adj, -abs(de_res$avg_log2FC)), ]
  print("B细胞相关基因差异表达结果：")
  print(head(de_res, 30))
  
  # 火山图
  de_res$significant <- de_res$p_val_adj < 0.05
  p_vol <- ggplot(de_res, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significant, label = gene)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red")) +
    theme_bw() + labs(title = "LIVE vs DEATH 在B_enriched spots的差异表达", x = "log2FC LIVE-DEATH", y = "-log10(adj p)")
  print(p_vol)
} else {
  message("B_enriched spots数量不足 跳过差异表达分析")
}

# 8 受体-配体与抗原呈递线索（简化版打分）
# B细胞与T细胞交互相关基因集
b_t_interaction <- c("CD40","CD40LG","ICOSLG","ICOS","CXCL13","CXCR5","IL21R","IL21")
b_t_interaction <- b_t_interaction[b_t_interaction %in% rownames(Spatial_integrated)]
if(length(b_t_interaction) >= 3) {
  Spatial_integrated <- AddModuleScore(Spatial_integrated, features = list(b_t_interaction), name = "BT_interact", assay = DefaultAssay(Spatial_integrated))
  md <- Spatial_integrated@meta.data
  p_bt <- ggplot(md %>% filter(B_enriched), aes(x = orig.ident, y = BT_interact1, fill = orig.ident)) +
    geom_violin(trim = TRUE) + geom_boxplot(width = 0.1, outlier.size = 0.5) +
    scale_fill_manual(values = c("LIVE" = "forestgreen", "DEATH" = "red3")) +
    theme_bw() + labs(title = "B-T交互打分在B_enriched spots中的差异", x = "", y = "BT_interact score")
  print(p_bt)
}

# 9 可选通路富集分析（基于差异结果rank）
if(exists("de_res") && nrow(de_res) > 0) {
  ranks <- de_res$avg_log2FC
  names(ranks) <- de_res$gene
  # 以MSigDB的免疫相关集合为例（需要msigdbr）
  imm_sets <- msigdbr(species = "Homo sapiens", category = "C7") %>% select(gs_name, gene_symbol)
  pathways <- split(imm_sets$gene_symbol, imm_sets$gs_name)
  pathways <- lapply(pathways, function(gs) gs[gs %in% rownames(Spatial_integrated)])
  set.seed(1)
  fgsea_res <- fgsea(pathways = pathways, stats = ranks, minSize = 10, maxSize = 300, nperm = 1000)
  fgsea_res <- fgsea_res %>% arrange(padj)
  print("B细胞相关差异的免疫通路富集结果（MSigDB C7）：")
  print(head(fgsea_res, 20))
}

# 10 汇总与保存
write.csv(b_comp, file = "Bcell_module_compare_LIVE_DEATH.csv", row.names = FALSE)
if(exists("de_res")) write.csv(de_res, file = "Bcell_DEG_LIVE_DEATH.csv", row.names = FALSE)
saveRDS(Spatial_integrated, file = "Spatial_integrated_Bcell_analysis.rds")

# 11 额外 可视化B细胞比例差异（若有celltypeprops）
if(!is.null(b_prop_name) && "celltypeprops" %in% names(Spatial_integrated@assays)) {
  DefaultAssay(Spatial_integrated) <- "celltypeprops"
  dfp <- FetchData(Spatial_integrated, vars = c(b_prop_name, "orig.ident"))
  p_bar <- dfp %>%
    group_by(orig.ident) %>%
    summarise(mean_B_prop = mean(.data[[b_prop_name]], na.rm = TRUE), .groups = "drop") %>%
    ggplot(aes(x = orig.ident, y = mean_B_prop, fill = orig.ident)) +
    geom_bar(stat = "identity", width = 0.6) +
    scale_fill_manual(values = c("LIVE" = "forestgreen", "DEATH" = "red3")) +
    theme_bw() + labs(title = "B细胞权重均值在LIVE与DEATH的对比", x = "", y = "Mean B weight")
  print(p_bar)
  DefaultAssay(Spatial_integrated) <- if("RNA" %in% names(Spatial_integrated@assays)) "RNA" else if("SCT" %in% names(Spatial_integrated@assays)) "SCT" else "Spatial"
}

# 12 复现实验性检查
sessionInfo()


#3 在每个region内，比较LIVE vs DEATH的B细胞功能
md <- Spatial_integrated@meta.data
md$sample <- md$orig.ident
stopifnot("regionLevel" %in% colnames(Spatial_integrated@meta.data)) regions <- sort(unique(Spatial_integrated$regionLevel)) print(regions)
region_compare_list <- lapply(regions, function(rg){
  d_sub <- md %>% filter(regionLevel == rg & B_enriched)
  if(nrow(d_sub) < 10) {
    return(data.frame(region = rg, module = NA, LIVE_mean = NA, DEATH_mean = NA, diff = NA, cohen_d = NA, pval = NA))
  }
  out <- lapply(score_cols, function(sc){
    x <- d_sub %>% filter(sample == "LIVE") %>% pull(sc)
    y <- d_sub %>% filter(sample == "DEATH") %>% pull(sc)
    LIVE_mean <- mean(x, na.rm=TRUE); DEATH_mean <- mean(y, na.rm=TRUE)
    if(length(na.omit(x)) >= 3 && length(na.omit(y)) >= 3) {
      tt <- tryCatch(t.test(x, y), error=function(e) NULL)
      pval <- if(!is.null(tt)) tt$p.value else NA
      cohen_d <- (LIVE_mean - DEATH_mean) / sqrt(((sd(x,na.rm=TRUE)^2 + sd(y,na.rm=TRUE)^2))/2)
    } else {
      pval <- NA; cohen_d <- NA
    }
    data.frame(region = rg,
               module = sub("_Score1","", sc),
               LIVE_mean = LIVE_mean,
               DEATH_mean = DEATH_mean,
               diff = LIVE_mean - DEATH_mean,
               cohen_d = cohen_d,
               pval = pval)
  })
  do.call(rbind, out)
})
region_compare <- do.call(rbind, region_compare_list)

#多重校正
region_compare$padj <- ave(region_compare$pval, region_compare$region, FUN = function(p) p.adjust(p, method="BH"))
print("各区域B细胞功能模块 LIVE vs DEATH 比较：")
print(head(region_compare[order(region_compare$region, region_compare$padj), ], 50))
score_cols <- paste0(c("B_GC","B_isotype","B_plasma","B_naive"), "_Score1")
#4 可视化：按区域展示模块分布（仅B_enriched）
plot_df <- md %>% filter(B_enriched) %>% select(sample, regionLevel, all_of(score_cols)) %>%
  melt(id.vars = c("sample","regionLevel"), variable.name = "module", value.name = "score")
plot_df$module <- sub("_Score1","", plot_df$module)

p_ridge_list <- lapply(unique(plot_df$regionLevel), function(rg){
  ggplot(plot_df %>% filter(regionLevel == rg), aes(x = module, y = score, fill = sample)) +
    geom_violin(aes(fill = sample), position = position_dodge(width = 0.8), trim = FALSE) +
    geom_boxplot(aes(fill = sample), width = 0.12, outlier.size = 0.4, position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = c("LIVE" = "forestgreen", "DEATH" = "red3")) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0("Region: ", rg), x = "Module", y = "Module score")
})
pdf("B细胞评分在各个区域的小提琴图.pdf",width=9,height=7)
wrap_plots(p_ridge_list, ncol = 3)
dev.off()
#5 区域内差异表达（可选）：在每个region的B_enriched spots中做LIVE vs DEATH
Idents(Spatial_integrated) <- "orig.ident"
DefaultAssay(Spatial_integrated) <- if("RNA" %in% names(Spatial_integrated@assays)) "RNA" else if("SCT" %in% names(Spatial_integrated@assays)) "SCT" else "Spatial"

de_gene_panel <- unique(unlist(b_sets))
de_gene_panel <- intersect(de_gene_panel, rownames(Spatial_integrated))

de_by_region <- list()
for(rg in regions) {
  cells_rg <- rownames(Spatial_integrated@meta.data)[Spatial_integrated$regionLevel == rg & Spatial_integrated$B_enriched %in% TRUE]
  if(length(cells_rg) >= 20 && length(intersect(cells_rg, WhichCells(Spatial_integrated, idents = "LIVE"))) >= 5 && length(intersect(cells_rg, WhichCells(Spatial_integrated, idents = "DEATH"))) >= 5) {
    de_res_rg <- tryCatch({
      FindMarkers(Spatial_integrated,
                  ident.1 = "LIVE",
                  ident.2 = "DEATH",
                  features = de_gene_panel,
                  subset.cells = cells_rg,
                  test.use = "wilcox",
                  min.pct = 0.1,
                  logfc.threshold = 0.1)
    }, error=function(e) NULL)
    if(!is.null(de_res_rg)) {
      de_res_rg$gene <- rownames(de_res_rg)
      de_res_rg$region <- rg
      de_by_region[[rg]] <- de_res_rg[order(de_res_rg$p_val_adj, -abs(de_res_rg$avg_log2FC)), ]
    }
  } else {
    message(paste("Region", rg, "B_enriched数据不足，跳过DE"))
  }
}

#合并并导出
if(length(de_by_region) > 0) {
  de_all <- do.call(rbind, de_by_region)
  write.csv(de_all, file = "Bcell_DE_by_region_LIVE_vs_DEATH.csv", row.names = FALSE)
  print(head(de_all, 30))
}

#6 区域层面的B细胞富集程度与关键功能热图
#计算每个region、每个样本下的均值
agg_fun <- md %>%
  group_by(regionLevel, sample) %>%
  summarise(B_weight_mean = mean(B_weight, na.rm = TRUE),
            across(all_of(score_cols), ~mean(.x, na.rm=TRUE), .names = "{.col}"),
            .groups = "drop")

#转成长表
agg_long <- agg_fun %>%
  melt(id.vars = c("regionLevel","sample"), variable.name = "feature", value.name = "mean_value")
head(agg_long)
#只显示模块分数（可选）
agg_long$feature_clean <- sub("_Score1","", agg_long$feature)
p_tile <- ggplot(
  agg_long %>% filter(feature == "B_weight_mean"),
  aes(x = sample, y = regionLevel, fill = mean_value)
) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "navy", mid = "white", high = "firebrick",
    midpoint = median(agg_long$mean_value[agg_long$feature == "B_weight_mean"], na.rm = TRUE)
  ) +
  theme_bw() +
  labs(title = "Regional-level B-cell weight (mean)", x = "", y = "", fill = "Mean")
print(p_tile)

#7 空间可视化：针对代表性区域绘制B权重/关键模块
#如果有celltypeprops的B权重列名 b_prop_name，可替换为它；否则使用B_weight或B_core模块分数
feat_for_spatial <- NULL
if("celltypeprops" %in% names(Spatial_integrated@assays)) {
  DefaultAssay(Spatial_integrated) <- "celltypeprops"
  props <- tryCatch(rownames(Spatial_integrated[["celltypeprops"]]@data), error=function(e) NULL)
  if(!is.null(props)) {
    cand <- grep("B|Bcell|B_cell|B.cells|Bcells", props, value = TRUE, ignore.case = TRUE)
    if(length(cand) > 0) feat_for_spatial <- cand[1]
  }
}
DefaultAssay(Spatial_integrated) <- if("RNA" %in% names(Spatial_integrated@assays)) "RNA" else if("SCT" %in% names(Spatial_integrated@assays)) "SCT" else "Spatial"

features_to_show <- c("B_naive_Score1","B_plasma_Score1","B_isotype_Score1","B_GC_Score1")
features_to_show <- intersect(features_to_show, colnames(Spatial_integrated@meta.data))

#仅可视化若干关键区域（例如 Bcellenrich、bronchial、Diffuse damage、macrophageenriched）
show_regions <- intersect(c("vessel","fibroblastexpanding","Bcellenrich","Bleeding","bronchial","Diffuse damage","macrophageenriched","alveolarenlarged"), regions)

rg_obj <- subset(Spatial_integrated, subset = regionLevel == rg)
DefaultAssay(rg_obj) <- "Spatial"
if (!is.null(feat_for_spatial)) {
  print(SpatialFeaturePlot(rg_obj, features = feat_for_spatial, pt.size.factor = 6500) + ggtitle(paste("B权重 -", rg)))
}
DefaultAssay(rg_obj) <- if ("RNA" %in% names(rg_obj@assays)) "RNA" else if ("SCT" %in% names(rg_obj@assays)) "SCT" else "Spatial"
for (fc in features_to_show) {
  print(SpatialFeaturePlot(Spatial_integrated, features = fc, pt.size.factor = 6500) + ggtitle(paste(fc, "-", rg)))
}

#8 区域内B细胞比例（若有B权重，用阈值>0或top阈值近似比例）
region_Bstats <- Spatial_integrated@meta.data %>%
  group_by(regionLevel, orig.ident) %>%
  summarise(
    n_spots = n(),
    n_B_enriched = sum(B_enriched %in% TRUE, na.rm = TRUE),
    B_enriched_ratio = n_B_enriched / n_spots,
    mean_B_weight = mean(B_weight, na.rm = TRUE),
    .groups = "drop"
  )
print(region_Bstats)
write.csv(region_Bstats, "Bcell_region_summary.csv", row.names = FALSE)

p_ratio <- ggplot(region_Bstats, aes(x = regionLevel, y = B_enriched_ratio, fill = orig.ident)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  scale_fill_manual(values = c("LIVE" = "forestgreen", "DEATH" = "red3")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "各区域B_enriched比例（LIVE vs DEATH）", x = "Region", y = "B_enriched ratio")
print(p_ratio)
pdf("各区域B_enriched比例（LIVE vs DEATH）.pdf",width = 7,height =4 )
p_ratio
dev.off()
#9 统计检验：每个region内 LIVE vs DEATH 的关键模块p值表
modules_key <- c("B_naive_Score1","B_plasma_Score1","B_isotype_Score1","B_GC_Score1")
modules_key <- intersect(modules_key, sub("_Score1","", score_cols))
pv_tbl <- do.call(rbind, lapply(regions, function(rg){
  d_sub <- md %>% filter(regionLevel == rg & B_enriched)
  if(nrow(d_sub) < 10) return(NULL)
  out <- lapply(modules_key, function(m){
    sc <- paste0(m,"_Score1")
    x <- d_sub %>% filter(sample == "LIVE") %>% pull(sc)
    y <- d_sub %>% filter(sample == "DEATH") %>% pull(sc)
    p <- tryCatch(t.test(x, y)$p.value, error=function(e) NA)
    data.frame(region = rg, module = m, pval = p,
               LIVE_mean = mean(x, na.rm=TRUE), DEATH_mean = mean(y, na.rm=TRUE))
  })
  do.call(rbind, out)
}))
if(!is.null(pv_tbl)) {
  pv_tbl$padj <- ave(pv_tbl$pval, pv_tbl$module, FUN = function(p) p.adjust(p, method="BH"))
  write.csv(pv_tbl, "Bcell_region_module_pvalues.csv", row.names = FALSE)
  print(head(pv_tbl[order(pv_tbl$padj), ], 30))
}

#10 保存对象
saveRDS(Spatial_integrated, file = "Spatial_integrated_Bcell_region_analysis.rds")
