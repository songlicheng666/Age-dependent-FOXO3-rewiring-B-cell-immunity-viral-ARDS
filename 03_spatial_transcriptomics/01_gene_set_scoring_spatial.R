# =============================================================================
# Script: 01_gene_set_scoring_spatial.R
# Description: AddModuleScore-based scoring of apoptosis, inflammation, IFN-alpha, and IFN-gamma
# gene sets at the Visium spot level; SpatialFeaturePlot visualization
# Author: Licheng Song
# Institution: Chinese PLA General Hospital (Eighth Medical Center)
# Study: Age-dependent B-cell immune dysregulation in virus-associated ARDS
# Year: 2025
# Input: Spatial_integrated (Seurat Visium object from prior publication GSE223793/PRJCA042363)
# Output: 空间IFNG打分.pdf, 空间IFNA打分.pdf, 空间Apoptosis打分.pdf, 空间Inflammation打分.pdf
# =============================================================================


#================================================================================
#                         10X空间转录组下游分析（2）
#================================================================================
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
setwd("C:/Users/dell/Desktop/新新分析/儿童空转")
#================================================================================
#                           1、基因集评分
#================================================================================
#AddModuleScore这些评分都可以
#评分有这么几个作用：
#第一，基因集是marker基因，可以根据评分看celltype
#第二，通路基因集合评分，可以看看通路活性变化
#第三，自定义基因集，例如与某些信号（如炎症、衰老等相关的基因）相关的变化


#这里我们用的是AddModuleScore，其他的Aucell，Ucell，VISION，ssGESEA等等都可以
#和转录组是一样的

markers <- list()
markers$apopotosis <- c("ADD1","AIFM3","ANKH","ANXA1","APP","ATF3","AVPR1A","BAX","BCAP31","BCL10","BCL2L1","BCL2L10","BCL2L11","BCL2L2","BID","BIK","BIRC3","BMF","BMP2","BNIP3L","BRCA1","BTG2","BTG3","CASP1","CASP2","CASP3","CASP4","CASP6","CASP7","CASP8","CASP9","CAV1","CCNA1","CCND1","CCND2","CD14","CD2","CD38","CD44","CD69","CDC25B","CDK2","CDKN1A","CDKN1B","CFLAR","CLU","CREBBP","CTH","CTNNB1","CYLD","DAP","DAP3","DCN","DDIT3","DFFA","DIABLO","DNAJA1","DNAJC3","DNM1L","DPYD","EBP","EGR3","EMP1","ENO2","ERBB2","ERBB3","EREG","ETF1","F2","F2R","FAS","FASLG","FDXR","FEZ1","GADD45A","GADD45B","GCH1","GNA15","GPX1","GPX3","GPX4","GSN","GSR","GSTM1","GUCY2D","HGF","HMGB2","HMOX1","HSPB1","IER3","IFITM3","IFNB1","IFNGR1","IGF2R", "IGFBP6", "IL18", "IL1A", "IL1B", "IL6", "IRF1", "ISG20", "JUN", "KRT18", "LEF1", "LGALS3", "LMNA", "PLPPR4", "LUM", "MADD", "MCL1", "MGMT", "MMP2", "NEDD9", "NEFH", "PAK1", "PDCD4", "PDGFRB", "PEA15", "PLAT", "PLCB2", "PMAIP1", "PPP2R5B", "PPP3R1", "PPT1", "PRF1", "PSEN1", "PSEN2", "PTK2", "RARA", "RELA", "RETSAT", "RHOB", "RHOT2", "RNASEL", "ROCK1", "SAT1", "SATB1", "SC5D", "SLC20A1", "SMAD7", "SOD1", "SOD2", "SPTAN1", "SQSTM1", "TAP1", "TGFB2", "TGFBR3", "TIMP1", "TIMP2", "TIMP3", "TNF", "TNFRSF12A", "TNFSF10", "TOP2A", "TSPO", "TXNIP", "VDAC2", "WEE1", "XIAP")
markers$Inflammation <- c("ABCA1","ABI1","ACVR1B","ACVR2A","ADM","ADORA2B","ADRM1","AHR","APLNR","AQP9","ATP2A2","ATP2B1","ATP2C1","AXL","BDKRB1","BEST1","BST2","BTG2","C3AR1","C5AR1","CALCRL","CCL17","CCL2","CCL20","CCL22","CCL24","CCL5","CCL7","CCR7","CCRL2","CD14","CD40","CD48","CD55","CD69","CD70","CD82","CDKN1A","CHST2","CLEC5A","CMKLR1","CSF1","CSF3","CSF3R","CX3CL1","CXCL10","CXCL11","CXCL6","CXCL9","CXCR6","CYBB","DCBLD2","EBI3","EDN1","EIF2AK2","EMP3","ADGRE1","EREG","F3","FFAR2","FPR1","FZD5","GABBR1","GCH1","GNA15","GNAI3","GP1BA","GPC3","GPR132","GPR183","HAS2","HBEGF","HIF1A","HPN","HRH1","ICAM1","ICAM4","ICOSLG","IFITM1","IFNAR1","IFNGR2","IL10","IL10RA","IL12B", "IL15", "IL15RA", "IL18", "IL18R1", "IL18RAP", "IL1A", "IL1B", "IL1R1", "IL2RB", "IL4R", "IL6", "IL7R", "CXCL8", "INHBA", "IRAK2", "IRF1", "IRF7", "ITGA5", "ITGB3", "ITGB8", "KCNA3", "KCNJ2", "KCNMB2", "KIF1B", "KLF6", "LAMP3", "LCK", "LCP2", "LDLR", "LIF", "LPAR1", "LTA", "LY6E", "LYN", "MARCO", "MEFV", "MEP1A", "MET", "MMP14", "MSR1", "MXD1", "MYC", "NAMPT", "NDP", "NFKB1", "NFKBIA", "NLRP3", "NMI", "NMUR1", "NOD2", "NPFFR2", "OLR1", "OPRK1", "OSM", "OSMR", "P2RX4", "P2RX7", "P2RY2", "PCDH7", "PDE4B", "PDPN", "PIK3R5", "PLAUR", "PROK2", "PSEN1", "PTAFR", "PTGER2", "PTGER4","PTGIR","PTPRE","PVR","RAF1","RASGRP1","RELA","RGS1","RGS16","RHOG","RIPK2","RNF144B","ROS1","RTP4","SCARF1","SCN1B","SELE","SELL","SELENOS","SEMA4D","SERPINE1","SGMS2","SLAMF1","SLC11A2","SLC1A2","SLC28A2","SLC31A1","SLC31A2","SLC4A4","SLC7A1","SLC7A2","SPHK1","SRI","STAB1","TACR1","TACR3","TAPBP","TIMP1","TLR1","TLR2","TLR3","TNFAIP6","TNFRSF1B","TNFRSF9","TNFSF10","TNFSF15","TNFSF9","TPBG","VIP")
markers$IFNA <- c("ADAR","B2M","BATF2","BST2","C1S","CASP1","CASP8","CCRL2","CD47","CD74","CMPK2","CNP","CSF1","CXCL10","CXCL11","DDX60","DHX58","EIF2AK2","ELF1","EPSTI1","MVB12A","TENT5A","CMTR1","GBP2","GBP4","GMPR","HERC6","HLA-C","IFI27","IFI30","IFI35","IFI44","IFI44L","IFIH1","IFIT2","IFIT3","IFITM1","IFITM2","IFITM3","IL15","IL4R","IL7","IRF1","IRF2","IRF7","IRF9","ISG15","ISG20","LAMP3","LAP3","LGALS3BP","LPAR6","LY6E","MOV10","MX1","NCOA7","NMI","NUB1","OAS1","OASL","OGFR","PARP12","PARP14","PARP9","PLSCR1","PNPT1","HELZ2","PROCR","PSMA3","PSMB8","PSMB9","PSME1","PSME2","RIPK2","RNF31","RSAD2","RTP4","SAMD9","SAMD9L","SELL","SLC25A28","SP110","STAT2","TAP1","TDRD7","TMEM140", "TRAFD1", "TRIM14", "TRIM21", "TRIM25", "TRIM26", "TRIM5", "TXNIP", "UBA7", "UBE2L6", "USP18", "WARS1"
)
markers$IFNG <- c("ADAR","APOL6","ARID5B","ARL4A","AUTS2","B2M","BANK1","BATF2","BPGM","BST2","BTG1","C1R","C1S","CASP1","CASP3","CASP4","CASP7","CASP8","CCL2","CCL5","CCL7","CD274","CD38","CD40","CD69","CD74","CD86","CDKN1A","CFB","CFH","CIITA","CMKLR1","CMPK2","CSF2RB","CXCL10","CXCL11","CXCL9","RIGI","DDX60","DHX58","EIF2AK2","EIF4E3","EPSTI1","FAS","FCGR1A","FGL2","FPR1","CMTR1","GBP4","GBP6","GCH1","GPR18","GZMA","HERC6","HIF1A","HLA-A","HLA-B","HLA-DMA","HLA-DQA1","HLA-DRB1","HLA-G","ICAM1","IDO1","IFI27","IFI30","IFI35","IFI44","IFI44L","IFIH1","IFIT1","IFIT2","IFIT3","IFITM2", "IFITM3", "IFNAR2", "IL10RA", "IL15", "IL15RA", "IL18BP", "IL2RB", "IL4R", "IL6", "IL7", "IRF1", "IRF2", "IRF4", "IRF5", "IRF7", "IRF8", "IRF9", "ISG15", "ISG20", "ISOC1", "ITGB7", "JAK2", "KLRK1", "LAP3", "LATS2", "LCP2", "LGALS3BP", "LY6E", "LYSMD2", "MARCHF1", "METTL7B", "MT2A", "MTHFD2", "MVP", "MX1", "MX2", "MYD88", "NAMPT", "NCOA3", "NFKB1", "NFKBIA", "NLRC5", "NMI", "NOD1", "NUP93", "OAS2", "OAS3", "OASL", "OGFR", "P2RY14", "PARP12", "PARP14", "PDE4B", "PELI1", "PFKP", "PIM1", "PLA2G4A", "PLSCR1", "PML", "PNP", "PNPT1", "HELZ2", "PSMA2", "PSMA3", "PSMB10", "PSMB2", "PSMB8", "PSMB9", "PSME1", "PSME2", "PTGS2", "PTPN1", "PTPN2", "PTPN6", "RAPGEF6", "RBCK1", "RIPK1", "RIPK2", "RNF213", "RNF31", "RSAD2", "RTP4",  "SAMD9L" ,"SAMHD1" ,"SECTM1" ,"SELP" ,"SERPING1" ,"SLAMF7" ,"SLC25A28" ,"SOCS1" ,"SOCS3" ,"SOD2" ,"SP110" ,"SPPL2A" ,"SRI" ,"SSPN" ,"ST3GAL5" ,"ST8SIA4" ,"STAT1" ,"STAT2" ,"STAT3" ,"STAT4" ,"TAP1" ,"TAPBP" ,"TDRD7" ,"TNFAIP2" ,"TNFAIP3" ,"TNFAIP6" ,"TNFSF10" ,"TOR1B" ,"TRAFD1" ,"TRIM14" ,"TRIM21" ,"TRIM25" ,"TRIM26" ,"TXNIP" ,"UBE2L6" ,"UPP1" ,"USP18" ,"VAMP5" ,"VAMP8" ,"VCAM1" ,"WARS1" ,"XAF1" ,"XCL1" ,"ZBP1" ,"ZNFX1")

library(UCell)

#计算评分
Spatial_integrated <- AddModuleScore(object = Spatial_integrated, 
                               features = markers)
Spatial_integrated$Apoptosis <- Spatial_integrated$Cluster1
Spatial_integrated$Inflammation <- Spatial_integrated$Cluster2
Spatial_integrated$IFNA <- Spatial_integrated$Cluster3
Spatial_integrated$IFNG <- Spatial_integrated$Cluster4

#plot==figure9
pdf("空间IFNG打分.pdf",width = 11,height = 4)
SpatialFeaturePlot(Spatial_integrated, 
                   features = c("IFNG"), 
                   alpha = c(0.1,1),
                   ncol = 2,pt.size.factor = 6000,image.alpha = 0.2,
                   min.cutoff = 0,
                   max.cutoff = 1)&
  theme_bw()&
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())
dev.off()


pdf("空间IFNA打分.pdf",width = 11,height = 4)
SpatialFeaturePlot(Spatial_integrated, 
                   features = c("IFNA"), 
                   alpha = c(0.1,1),
                   ncol = 2,pt.size.factor = 6000,image.alpha = 0.2,
                   min.cutoff = 0,
                   max.cutoff = 1)&
  theme_bw()&
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())
dev.off()
pdf("空间Apoptosis打分.pdf",width = 11,height = 4)
SpatialFeaturePlot(Spatial_integrated, 
                   features = c("Apoptosis"), 
                   alpha = c(0.1,1),
                   ncol = 2,pt.size.factor = 6000,image.alpha = 0.2,
                   min.cutoff = 0,
                   max.cutoff = 1)&
  theme_bw()&
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())
dev.off()

pdf("空间Inflammation打分.pdf",width = 11,height = 4)
SpatialFeaturePlot(Spatial_integrated, 
                   features = c("Inflammation"), 
                   alpha = c(0.1,1),
                   ncol = 2,pt.size.factor = 6000,image.alpha = 0.2,
                   min.cutoff = 0,
                   max.cutoff = 1)&
  theme_bw()&
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())
dev.off()


Idents(Spatial_integrated) <- "regionLevel"
# 获取 regionLevel 的因子水平顺序
region_levels <- levels(Spatial_integrated$regionLevel)
table(Spatial_integrated$regionLevel)
# 创建颜色向量（示例用6个颜色，需按实际分组数调整）
custom_colors <- c("#F8766D", "#00BFC4", "#7CAE00", "#C77CFF", "#FF61CC", 
                   "#00B0F6", "#FF9A00", "#A3A500", "#00C19F", "#DB72FB", "#FF5E5B")
names(custom_colors) <- levels(Spatial_integrated$regionLevel)  # 按因子命名
pdf("两个样本dimplot.pdf",width = 9,height = 4.5)
DimPlot(Spatial_integrated, group.by = "regionLevel", 
              split.by = "orig.ident", cols = custom_colors)
dev.off()
pdf("两个样本空间dimplot.pdf",width = 9,height = 4.5)
SpatialDimPlot(Spatial_integrated, 
                     stroke = 0, label = F,image.alpha = 0.2, 
                     pt.size.factor = 6500, label.box = FALSE, 
                     label.color = 'black', repel = TRUE, 
                     cols = custom_colors) #+  # 关键：统一颜色
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")
dev.off()
save(Spatial_integrated,file = "带评分的空间seurat对象.Rdata")
library(ggplot2)
library(dplyr)
library(Seurat)
library(tidyr)
# 提取评分数据
scores <- Spatial_integrated@meta.data %>%
  select(orig.ident, regionLevel, Apoptosis, Inflammation, IFNA, IFNG)
# 将 orig.ident 转换为因子，并设置顺序
scores <- scores %>%
  mutate(orig.ident = factor(orig.ident, levels = c("LIVE", "DEATH")))
# 将数据转换为长格式，方便 ggplot2 作图
scores_long <- scores %>%
  pivot_longer(cols = c(Apoptosis, Inflammation, IFNA, IFNG),
               names_to = "ScoreType", values_to = "Score")

# 绘制小提琴图
pdf("各样本regionLevel评分小提琴图.pdf", width = 10, height = 5)
ggplot(scores_long, aes(x = orig.ident, y = Score, fill = orig.ident)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.7) +
  facet_grid(ScoreType ~ regionLevel, scales = "free_y") +  # 使用 facet_grid
  scale_fill_manual(values = cols) +  # 使用自定义颜色
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  labs(y = "Score", title = "scores across region")
dev.off()
#================================================================================
#           2、marker基因作图（一些和scRNA一样的常规操作）
#================================================================================
#特征（marker）基因可视化
#-----plot marker genes average expression heatmap
#figure11
library(dplyr)
library(data.table)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
#marker genes
Idents(Spatial_integrated) <- "clusterLevel"
sp_markers <- FindAllMarkers(object = Spatial_integrated, only.pos = TRUE, min.pct = 0.25, 
                          thresh.use = 0.25, verbose = T, assay = 'SCT')
#write.csv(sp_markers,file = "自定义差异基因.csv")
#sp_markers <- read.csv("自定义差异基因.csv")
# 定义一个不在的操作符
'%!in%' <- function(x, y) {!('%in%'(x, y))}

# 筛选每个 cluster 的 top 5 基因
sp_markers_top_5 <- sp_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(
    avg_log2FC > 0 & 
      p_val_adj < 0.05 & 
      gene %!in% grep('^RPL|^RPS|^MT-|^KHD|^ALD|^MCEM|^ECH|^DNA|^KCNA|^BPI|^CFAP|^SNX|^ANK', sp_markers$gene, value = TRUE) &  # 排除以 "Rpl"、"Rps" 或 "Mt-" 开头的基因
      !grepl('D-E', gene) &  # 排除包含 "D-E" 的基因
      !grepl('[a-z]', gene)  # 排除包含小写字母的基因
  ) %>%
  top_n(n = 5, wt = avg_log2FC)


# 求每个 cluster 的平均表达值
sp_markers_top_5_avg <- AverageExpression(Spatial_integrated, assays = 'SCT', return.seurat = TRUE)

cols= c("#EDB931","#eb6841","#cc2a36","#00a0b0","#7A989A", "#849271", "#CF9546", "#C67052", "#C1AE8D",
                  "#3F6F76", "#C65840", "#62496F", "#69B7CE","#91323A", "#3A4960", "#6D7345", "#D7C969",
                  "#C1395E", "#AEC17B", "#E07B42", "#89A7C2", "#F0CA50","#a53e1f", "#457277", "#8f657d", "#8dcee2",
                  "#E69253", "#EDB931", "#E4502E", "#4378A0", "#272A2A","#3F6148", "#A4804C", "#4B5F80", "#DBD3A4")
pdf("每个区域的差异基因.pdf",width = 4,height = 8)
 DoHeatmap(sp_markers_top_5_avg, 
                      features=sp_markers_top_5$gene, 
                      size=3, 
                      angle=90,
                      group.bar.height=0.01, 
                     group.colors=cols, 
                     disp.min = -2, 
                       disp.max = 2, 
                      draw.lines=F) + 
    scale_fill_gradientn(colours = rev(brewer.pal(10, 'RdBu')))+
   theme(text = element_text(size = 8, color='black'))
dev.off()



library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr) # 用于添加统计学显著性 P值

# 1. 提取 B 细胞富集区域
# 假设你的 regionLevel 中确实有 "Bcellenrich" 这一项
b_cell_subset <- subset(Spatial_integrated, subset = regionLevel == "Bcellenrich")

# 检查一下提取出来的细胞数
print(table(b_cell_subset$orig.ident))

# 2. 确保分组因子顺序 (LIVE 在前，DEATH 在后，或者反过来，看你想怎么比)
b_cell_subset$orig.ident <- factor(b_cell_subset$orig.ident, levels = c("LIVE", "DEATH"))
# 提取 FOXO3 的表达数据
# 注意：使用 SCT 或 RNA assay 数据，建议用 data slot
foxo3_data <- FetchData(b_cell_subset, vars = c("FOXO3", "orig.ident", "Apoptosis", "Inflammation"))

# --- 绘制箱线图/小提琴图并添加 P 值 ---
pdf("B细胞区域_FOXO3_表达差异.pdf", width = 6, height = 5)
ggplot(foxo3_data, aes(x = orig.ident, y = FOXO3, fill = orig.ident)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  scale_fill_manual(values = c("LIVE" = "#00BFC4", "DEATH" = "#F8766D")) + # 呼应你之前的配色
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 6, label.y.npc = 0.95) + # 添加星号
  stat_compare_means(method = "wilcox.test", label.y.npc = 0.9) + # 添加具体数值
  theme_bw() +
  labs(title = "FOXO3 Expression in B-cell Enriched Regions", 
       y = "FOXO3 Expression Level", x = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()
# 定义 B 细胞经典 Marker，用于对比
b_marker <- "IGHD" # 或者 CD79A

# --- 绘制空间对比图 ---
# 我们对比展示：B细胞Marker分布 vs FOXO3分布
pdf("空间共定位_FOXO3_vs_Bcell.pdf", width = 9, height =7)
SpatialFeaturePlot(Spatial_integrated, 
                   features = c(b_marker, "FOXO3"), 
                   #split.by = "orig.ident", # 按样本分开画
                   ncol = 2, 
                   pt.size.factor = 6000, # 根据你的切片调整点的大小
                   image.alpha = 0.3,
                   stroke = 0,
                   min.cutoff = "q10", 
                   max.cutoff = "q95") &
  scale_fill_gradientn(colours = c("grey90", "yellow", "red")) & # 热图色：灰-黄-红
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"))
dev.off()
# 使用之前提取的 foxo3_data (仅包含 Bcellenrich 区域)

# --- 1. FOXO3 vs Apoptosis (凋亡) ---
p_apo <- ggplot(foxo3_data, aes(x = FOXO3, y = Apoptosis)) +
  geom_point(aes(color = orig.ident), size = 1, alpha = 0.6) +
  geom_smooth(method = "lm", color = "black", se = TRUE) + # 添加线性拟合线
  stat_cor(method = "pearson", label.x = 0.5) + # 添加相关性系数 R 和 p 值
  scale_color_manual(values = c("LIVE" = "#00BFC4", "DEATH" = "#F8766D")) +
  theme_bw() +
  labs(title = "Correlation: FOXO3 vs Apoptosis (B-cell Region)",
       x = "FOXO3 Expression", y = "Apoptosis Score")

# --- 2. FOXO3 vs Inflammation (炎症) ---
p_inf <- ggplot(foxo3_data, aes(x = FOXO3, y = Inflammation)) +
  geom_point(aes(color = orig.ident), size = 1, alpha = 0.6) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  stat_cor(method = "pearson", label.x = 0.5) +
  scale_color_manual(values = c("LIVE" = "#00BFC4", "DEATH" = "#F8766D")) +
  theme_bw() +
  labs(title = "Correlation: FOXO3 vs Inflammation (B-cell Region)",
       x = "FOXO3 Expression", y = "Inflammation Score")

# 拼图输出
pdf("B细胞区域_FOXO3_功能相关性.pdf", width = 10, height = 5)
p_apo + p_inf
dev.off()



DefaultAssay(b_cell_subset) <- "integrated"
# 切换 Idents 到分组
Idents(b_cell_subset) <- "orig.ident"

# 寻找 LIVE 相对于 DEATH 上调的基因
de_b_cell <- FindMarkers(b_cell_subset, ident.1 = "LIVE", ident.2 = "DEATH", 
                         min.pct = 0.1, logfc.threshold = 0.25)

# 添加基因名列
de_b_cell$gene <- rownames(de_b_cell)

# 看看 FOXO3 是否在差异基因列表里
print(de_b_cell["SFTPC", ])

# --- 绘制火山图 ---
library(ggrepel)

# 标记显著基因
de_b_cell$diff <- "NO"
de_b_cell$diff[de_b_cell$avg_log2FC > 0.5 & de_b_cell$p_val_adj < 0.05] <- "UP (Live)"
de_b_cell$diff[de_b_cell$avg_log2FC < -0.5 & de_b_cell$p_val_adj < 0.05] <- "DOWN (Death)"

# 标记 FOXO3 和一些 B 细胞相关基因
genes_to_label <- c("KRT17", "MS4A1", "IGHG3", "IGHG1", "BCL6", "JCHAIN", "IL10") # IL10 是 Breg 标志

pdf("B细胞区域_Live_vs_Death_火山图.pdf", width = 6, height = 5)
ggplot(de_b_cell, aes(x = avg_log2FC, y = -log10(p_val_adj), color = diff)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("DOWN (Death)" = "#F8766D", "NO" = "grey", "UP (Live)" = "#00BFC4")) +
  geom_text_repel(data = subset(de_b_cell, gene %in% genes_to_label),
                  aes(label = gene), color = "black", box.padding = 0.5, max.overlaps = Inf) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw() +
  labs(title = "DE Genes in B-cell Region: LIVE vs DEATH",
       subtitle = "Highlighting FOXO3 and B-cell markers")
dev.off()
