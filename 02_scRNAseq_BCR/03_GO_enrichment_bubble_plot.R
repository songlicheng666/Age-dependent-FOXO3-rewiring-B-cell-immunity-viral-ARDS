# =============================================================================
# Script: 03_GO_enrichment_bubble_plot.R
# Description: Multi-cluster GO enrichment bubble plot from Metascape/clusterProfiler output
# (gene ratio x -log10 p-value x gene count)
# Author: Licheng Song
# Institution: Chinese PLA General Hospital (Eighth Medical Center)
# Study: Age-dependent B-cell immune dysregulation in virus-associated ARDS
# Year: 2025
# Input: GO_bubble_4.xlsx (multi-sheet Excel: one sheet per cluster)
# Output: GO_enrichment_bubble_plot.pdf
# =============================================================================

library(dplyr)
library(ggplot2)  
library(tidyverse)
library(openxlsx)
library(readxl)
setwd("D:/博士毕业/毕业文章/F4ARDS组比较")

# cluster0数值导入
file_path <- "GO_bubble_4.xlsx"
sheet_names <- excel_sheets(file_path)  
# 创建一个空的列表来存储每个sheet的数据
all_data <- list()

# 改进的分数转换函数
fraction_to_decimal <- function(x) {
  x <- as.character(x)
  result <- sapply(x, function(val) {
    if(grepl("/", val)) {
      nums <- strsplit(val, "/")[[1]]
      as.numeric(nums[1]) / as.numeric(nums[2])
    } else {
      as.numeric(val)
    }
  })
  return(as.numeric(result))
}

# 读取每个sheet并添加cluster信息
for(i in seq_along(sheet_names)) {
  temp_data <- read_excel(file_path, sheet = sheet_names[i])
  temp_data$gene_ratio <- as.character(temp_data$gene_ratio)
  temp_data$gene_ratio <- fraction_to_decimal(temp_data$gene_ratio)
  temp_data$Cluster <- paste0("cluster_", sheet_names[i])
  all_data[[i]] <- temp_data
}

# 合并所有数据
combined_data <- bind_rows(all_data)

# **添加基因计数列**
# 从Symbols列计算基因数量
combined_data$Gene_Count <- sapply(combined_data$Symbols, function(x) {
  if(is.na(x) || x == "") return(0)
  length(strsplit(x, ",")[[1]])
})

# 文本换行函数
wrap_text <- function(text, width = 50) {
  sapply(text, function(x) {
    paste(strwrap(x, width = width), collapse = "\n")
  })
}

# 对Description应用换行
combined_data$Description <- wrap_text(combined_data$Description)

# **关键修改1：按gene_ratio排序GO通路**
# 计算每个通路的平均gene_ratio，然后排序
pathway_order <- combined_data %>%
  group_by(Description) %>%
  summarise(mean_gene_ratio = mean(gene_ratio, na.rm = TRUE)) %>%
  arrange(mean_gene_ratio) %>%
  pull(Description)

# 将Description转换为有序因子
combined_data$Description <- factor(combined_data$Description, levels = pathway_order)

# 创建气泡图
p <- ggplot(combined_data, 
            aes(x = gene_ratio, 
                y = Description,
                size = Gene_Count,  # 使用计算得到的基因数量
                color = LogP)) +
  # 添加带黑色边框的点
  geom_point(alpha = 0.7, shape = 21, stroke = 0.8, 
             aes(fill = LogP),  # 填充色基于LogP
             color = "black") +  # 边框颜色设为黑色
  scale_size_continuous(range = c(3, 15), name = "Gene Count") +
  scale_fill_gradient(low = "lightblue", high = "red", name = "-log10(p-value)") +
  scale_color_gradient(low = "lightblue", high = "red", name = "-log10(p-value)") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), 
                     limits = c(0, max(combined_data$gene_ratio) * 1.1)) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 10, color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    axis.line = element_line(color = "black", linewidth = 0.5)  # 修复：使用linewidth替代size
  ) +
  labs(title = "GO Pathway Enrichment Analysis of Cluster 8",
       x = "Gene Ratio",
       y = "GO Pathway Description")

print(p)

# 获取唯一的Description数量
n_unique_descriptions <- length(unique(combined_data$Description))
# 根据唯一Description的数量调整高度
height <- max(8, n_unique_descriptions * 0.4)

# 保存图片
ggsave("1_GO_enrichment_bubble_plot_sorted.pdf", plot = p, width = 7.5, height = 7.5)

# 查看数据结构（用于调试）
print("Data structure:")
print(str(combined_data))
print("Gene count range:")
print(range(combined_data$Gene_Count))








library(dplyr)
library(ggplot2)  
library(tidyverse)
library(openxlsx)
library(readxl)
setwd("D:/空转合并区域/FIB分析/myo_cluster差异基因heatmap")

# cluster0数值导入
file_path <- "GO_bubble_0.xlsx"
sheet_names <- excel_sheets(file_path)  
# 创建一个空的列表来存储每个sheet的数据
all_data <- list()

# 改进的分数转换函数
fraction_to_decimal <- function(x) {
  x <- as.character(x)
  result <- sapply(x, function(val) {
    if(grepl("/", val)) {
      nums <- strsplit(val, "/")[[1]]
      as.numeric(nums[1]) / as.numeric(nums[2])
    } else {
      as.numeric(val)
    }
  })
  return(as.numeric(result))
}

# 读取每个sheet并添加cluster信息
for(i in seq_along(sheet_names)) {
  temp_data <- read_excel(file_path, sheet = sheet_names[i])
  temp_data$gene_ratio <- as.character(temp_data$gene_ratio)
  temp_data$gene_ratio <- fraction_to_decimal(temp_data$gene_ratio)
  temp_data$Cluster <- paste0("cluster_", sheet_names[i])
  all_data[[i]] <- temp_data
}

# 合并所有数据
combined_data <- bind_rows(all_data)

# **添加基因计数列**
# 从Symbols列计算基因数量
combined_data$Gene_Count <- sapply(combined_data$Symbols, function(x) {
  if(is.na(x) || x == "") return(0)
  length(strsplit(x, ",")[[1]])
})

# 文本换行函数
wrap_text <- function(text, width = 50) {
  sapply(text, function(x) {
    paste(strwrap(x, width = width), collapse = "\n")
  })
}

# 对Description应用换行
combined_data$Description <- wrap_text(combined_data$Description)

# **关键修改1：按gene_ratio排序GO通路**
# 计算每个通路的平均gene_ratio，然后排序
pathway_order <- combined_data %>%
  group_by(Description) %>%
  summarise(mean_gene_ratio = mean(gene_ratio, na.rm = TRUE)) %>%
  arrange(mean_gene_ratio) %>%
  pull(Description)

# 将Description转换为有序因子
combined_data$Description <- factor(combined_data$Description, levels = pathway_order)

# 创建气泡图
p <- ggplot(combined_data, 
            aes(x = gene_ratio, 
                y = Description,
                size = Gene_Count,  # 使用计算得到的基因数量
                color = LogP)) +
  # 添加带黑色边框的点
  geom_point(alpha = 0.7, shape = 21, stroke = 0.8, 
             aes(fill = LogP),  # 填充色基于LogP
             color = "black") +  # 边框颜色设为黑色
  scale_size_continuous(range = c(3, 15), name = "Gene Count") +
  scale_fill_gradient(low = "lightblue", high = "red", name = "-log10(p-value)") +
  scale_color_gradient(low = "lightblue", high = "red", name = "-log10(p-value)") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), 
                     limits = c(0, max(combined_data$gene_ratio) * 1.1)) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 10, color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    axis.line = element_line(color = "black", linewidth = 0.5)  # 修复：使用linewidth替代size
  ) +
  labs(title = "GO Pathway Enrichment Analysis of Cluster 0",
       x = "Gene Ratio",
       y = "GO Pathway Description")

print(p)

# 获取唯一的Description数量
n_unique_descriptions <- length(unique(combined_data$Description))
# 根据唯一Description的数量调整高度
height <- max(8, n_unique_descriptions * 0.4)

# 保存图片
ggsave("0_GO_enrichment_bubble_plot_sorted.pdf", plot = p, width = 7.5, height = 7.5)

# 查看数据结构（用于调试）
print("Data structure:")
print(str(combined_data))
print("Gene count range:")
print(range(combined_data$Gene_Count))








library(dplyr)
library(ggplot2)  
library(tidyverse)
library(openxlsx)
library(readxl)
setwd("D:/空转合并区域/FIB分析/myo_cluster差异基因heatmap")

# cluster0数值导入
file_path <- "GO_bubble_9.xlsx"
sheet_names <- excel_sheets(file_path)  
# 创建一个空的列表来存储每个sheet的数据
all_data <- list()

# 改进的分数转换函数
fraction_to_decimal <- function(x) {
  x <- as.character(x)
  result <- sapply(x, function(val) {
    if(grepl("/", val)) {
      nums <- strsplit(val, "/")[[1]]
      as.numeric(nums[1]) / as.numeric(nums[2])
    } else {
      as.numeric(val)
    }
  })
  return(as.numeric(result))
}

# 读取每个sheet并添加cluster信息
for(i in seq_along(sheet_names)) {
  temp_data <- read_excel(file_path, sheet = sheet_names[i])
  temp_data$gene_ratio <- as.character(temp_data$gene_ratio)
  temp_data$gene_ratio <- fraction_to_decimal(temp_data$gene_ratio)
  temp_data$Cluster <- paste0("cluster_", sheet_names[i])
  all_data[[i]] <- temp_data
}

# 合并所有数据
combined_data <- bind_rows(all_data)

# **添加基因计数列**
# 从Symbols列计算基因数量
combined_data$Gene_Count <- sapply(combined_data$Symbols, function(x) {
  if(is.na(x) || x == "") return(0)
  length(strsplit(x, ",")[[1]])
})

# 文本换行函数
wrap_text <- function(text, width = 50) {
  sapply(text, function(x) {
    paste(strwrap(x, width = width), collapse = "\n")
  })
}

# 对Description应用换行
combined_data$Description <- wrap_text(combined_data$Description)

# **关键修改1：按gene_ratio排序GO通路**
# 计算每个通路的平均gene_ratio，然后排序
pathway_order <- combined_data %>%
  group_by(Description) %>%
  summarise(mean_gene_ratio = mean(gene_ratio, na.rm = TRUE)) %>%
  arrange(mean_gene_ratio) %>%
  pull(Description)

# 将Description转换为有序因子
combined_data$Description <- factor(combined_data$Description, levels = pathway_order)

# 创建气泡图
p <- ggplot(combined_data, 
            aes(x = gene_ratio, 
                y = Description,
                size = Gene_Count,  # 使用计算得到的基因数量
                color = LogP)) +
  # 添加带黑色边框的点
  geom_point(alpha = 0.7, shape = 21, stroke = 0.8, 
             aes(fill = LogP),  # 填充色基于LogP
             color = "black") +  # 边框颜色设为黑色
  scale_size_continuous(range = c(3, 15), name = "Gene Count") +
  scale_fill_gradient(low = "lightblue", high = "red", name = "-log10(p-value)") +
  scale_color_gradient(low = "lightblue", high = "red", name = "-log10(p-value)") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), 
                     limits = c(0, max(combined_data$gene_ratio) * 1.1)) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 10, color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    axis.line = element_line(color = "black", linewidth = 0.5)  # 修复：使用linewidth替代size
  ) +
  labs(title = "GO Pathway Enrichment Analysis of Cluster 9",
       x = "Gene Ratio",
       y = "GO Pathway Description")

print(p)

# 获取唯一的Description数量
n_unique_descriptions <- length(unique(combined_data$Description))
# 根据唯一Description的数量调整高度
height <- max(8, n_unique_descriptions * 0.4)

# 保存图片
ggsave("9_GO_enrichment_bubble_plot_sorted.pdf", plot = p, width = 7.5, height = 7.5)

# 查看数据结构（用于调试）
print("Data structure:")
print(str(combined_data))
print("Gene count range:")
print(range(combined_data$Gene_Count))








library(dplyr)
library(ggplot2)  
library(tidyverse)
library(openxlsx)
library(readxl)
setwd("D:/空转合并区域/FIB分析/myo_cluster差异基因heatmap")

# cluster0数值导入
file_path <- "GO_bubble_2.xlsx"
sheet_names <- excel_sheets(file_path)  
# 创建一个空的列表来存储每个sheet的数据
all_data <- list()

# 改进的分数转换函数
fraction_to_decimal <- function(x) {
  x <- as.character(x)
  result <- sapply(x, function(val) {
    if(grepl("/", val)) {
      nums <- strsplit(val, "/")[[1]]
      as.numeric(nums[1]) / as.numeric(nums[2])
    } else {
      as.numeric(val)
    }
  })
  return(as.numeric(result))
}

# 读取每个sheet并添加cluster信息
for(i in seq_along(sheet_names)) {
  temp_data <- read_excel(file_path, sheet = sheet_names[i])
  temp_data$gene_ratio <- as.character(temp_data$gene_ratio)
  temp_data$gene_ratio <- fraction_to_decimal(temp_data$gene_ratio)
  temp_data$Cluster <- paste0("cluster_", sheet_names[i])
  all_data[[i]] <- temp_data
}

# 合并所有数据
combined_data <- bind_rows(all_data)

# **添加基因计数列**
# 从Symbols列计算基因数量
combined_data$Gene_Count <- sapply(combined_data$Symbols, function(x) {
  if(is.na(x) || x == "") return(0)
  length(strsplit(x, ",")[[1]])
})

# 文本换行函数
wrap_text <- function(text, width = 50) {
  sapply(text, function(x) {
    paste(strwrap(x, width = width), collapse = "\n")
  })
}

# 对Description应用换行
combined_data$Description <- wrap_text(combined_data$Description)

# **关键修改1：按gene_ratio排序GO通路**
# 计算每个通路的平均gene_ratio，然后排序
pathway_order <- combined_data %>%
  group_by(Description) %>%
  summarise(mean_gene_ratio = mean(gene_ratio, na.rm = TRUE)) %>%
  arrange(mean_gene_ratio) %>%
  pull(Description)

# 将Description转换为有序因子
combined_data$Description <- factor(combined_data$Description, levels = pathway_order)

# 创建气泡图
p <- ggplot(combined_data, 
            aes(x = gene_ratio, 
                y = Description,
                size = Gene_Count,  # 使用计算得到的基因数量
                color = LogP)) +
  # 添加带黑色边框的点
  geom_point(alpha = 0.7, shape = 21, stroke = 0.8, 
             aes(fill = LogP),  # 填充色基于LogP
             color = "black") +  # 边框颜色设为黑色
  scale_size_continuous(range = c(3, 15), name = "Gene Count") +
  scale_fill_gradient(low = "lightblue", high = "red", name = "-log10(p-value)") +
  scale_color_gradient(low = "lightblue", high = "red", name = "-log10(p-value)") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), 
                     limits = c(0, max(combined_data$gene_ratio) * 1.1)) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 10, color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    axis.line = element_line(color = "black", linewidth = 0.5)  # 修复：使用linewidth替代size
  ) +
  labs(title = "GO Pathway Enrichment Analysis of Cluster 2",
       x = "Gene Ratio",
       y = "GO Pathway Description")

print(p)

# 获取唯一的Description数量
n_unique_descriptions <- length(unique(combined_data$Description))
# 根据唯一Description的数量调整高度
height <- max(8, n_unique_descriptions * 0.4)

# 保存图片
ggsave("2_GO_enrichment_bubble_plot_sorted.pdf", plot = p, width = 7.5, height = 7.5)

# 查看数据结构（用于调试）
print("Data structure:")
print(str(combined_data))
print("Gene count range:")
print(range(combined_data$Gene_Count))








library(dplyr)
library(ggplot2)  
library(tidyverse)
library(openxlsx)
library(readxl)
setwd("D:/空转合并区域/FIB分析/myo_cluster差异基因heatmap")

# cluster0数值导入
file_path <- "GO_bubble_1.xlsx"
sheet_names <- excel_sheets(file_path)  
# 创建一个空的列表来存储每个sheet的数据
all_data <- list()

# 改进的分数转换函数
fraction_to_decimal <- function(x) {
  x <- as.character(x)
  result <- sapply(x, function(val) {
    if(grepl("/", val)) {
      nums <- strsplit(val, "/")[[1]]
      as.numeric(nums[1]) / as.numeric(nums[2])
    } else {
      as.numeric(val)
    }
  })
  return(as.numeric(result))
}

# 读取每个sheet并添加cluster信息
for(i in seq_along(sheet_names)) {
  temp_data <- read_excel(file_path, sheet = sheet_names[i])
  temp_data$gene_ratio <- as.character(temp_data$gene_ratio)
  temp_data$gene_ratio <- fraction_to_decimal(temp_data$gene_ratio)
  temp_data$Cluster <- paste0("cluster_", sheet_names[i])
  all_data[[i]] <- temp_data
}

# 合并所有数据
combined_data <- bind_rows(all_data)

# **添加基因计数列**
# 从Symbols列计算基因数量
combined_data$Gene_Count <- sapply(combined_data$Symbols, function(x) {
  if(is.na(x) || x == "") return(0)
  length(strsplit(x, ",")[[1]])
})

# 文本换行函数
wrap_text <- function(text, width = 50) {
  sapply(text, function(x) {
    paste(strwrap(x, width = width), collapse = "\n")
  })
}

# 对Description应用换行
combined_data$Description <- wrap_text(combined_data$Description)

# **关键修改1：按gene_ratio排序GO通路**
# 计算每个通路的平均gene_ratio，然后排序
pathway_order <- combined_data %>%
  group_by(Description) %>%
  summarise(mean_gene_ratio = mean(gene_ratio, na.rm = TRUE)) %>%
  arrange(mean_gene_ratio) %>%
  pull(Description)

# 将Description转换为有序因子
combined_data$Description <- factor(combined_data$Description, levels = pathway_order)

# 创建气泡图
p <- ggplot(combined_data, 
            aes(x = gene_ratio, 
                y = Description,
                size = Gene_Count,  # 使用计算得到的基因数量
                color = LogP)) +
  # 添加带黑色边框的点
  geom_point(alpha = 0.7, shape = 21, stroke = 0.8, 
             aes(fill = LogP),  # 填充色基于LogP
             color = "black") +  # 边框颜色设为黑色
  scale_size_continuous(range = c(3, 15), name = "Gene Count") +
  scale_fill_gradient(low = "lightblue", high = "red", name = "-log10(p-value)") +
  scale_color_gradient(low = "lightblue", high = "red", name = "-log10(p-value)") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), 
                     limits = c(0, max(combined_data$gene_ratio) * 1.1)) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 10, color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    axis.line = element_line(color = "black", linewidth = 0.5)  # 修复：使用linewidth替代size
  ) +
  labs(title = "GO Pathway Enrichment Analysis of Cluster 1",
       x = "Gene Ratio",
       y = "GO Pathway Description")

print(p)

# 获取唯一的Description数量
n_unique_descriptions <- length(unique(combined_data$Description))
# 根据唯一Description的数量调整高度
height <- max(8, n_unique_descriptions * 0.4)

# 保存图片
ggsave("1_GO_enrichment_bubble_plot_sorted.pdf", plot = p, width = 7.5, height = 7.5)

# 查看数据结构（用于调试）
print("Data structure:")
print(str(combined_data))
print("Gene count range:")
print(range(combined_data$Gene_Count))








library(dplyr)
library(ggplot2)  
library(tidyverse)
library(openxlsx)
library(readxl)
setwd("D:/空转合并区域/FIB分析/myo_cluster差异基因heatmap")

# cluster0数值导入
file_path <- "GO_bubble_6.xlsx"
sheet_names <- excel_sheets(file_path)  
# 创建一个空的列表来存储每个sheet的数据
all_data <- list()

# 改进的分数转换函数
fraction_to_decimal <- function(x) {
  x <- as.character(x)
  result <- sapply(x, function(val) {
    if(grepl("/", val)) {
      nums <- strsplit(val, "/")[[1]]
      as.numeric(nums[1]) / as.numeric(nums[2])
    } else {
      as.numeric(val)
    }
  })
  return(as.numeric(result))
}

# 读取每个sheet并添加cluster信息
for(i in seq_along(sheet_names)) {
  temp_data <- read_excel(file_path, sheet = sheet_names[i])
  temp_data$gene_ratio <- as.character(temp_data$gene_ratio)
  temp_data$gene_ratio <- fraction_to_decimal(temp_data$gene_ratio)
  temp_data$Cluster <- paste0("cluster_", sheet_names[i])
  all_data[[i]] <- temp_data
}

# 合并所有数据
combined_data <- bind_rows(all_data)

# **添加基因计数列**
# 从Symbols列计算基因数量
combined_data$Gene_Count <- sapply(combined_data$Symbols, function(x) {
  if(is.na(x) || x == "") return(0)
  length(strsplit(x, ",")[[1]])
})

# 文本换行函数
wrap_text <- function(text, width = 50) {
  sapply(text, function(x) {
    paste(strwrap(x, width = width), collapse = "\n")
  })
}

# 对Description应用换行
combined_data$Description <- wrap_text(combined_data$Description)

# **关键修改1：按gene_ratio排序GO通路**
# 计算每个通路的平均gene_ratio，然后排序
pathway_order <- combined_data %>%
  group_by(Description) %>%
  summarise(mean_gene_ratio = mean(gene_ratio, na.rm = TRUE)) %>%
  arrange(mean_gene_ratio) %>%
  pull(Description)

# 将Description转换为有序因子
combined_data$Description <- factor(combined_data$Description, levels = pathway_order)

# 创建气泡图
p <- ggplot(combined_data, 
            aes(x = gene_ratio, 
                y = Description,
                size = Gene_Count,  # 使用计算得到的基因数量
                color = LogP)) +
  # 添加带黑色边框的点
  geom_point(alpha = 0.7, shape = 21, stroke = 0.8, 
             aes(fill = LogP),  # 填充色基于LogP
             color = "black") +  # 边框颜色设为黑色
  scale_size_continuous(range = c(3, 15), name = "Gene Count") +
  scale_fill_gradient(low = "lightblue", high = "red", name = "-log10(p-value)") +
  scale_color_gradient(low = "lightblue", high = "red", name = "-log10(p-value)") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), 
                     limits = c(0, max(combined_data$gene_ratio) * 1.1)) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 10, color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    axis.line = element_line(color = "black", linewidth = 0.5)  # 修复：使用linewidth替代size
  ) +
  labs(title = "GO Pathway Enrichment Analysis of Cluster 6",
       x = "Gene Ratio",
       y = "GO Pathway Description")

print(p)

# 获取唯一的Description数量
n_unique_descriptions <- length(unique(combined_data$Description))
# 根据唯一Description的数量调整高度
height <- max(8, n_unique_descriptions * 0.4)

# 保存图片
ggsave("6_GO_enrichment_bubble_plot_sorted.pdf", plot = p, width = 7.5, height = 7.5)

# 查看数据结构（用于调试）
print("Data structure:")
print(str(combined_data))
print("Gene count range:")
print(range(combined_data$Gene_Count))








# cluster0数值导入
file_path <- "GO_bubble_7.xlsx"
sheet_names <- excel_sheets(file_path)  
# 创建一个空的列表来存储每个sheet的数据
all_data <- list()

# 改进的分数转换函数
fraction_to_decimal <- function(x) {
  x <- as.character(x)
  result <- sapply(x, function(val) {
    if(grepl("/", val)) {
      nums <- strsplit(val, "/")[[1]]
      as.numeric(nums[1]) / as.numeric(nums[2])
    } else {
      as.numeric(val)
    }
  })
  return(as.numeric(result))
}

# 读取每个sheet并添加cluster信息
for(i in seq_along(sheet_names)) {
  temp_data <- read_excel(file_path, sheet = sheet_names[i])
  temp_data$gene_ratio <- as.character(temp_data$gene_ratio)
  temp_data$gene_ratio <- fraction_to_decimal(temp_data$gene_ratio)
  temp_data$Cluster <- paste0("cluster_", sheet_names[i])
  all_data[[i]] <- temp_data
}

# 合并所有数据
combined_data <- bind_rows(all_data)

# **添加基因计数列**
# 从Symbols列计算基因数量
combined_data$Gene_Count <- sapply(combined_data$Symbols, function(x) {
  if(is.na(x) || x == "") return(0)
  length(strsplit(x, ",")[[1]])
})

# 文本换行函数
wrap_text <- function(text, width = 50) {
  sapply(text, function(x) {
    paste(strwrap(x, width = width), collapse = "\n")
  })
}

# 对Description应用换行
combined_data$Description <- wrap_text(combined_data$Description)

# **关键修改1：按gene_ratio排序GO通路**
# 计算每个通路的平均gene_ratio，然后排序
pathway_order <- combined_data %>%
  group_by(Description) %>%
  summarise(mean_gene_ratio = mean(gene_ratio, na.rm = TRUE)) %>%
  arrange(mean_gene_ratio) %>%
  pull(Description)

# 将Description转换为有序因子
combined_data$Description <- factor(combined_data$Description, levels = pathway_order)

# 创建气泡图
p <- ggplot(combined_data, 
            aes(x = gene_ratio, 
                y = Description,
                size = Gene_Count,  # 使用计算得到的基因数量
                color = LogP)) +
  # 添加带黑色边框的点
  geom_point(alpha = 0.7, shape = 21, stroke = 0.8, 
             aes(fill = LogP),  # 填充色基于LogP
             color = "black") +  # 边框颜色设为黑色
  scale_size_continuous(range = c(3, 15), name = "Gene Count") +
  scale_fill_gradient(low = "lightblue", high = "red", name = "-log10(p-value)") +
  scale_color_gradient(low = "lightblue", high = "red", name = "-log10(p-value)") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), 
                     limits = c(0, max(combined_data$gene_ratio) * 1.1)) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 10, color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 14, hjust = 0.5, color = "black"),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    axis.line = element_line(color = "black", linewidth = 0.5)  # 修复：使用linewidth替代size
  ) +
  labs(title = "GO Pathway Enrichment Analysis of Cluster 7",
       x = "Gene Ratio",
       y = "GO Pathway Description")

print(p)

# 获取唯一的Description数量
n_unique_descriptions <- length(unique(combined_data$Description))
# 根据唯一Description的数量调整高度
height <- max(8, n_unique_descriptions * 0.4)

# 保存图片
ggsave("7_GO_enrichment_bubble_plot_sorted.pdf", plot = p, width = 7.5, height = 7.5)

# 查看数据结构（用于调试）
print("Data structure:")
print(str(combined_data))
print("Gene count range:")
print(range(combined_data$Gene_Count))