library(Seurat)
library(ggplot2)
library(cowplot)
library(BiocParallel)
library(BiocNeighbors)
library(data.table)
library(dplyr)
library(Matrix)
library(readr)

setwd("~/2026sc/important_data/rds")
total<-readRDS(file="nk.rds")
DimPlot(total, reduction = "umap", label = TRUE, pt.size = 0.5) 
FeaturePlot(total,features=c("Lgmn","Ctsb","Ctsd","Ctss","H2-Ab1"))
# 过滤：5个基因 同时 大于0 的细胞
total_filtered <- subset(
  total,
  Lgmn > 0 & Ctsb > 0 & Ctsd > 0 & Ctss > 0 & 'H2-Ab1' > 0
)

# 画小提琴图（只显示 >0 的细胞）
VlnPlot(
  total_filtered,
  features = c("Lgmn", "Ctsb", "Ctsd", "Ctss", "H2-Ab1"),
  pt.size = 0.001,
  group.by = "tech"
)

# 1. 先提取5个基因（H2-Ab1 正常写）
expr <- FetchData(total, vars = c("Lgmn","Ctsb","Ctsd","Ctss","H2-Ab1","tech"))

# 2. 批量计算非零中位数（重点在这里！）
MedianExp <- expr %>%
  group_by(tech) %>%
  summarise(
    across(
      # 带 - 的基因必须用 反引号 ` 包起来！！
      .cols = c(Lgmn, Ctsb, Ctsd, Ctss, `H2-Ab1`),
      .fns = ~ median(.x[.x > 2], na.rm = TRUE)
    ),
    .groups = "drop"
  )

# 查看结果
MedianExp
# 直接导出 CSV（中文不乱码）
write.csv(MedianExp, "五个基因非零中位数.csv", row.names = FALSE, fileEncoding = "UTF-8")


