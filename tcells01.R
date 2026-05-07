library(Seurat)
library(ggplot2)
library(cowplot)
library(BiocParallel)
library(BiocNeighbors)
library(data.table)
library(dplyr)
library(Matrix)

setwd("~/2026sc")
hms_cluster<-readRDS(file="t_test_10.rds")
DimPlot(hms_cluster, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster, reduction = "umap", label = FALSE, pt.size = 0.5) 
new.cluster.ids <- c("stromal","stromal","stromal","activated cd4","Double-negative (DN) γδ T","exhausted cd8","nk","nk","effector cd8","stromal","treg",
                     "exhausted cd8","stromal","stromal","exhausted cd8","Double-negative (DN) γδ T","effector cd8","nk","Double-negative (DN) γδ T","nk","stromal",
                     "nk","Double-negative (DN) γδ T","stromal","stromal","nk","stromal","proliferating cd8","macrophage","nk","Double-negative (DN) γδ T",
                     "cd4 Th2","nk","Double-negative (DN) γδ T","b","cd4 Th2","stromal","treg","stromal","nk","Double-negative (DN) γδ T",
                     "Double-negative (DN) γδ T","macrophage","Double-negative (DN) γδ T","exhausted cd8","exhausted cd8","exhausted cd8","macrophage","macrophage","Double-negative (DN) γδ T","stromal",
                     "nk","exhausted cd8","stromal","endothelial","stromal","activated cd8","activated cd8","naive cd4","nk","mast",
                     "stromal","stromal","activated cd4","endothelial","stromal","activated cd4","nk","treg","nk","nk",
                     "nk","stromal","activated cd4","stromal","stromal","macrophage","stromal","nk","stromal","Double-negative (DN) γδ T",
                     "stromal","treg","treg","activated cd8","Double-negative (DN) γδ T","treg","treg","macrophage","activated cd8","Double-negative (DN) γδ T",
                     "Double-negative (DN) γδ T","nk","b","Double-negative (DN) γδ T","stromal","macrophage") 

names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(hms_cluster_id, file = "hms_cluster_id_only_t.rds")
Only_lymphoid<-subset(hms_cluster_id, idents=c("activated cd4", "activated cd8","Double-negative (DN) γδ T","exhausted cd8","effector cd8","treg",
                                               "cd4 Th2","naive cd4","proliferating cd8"))
DimPlot(Only_lymphoid, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(Only_lymphoid, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(Only_lymphoid, file = "Only_t_cluster_id_test.rds")

table(Idents(Only_lymphoid), Only_lymphoid$orig.ident)

# 分组：自动把样本归类到 G1、G2、G3、G4
Only_lymphoid$group <- ifelse(
  Only_lymphoid$orig.ident %in% c("G1_71", "G1_84", "G1_87"), "G1",
  ifelse(Only_lymphoid$orig.ident %in% c("G2_63", "G2_73", "G2_79"), "G2",
         ifelse(Only_lymphoid$orig.ident %in% c("G3_68", "G3_80", "G3_88"), "G3",
                ifelse(Only_lymphoid$orig.ident %in% c("G4_83", "G4_90", "G4_92"), "G4", NA)
         )
  )
)

# 检查分组是否正确（运行看看结果对不对）
table(Only_lymphoid$orig.ident, Only_lymphoid$group)

DimPlot(
  Only_lymphoid,
  group.by = "group",      # 按我们分好的4组上色
  label = TRUE,            # 显示组名
  repel = TRUE,            # 标签不重叠
  cols = c("G1"="#E53935", "G2"="#43A047", "G3"="#1E88E5", "G4"="#FDD835") # 红/绿/蓝/黄
)
DimPlot(
  Only_lymphoid,
  group.by = "ident",    # 按细胞类型上色
  split.by = "group",    # 按 G1/G2/G3/G4 分成4张图
  ncol = 2,              # 2列排列
  label = TRUE,
  repel = TRUE
)
DimPlot(
  Only_lymphoid,
  group.by = "ident",    # 按细胞类型上色
  split.by = "group",    # 按 G1/G2/G3/G4 分成4张图
  ncol = 2,              # 2列排列
  label = FALSE,
  repel = TRUE
)
