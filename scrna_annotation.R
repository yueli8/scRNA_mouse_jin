# ============================================================
# 118926 个细胞一起分群 + 注释
# 输入：total_annotated00.rds
# 输出：total_reclustered_allcells.rds
#       recluster_top20_genes_per_cluster_allcells.csv
#       recluster_top20_genes_per_cluster_allcells_transposed.csv
# ============================================================
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/2026sc/important_data/rds")
combined <- readRDS("total_annotated00.rds")
print(combined)

# ============================================================
# 1. 重新降维（RNA assay）
# ============================================================
DefaultAssay(combined) <- "RNA"

combined <- NormalizeData(combined, verbose = FALSE)
combined <- FindVariableFeatures(combined, nfeatures = 2000, verbose = FALSE)
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 50, verbose = FALSE,
                   reduction.name = "pca_reclust")

pdf("elbowplot.pdf", width = 8, height = 5)
print(ElbowPlot(combined, ndims = 50))
dev.off()

# ============================================================
# 2. 全部 118926 个细胞一起聚类
# ============================================================
use_dims <- 1:30
res_use  <- 15.0   # 越大 cluster 越多；建议 2.0-5.0；15.0 会出几百个

combined <- FindNeighbors(combined,
                          reduction = "pca_reclust",
                          dims = use_dims,
                          graph.name = "reclust_snn")

combined <- FindClusters(combined,
                         graph.name = "reclust_snn",
                         resolution = res_use)

cat("\ncluster 数：", length(unique(Idents(combined))), "\n")
print(table(Idents(combined)))

combined <- RunUMAP(combined,
                    reduction = "pca_reclust",
                    dims = use_dims,
                    reduction.name = "umap_reclust")

p1 <- DimPlot(combined, reduction = "umap_reclust", label = TRUE) +
  ggtitle(paste0("Re-clustered (res=", res_use, ")"))
p2 <- DimPlot(combined, reduction = "umap_reclust", group.by = "tech") +
  ggtitle("By group (G1-G4)")
print(p1 | p2)
ggsave("umap_reclust_allcells.pdf", p1 | p2, width = 16, height = 7)

# ============================================================
# 3. 按 tech 拆开跑 FindMarkers（避免整体 JoinLayers 爆内存）
# ============================================================
techs <- c("G1", "G2", "G3", "G4")
markers_all_tech <- list()

for (g in techs) {
  cat("\n>>> 处理：", g, "\n")
  
  cells_g <- rownames(combined@meta.data)[combined$tech == g]
  cat("  细胞数：", length(cells_g), "\n")
  
  sub_obj <- subset(combined, cells = cells_g)
  sub_obj <- DietSeurat(sub_obj, assays = "RNA")
  sub_obj[["RNA"]] <- JoinLayers(sub_obj[["RNA"]])
  
  # 用 colnames(sub_obj) 保证顺序一致
  sub_obj$seurat_clusters <- combined$seurat_clusters[colnames(sub_obj)]
  
  DefaultAssay(sub_obj) <- "RNA"
  Idents(sub_obj)       <- "seurat_clusters"
  
  m_list <- list()
  for (cl in levels(Idents(sub_obj))) {
    n_cl <- sum(Idents(sub_obj) == cl)
    if (n_cl < 3) next   # 跳过细胞数 < 3 的 cluster
    
    m <- FindMarkers(sub_obj,
                     ident.1         = cl,
                     only.pos        = TRUE,
                     min.pct         = 0.1,
                     logfc.threshold = 0.1,
                     verbose         = FALSE)
    if (nrow(m) > 0) {
      m$cluster <- cl
      m$gene    <- rownames(m)
      m$tech    <- g
      m_list[[cl]] <- m
    }
  }
  
  markers_all_tech[[g]] <- do.call(rbind, m_list)
  cat("  该组 marker 数：", nrow(markers_all_tech[[g]]), "\n")
  rm(sub_obj); gc()
}

markers_all <- do.call(rbind, markers_all_tech)
cat("\nmarker 总行数：", nrow(markers_all), "\n")

# ============================================================
# 4. 合并四组：每个 cluster + gene 取平均 log2FC
# ============================================================
markers_avg <- markers_all %>%
  group_by(cluster, gene) %>%
  summarise(avg_log2FC = mean(avg_log2FC, na.rm = TRUE),
            n_tech = n(),
            .groups = "drop")

# 每个 cluster top 20
markers_top <- markers_avg %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# 宽表：行=cluster，列=Top1..Top20
top20_wide <- markers_top %>%
  select(cluster, rank, gene) %>%
  pivot_wider(names_from = rank, values_from = gene,
              names_prefix = "Top") %>%
  arrange(as.numeric(as.character(cluster)))

write.csv(top20_wide,
          "recluster_top20_genes_per_cluster_allcells.csv",
          row.names = FALSE)
print(as.data.frame(top20_wide), right = FALSE)

# ============================================================
# 5. 转置版：行=Top1..Top20，列=cluster
# ============================================================
top20_t <- as.data.frame(t(top20_wide[, -1]), stringsAsFactors = FALSE)
colnames(top20_t) <- top20_wide$cluster
top20_t <- cbind(Top = rownames(top20_t), top20_t)
rownames(top20_t) <- NULL

write.csv(top20_t,
          "recluster_top20_genes_per_cluster_allcells_transposed.csv",
          row.names = FALSE)
print(top20_t, right = FALSE)

# ============================================================
# 6. 人工注释：把 cluster 映射成细胞类型
# ============================================================
# 看 top20_t 后，按你的判断填这张表：
cluster2celltype <- c(
  # "0"  = "Fibroblast",
  # "1"  = "DC",
  # "2"  = "NK",
  # ... 按实际填
)

if (length(cluster2celltype) > 0) {
  combined$celltype_new <- cluster2celltype[as.character(combined$seurat_clusters)]
  print(table(combined$celltype_new))
  
  p3 <- DimPlot(combined, reduction = "umap_reclust",
                group.by = "celltype_new", label = TRUE) +
    ggtitle("Annotated celltypes")
  print(p3)
}

# ============================================================
# 7. 保存
# ============================================================
saveRDS(combined, "total_reclustered_allcells.rds")
cat("\n已保存：total_reclustered_allcells.rds\n")
cat("marker 表：recluster_top20_genes_per_cluster_allcells.csv\n")
cat("转置表：  recluster_top20_genes_per_cluster_allcells_transposed.csv\n")



# =========================================================
# 加载包
# =========================================================
library(Seurat)
library(ggplot2)

# =========================================================
# 细胞类型注释向量（按 reclust_snn_res.15 的 cluster 顺序）
# =========================================================
celltype_annotations <- c(
  "n","T","T","T","T","T","macrophage","DC","neutrophil","B",
  "neutrophil","DC","B","macrophage","macrophage","macrophage","T","NK","DC","ILC",
  "B","macrophage","T","macrophage","T","DC","NK","T","DC","macrophage",
  "DC","NK","DC","T","T","T","T","monocyte","T","neutrophil",
  "DC","monocyte","NK","ILC","T","ILC","mast","NK","monocyte","monocyte",
  "granulocyte","macrophage","T","monocyte","macrophage","neutrophil","granulocyte","T","neutrophil","monocyte",
  "macrophage","macrophage","monocyte","DC","monocyte","monocyte","monocyte","granulocyte","macrophage","T",
  "monocyte","monocyte","T","T","granulocyte","macrophage","DC","granulocyte","DC","T",
  "T","T","macrophage","NK","T","macrophage","DC","granulocyte","DC","T",
  "macrophage","macrophage","B","T","T","macrophage","DC","macrophage","neutrophil","DC",
  "n","monocyte","n","n","B","DC","NK","n","monocyte","n",
  "n","n","DC","T","n","DC","n","NK","n","NK",
  "macrophage","DC","DC","n","NK","n","n","macrophage","n","B",
  "monocyte","n","T","B","T","n","n","macrophage","T","n",
  "neutrophil","DC","T","neutrophil","macrophage","T","NK","n","NK","n",
  "T","n","n","T","T","n","n","n","neutrophil","T",
  "n","n","n","B","T","n","n","granulocyte","n","DC",
  "macrophage","n","DC","DC","n","macrophage","T","NK","DC","macrophage",
  "n","n","n","B","n","n","n","macrophage","n","T",
  "macrophage","n","n","n","NK","DC","T","macrophage","neutrophil","n",
  "n","n","n","n","macrophage","n","T","n","DC","n",
  "T","T","n","n","n","n","n","n","macrophage","n",
  "n","DC","n","granulocyte","B","n"
)

# =========================================================
# 把注释写入 combined
# =========================================================
cluster_col <- "reclust_snn_res.15"

cluster_levels <- levels(combined@meta.data[[cluster_col]])

stopifnot(length(cluster_levels) == length(celltype_annotations))

names(celltype_annotations) <- cluster_levels

new_celltype <- celltype_annotations[
  as.character(combined@meta.data[[cluster_col]])
]
new_celltype <- unname(new_celltype)

combined$celltype <- new_celltype

Idents(combined) <- "celltype"

# =========================================================
# 画注释好细胞类型的 UMAP 图
# =========================================================
if (!"umap" %in% Reductions(combined)) {
  combined <- RunUMAP(combined, dims = 1:30)
}

# 带 label
DimPlot(
  combined,
  reduction = "umap",
  group.by  = "celltype",
  label     = TRUE,
  repel     = TRUE
) + ggtitle("UMAP by celltype")

# 不带 label
DimPlot(
  combined,
  reduction = "umap",
  group.by  = "celltype"
) + ggtitle("UMAP by celltype")


# =========================================================
# 加载包
# =========================================================
library(Seurat)
library(ggplot2)

# =========================================================
# 先确认 Idents 是 celltype
# =========================================================
Idents(combined) <- "celltype"

# 查看当前 celltype 种类
print(table(combined$celltype))

# =========================================================
# 去掉 n，保留其他 10 种细胞
# =========================================================
keep_celltypes <- c(
  "B", "DC", "granulocyte", "ILC", "macrophage",
  "mast", "monocyte", "neutrophil", "NK", "T"
)

combined_sub <- subset(
  combined,
  subset = celltype %in% keep_celltypes
)

# 检查
print(table(combined_sub$celltype))

# =========================================================
# 重新画图
# =========================================================

# 带 label
DimPlot(
  combined_sub,
  reduction = "umap",
  group.by  = "celltype",
  label     = TRUE,
  repel     = TRUE
) + ggtitle("UMAP by celltype (n removed)")

# 不带 label
DimPlot(
  combined_sub,
  reduction = "umap",
  group.by  = "celltype"
) + ggtitle("UMAP by celltype (n removed)")

# =========================================================
# 如需保存
# =========================================================
saveRDS(combined_sub, "total_reclustered_allcells_no_n.rds")


DimPlot(
  combined_sub,
  reduction = "umap",
  group.by  = "celltype",
  label     = FALSE,
  repel     = TRUE,
  pt.size   = 0.00001     # 默认是 auto，改小即可
) + ggtitle("UMAP by celltype (n removed)")