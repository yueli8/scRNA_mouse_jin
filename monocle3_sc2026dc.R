# 完整分析流程（Dendritic Cells, DC）：
# 1. 加载 dc_annotated01.rds
# 2. 提取 DC 亚群 (moDC, pDC, tolerogenic DC, mregDC, cDC1, cDC2)
# 3. 拆分为 G1~G4 四个组
# 4. 对每组进行 Monocle3 轨迹分析
# 5. 绘制 40 个基因沿拟时序表达图
# 6. 绘制拟时序分布箱线图
# =====================================================================

# ==================== 0. 加载必要的包 ====================
library(Seurat)
library(monocle3)
library(ggplot2)
library(tidyr)
library(SingleCellExperiment)

# ==================== 1. 加载原始数据 ====================
setwd("~/2026sc/important_data/rds")
hms_cluster_id <- readRDS("dc_annotated01.rds")

DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.5)

# ==================== 2. 提取 DC 亚群 ====================
dc <- subset(
  hms_cluster_id, 
  idents = c("moDC", "pDC", "tolerogenic DC", "mregDC", "cDC1", "cDC2")
)

# 可视化 subset 后的数据
DimPlot(dc, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(dc, reduction = "umap", label = FALSE, pt.size = 0.5)

# 查看细胞类型分布
print("DC 亚群细胞类型分布：")
print(table(Idents(dc)))

# ==================== 3. 按样本拆分为 G1~G4 四组 ====================
g1_samples <- c("G1_71", "G1_84", "G1_87")
g2_samples <- c("G2_63", "G2_73", "G2_79")
g3_samples <- c("G3_68", "G3_80", "G3_88")
g4_samples <- c("G4_83", "G4_90", "G4_92")

# 检查样本是否存在
all_samples <- unique(dc@meta.data$orig.ident)
print("当前数据中包含的样本：")
print(all_samples)

# 拆分
dc_G1 <- subset(dc, subset = orig.ident %in% g1_samples)
dc_G2 <- subset(dc, subset = orig.ident %in% g2_samples)
dc_G3 <- subset(dc, subset = orig.ident %in% g3_samples)
dc_G4 <- subset(dc, subset = orig.ident %in% g4_samples)

# 打印各组细胞数
print(paste("G1 细胞数:", ncol(dc_G1)))
print(paste("G2 细胞数:", ncol(dc_G2)))
print(paste("G3 细胞数:", ncol(dc_G3)))
print(paste("G4 细胞数:", ncol(dc_G4)))

# 保存拆分后的 Seurat 对象（可选）
saveRDS(dc_G1, "dc_G1_subset.rds")
saveRDS(dc_G2, "dc_G2_subset.rds")
saveRDS(dc_G3, "dc_G3_subset.rds")
saveRDS(dc_G4, "dc_G4_subset.rds")

# =====================================================================
# 4. 对 G1~G4 分别进行 Monocle3 轨迹分析（封装为函数）
# =====================================================================
run_monocle_analysis <- function(seurat_obj, group_name) {
  
  cat("\n========================================\n")
  cat("开始 Monocle3 分析:", group_name, "\n")
  cat("========================================\n")
  
  # --- 4.1 准备表达矩阵 ---
  hvg_genes <- rownames(seurat_obj)
  expression_data <- FetchData(seurat_obj, vars = hvg_genes, layer = "data")
  expression_matrix <- t(as.matrix(expression_data))
  expression_matrix <- as(expression_matrix, "sparseMatrix")
  
  # --- 4.2 构建 Monocle3 对象 ---
  cell_metadata <- seurat_obj@meta.data
  cell_metadata$celltype <- as.character(Idents(seurat_obj))
  
  gene_metadata <- data.frame(
    gene_short_name = rownames(expression_matrix),
    row.names = rownames(expression_matrix)
  )
  
  cds <- new_cell_data_set(
    expression_data = expression_matrix,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata
  )
  
  # --- 4.3 预处理 ---
  cds <- preprocess_cds(cds, num_dim = 50)
  
  # --- 4.4 降维 ---
  cds <- reduce_dimension(cds)
  
  # --- 4.5 聚类（数据量小，分辨率低） ---
  cds <- cluster_cells(cds, resolution = 1e-3)
  
  # --- 4.6 学习轨迹 ---
  cds <- learn_graph(cds)
  
  # --- 4.7 排序细胞（拟时序） ---
  # DC 起始：优先使用 moDC 或 cDC1
  root_cells <- colnames(cds)[colData(cds)$celltype == "moDC"]
  if (length(root_cells) == 0) {
    root_cells <- colnames(cds)[colData(cds)$celltype == "cDC1"]
  }
  if (length(root_cells) > 0) {
    cds <- order_cells(cds, root_cells = root_cells[1])
    print(paste("使用根节点:", root_cells[1]))
  } else {
    cds <- order_cells(cds)
    print("自动选择根节点")
  }
  
  # --- 4.8 保存 Monocle3 对象 ---
  saveRDS(cds, paste0(group_name, "_dc_monocle3_final.rds"))
  
  # --- 4.9 可视化（拟时序 + 细胞类型） ---
  p1 <- plot_cells(cds, color_cells_by = "pseudotime",
                   label_cell_groups = FALSE, label_leaves = FALSE,
                   label_branch_points = FALSE) +
    ggtitle(paste(group_name, "DC Pseudotime")) + theme_minimal() +
    scale_color_gradientn(colors = c("blue", "cyan", "green", "yellow", "red"))
  print(p1)
  
  p2 <- plot_cells(cds, color_cells_by = "celltype",
                   label_cell_groups = TRUE, label_leaves = TRUE,
                   label_branch_points = TRUE) +
    ggtitle(paste(group_name, "DC Cell Types in Trajectory")) + theme_minimal()
  print(p2)
  
  ggsave(paste0(group_name, "_dc_pseudotime.pdf"), p1, width = 10, height = 8)
  ggsave(paste0(group_name, "_dc_celltypes.pdf"), p2, width = 10, height = 8)
  
  return(cds)
}

# 分别对四组运行 Monocle3 分析
cds_G1 <- run_monocle_analysis(dc_G1, "G1")
cds_G2 <- run_monocle_analysis(dc_G2, "G2")
cds_G3 <- run_monocle_analysis(dc_G3, "G3")
cds_G4 <- run_monocle_analysis(dc_G4, "G4")

# =====================================================================
# 5. 对 G1~G4 分别绘制 40 基因拟时序表达图（显示 + 保存）
# =====================================================================

# 40 个基因列表（小鼠基因，DC 相关）
target_genes <- c(
  # --- DC 身份与谱系 ---
  "Itgax", "Cd11c", "Cd11b", "Itgam", "Mhc2", "H2-Ab1", "Cd74",
  "Zbtb46", "Batf3", "Irf4", "Irf8", "Irf7", "Spib", "Tcf4", "Bcl11a",
  
  # --- cDC1 相关 ---
  "Xcr1", "Clec9a", "Batf3", "Irf8", "Cadm1", "Itgae",
  
  # --- cDC2 相关 ---
  "Cd1c", "Clec10a", "Cd209a", "Sirpa", "Itgb7", "Irf4", "Klf4",
  
  # --- pDC 相关 ---
  "Siglech", "Bst2", "Gzmb", "Irf7", "Tcf4", "Ly6c2",
  
  # --- moDC 相关 ---
  "Ly6c2", "Ccr2", "Fcgr1", "Fcgr3", "Cd14", "Itgam", "Nos2", "Cxcl10",
  
  # --- mregDC / tolerogenic DC 相关 ---
  "Fscn1", "Ccr7", "Cd83", "Cd40", "Il12b", "Aldh1a2", "Ido1", "Tgfb1",
  "Socs2", "Socs3", "Il10", "Pdcd1lg2", "Cd274",
  
  # --- 共刺激与抗原提呈 ---
  "Cd80", "Cd86", "Cd40", "H2-Aa", "H2-Eb1", "Tap1", "Tap2", "B2m",
  
  # --- 趋化因子与受体 ---
  "Ccr7", "Ccr2", "Ccr5", "Cxcr4", "Ccr6", "Ccl5", "Ccl19", "Ccl22",
  
  # --- 激活与炎症 ---
  "Tnfsf9", "Tnf", "Il6", "Il1b", "Il12a", "Il23a", "Cxcl9", "Cxcl10"
)

group_list <- list(
  G1 = cds_G1,
  G2 = cds_G2,
  G3 = cds_G3,
  G4 = cds_G4
)

for (group_name in names(group_list)) {
  
  cat("\n========================================\n")
  cat("绘制 40 基因拟时序图:", group_name, "\n")
  cat("========================================\n")
  
  cds_obj <- group_list[[group_name]]
  
  # --- 5.1 过滤基因 ---
  valid_genes <- target_genes[target_genes %in% rownames(cds_obj)]
  missing_genes <- target_genes[!target_genes %in% rownames(cds_obj)]
  if (length(missing_genes) > 0) {
    cat(paste(group_name, "缺失基因：", paste(missing_genes, collapse = ", "), "\n"))
  }
  cat(paste(group_name, "即将绘制", length(valid_genes), "个基因\n"))
  if (length(valid_genes) == 0) next
  
  # --- 5.2 自动选择 assay ---
  if ("logcounts" %in% assayNames(cds_obj)) {
    target_assay <- "logcounts"
  } else if ("counts" %in% assayNames(cds_obj)) {
    target_assay <- "counts"
  } else {
    target_assay <- assayNames(cds_obj)[1]
  }
  
  # --- 5.3 提取数据 ---
  expr_matrix <- assay(cds_obj, target_assay)
  expr_subset <- t(as.matrix(expr_matrix[valid_genes, , drop = FALSE]))
  pseudotime_vals <- pseudotime(cds_obj)
  celltype_vals <- colData(cds_obj)$celltype
  
  expression_data <- data.frame(
    expr_subset,
    pseudotime = pseudotime_vals,
    celltype = celltype_vals,
    check.names = FALSE
  )
  
  # --- 5.4 转长数据 ---
  plot_data <- pivot_longer(
    data = expression_data,
    cols = all_of(valid_genes),
    names_to = "gene",
    values_to = "expression"
  )
  plot_data$gene <- factor(plot_data$gene, levels = valid_genes)
  
  # --- 5.5 绘制基因图 ---
  p_genes <- ggplot(plot_data, aes(x = pseudotime, y = expression)) +
    geom_point(aes(color = celltype), alpha = 0.4, size = 0.6) +
    geom_smooth(method = "loess", color = "black", size = 1, se = FALSE, span = 0.8) +
    facet_wrap(~ gene, scales = "free_y", ncol = 8) +
    theme_minimal() +
    labs(
      title = paste0(group_name, " Group: ", length(valid_genes),
                     " Genes Expression Along Pseudotime (Dendritic Cells)"),
      x = "Pseudotime", y = "Expression", color = "Cell Type"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      strip.text = element_text(face = "bold", size = 9),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 12),
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      panel.grid.major = element_line(color = "grey90", size = 0.2),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey80", fill = NA, size = 0.4)
    ) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))
  
  print(p_genes)
  
  output_file <- paste0(group_name, "_40_Genes_Pseudotime_DC.pdf")
  ggsave(output_file, p_genes, width = 24, height = 16)
  cat(paste(group_name, "基因图已保存为:", output_file, "\n"))
  
  # --- 5.6 拟时序分布箱线图 ---
  cat(paste("绘制", group_name, "拟时序箱线图...\n"))
  
  pseudotime_data <- data.frame(
    celltype = colData(cds_obj)$celltype,
    pseudotime = pseudotime(cds_obj)
  )
  pseudotime_data <- pseudotime_data[!is.na(pseudotime_data$pseudotime), ]
  
  p_box <- ggplot(pseudotime_data, aes(x = celltype, y = pseudotime, fill = celltype)) +
    geom_boxplot(outlier.size = 0.8, outlier.color = "black", color = "black") +
    theme_minimal() +
    labs(
      title = paste0(group_name, " Group: Pseudotime Distribution by Cell Type"),
      x = "Cell Type", y = "Pseudotime", fill = "celltype"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      legend.position = "right",
      panel.grid.major = element_line(color = "grey90", size = 0.2),
      panel.grid.minor = element_blank(),
      panel.border = element_blank()
    )
  
  print(p_box)
  
  box_file <- paste0(group_name, "_dc_Pseudotime_Boxplot.pdf")
  ggsave(box_file, p_box, width = 10, height = 8)
  cat(paste(group_name, "箱线图已保存为:", box_file, "\n"))
}

cat("\n========================================\n")
cat("全部流程（G1~G4 DC）分析完成！\n")
cat("========================================\n")
