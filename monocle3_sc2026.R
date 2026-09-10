# ==================== 1. 加载必要的库 ====================
library(Seurat)
library(monocle3)
library(Matrix)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(patchwork)
# ==================== 2. 设置工作目录并读取数据 ====================
setwd("~/2026sc/important_data/rds")
hms_cluster_id <- readRDS("t_annotated01.rds")

# 查看原始数据
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.5)

# ==================== 3. 提取 CD+ T 细胞亚群 ====================
Only_cd4T <- subset(hms_cluster_id, idents = c("treg", "naive cd4", 
                                               "cd4 Th2", "activated cd4"))

# 可视化 subset 后的数据
DimPlot(Only_cd4T, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(Only_cd4T, reduction = "umap", label = FALSE, pt.size = 0.5)

# 查看细胞类型分布
print("细胞类型分布：")
print(table(Only_cd4T@active.ident))

# 保存 subset 后的数据
saveRDS(Only_cd4T, file = "cd4_id_test.rds")

# ==================== 4. 准备 Monocle3 数据 ====================
# 获取高变基因（HVGs）
hvg_genes <- rownames(Only_cd4T)
print(paste("高变基因数量:", length(hvg_genes)))

# 获取表达数据（使用 normalized data）
expression_data <- FetchData(Only_cd4T, vars = hvg_genes, layer = "data")

# 转置为基因×细胞矩阵
expression_matrix <- t(as.matrix(expression_data))
expression_matrix <- as(expression_matrix, "sparseMatrix")

# 验证矩阵维度
print(paste("表达矩阵维度:", nrow(expression_matrix), "基因 x", ncol(expression_matrix), "细胞"))

# ==================== 5. 创建 Monocle3 对象 ====================
# 准备细胞元数据
cell_metadata <- Only_cd4T@meta.data
cell_metadata$celltype <- as.character(Only_cd4T@active.ident)

# 准备基因元数据
gene_metadata <- data.frame(
  gene_short_name = rownames(expression_matrix),
  row.names = rownames(expression_matrix)
)

# 创建 Monocle3 对象
cds <- new_cell_data_set(
  expression_data = expression_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

print("Monocle3 对象创建成功！")
print(paste("维度:", nrow(cds), "基因 x", ncol(cds), "细胞"))

# ==================== 6. 预处理 ====================
print("开始预处理...")
cds <- preprocess_cds(cds, num_dim = 50)
print("预处理完成！")

# 查看方差解释
plot_pc_variance_explained(cds)

# ==================== 7. 降维 ====================
print("开始降维...")
cds <- reduce_dimension(cds)
print("降维完成！")

# 查看降维结果
plot_cells(cds)

# ==================== 8. 聚类 ====================
print("开始聚类...")
cds <- cluster_cells(cds, resolution = 1e-4)
print("聚类完成！")

# 查看聚类结果
plot_cells(cds, color_cells_by = "partition")

# ==================== 9. 学习轨迹 ====================
print("学习轨迹...")
cds <- learn_graph(cds)
print("轨迹学习完成！")

# 查看轨迹
plot_cells(cds, color_cells_by = "celltype")

# ==================== 10. 排序细胞（拟时序分析） ====================
print("排序细胞...")

# 选择根节点（naive cd4 作为起始）
root_cells <- colnames(cds)[colData(cds)$celltype == "naive cd4"]
print(paste("naive cd4 细胞数量:", length(root_cells)))

if(length(root_cells) > 0) {
  cds <- order_cells(cds, root_cells = root_cells[1])
  print(paste("使用根节点:", root_cells[1]))
} else {
  cds <- order_cells(cds)
  print("自动选择根节点")
}

print("排序完成！")

# ==================== 11. 可视化结果 ====================
# 11.1 按拟时序着色
p1 <- plot_cells(cds, 
                 color_cells_by = "pseudotime", 
                 label_cell_groups = FALSE,
                 label_leaves = FALSE,
                 label_branch_points = FALSE)
p1 <- p1 + ggtitle("CD4+ T Cell Pseudotime") + 
  theme_minimal() + 
  scale_color_gradientn(colors = c("blue", "cyan", "green", "yellow", "red"))
print(p1)

# 11.2 按细胞类型着色
p2 <- plot_cells(cds, 
                 color_cells_by = "celltype", 
                 label_cell_groups = TRUE,
                 label_leaves = TRUE,
                 label_branch_points = TRUE)
p2 <- p2 + ggtitle("CD4+ T Cell Types in Trajectory") + theme_minimal()
print(p2)

# 11.3 按分区着色
p3 <- plot_cells(cds, 
                 color_cells_by = "partition", 
                 label_cell_groups = TRUE)
p3 <- p3 + ggtitle("CD4+ T Cell Partitions") + theme_minimal()
print(p3)

# 11.4 按样本着色
p4 <- plot_cells(cds, 
                 color_cells_by = "orig.ident", 
                 label_cell_groups = FALSE)
p4 <- p4 + ggtitle("CD4+ T Cells by Sample") + theme_minimal()
print(p4)

# ==================== 12. 保存可视化结果 ====================
ggsave("cd4_monocle3_pseudotime.pdf", p1, width = 10, height = 8)
ggsave("cd4_monocle3_celltypes.pdf", p2, width = 10, height = 8)
ggsave("cd4_monocle3_partitions.pdf", p3, width = 10, height = 8)
ggsave("cd4_monocle3_samples.pdf", p4, width = 10, height = 8)

# ==================== 13. 保存 Monocle3 对象 ====================
saveRDS(cds, file = "cd4_monocle3_final.rds")

# ==================== 14. 差异基因分析 ====================
print("进行差异基因分析...")
pr_test_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)

# 筛选显著基因
sig_genes <- subset(pr_test_res, q_value < 0.05)
print(paste("显著基因数量:", nrow(sig_genes)))

# 按 morans_I 排序
top_genes <- row.names(sig_genes)[order(sig_genes$morans_I, decreasing = TRUE)][1:20]
print("Top 20 变化基因:")
print(top_genes)

# 保存差异基因结果
write.csv(pr_test_res, file = "cd4_monocle3_graph_test.csv")
write.csv(sig_genes, file = "cd4_monocle3_significant_genes.csv")

# ==================== 15. 基因沿拟时序的表达 ====================
# 使用 CD4+ T 细胞关键基因
cd4_key_genes <- c("Gzmb", "Pdcd1", "Ctla4", "Havcr2", "Tox", 
                   "Tcf7", "Mki67", "Cd69", "Ifng", "Lag3",
                   "Tigit", "Eomes", "Tbx21", "Cxcr3", "Il7r")
cd4_key_genes <- cd4_key_genes[cd4_key_genes %in% rownames(cds)]
print(paste("可用的关键基因:", length(cd4_key_genes)))

if(length(cd4_key_genes) > 0) {
  # 伪时序基因表达热图
  p5 <- plot_genes_in_pseudotime(cds[cd4_key_genes, ], 
                                 color_cells_by = "celltype",
                                 min_expr = 0.1)
  ggsave("cd4_key_genes_pseudotime.pdf", p5, width = 14, height = 10)
}

# ==================== 16. 拟时序分布分析 ====================
pseudotime_values <- pseudotime(cds)
summary(pseudotime_values)

# 按细胞类型的拟时序分布
pseudotime_data <- data.frame(
  celltype = colData(cds)$celltype,
  pseudotime = pseudotime(cds)
)

p_box <- ggplot(pseudotime_data, aes(x = celltype, y = pseudotime, fill = celltype)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("Pseudotime Distribution by Cell Type") +
  xlab("Cell Type") + ylab("Pseudotime") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_box)
ggsave("cd4_pseudotime_boxplot.pdf", p_box, width = 10, height = 6)

print("Monocle3 分析完成！所有结果已保存。")


# 1. 获取表达矩阵
expr_matrix <- exprs(cds)

# 2. 计算每个基因的统计信息
gene_stats <- data.frame(
  gene = rownames(expr_matrix),
  mean_expr = rowMeans(expr_matrix),
  pct_expr = rowMeans(expr_matrix > 0) * 100
)

# 3. 添加 morans_I（如果有）
if(exists("pr_test_res")) {
  gene_stats$morans_I <- pr_test_res$morans_I[match(gene_stats$gene, rownames(pr_test_res))]
} else {
  # 如果没有 pr_test_res，运行 graph_test
  pr_test_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
  gene_stats$morans_I <- pr_test_res$morans_I[match(gene_stats$gene, rownames(pr_test_res))]
}

# 4. 筛选可靠基因
reliable_genes <- gene_stats %>%
  filter(pct_expr > 20,      # 至少在20%的细胞中表达
         mean_expr > 0.1) %>% # 平均表达量 > 0.1
  arrange(desc(morans_I))     # 按 morans_I 降序

# 5. 查看前50个
head(reliable_genes, 50)

# 6. 提取基因名
reliable_genes_50 <- reliable_genes$gene[1:50]
print(reliable_genes_50)

reliable_genes_200 <- reliable_genes$gene[1:200]
print(reliable_genes_200)


# 提取拟时序和细胞类型
pseudotime_values <- pseudotime(cds)
cell_types <- colData(cds)$celltype

# 定义核心功能基因分组
core_functional_groups <- list(
  # 细胞毒性
  Cytotoxicity = c("Gzmb", "Prf1", "Ifng"),
  
  # T细胞受体
  TCR = c("Cd3g", "Cd8a", "Trac"),
  
  # 趋化因子
  Chemokines = c("Cxcl10", "Ccl5", "Cxcr6"),
  
  # 免疫检查点
  Checkpoints = c("Pdcd1", "Ctla4", "Lag3", "Tox"),
  
  # 增殖
  Proliferation = c("Mki67", "Top2a", "Stmn1"),
  
  # 转录因子
  TFs = c("Eomes", "Ikzf2", "Tox"),
  
  # 活化标记
  Activation = c("Cd69", "Il2ra", "Icos")
)

# 为每个功能组绘制（每组一页）
for(group_name in names(core_functional_groups)) {
  genes <- core_functional_groups[[group_name]]
  genes <- genes[genes %in% rownames(cds)]
  
  print(paste("绘制", group_name, "组，共", length(genes), "个基因"))
  print(genes)
  
  if(length(genes) > 0) {
    # 提取数据
    plot_data <- data.frame()
    for(gene in genes) {
      gene_expr <- exprs(cds)[gene, ]
      temp_data <- data.frame(
        gene = gene,
        pseudotime = pseudotime_values,
        expression = gene_expr,
        celltype = cell_types
      )
      plot_data <- rbind(plot_data, temp_data)
    }
    
    plot_data$gene <- factor(plot_data$gene, levels = genes)
    
    # 绘制
    p <- ggplot(plot_data, aes(x = pseudotime, y = expression)) +
      geom_point(aes(color = celltype), size = 1, alpha = 0.4) +
      geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 1.5) +
      facet_wrap(~gene, scales = "free_y", ncol = length(genes)) +
      theme_minimal() +
      ggtitle(paste(group_name, "Genes Along Pseudotime")) +
      xlab("Pseudotime") +
      ylab("Expression") +
      theme(
        legend.position = "bottom",
        axis.text = element_text(size = 8),
        strip.text = element_text(size = 10, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
      ) +
      guides(color = guide_legend(title = "Cell Type", nrow = 2))
    
    print(p)
    ggsave(paste0("cd4_core_", group_name, ".pdf"), p, 
           width = 12, height = 6, limitsize = FALSE)
  }
}





