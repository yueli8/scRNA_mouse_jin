library(Seurat)
library(CellChat)
library(patchwork)
library(future)
library(ggplot2)
library(dplyr)

setwd("~/2026sc/important_data/rds")
plan("sequential")

# ---------- 1. 读入两个对象 ----------
myeloid <- readRDS("myeloid_annotated01.rds")
tcell   <- readRDS("t_annotated01.rds")

# ---------- 2. 把 Idents 写回 meta.data ----------
myeloid$celltype_real <- as.character(Idents(myeloid))
tcell$celltype_real   <- as.character(Idents(tcell))

# ---------- 3. 合并 ----------
combined <- merge(myeloid, y = tcell)
ncol(combined)                 # 51889
table(combined$celltype_real)  # 15 种细胞类型

# ---------- 3b. 确认分组列（你的对象里是 tech）----------
table(combined$tech)           # G1 G2 G3 G4

# 统一命名成 group（后续循环用）
combined$group <- combined$tech
table(combined$group)

# ---------- 4. 标准化 ----------
DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(combined)

# ---------- 5. 合并 layer ----------
combined <- JoinLayers(combined)

# ============================================================
# 分组循环：G1 / G2 / G3 / G4
# ============================================================

groups <- c("G1", "G2", "G3", "G4")
outdir <- "cellchat_results"
dir.create(outdir, showWarnings = FALSE)

for (g in groups) {
  
  cat("\n============================================\n")
  cat("正在处理分组：", g, "\n")
  cat("============================================\n")
  
  # ---------- A. 取子集 ----------
  sub_obj <- subset(combined, subset = group == g)
  cat("细胞数：", ncol(sub_obj), "\n")
  
  if (ncol(sub_obj) < 50) {
    cat("警告：", g, "细胞数过少，跳过\n"); next
  }
  
  # ---------- B. 提取表达矩阵 ----------
  data.input <- sub_obj[["RNA"]]$data
  
  # ---------- C. 构建 meta ----------
  meta <- data.frame(
    labels = sub_obj$celltype_real,
    row.names = colnames(data.input)
  )
  
  # 去掉该组中细胞数 < 10 的细胞类型（CellChat 建议）
  tb <- table(meta$labels)
  keep_labels <- names(tb)[tb >= 10]
  meta <- meta[meta$labels %in% keep_labels, , drop = FALSE]
  data.input <- data.input[, rownames(meta), drop = FALSE]
  
  cat("保留的细胞类型数：", length(keep_labels), "\n")
  print(table(meta$labels))
  
  if (length(keep_labels) < 2) {
    cat("警告：", g, "有效细胞类型 < 2，跳过\n"); next
  }
  
  # ---------- D. 创建 CellChat 对象 ----------
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
  
  # ---------- E. 小鼠数据库 ----------
  cellchat@DB <- CellChatDB.mouse
  
  # ---------- F. 检查基因匹配率 ----------
  db_genes <- unique(c(CellChatDB.mouse$interaction$ligand,
                       CellChatDB.mouse$interaction$receptor))
  matched <- intersect(rownames(data.input), db_genes)
  cat("匹配到的 LR 基因数：", length(matched), "\n")
  
  # ---------- G. 预处理 ----------
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # ---------- H. 开并行 + 计算通讯概率 ----------
  options(future.globals.maxSize = 8 * 1024^3)
  plan("multisession", workers = 4, maxSizeOfObjects = 8 * 1024^3)
  
  cellchat <- computeCommunProb(cellchat, type = "triMean", nboot = 100)
  
  # ---------- I. 立刻关并行 ----------
  plan("sequential")
  
  # ---------- J. 过滤 + 汇总 ----------
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  # ---------- K. 查看显著通路 ----------
  cat("显著通路数：", length(cellchat@netP$pathways), "\n")
  print(cellchat@netP$pathways)
  
  # ---------- L. 保存 ----------
  saveRDS(cellchat, file = file.path(outdir, paste0("cellchat_", g, ".rds")))
  cat("已保存：", file.path(outdir, paste0("cellchat_", g, ".rds")), "\n")
}

cat("\n全部分组处理完成！\n")

# ============================================================
# 四组 CellChat 图：全部保存到文件夹（不含 pathway 图）
# 直接复制整段运行
# ============================================================

#保存到一个文件夹里（不再屏幕显示）

# 结果总目录
outdir <- "cellchat_plots"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

groups <- c("G1", "G2", "G3", "G4")

myeloid_types <- c("macrophage", "dc", "neutrophil",
                   "monocyte", "granulocyte", "mast")
tcell_types   <- c("activated cd4", "Double-negative (DN) γδ T",
                   "exhausted cd8", "effector cd8", "treg",
                   "proliferating cd8", "cd4 Th2",
                   "activated cd8", "naive cd4")

# rds 所在目录（你前面跑计算时的输出目录）
rds_dir <- "cellchat_results"
#从这里开始调用已经保存的rds文件

setwd("~/2026sc/important_data/rds/cellchat_results")

library(CellChat)
library(ggplot2)
library(dplyr)

# ============================================================
# 路径设置
# ============================================================
rds_dir <- "~/2026sc/important_data/rds/cellchat_results"
outdir  <- "~/2026sc/important_data/rds/cellchat_plots"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

groups <- c("G1", "G2", "G3", "G4")

myeloid_types <- c("macrophage", "dc", "neutrophil",
                   "monocyte", "granulocyte", "mast")
tcell_types   <- c("activated cd4", "Double-negative (DN) γδ T",
                   "exhausted cd8", "effector cd8", "treg",
                   "proliferating cd8", "cd4 Th2",
                   "activated cd8", "naive cd4")

# ============================================================
# 循环：导入 rds + 画图
# ============================================================
for (g in groups) {
  
  cat("\n============================================\n")
  cat("  分组：", g, "\n")
  cat("============================================\n")
  
  # ---------- 找 rds ----------
  cc_file <- file.path(rds_dir, g, paste0("cellchat_", g, ".rds"))
  if (!file.exists(cc_file)) {
    cc_file <- file.path(rds_dir, paste0("cellchat_", g, ".rds"))
  }
  if (!file.exists(cc_file)) {
    cat("找不到 ", g, " 的 rds，跳过\n"); next
  }
  
  cellchat <- readRDS(cc_file)
  
  # ---------- 输出目录 ----------
  gdir <- file.path(outdir, g)
  dir.create(gdir, showWarnings = FALSE, recursive = TRUE)
  
  # ---------- 变量 ----------
  present_types <- levels(cellchat@idents)
  my_types_in   <- intersect(myeloid_types, present_types)
  t_types_in    <- intersect(tcell_types,   present_types)
  groupSize     <- as.numeric(table(cellchat@idents))
  
  # ==========================================================
  # 1. 圆图（count + weight）
  # ==========================================================
  pdf(file.path(gdir, "01_circle_overall.pdf"), width = 14, height = 7)
  par(mfrow = c(1, 2), xpd = TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                   weight.scale = TRUE, label.edge = FALSE,
                   vertex.label.cex = 0.9,
                   title.name = paste0(g, " - Number of interactions"))
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                   weight.scale = TRUE, label.edge = FALSE,
                   vertex.label.cex = 0.9,
                   title.name = paste0(g, " - Interaction weights"))
  par(mfrow = c(1, 1))
  dev.off()
  
  # ==========================================================
  # 2. 和弦图 count
  # ==========================================================
  pdf(file.path(gdir, "02_chord_count.pdf"), width = 10, height = 10)
  netVisual_chord_cell(cellchat,
                       net = cellchat@net$count,
                       title.name = paste0(g, " - Chord (count)"),
                       legend.pos.x = 8)
  dev.off()
  
  # ==========================================================
  # 3. 和弦图 weight
  # ==========================================================
  pdf(file.path(gdir, "03_chord_weight.pdf"), width = 10, height = 10)
  netVisual_chord_cell(cellchat,
                       net = cellchat@net$weight,
                       title.name = paste0(g, " - Chord (weight)"),
                       legend.pos.x = 8)
  dev.off()
  
  # ==========================================================
  # 4. 热图 count
  # ==========================================================
  pdf(file.path(gdir, "04_heatmap_count.pdf"), width = 10, height = 8)
  print(netVisual_heatmap(cellchat, measure = "count",
                          color.heatmap = "Reds",
                          title.name = paste0(g, " - Interaction count")))
  dev.off()
  
  # ==========================================================
  # 5. 热图 weight
  # ==========================================================
  pdf(file.path(gdir, "05_heatmap_weight.pdf"), width = 10, height = 8)
  print(netVisual_heatmap(cellchat, measure = "weight",
                          color.heatmap = "Reds",
                          title.name = paste0(g, " - Interaction weight")))
  dev.off()
  
  # ==========================================================
  # 6. 原生 bubble 图：myeloid → T
  # ==========================================================
  if (length(my_types_in) >= 1 && length(t_types_in) >= 1) {
    p_bubble <- tryCatch(
      netVisual_bubble(cellchat,
                       sources.use = my_types_in,
                       targets.use = t_types_in,
                       font.size = 6,
                       remove.isolate = TRUE),
      error = function(e) NULL)
    if (!is.null(p_bubble)) {
      ggsave(file.path(gdir, "06_bubble_native.pdf"),
             plot = p_bubble, width = 16, height = 22, limitsize = FALSE)
    }
  }
  
  # ==========================================================
  # 7. 自定义 bubble（X = LR / X = cellpair）
  # ==========================================================
  df <- subsetCommunication(cellchat)
  if (nrow(df) > 0) {
    df$cellpair <- paste0(df$source, " → ", df$target)
    df$LR       <- paste0(df$ligand, " – ", df$receptor)
    
    df$pval_adj <- df$pval
    df$pval_adj[df$pval_adj == 0] <- min(df$pval[df$pval > 0], na.rm = TRUE) / 10
    df$size_val <- -log10(df$pval_adj)
    
    df_sub <- df %>% filter(source %in% my_types_in, target %in% t_types_in)
    
    if (nrow(df_sub) > 0) {
      df_top <- df_sub %>%
        group_by(source, target) %>%
        arrange(desc(prob)) %>%
        slice_head(n = 10) %>%
        ungroup()
      
      # X = LR
      p_custom <- ggplot(df_top, aes(x = LR, y = cellpair,
                                     size = size_val, color = prob)) +
        geom_point() +
        scale_color_gradientn(colors = c("#4575b4","#91bfdb","#ffffbf",
                                         "#fc8d59","#d73027")) +
        scale_size_continuous(range = c(0.5, 5)) +
        labs(title = paste0(g, ": myeloid → T"),
             x = "Ligand – Receptor", y = "Cell pair") +
        theme_bw(base_size = 14) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
              axis.text.y = element_text(size = 10))
      ggsave(file.path(gdir, "07_bubble_custom.pdf"),
             plot = p_custom, width = 24, height = 16, limitsize = FALSE)
      
      # X = cellpair
      p_swap <- ggplot(df_top, aes(x = cellpair, y = LR,
                                   size = size_val, color = prob)) +
        geom_point() +
        scale_color_gradientn(colors = c("#4575b4","#91bfdb","#ffffbf",
                                         "#fc8d59","#d73027")) +
        scale_size_continuous(range = c(0.5, 5)) +
        labs(title = paste0(g, ": myeloid → T (swapped)")) +
        theme_bw(base_size = 14) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
              axis.text.y = element_text(size = 8))
      ggsave(file.path(gdir, "08_bubble_custom_swapped.pdf"),
             plot = p_swap, width = 32, height = 22, limitsize = FALSE)
    }
  }
  
  # ==========================================================
  # 8. 信号通路信息流排名
  # ==========================================================
  pdf(file.path(gdir, "09_rankNet.pdf"), width = 10, height = 6)
  print(tryCatch(
    rankNet(cellchat, mode = "single",
            measure = "weight", do.stat = FALSE),
    error = function(e) NULL))
  dev.off()
  
  cat("  ", g, " 所有图已保存到：", normalizePath(gdir), "\n")
}

cat("\n\n四组图全部保存完成！\n")
cat("总目录：", normalizePath(outdir), "\n")

# ============================================================
# 合并四组 CellChat 对象
# ============================================================
cc_list <- lapply(groups, function(g) {
  f <- file.path(rds_dir, g, paste0("cellchat_", g, ".rds"))
  if (!file.exists(f)) f <- file.path(rds_dir, paste0("cellchat_", g, ".rds"))
  readRDS(f)
})
names(cc_list) <- groups

cellchat_merged <- mergeCellChat(cc_list, add.names = groups)

# ============================================================
# 汇总每组每个 celltype 的总交互数
# ============================================================
df_all <- do.call(rbind, lapply(groups, function(g) {
  m <- cellchat_merged@net[[g]]$count
  data.frame(group    = g,
             celltype = rownames(m),
             total    = rowSums(m))
}))

# ---------- 图 A：包含 dc ----------
target_A <- c("macrophage", "dc",
              "activated cd8", "activated cd4",
              "naive cd4", "mast", "treg")
df_A <- df_all %>% filter(celltype %in% target_A)
df_A$celltype <- factor(df_A$celltype, levels = target_A)

pA <- ggplot(df_A, aes(x = celltype, y = total, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.75, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("G1" = "#B0C4DE", "G2" = "#9ACD9A",
                               "G3" = "#F4A582", "G4" = "#D73027")) +
  labs(title = "Interaction counts of key cell types across G1–G4 (with dc)",
       x = "", y = "Total interactions (count)") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 13),
        plot.title = element_text(hjust = 0.5, face = "bold"))

print(pA)

# ---------- 图 B：不含 dc ----------
target_B <- c("macrophage",
              "activated cd8", "activated cd4",
              "naive cd4", "mast", "treg")
df_B <- df_all %>% filter(celltype %in% target_B)
df_B$celltype <- factor(df_B$celltype, levels = target_B)

pB <- ggplot(df_B, aes(x = celltype, y = total, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.75, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("G1" = "#B0C4DE", "G2" = "#9ACD9A",
                               "G3" = "#F4A582", "G4" = "#D73027")) +
  labs(title = "Interaction counts of key cell types across G1–G4 (without dc)",
       x = "", y = "Total interactions (count)") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 13),
        plot.title = element_text(hjust = 0.5, face = "bold"))

print(pB)


# ============================================================
# 提取 merged 对象的通讯结果
# ============================================================
df_list <- subsetCommunication(cellchat_merged)

df <- do.call(rbind, lapply(names(df_list), function(g) {
  d <- df_list[[g]]
  if (is.null(d) || nrow(d) == 0) return(NULL)
  d$group <- g
  d
}))

# ---------- 改名 ----------
df <- df %>%
  rename(
    Sender   = source,
    Receiver = target,
    LR_pair  = interaction_name_2
  ) %>%
  mutate(
    Cell_pair = paste(Sender, "→", Receiver)
  )

# ---------- 四组并排 ----------
df_wide <- df %>%
  group_by(Sender, Receiver, Cell_pair, LR_pair) %>%
  summarise(
    G1 = ifelse(any(group == "G1"), prob[group == "G1"][1], NA_real_),
    G2 = ifelse(any(group == "G2"), prob[group == "G2"][1], NA_real_),
    G3 = ifelse(any(group == "G3"), prob[group == "G3"][1], NA_real_),
    G4 = ifelse(any(group == "G4"), prob[group == "G4"][1], NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    mean_other = (G1 + G2 + G3) / 3,
    ratio_G4   = G4 / mean_other,
    G4_strict  = (G4 > G1) & (G4 > G2) & (G4 > G3)
  )

# ---------- 筛 G4 优势 LR ----------
g4_strict <- df_wide %>%
  filter(!is.na(G4), !is.na(G1), !is.na(G2), !is.na(G3)) %>%
  filter(G4_strict, ratio_G4 >= 1.2) %>%
  arrange(desc(ratio_G4))

cat("G4 优势 LR 数：", nrow(g4_strict), "\n")

# ---------- 取 Top 15 ----------
N <- 15
g4_top <- g4_strict %>% slice_head(n = N)

# ---------- 匹配回原始 df ----------
df_plot <- df %>%
  filter(LR_pair %in% g4_top$LR_pair,
         Cell_pair %in% g4_top$Cell_pair)

cat("df_plot 行数：", nrow(df_plot), "\n")

# ---------- 固定顺序 ----------
df_plot$LR_pair <- factor(df_plot$LR_pair,
                          levels = rev(unique(g4_top$LR_pair)))
df_plot$Cell_pair <- factor(df_plot$Cell_pair,
                            levels = unique(g4_top$Cell_pair))

# ---------- pval = 0 统一处理 ----------
df_plot <- df_plot %>%
  mutate(
    pval_adj = ifelse(is.na(pval) | pval == 0, 1e-5, pval),
    size_val = -log10(pval_adj)
  )

# ============================================================
# 四张单组图
# ============================================================
groups <- c("G1", "G2", "G3", "G4")

for (g in groups) {
  cat("\n>>> 画图：", g, "\n")
  
  d_sub <- df_plot %>% filter(group == g)
  cat("  行数：", nrow(d_sub), "\n")
  
  if (nrow(d_sub) == 0) {
    cat("  ", g, " 没有数据，跳过\n")
    next
  }
  
  p <- ggplot(d_sub, aes(x = Cell_pair, y = LR_pair,
                         size = size_val, color = prob)) +
    geom_point() +
    scale_color_gradientn(colors = c("#4575b4", "#91bfdb", "#ffffbf",
                                     "#fc8d59", "#d73027"),
                          name = "Commun. Prob.") +
    scale_size_continuous(name = "-log10(p)", range = c(2, 8)) +
    labs(title = g,
         x = "Cell pair (Sender → Receiver)",
         y = "Ligand–Receptor pair") +
    theme_bw(base_size = 13) +
    theme(
      axis.text.x  = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5,
                                  family = "Arial", face = "plain", color = "black"),
      axis.text.y  = element_text(size = 11, family = "Arial",
                                  face = "plain", color = "black"),
      axis.title.x = element_text(size = 12, family = "Arial", face = "plain"),
      axis.title.y = element_text(size = 12, family = "Arial", face = "plain"),
      plot.title   = element_text(size = 14, family = "Arial",
                                  face = "bold", hjust = 0.5),
      legend.text  = element_text(size = 10, family = "Arial"),
      legend.title = element_text(size = 10, family = "Arial"),
      panel.grid.minor = element_blank()
    )
  
  print(p)
  
  # ============================================================
  # 只画 dc → T 细胞的 bubble 图
  # ============================================================
  setwd("~/2026sc/important_data/rds/cellchat_results")
  
  library(CellChat)
  library(ggplot2)
  library(dplyr)
  
  rds_dir <- "~/2026sc/important_data/rds/cellchat_results"
  outdir  <- "~/2026sc/important_data/rds/cellchat_plots"
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  groups <- c("G1", "G2", "G3", "G4")
  
  # 只保留 dc 作为 source，T 细胞作为 target
  dc_type     <- "dc"
  tcell_types <- c("activated cd4", "Double-negative (DN) γδ T",
                   "exhausted cd8", "effector cd8", "treg",
                   "proliferating cd8", "cd4 Th2",
                   "activated cd8", "naive cd4")
  
  for (g in groups) {
    
    cat("\n============================================\n")
    cat("  分组：", g, "\n")
    cat("============================================\n")
    
    # ---------- 找 rds ----------
    cc_file <- file.path(rds_dir, g, paste0("cellchat_", g, ".rds"))
    if (!file.exists(cc_file)) {
      cc_file <- file.path(rds_dir, paste0("cellchat_", g, ".rds"))
    }
    if (!file.exists(cc_file)) {
      cat("找不到 ", g, " 的 rds，跳过\n"); next
    }
    
    cellchat <- readRDS(cc_file)
    
    # ---------- 输出目录 ----------
    gdir <- file.path(outdir, g)
    dir.create(gdir, showWarnings = FALSE, recursive = TRUE)
    
    # ---------- 实际存在的 T 细胞类型 ----------
    present_types <- levels(cellchat@idents)
    t_types_in    <- intersect(tcell_types, present_types)
    
    if (!(dc_type %in% present_types)) {
      cat("  ", g, " 中没有 dc，跳过\n"); next
    }
    if (length(t_types_in) == 0) {
      cat("  ", g, " 中没有 T 细胞，跳过\n"); next
    }
    
    # ---------- 画 dc → T 的 bubble 图 ----------
    p_bubble <- tryCatch(
      netVisual_bubble(cellchat,
                       sources.use = dc_type,
                       targets.use = t_types_in,
                       font.size   = 7,
                       remove.isolate = TRUE),
      error = function(e) NULL)
    
    if (is.null(p_bubble)) {
      cat("  ", g, " dc → T bubble 图失败，跳过\n"); next
    }
    
    ggsave(file.path(gdir, "dc_to_T_bubble.pdf"),
           plot = p_bubble, width = 10, height = 14, limitsize = FALSE)
    
    cat("  ", g, " dc → T bubble 图已保存：", gdir, "\n")
  }
  
  cat("\n全部完成！\n")
  
  # 合并四组
  cc_list <- lapply(groups, function(g) {
    f <- file.path(rds_dir, g, paste0("cellchat_", g, ".rds"))
    if (!file.exists(f)) f <- file.path(rds_dir, paste0("cellchat_", g, ".rds"))
    readRDS(f)
  })
  names(cc_list) <- groups
  cellchat_merged <- mergeCellChat(cc_list, add.names = groups)
  
  # 取实际存在的 T 细胞类型
  present_types <- levels(cc_list[[1]]@idents)
  t_types_in    <- intersect(tcell_types, present_types)
  
  # 四组对比 bubble（一张图）
  p <- netVisual_bubble(cellchat_merged,
                        sources.use = dc_type,
                        targets.use = t_types_in,
                        comparison    = c(1, 2, 3, 4),   # G1-G4
                        angle.x       = 45,
                        remove.isolate = TRUE)
  
  ggsave("dc_to_T_bubble_4groups.pdf",
         plot = p, width = 16, height = 18, limitsize = FALSE)
  
  library(CellChat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  
  # ============================================================
  # 前提：cellchat_merged 已存在
  # ============================================================
  
  rds_dir <- "~/2026sc/important_data/rds/cellchat_results"
  outdir  <- "~/2026sc/important_data/rds/cellchat_plots"
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  groups  <- c("G1", "G2", "G3", "G4")
  dc_type <- "dc"
  tcell_types <- c("activated cd4", "Double-negative (DN) γδ T",
                   "exhausted cd8", "effector cd8", "treg",
                   "proliferating cd8", "cd4 Th2",
                   "activated cd8", "naive cd4")
  
  # ============================================================
  # 1. 提取数据
  # ============================================================
  df_list <- subsetCommunication(cellchat_merged)
  
  df <- do.call(rbind, lapply(names(df_list), function(g) {
    d <- df_list[[g]]
    if (is.null(d) || nrow(d) == 0) return(NULL)
    d$group <- g
    d
  }))
  
  df_dcT <- df %>%
    filter(source == dc_type, target %in% tcell_types)
  
  # ============================================================
  # 2. 筛 G4 优势 LR
  # ============================================================
  df_wide <- df_dcT %>%
    group_by(target, ligand, receptor, interaction_name_2) %>%
    summarise(
      G3 = ifelse(any(group == "G3"), prob[group == "G3"][1], NA_real_),
      G4 = ifelse(any(group == "G4"), prob[group == "G4"][1], NA_real_),
      .groups = "drop"
    ) %>%
    filter(!is.na(G3), !is.na(G4)) %>%
    mutate(
      diff_G4_G3 = G4 - G3,
      ratio      = G4 / G3
    ) %>%
    filter(diff_G4_G3 > 0, ratio > 1.2) %>%
    arrange(desc(diff_G4_G3))
  
  N <- 15
  top_lr  <- df_wide %>% slice_head(n = N)
  
  # 去重 + 排除指定 LR + 强制保留 Ccl6-Ccr2
  lr_keep <- unique(paste0(top_lr$ligand, " - ", top_lr$receptor))
  lr_keep <- setdiff(lr_keep, c("Cd86 - Ctla4",
                                "Cd80 - Ctla4",
                                "Cd274 - Pdcd1"))
  lr_keep <- unique(c(lr_keep, "Ccl6 - Ccr2"))
  
  cat("保留的 LR：\n"); print(lr_keep)
  
  # ============================================================
  # 3. 四组各画一张，补全缺失 LR（空行）
  # ============================================================
  p_list <- list()
  
  # 所有可能的 cellpair（4 组共有的）
  all_cellpairs <- unique(paste0(df_dcT$source, " → ", df_dcT$target))
  
  for (g in groups) {
    
    # 该组实际数据
    d_sub <- df_dcT %>%
      filter(group == g,
             paste0(ligand, " - ", receptor) %in% lr_keep) %>%
      mutate(
        LR       = paste0(ligand, " - ", receptor),
        cellpair = paste0(source, " → ", target),
        pval_adj = ifelse(is.na(pval) | pval == 0, 1e-5, pval),
        size_val = -log10(pval_adj)
      ) %>%
      select(LR, cellpair, prob, size_val)
    
    # ★ 补全所有 LR × cellpair 组合（缺失的为 NA）
    d_full <- expand.grid(
      LR       = lr_keep,
      cellpair = all_cellpairs,
      stringsAsFactors = FALSE
    )
    
    d_plot <- d_full %>%
      left_join(d_sub, by = c("LR", "cellpair"))
    
    # 固定 Y 轴顺序
    d_plot$LR <- factor(d_plot$LR, levels = rev(lr_keep))
    
    # 画图
    p <- ggplot(d_plot, aes(x = cellpair, y = LR,
                            size = size_val, color = prob)) +
      geom_point(na.rm = TRUE) +   # ★ NA 不画点 → 空行
      scale_color_gradientn(
        colors = c("#4575b4", "#91bfdb", "#ffffbf", "#fc8d59", "#d73027"),
        name   = "Commun. Prob.",
        na.value = "transparent"
      ) +
      scale_size_continuous(name = "-log10(p)", range = c(2, 8)) +
      labs(title = g,
           x = "Cell pair", y = "Ligand – Receptor") +
      theme_bw(base_size = 13) +
      theme(
        axis.text.x  = element_text(size = 11, angle = 45, hjust = 1, color = "black"),
        axis.text.y  = element_text(size = 11, color = "black"),
        plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank()
      )
    
    # G1/G3 隐藏图例
    if (g %in% c("G1", "G3")) {
      p <- p + theme(legend.position = "none")
    }
    
    p_list[[g]] <- p
  }
  
  # ============================================================
  # 4. 拼成 2×2 大图
  # ============================================================
  combined <- (p_list$G1 | p_list$G2) / (p_list$G3 | p_list$G4) +
    plot_annotation(
      title = "dc → T : G4-advantage LR (G1–G4)",
      theme = theme(plot.title = element_text(size = 16,
                                              face = "bold",
                                              hjust = 0.5))
    )
  
  ggsave(file.path(outdir, "dc_to_T_G4advantage_4groups.pdf"),
         plot = combined, width = 20, height = 14, limitsize = FALSE)
  
  print(combined)
  

