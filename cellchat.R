# ============================================================
# CellChat 分析完整流程：myeloid + T 细胞（小鼠）
# 从读入数据到可视化，一个脚本跑到底
# ============================================================

# ============================================================
# 第 0 部分：加载包 + 全局设置
# ============================================================

library(Seurat)
library(CellChat)
library(patchwork)
library(future)
library(ggplot2)
library(dplyr)

setwd("~/2026sc/important_data/rds")

# Seurat 阶段关并行（避免 future 复制大对象报错）
plan("sequential")

# ============================================================
# 第 1 部分：数据准备（Seurat）
# ============================================================

# ---------- 1. 读入两个对象 ----------
myeloid <- readRDS("myeloid_annotated01.rds")
tcell   <- readRDS("t_annotated01.rds")

# ---------- 2. 把 Idents 写回 meta.data ----------
myeloid$celltype_real <- as.character(Idents(myeloid))
tcell$celltype_real   <- as.character(Idents(tcell))

table(myeloid$celltype_real)
table(tcell$celltype_real)

# ---------- 3. 合并 ----------
combined <- merge(myeloid, y = tcell)
ncol(combined)                 # 应为 51889
table(combined$celltype_real)  # 应为 15 种细胞类型

# ---------- 4. 标准化 ----------
DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(combined)

# ---------- 5. 合并 layer，提取表达矩阵 ----------
combined <- JoinLayers(combined)
data.input <- combined[["RNA"]]$data
ncol(data.input)               # 应为 51889

# ---------- 6. 构建 meta ----------
meta <- data.frame(
  labels = combined$celltype_real[match(colnames(data.input), colnames(combined))],
  row.names = colnames(data.input)
)
sum(is.na(meta$labels))        # 应为 0
table(meta$labels)             # 应为 15 种细胞类型

# ============================================================
# 第 2 部分：创建 CellChat 对象 + 小鼠数据库
# ============================================================

# ---------- 7. 创建 CellChat 对象 ----------
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

# ---------- 8. 设置小鼠数据库（关键！）----------
cellchat@DB <- CellChatDB.mouse

# ---------- 9. 检查基因名匹配率（非常重要）----------
db_genes <- unique(c(CellChatDB.mouse$interaction$ligand,
                     CellChatDB.mouse$interaction$receptor))
matched <- intersect(rownames(data.input), db_genes)
cat("匹配到的 LR 基因数：", length(matched), "\n")
# 期望：几千个。如果只有几十个，说明基因名格式不对，需要 str_to_title()

# ============================================================
# 第 3 部分：CellChat 预处理 + 计算通讯
# ============================================================

# ---------- 10. 预处理 ----------
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# ---------- 11. 开并行 + 计算通讯概率 ----------
options(future.globals.maxSize = 8 * 1024^3)
plan("multisession", workers = 4, maxSizeOfObjects = 8 * 1024^3)

cellchat <- computeCommunProb(cellchat, type = "triMean", nboot = 100)

# ---------- 12. 计算完立刻关并行 ----------
plan("sequential")

# ---------- 13. 过滤低表达细胞 ----------
cellchat <- filterCommunication(cellchat, min.cells = 10)

# ---------- 14. 按通路汇总 ----------
cellchat <- computeCommunProbPathway(cellchat)

# ---------- 15. 聚合网络 ----------
cellchat <- aggregateNet(cellchat)

# ---------- 16. 查看显著通路 ----------
cellchat@netP$pathways

# ---------- 17. 保存（关键！先存再画）----------
saveRDS(cellchat, file = "cellchat_myeloid_tcell_mouse.rds")

# ============================================================
# 第 4 部分：可视化
# ============================================================

# ============================================================
# 4.1 整体互作圆图
# ============================================================

groupSize <- as.numeric(table(cellchat@idents))

par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = TRUE, label.edge = FALSE,
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                 weight.scale = TRUE, label.edge = FALSE,
                 title.name = "Interaction weights/strength")
par(mfrow = c(1, 1))

# ============================================================
# 4.2 髓系 → T 气泡图（原生 netVisual_bubble）
# ============================================================

x11(width = 16, height = 22)   # Linux 用 x11；Windows 用 windows()；Mac 用 quartz()

netVisual_bubble(cellchat,
                 sources.use = c("macrophage", "dc", "neutrophil",
                                 "monocyte", "granulocyte", "mast"),
                 targets.use = c("activated cd4", "Double-negative (DN) γδ T",
                                 "exhausted cd8", "effector cd8", "treg",
                                 "proliferating cd8", "cd4 Th2",
                                 "activated cd8", "naive cd4"),
                 font.size = 6,
                 remove.isolate = TRUE)

# ============================================================
# 4.3 自定义 bubble 图（横纵坐标互换）
# ============================================================

# ---------- 18. 提取通讯数据 ----------
df <- subsetCommunication(cellchat)

df$cellpair <- paste0(df$source, " → ", df$target)
df$LR       <- paste0(df$ligand, " – ", df$receptor)

# ---------- 19. 处理 pval = 0（关键！否则点消失）----------
df$pval_adj <- df$pval
df$pval_adj[df$pval_adj == 0] <- min(df$pval[df$pval > 0], na.rm = TRUE) / 10
df$size_val <- -log10(df$pval_adj)

# ---------- 20. 筛选 myeloid → T ----------
df_sub <- df %>%
  filter(source %in% c("macrophage", "dc", "neutrophil",
                       "monocyte", "granulocyte", "mast"),
         target %in% c("activated cd4", "Double-negative (DN) γδ T",
                       "exhausted cd8", "effector cd8", "treg",
                       "proliferating cd8", "cd4 Th2",
                       "activated cd8", "naive cd4"))

cat("筛选后行数：", nrow(df_sub), "\n")

# 可选：每个细胞对只保留 Top 10 LR，避免太挤
df_top <- df_sub %>%
  group_by(source, target) %>%
  arrange(desc(prob)) %>%
  slice_head(n = 10) %>%
  ungroup()

# ---------- 21. 画图：X = LR，Y = 细胞对 ----------
x11(width = 20, height = 16)

ggplot(df_top, aes(x = LR, y = cellpair,
                   size = size_val, color = prob)) +
  geom_point() +
  scale_color_gradientn(colors = c("#4575b4", "#91bfdb", "#ffffbf",
                                   "#fc8d59", "#d73027"),
                        name = "Commun. Prob.") +
  scale_size_continuous(name = "-log10(p-value)", range = c(1, 6)) +
  labs(x = "Ligand – Receptor", y = "Cell pair (source → target)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        axis.text.y = element_text(size = 8),
        panel.grid.minor = element_blank())

# ============================================================
# 完成
# ============================================================
cat("\n✅ CellChat 全流程完成！\n")
cat("结果已保存到：cellchat_myeloid_tcell_mouse.rds\n")


# 不要 str_wrap，直接用 df_top$LR
pdf("bubble_custom_big.pdf", width = 24, height = 16)
ggplot(df_top, aes(x = LR, y = cellpair,
                   size = size_val, color = prob)) +
  geom_point() +
  scale_color_gradientn(colors = c("#4575b4", "#91bfdb", "#ffffbf",
                                   "#fc8d59", "#d73027")) +
  scale_size_continuous(range = c(0.5, 5)) +
  theme_bw(base_size = 19) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 11),   # 字体小一点
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 17))
dev.off()

#bubble_custom_swapped
p <- ggplot(df_top, aes(x = cellpair, y = LR,
                        size = size_val, color = prob)) +
  geom_point() +
  scale_color_gradientn(colors = c("#4575b4", "#91bfdb", "#ffffbf",
                                   "#fc8d59", "#d73027")) +
  scale_size_continuous(range = c(0.5, 5)) +
  theme_bw(base_size = 19) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 13),
        axis.text.y = element_text(size = 13))

ggsave("bubble_custom_swapped.pdf", plot = p,
       width = 32, height = 22, units = "in")



# ============================================================
# CellChat 圆图：全名标签 + 字体大一点（RStudio 直接显示）
# ============================================================

library(CellChat)

# ---------- 1. 确认 idents 和 net 都正常 ----------
levels(cellchat@idents)
dimnames(cellchat@net$count)[[1]]
table(cellchat@idents)          # 必须非零！如果全 0 先修复 idents

# ---------- 2. 节点大小 ----------
groupSize <- as.numeric(table(cellchat@idents))

# ---------- 3. 画圆图（RStudio 直接显示）----------
par(mfrow = c(1, 2), xpd = TRUE)

netVisual_circle(cellchat@net$count,
                 vertex.weight = groupSize,
                 weight.scale = TRUE,
                 label.edge = FALSE,
                 vertex.label.cex = 0.9,       # ← 从 0.7 调到 0.9
                 title.name = "Number of interactions")

netVisual_circle(cellchat@net$weight,
                 vertex.weight = groupSize,
                 weight.scale = TRUE,
                 label.edge = FALSE,
                 vertex.label.cex = 0.9,       # ← 同上
                 title.name = "Interaction weights/strength")

par(mfrow = c(1, 1))


#pathway
# 圆图：某个通路在所有细胞间的通讯
netVisual_aggregate(cellchat, signaling = "TGFb", layout = "circle")

# 和弦图：更清晰展示细胞对
netVisual_aggregate(cellchat, signaling = "TGFb", layout = "chord")

# 只画 myeloid → T 子网络
netVisual_aggregate(cellchat, signaling = "TGFb", layout = "circle",
                    sources.use = c("macrophage", "dc", "neutrophil",
                                    "monocyte", "granulocyte", "mast"),
                    targets.use = c("activated cd4", "Double-negative (DN) γδ T",
                                    "exhausted cd8", "effector cd8", "treg",
                                    "proliferating cd8", "cd4 Th2",
                                    "activated cd8", "naive cd4"))



netVisual_individual(cellchat, signaling = "TGFb", layout = "circle")

netVisual_heatmap(cellchat, signaling = "TGFb")