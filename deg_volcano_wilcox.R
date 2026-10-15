library(Seurat)
library(ggplot2)
library(ggrepel)

Dendritic <- readRDS("dc_annotated01.rds")
dc_barcodes <- names(group)
Dendritic_tmp <- Dendritic[, dc_barcodes]
Dendritic_tmp$condition <- group

de_wilcox <- FindMarkers(
  object = Dendritic_tmp,
  group.by = "condition",
  ident.1 = "G4",
  ident.2 = "G3",
  test.use = "wilcox",
  min.pct = 0.05,
  logfc.threshold = 0
)

de_up_G4 <- subset(de_wilcox, p_val_adj < 0.05 & avg_log2FC > 0)
write.csv(de_up_G4, "G4_vs_G3_up_genes_wilcox.csv", row.names = TRUE)

mark_gene_lower <- c("il12a","ccl4","cxcl10")

de_wilcox$sig <- ifelse(de_wilcox$p_val_adj < 0.05 & de_wilcox$avg_log2FC > 0, "G4 up", "NS")
de_wilcox$gene_lower <- tolower(rownames(de_wilcox))
de_wilcox$label <- ifelse(de_wilcox$gene_lower %in% mark_gene_lower, toupper(de_wilcox$gene_lower), NA)

ggplot(de_wilcox, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(data = subset(de_wilcox, sig == "NS"), color = "grey70", alpha = 0.5, size = 0.8) +
  geom_point(data = subset(de_wilcox, sig == "G4 up"), color = "#E41A1C", alpha = 0.7, size = 1.2) +
  geom_point(data = subset(de_wilcox, gene_lower %in% mark_gene_lower), color = "black", size = 2.8) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_text_repel(
    data = subset(de_wilcox, gene_lower %in% mark_gene_lower),
    aes(label = label),
    size = 3.8,
    max.overlaps = 50,
    segment.size = 0.3,
    fontface = "bold"
  ) +
  # 固定坐标轴范围：x -4 ~ 4；y 0 ~ 50
  xlim(-4, 4) +
  ylim(0, 50) +
  labs(x = "log2FC (G4 / G3)", y = "-log10(adj.pval)", title = "G4 vs G3") +
  theme_bw()

ggsave("volcano_G4vsG3_right_only_x-4_4_y0_50.pdf", width = 7, height = 6, dpi = 300)

