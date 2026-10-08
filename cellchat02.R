# ============================================================
# 完整脚本：DC → T 通路分析（读数据 → 打分 → 画图 → IL12 伪通路）
# 包含 Part 0 - Part 5
# ============================================================

# ============================================================
# Part 0：环境与路径
# ============================================================
library(CellChat)
library(ggplot2)
library(dplyr)
library(patchwork)

rds_dir <- "~/2026sc/important_data/rds/cellchat_results"
outdir  <- "~/2026sc/important_data/rds/cellchat_plots"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

groups <- c("G1", "G2", "G3", "G4")

tcell_types <- c("activated cd4", "Double-negative (DN) γδ T",
                 "exhausted cd8", "effector cd8", "treg",
                 "proliferating cd8", "cd4 Th2",
                 "activated cd8", "naive cd4")
dc_type <- "dc"

# ============================================================
# Part 1：读入 4 组 rds + 构建 cellchat_merged
# ============================================================
cc_list <- list()
for (g in groups) {
  f1 <- file.path(rds_dir, g, paste0("cellchat_", g, ".rds"))
  f2 <- file.path(rds_dir, paste0("cellchat_", g, ".rds"))
  f  <- if (file.exists(f1)) f1 else if (file.exists(f2)) f2 else NA
  if (is.na(f)) { cat("找不到 ", g, " 的 rds\n"); next }
  cc_list[[g]] <- readRDS(f)
  cat("已读入：", g, "\n")
}

print(sapply(cc_list, function(x) length(levels(x@idents))))

cellchat_merged <- mergeCellChat(cc_list, add.names = groups)
cat("mergeCellChat 完成\n")
saveRDS(cellchat_merged, file.path(rds_dir, "cellchat_merged.rds"))

for (g in groups) {
  cat("\n>>> ", g, " 显著通路：\n", sep = "")
  print(cc_list[[g]]@netP$pathways)
}

# ============================================================
# Part 2：所有通路 dc → T 的打分表
# 输出：pathway_dc_to_T_all_groups.csv
#       pathway_dc_to_T_G4_advantage.csv
# ============================================================
all_pw <- sort(unique(unlist(lapply(cc_list, function(x) x@netP$pathways))))
cat("\n四组所有显著通路（共", length(all_pw), "个）：\n")
print(all_pw)

score_pw <- function(cellchat_obj, pw, dc, tcell) {
  if (!(pw %in% cellchat_obj@netP$pathways)) return(NA_real_)
  m <- cellchat_obj@netP$prob[, , pw]
  present <- rownames(m)
  src <- intersect(dc, present)
  tgt <- intersect(tcell, present)
  if (length(src) == 0 || length(tgt) == 0) return(0)
  sum(m[src, tgt, drop = FALSE])
}

res <- do.call(rbind, lapply(all_pw, function(pw) {
  vals <- sapply(groups, function(g) score_pw(cc_list[[g]], pw, dc_type, tcell_types))
  vals_num <- vals[!is.na(vals)]
  if (length(vals_num) == 0) return(NULL)
  data.frame(
    pathway    = pw,
    G1 = vals["G1"], G2 = vals["G2"], G3 = vals["G3"], G4 = vals["G4"],
    mean_other = mean(vals[c("G1","G2","G3")], na.rm = TRUE),
    ratio_G4   = vals["G4"] / mean(vals[c("G1","G2","G3")], na.rm = TRUE),
    G4_max     = vals["G4"] == max(vals_num, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

res <- res %>%
  filter(rowSums(cbind(G1, G2, G3, G4), na.rm = TRUE) > 0) %>%
  arrange(desc(ratio_G4))

g4_advantage <- res %>%
  filter(!is.na(ratio_G4), G4_max, ratio_G4 >= 1.2) %>%
  arrange(desc(ratio_G4))

cat("\n>>> 四组所有通路打分（dc → T）：\n")
print(res, row.names = FALSE)

cat("\n>>> G4 优势通路：\n")
print(g4_advantage, row.names = FALSE)

write.csv(res,
          file.path(rds_dir, "pathway_dc_to_T_all_groups.csv"),
          row.names = FALSE)
write.csv(g4_advantage,
          file.path(rds_dir, "pathway_dc_to_T_G4_advantage.csv"),
          row.names = FALSE)

# ============================================================
# Part 3：所有显著通路的 dc → T 圆图
# 每个通路一张 PDF，4 组 2×2
# ============================================================
pathways_show <- all_pw   # 若只画部分，可手动写 c("SPP1", ...)

for (pw in pathways_show) {
  
  cat("\n>>> 正在画通路：", pw, "\n")
  
  ok_groups <- groups[sapply(groups, function(g) pw %in% cc_list[[g]]@netP$pathways)]
  if (length(ok_groups) == 0) {
    cat("  ", pw, " 在四组中都不显著，跳过\n"); next
  }
  cat("  显著组：", paste(ok_groups, collapse = ", "), "\n")
  
  pdf(file.path(outdir, paste0("dc_to_T_pathway_", pw, "_4groups.pdf")),
      width = 20, height = 20)
  par(mfrow = c(2, 2), xpd = TRUE)
  
  for (g in groups) {
    
    if (!(pw %in% cc_list[[g]]@netP$pathways)) {
      plot.new()
      title(main = paste0(g, " : ", pw, " (not significant)"), cex.main = 1.5)
      next
    }
    
    prob_mat <- cc_list[[g]]@netP$prob[, , pw]
    present  <- rownames(prob_mat)
    keep_src <- intersect(dc_type, present)
    keep_tgt <- intersect(tcell_types, present)
    
    if (length(keep_src) == 0 || length(keep_tgt) == 0) {
      plot.new()
      title(main = paste0(g, " : ", pw, " (no dc/T)"), cex.main = 1.5)
      next
    }
    
    prob_sub <- prob_mat * 0
    prob_sub[keep_src, keep_tgt] <- prob_mat[keep_src, keep_tgt, drop = FALSE]
    
    groupSize <- as.numeric(table(cc_list[[g]]@idents))
    names(groupSize) <- levels(cc_list[[g]]@idents)
    vsize <- groupSize[present]
    
    netVisual_circle(prob_sub,
                     vertex.weight    = vsize,
                     weight.scale     = TRUE,
                     label.edge       = FALSE,
                     vertex.label.cex = 0.9,
                     title.name       = paste0(g, " : dc → T : ", pw))
  }
  
  par(mfrow = c(1, 1))
  dev.off()
  cat("  已保存：",
      file.path(outdir, paste0("dc_to_T_pathway_", pw, "_4groups.pdf")), "\n")
}

cat("\n全部 dc → T 通路图已保存到：", normalizePath(outdir), "\n")

# ============================================================
# Part 4：指定通路圆图（CD86、IL12、MHC-I、MHC-II、SPP1 等）
# ============================================================
pathways_show <- c(
  "CD86", "IL12", "MHC-I", "MHC-II", "SPP1",
  "ALCAM", "ITGAL-ITGB2", "CD6", "ICAM", "VCAM", "FN1",
  "CCL", "CXCL", "MIF",
  "TGFb", "PD-L1", "PDL2", "TIGIT", "IL10", "GALECTIN", "CD39"
)

for (pw in pathways_show) {
  
  cat("\n>>> 正在画通路：", pw, "\n")
  
  ok_groups <- groups[sapply(groups, function(g) pw %in% cc_list[[g]]@netP$pathways)]
  if (length(ok_groups) == 0) {
    cat("  ", pw, " 在四组中都不显著，跳过\n"); next
  }
  
  pdf(file.path(outdir, paste0("dc_to_T_pathway_", pw, "_4groups.pdf")),
      width = 20, height = 20)
  par(mfrow = c(2, 2), xpd = TRUE)
  
  for (g in groups) {
    if (!(pw %in% cc_list[[g]]@netP$pathways)) {
      plot.new()
      title(main = paste0(g, " : ", pw, " (not significant)"), cex.main = 1.5)
      next
    }
    prob_mat <- cc_list[[g]]@netP$prob[, , pw]
    present  <- rownames(prob_mat)
    keep_src <- intersect(dc_type, present)
    keep_tgt <- intersect(tcell_types, present)
    if (length(keep_src) == 0 || length(keep_tgt) == 0) {
      plot.new()
      title(main = paste0(g, " : ", pw, " (no dc/T)"), cex.main = 1.5)
      next
    }
    prob_sub <- prob_mat * 0
    prob_sub[keep_src, keep_tgt] <- prob_mat[keep_src, keep_tgt, drop = FALSE]
    groupSize <- as.numeric(table(cc_list[[g]]@idents))
    names(groupSize) <- levels(cc_list[[g]]@idents)
    vsize <- groupSize[present]
    netVisual_circle(prob_sub,
                     vertex.weight = vsize,
                     weight.scale = TRUE,
                     label.edge = FALSE,
                     vertex.label.cex = 0.9,
                     title.name = paste0(g, " : dc → T : ", pw))
  }
  par(mfrow = c(1, 1))
  dev.off()
  cat("  已保存：",
      file.path(outdir, paste0("dc_to_T_pathway_", pw, "_4groups.pdf")), "\n")
}

# ============================================================
# Part 5：IL12 伪通路图（CellChat 里不显著，手动构造）
# ============================================================
il12_pairs <- list(
  c("Il12a", "Il12rb1"),
  c("Il12b", "Il12rb1"),
  c("Il12b", "Il12rb2")
)

build_il12_matrix <- function(cellchat_obj, pairs, dc, tcell_types) {
  dat    <- cellchat_obj@data.signaling
  idents <- cellchat_obj@idents
  cts    <- levels(idents)
  pseudo <- matrix(0, nrow = length(cts), ncol = length(cts),
                   dimnames = list(cts, cts))
  for (p in pairs) {
    lig <- p[1]; rec <- p[2]
    if (!(lig %in% rownames(dat)) || !(rec %in% rownames(dat))) next
    avg_lig <- sapply(cts, function(ct) mean(dat[lig, idents == ct]))
    avg_rec <- sapply(cts, function(ct) mean(dat[rec, idents == ct]))
    m <- outer(avg_lig, avg_rec, "*")
    pseudo <- pmax(pseudo, m)
  }
  keep_src <- intersect(dc, cts)
  keep_tgt <- intersect(tcell_types, cts)
  out <- pseudo * 0
  out[keep_src, keep_tgt] <- pseudo[keep_src, keep_tgt, drop = FALSE]
  out
}

pdf(file.path(outdir, "dc_to_T_pathway_IL12_4groups.pdf"),
    width = 20, height = 20)
par(mfrow = c(2, 2), xpd = TRUE)

for (g in groups) {
  pseudo_sub <- build_il12_matrix(cc_list[[g]], il12_pairs,
                                  dc_type, tcell_types)
  if (sum(pseudo_sub) == 0) {
    plot.new()
    title(main = paste0(g, " : dc → T : IL12 (no expression)"),
          cex.main = 1.5)
    next
  }
  groupSize <- as.numeric(table(cc_list[[g]]@idents))
  names(groupSize) <- levels(cc_list[[g]]@idents)
  vsize <- groupSize[rownames(pseudo_sub)]
  netVisual_circle(pseudo_sub,
                   vertex.weight = vsize,
                   weight.scale = TRUE,
                   label.edge = FALSE,
                   vertex.label.cex = 0.9,
                   title.name = paste0(g, " : dc → T : IL12"))
}
par(mfrow = c(1, 1))
dev.off()

cat("\nIL12 图已保存：",
    file.path(outdir, "dc_to_T_pathway_IL12_4groups.pdf"), "\n")

# ============================================================
# 全部完成
# ============================================================
cat("\n===== 全部完成 =====\n")
cat("CSV 输出目录：", normalizePath(rds_dir), "\n")
cat("PDF 输出目录：", normalizePath(outdir), "\n")