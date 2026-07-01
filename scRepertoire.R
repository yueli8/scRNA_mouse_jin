library(tibble)
library(Seurat)
library(scRepertoire)
library(dplyr)
library(tidyr)

setwd("~/2026sc/scRepertoire")

# ===============================
# 1. 读取数据
# ===============================
G2_63 <- read.csv("63_filtered_contig_annotations.csv")
G3_68 <- read.csv("68_filtered_contig_annotations.csv")
G1_71 <- read.csv("71_filtered_contig_annotations.csv")
G2_73 <- read.csv("73_filtered_contig_annotations.csv")
G2_79 <- read.csv("79_filtered_contig_annotations.csv")
G3_80 <- read.csv("80_filtered_contig_annotations.csv")
G4_83 <- read.csv("83_filtered_contig_annotations.csv")
G1_84 <- read.csv("84_filtered_contig_annotations.csv")
G1_87 <- read.csv("87_filtered_contig_annotations.csv")
G3_88 <- read.csv("88_filtered_contig_annotations.csv")
G4_90 <- read.csv("90_filtered_contig_annotations.csv")
G4_92 <- read.csv("92_filtered_contig_annotations.csv")

contig_list <- list(G1_71, G1_84, G1_87, G2_63, G2_73, G2_79, 
                    G3_68, G3_80, G3_88, G4_83, G4_90, G4_92)

combined <- combineTCR(contig_list, 
                       samples = c("G1_71","G1_84","G1_87","G2_63","G2_73","G2_79",
                                   "G3_68","G3_80","G3_88","G4_83","G4_90","G4_92"))

clonalQuant(combined, cloneCall = "strict", scale = TRUE)

# ===============================
# 2. 合并原始 contig，添加 Sample 标签和 barcode_full
# ===============================
samp_order <- c("G1_71","G1_84","G1_87","G2_63","G2_73","G2_79",
                "G3_68","G3_80","G3_88","G4_83","G4_90","G4_92")

raw_list <- list(G1_71, G1_84, G1_87, G2_63, G2_73, G2_79,
                 G3_68, G3_80, G3_88, G4_83, G4_90, G4_92)

for(i in seq_along(raw_list)) {
  raw_list[[i]]$Sample <- samp_order[i]
}

contig_raw <- bind_rows(raw_list) %>%
  mutate(barcode_full = paste0(Sample, "_", barcode))

# ===============================
# 3. 构建克隆信息表（包含 FWR/CDR 氨基酸序列）
# ===============================

# 第一步：过滤 productive 链，保留所有需要的列（包括框架区和CDR区氨基酸）
clone_info <- contig_raw %>%
  filter(productive == "true", chain %in% c("TRA","TRB")) %>%
  select(Sample, barcode, barcode_full, chain, 
         v_gene, j_gene, c_gene, cdr3,  # 原有信息
         fwr1, cdr1, fwr2, cdr2, fwr3, cdr3, fwr4,  # 氨基酸序列（注意cdr3重复，保留一个）
         reads)

# 第二步：处理可能的多条链（每个细胞每个链只保留 reads 最高的一条）
clone_info <- clone_info %>%
  group_by(Sample, barcode, chain) %>%
  arrange(desc(reads)) %>%
  slice(1) %>%
  ungroup()

# 第三步：宽表转换（每个细胞一行，TRA和TRB分开列）
clone_info_wide <- clone_info %>%
  select(Sample, barcode, barcode_full, chain,
         v_gene, j_gene, c_gene, cdr3,
         fwr1, cdr1, fwr2, cdr2, fwr3, fwr4) %>%
  pivot_wider(
    names_from = chain,
    values_from = c(v_gene, j_gene, c_gene, cdr3,
                    fwr1, cdr1, fwr2, cdr2, fwr3, fwr4),
    values_fill = list(v_gene = NA, j_gene = NA, c_gene = NA, cdr3 = NA,
                       fwr1 = NA, cdr1 = NA, fwr2 = NA, cdr2 = NA, fwr3 = NA, fwr4 = NA)
  ) %>%
  rename(
    # V/J/C 基因
    TRAV = v_gene_TRA, TRBV = v_gene_TRB,
    Ja = j_gene_TRA, Jb = j_gene_TRB,
    TRAC = c_gene_TRA, TRBC = c_gene_TRB,
    # CDR3（来自原始列）
    CDR3a = cdr3_TRA, CDR3b = cdr3_TRB,
    # TRA 框架区和CDR区
    TRA_fwr1 = fwr1_TRA, TRA_cdr1 = cdr1_TRA,
    TRA_fwr2 = fwr2_TRA, TRA_cdr2 = cdr2_TRA,
    TRA_fwr3 = fwr3_TRA, TRA_fwr4 = fwr4_TRA,
    # TRB 框架区和CDR区
    TRB_fwr1 = fwr1_TRB, TRB_cdr1 = cdr1_TRB,
    TRB_fwr2 = fwr2_TRB, TRB_cdr2 = cdr2_TRB,
    TRB_fwr3 = fwr3_TRB, TRB_fwr4 = fwr4_TRB
  ) %>%
  # 只保留同时有 TRA 和 TRB 的细胞
  filter(!is.na(TRAV) & !is.na(TRBV)) %>%
  mutate(CDR3_comb = paste(CDR3a, CDR3b, sep = "_"))

# 第四步：统计每个克隆型的细胞数（Count）
clone_count <- clone_info_wide %>%
  group_by(Sample, CDR3_comb) %>%
  summarise(Count = n(), .groups = "drop")

# 第五步：合并 Count 信息
clone_info_final <- left_join(clone_info_wide, clone_count, by = c("Sample", "CDR3_comb"))

# 第六步：处理多个 V/J 基因的情况（用分号合并），包括框架区
clone_info_fix <- clone_info_final %>%
  mutate(across(
    c(TRAV, TRBV, Ja, Jb, TRAC, TRBC,
      CDR3a, CDR3b, CDR3_comb,
      TRA_fwr1, TRA_cdr1, TRA_fwr2, TRA_cdr2, TRA_fwr3, TRA_fwr4,
      TRB_fwr1, TRB_cdr1, TRB_fwr2, TRB_cdr2, TRB_fwr3, TRB_fwr4),
    ~sapply(.x, function(z) paste0(z, collapse = ";"))
  ))

# 第七步：导出（包含所有框架区和CDR区）
write.csv(clone_info_fix, "~/2026sc/clone_info_final_fix.csv", row.names = FALSE)

# 验证：检查 G1_71 的 577 克隆是否包含框架区信息
clone_info_fix %>%
  filter(Sample == "G1_71", Count == 577) %>%
  select(Sample, barcode_full, CDR3_comb, Count, 
         TRA_fwr1, TRA_cdr1, TRA_fwr2, TRA_cdr2, TRA_fwr3, TRA_fwr4,
         TRB_fwr1, TRB_cdr1, TRB_fwr2, TRB_cdr2, TRB_fwr3, TRB_fwr4) %>%
  head()



# ===============================
# 从 clone_info_final_fix.csv 提取湿实验候选
# ===============================

# 读取数据（如果还没有在内存中）
clone_info_fix <- read.csv("~/2026sc/clone_info_final_fix.csv")


# ===============================
# 4. 导出克隆匹配键（细胞级 + 克隆型级）
# ===============================
clone_key <- clone_info_fix %>%
  select(Sample, CDR3_comb, TRAV, TRBV, CDR3a, CDR3b, Count, Ja, Jb)
write.csv(clone_key, "~/2026sc/clone_matching_key.csv", row.names = FALSE)

clone_key_unique <- clone_info_fix %>%
  distinct(Sample, CDR3_comb, .keep_all = TRUE) %>%
  select(Sample, CDR3_comb, TRAV, TRBV, CDR3a, CDR3b, Count, Ja, Jb)
write.csv(clone_key_unique, "~/2026sc/clone_matching_key_unique.csv", row.names = FALSE)

# ===============================
# 5. vizGenes (TRAV)
# ===============================
viz_prop <- vizGenes(combined, x.axis = "TRAV", summary.fun = "proportion", 
                     order = "variance", exportTable = TRUE)
write.csv(viz_prop, "~/2026sc/vizTRAV_prop.csv", row.names = FALSE)

viz_count <- vizGenes(combined, x.axis = "TRAV", summary.fun = "count", 
                      order = "variance", exportTable = TRUE)
write.csv(viz_count, "~/2026sc/vizTRAV_count.csv", row.names = FALSE)

# ===============================
# 6. TRAV 聚合统计 + 注释表
# ===============================
trav_aggregated <- clone_info_fix %>%
  separate_rows(TRAV, sep = ";") %>%
  group_by(Sample, TRAV) %>%
  summarise(
    n_cells = n(),
    n_clones = n_distinct(CDR3_comb),
    avg_clone_size = mean(Count, na.rm = TRUE),
    max_clone_size = max(Count, na.rm = TRUE),
    top_CDR3a = paste(head(sort(table(CDR3a), decreasing = TRUE), 3) %>% names(), collapse = ";"),
    top_CDR3b = paste(head(sort(table(CDR3b), decreasing = TRUE), 3) %>% names(), collapse = ";"),
    top_CDR3_comb = paste(head(sort(table(CDR3_comb), decreasing = TRUE), 3) %>% names(), collapse = ";"),
    top_Counts = paste(head(sort(table(CDR3_comb), decreasing = TRUE), 3), collapse = ";"),
    .groups = "drop"
  )

viz_prop_long <- as.data.frame(viz_prop) %>%
  mutate(TRAV = rownames(.)) %>%
  pivot_longer(cols = -TRAV, names_to = "Sample", values_to = "proportion")

viz_prop_anno <- left_join(viz_prop_long, trav_aggregated, by = c("TRAV", "Sample"))
write.csv(viz_prop_anno, "~/2026sc/vizTRAV_prop_anno.csv", row.names = FALSE)

viz_count_long <- as.data.frame(viz_count) %>%
  mutate(TRAV = rownames(.)) %>%
  pivot_longer(cols = -TRAV, names_to = "Sample", values_to = "count")

viz_count_anno <- left_join(viz_count_long, trav_aggregated, by = c("TRAV", "Sample"))
write.csv(viz_count_anno, "~/2026sc/vizTRAV_count_anno.csv", row.names = FALSE)

# ===============================
# 7. TRBV 聚合统计 + 注释表
# ===============================
trbv_aggregated <- clone_info_fix %>%
  separate_rows(TRBV, sep = ";") %>%
  group_by(Sample, TRBV) %>%
  summarise(
    n_cells = n(),
    n_clones = n_distinct(CDR3_comb),
    avg_clone_size = mean(Count, na.rm = TRUE),
    max_clone_size = max(Count, na.rm = TRUE),
    top_CDR3a = paste(head(sort(table(CDR3a), decreasing = TRUE), 3) %>% names(), collapse = ";"),
    top_CDR3b = paste(head(sort(table(CDR3b), decreasing = TRUE), 3) %>% names(), collapse = ";"),
    top_CDR3_comb = paste(head(sort(table(CDR3_comb), decreasing = TRUE), 3) %>% names(), collapse = ";"),
    top_Counts = paste(head(sort(table(CDR3_comb), decreasing = TRUE), 3), collapse = ";"),
    .groups = "drop"
  )

viz_trbv_prop <- vizGenes(combined, x.axis = "TRBV", summary.fun = "proportion", exportTable = TRUE)
viz_trbv_prop_long <- as.data.frame(viz_trbv_prop) %>%
  mutate(TRBV = rownames(.)) %>%
  pivot_longer(cols = -TRBV, names_to = "Sample", values_to = "proportion")
viz_trbv_prop_anno <- left_join(viz_trbv_prop_long, trbv_aggregated, by = c("TRBV", "Sample"))
write.csv(viz_trbv_prop_anno, "~/2026sc/vizTRBV_prop_anno.csv", row.names = FALSE)

viz_trbv_count <- vizGenes(combined, x.axis = "TRBV", summary.fun = "count", exportTable = TRUE)
viz_trbv_count_long <- as.data.frame(viz_trbv_count) %>%
  mutate(TRBV = rownames(.)) %>%
  pivot_longer(cols = -TRBV, names_to = "Sample", values_to = "count")
viz_trbv_count_anno <- left_join(viz_trbv_count_long, trbv_aggregated, by = c("TRBV", "Sample"))
write.csv(viz_trbv_count_anno, "~/2026sc/vizTRBV_count_anno.csv", row.names = FALSE)

# ===============================
# 8. 克隆型级别完整表（最终推荐）
# ===============================
clone_type_summary <- clone_info_fix %>%
  group_by(Sample, CDR3_comb) %>%
  summarise(
    TRAV = first(TRAV),
    TRBV = first(TRBV),
    CDR3a = first(CDR3a),
    CDR3b = first(CDR3b),
    Ja = first(Ja),
    Jb = first(Jb),
    Count = first(Count),
    .groups = "drop"
  ) %>%
  arrange(Sample, desc(Count))

write.csv(clone_type_summary, "~/2026sc/clone_type_summary.csv", row.names = FALSE)

# ===============================
# 9. （可选）Top 克隆 per TRAV / TRBV
# ===============================
trav_top_clones <- clone_info_fix %>%
  separate_rows(TRAV, sep = ";") %>%
  group_by(Sample, TRAV) %>%
  arrange(desc(Count), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  select(Sample, TRAV, barcode, CDR3_comb, CDR3a, CDR3b, Count)
write.csv(trav_top_clones, "~/2026sc/TRAV_top_clones.csv", row.names = FALSE)

trbv_top_clones <- clone_info_fix %>%
  separate_rows(TRBV, sep = ";") %>%
  group_by(Sample, TRBV) %>%
  arrange(desc(Count), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  select(Sample, TRBV, barcode, CDR3_comb, CDR3a, CDR3b, Count)
write.csv(trbv_top_clones, "~/2026sc/TRBV_top_clones.csv", row.names = FALSE)

# ===============================
# 1. 克隆型级别汇总（每个克隆型一行）
# ===============================
clone_type_for_wetlab <- clone_info_fix %>%
  # 按克隆型分组
  group_by(Sample, CDR3_comb) %>%
  summarise(
    # 基本信息（取第一条，同一克隆型下相同）
    TRAV = first(TRAV),
    TRBV = first(TRBV),
    CDR3a = first(CDR3a),
    CDR3b = first(CDR3b),
    Ja = first(Ja),
    Jb = first(Jb),
    Count = first(Count),
    
    # 框架区和CDR区（取第一条）
    TRA_fwr1 = first(TRA_fwr1),
    TRA_cdr1 = first(TRA_cdr1),
    TRA_fwr2 = first(TRA_fwr2),
    TRA_cdr2 = first(TRA_cdr2),
    TRA_fwr3 = first(TRA_fwr3),
    TRA_fwr4 = first(TRA_fwr4),
    TRB_fwr1 = first(TRB_fwr1),
    TRB_cdr1 = first(TRB_cdr1),
    TRB_fwr2 = first(TRB_fwr2),
    TRB_cdr2 = first(TRB_cdr2),
    TRB_fwr3 = first(TRB_fwr3),
    TRB_fwr4 = first(TRB_fwr4),
    
    # 细胞列表（用于追溯）
    cell_barcodes = paste(barcode_full, collapse = ";"),
    n_cells = n(),
    
    .groups = "drop"
  ) %>%
  # 按样本和克隆大小排序
  arrange(Sample, desc(Count))

# ===============================
# 2. 筛选条件
# ===============================

# 条件说明：
# - Count > 1：有扩增的克隆（排除单细胞克隆）
# - 可以根据需要调整阈值（如 Count >= 5, >= 10 等）

wetlab_candidates <- clone_type_for_wetlab %>%
  filter(Count > 1) %>%
  # 按 Count 降序排列（最大的克隆在前）
  arrange(desc(Count))

# ===============================
# 3. 添加优先级分级
# ===============================

wetlab_candidates <- wetlab_candidates %>%
  mutate(
    Priority = case_when(
      Count >= 50 ~ "A_最高优先（显著扩增）",
      Count >= 20 ~ "B_高优先（中等扩增）",
      Count >= 10 ~ "C_中优先（低度扩增）",
      Count >= 5  ~ "D_低优先（轻微扩增）",
      Count > 1   ~ "E_备选（极小扩增）",
      TRUE        ~ "F_单细胞（不推荐）"
    ),
    
    # 湿实验可行性标记
    WetLab_Ready = case_when(
      # 同时有完整的 TRA 和 TRB 信息
      !is.na(TRAV) & !is.na(TRBV) & 
        !is.na(CDR3a) & !is.na(CDR3b) & 
        nchar(CDR3a) > 5 & nchar(CDR3b) > 5 ~ "Ready",
      
      # 缺少某些信息
      is.na(TRAV) | is.na(TRBV) ~ "Missing_V_gene",
      is.na(CDR3a) | is.na(CDR3b) ~ "Missing_CDR3",
      TRUE ~ "Check_Manually"
    ),
    
    # TRBV 复杂度标记（处理双V基因）
    TRBV_Complexity = case_when(
      grepl("\\+", TRBV) ~ "Multi_Allele（需选择主要V基因）",
      TRUE ~ "Single"
    )
  )

# ===============================
# 4. 导出结果
# ===============================

# 4.1 所有 Count > 1 的克隆型（完整版）
write.csv(wetlab_candidates, "~/2026sc/WetLab_Candidates_All_CountGT1.csv", row.names = FALSE)

# 4.2 按优先级分文件导出
for(priority in c("A_最高优先（显著扩增）", "B_高优先（中等扩增）", 
                  "C_中优先（低度扩增）", "D_低优先（轻微扩增）", 
                  "E_备选（极小扩增）")) {
  
  subset_data <- wetlab_candidates %>% filter(Priority == priority)
  if(nrow(subset_data) > 0) {
    filename <- paste0("~/2026sc/WetLab_", gsub(" ", "_", priority), ".csv")
    write.csv(subset_data, filename, row.names = FALSE)
  }
}

# 4.3 仅湿实验就绪的克隆（完整信息）
wetlab_ready <- wetlab_candidates %>%
  filter(WetLab_Ready == "Ready")

write.csv(wetlab_ready, "~/2026sc/WetLab_Ready_Clones.csv", row.names = FALSE)

# 4.4 处理双TRBV的克隆（需要额外选择）
multi_trbv <- wetlab_candidates %>%
  filter(TRBV_Complexity == "Multi_Allele")

write.csv(multi_trbv, "~/2026sc/WetLab_Multi_TRBV_NeedsReview.csv", row.names = FALSE)

# 4.5 按样本分组导出
for(sample in unique(wetlab_candidates$Sample)) {
  sample_data <- wetlab_candidates %>%
    filter(Sample == sample) %>%
    arrange(desc(Count))
  
  if(nrow(sample_data) > 0) {
    filename <- paste0("~/2026sc/WetLab_", sample, ".csv")
    write.csv(sample_data, filename, row.names = FALSE)
  }
}

# ===============================
# 5. 统计报告
# ===============================

cat("\n========== 湿实验候选统计 ==========\n")
cat("总克隆型数:", nrow(clone_type_for_wetlab), "\n")
cat("Count > 1 的克隆型数:", nrow(wetlab_candidates), "\n")
cat("湿实验就绪克隆数:", nrow(wetlab_ready), "\n\n")

cat("按优先级分布:\n")
print(table(wetlab_candidates$Priority))

cat("\n按样本分布:\n")
print(table(wetlab_candidates$Sample))

cat("\n按 TRBV 复杂度分布:\n")
print(table(wetlab_candidates$TRBV_Complexity))

# 预览前20个最高优先级的克隆
cat("\n========== Top 20 最高优先级克隆 ==========\n")
print(wetlab_candidates %>%
        select(Sample, CDR3_comb, Count, Priority, WetLab_Ready, TRBV_Complexity) %>%
        head(20))


# 最小化版本：只导出核心表格
# 1-4: G1-G4 就绪克隆（湿实验直接用）
# 5: 组间汇总统计

sample_groups <- split(
  unique(wetlab_candidates$Sample),
  sub("_.*$", "", unique(wetlab_candidates$Sample))
)
group_name <-
for(group_name in c("G1", "G2", "G3", "G4")) {
  group_samples <- sample_groups[[group_name]]
  group_ready <- wetlab_candidates %>%
    filter(Sample %in% group_samples, WetLab_Ready == "Ready") %>%
    arrange(desc(Count))
  
  write.csv(group_ready, paste0("~/2026sc/WetLab_", group_name, "_Ready.csv"), row.names = FALSE)
}


group_summary <- wetlab_candidates %>%
  mutate(Group = sub("_.*$", "", Sample)) %>%
  group_by(Group) %>%
  summarise(
    Total_Clonotypes  = n(),
    Ready_For_WetLab  = sum(WetLab_Ready == "Ready", na.rm = TRUE),
    Max_Clone_Size    = max(n_cells, na.rm = TRUE),
    Mean_Clone_Size   = mean(n_cells, na.rm = TRUE),
    Median_Clone_Size = median(n_cells, na.rm = TRUE),
    Clones_GE50       = sum(n_cells >= 50, na.rm = TRUE),
    Clones_GE20       = sum(n_cells >= 20, na.rm = TRUE),
    Clones_GE10       = sum(n_cells >= 10, na.rm = TRUE),
    Multi_TRBV_Count  = sum(TRBV_Complexity != "Single", na.rm = TRUE),
    .groups = "drop"
  )

write.csv(group_summary, "~/2026sc/WetLab_Group_Summary.csv", row.names = FALSE)


seurat <- readRDS("~/2026sc/total_105cluster.rds")

# 直接从细胞级别数据提取 577 克隆的 barcode
cell_barcodes_all <- clone_info_fix %>%
  filter(Sample == "G1_71", Count == 577) %>%
  pull(barcode_full) %>%
  unique()

# 查看结果
length(cell_barcodes_all)
head(cell_barcodes_all, 10)

# 查询这些细胞在 Seurat 中的细胞类型分布
celltype_distribution <- seurat@meta.data %>%
  rownames_to_column("barcode") %>%
  filter(barcode %in% cell_barcodes_all) %>%
  group_by(celltype) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

# 查看结果
print(celltype_distribution)

# 计算匹配率
matched_count <- sum(cell_barcodes_all %in% rownames(seurat@meta.data))
cat("\n匹配率:", matched_count, "/", length(cell_barcodes_all), 
    "(", round(matched_count/length(cell_barcodes_all)*100, 2), "%)\n")


# 查看分布情况
seurat@meta.data %>%
  rownames_to_column("barcode") %>%
  filter(barcode %in% cell_barcodes_all) %>%
  summarise(
    nCount_RNA_min = min(nCount_RNA),
    nCount_RNA_max = max(nCount_RNA),
    nCount_RNA_median = median(nCount_RNA),
    nFeature_RNA_min = min(nFeature_RNA),
    nFeature_RNA_max = max(nFeature_RNA),
    nFeature_RNA_median = median(nFeature_RNA),
    percent_mito_min = min(percent.mt),
    percent_mito_max = max(percent.mt),
    percent_mito_median = median(percent.mt)
  )

# 方法2：直接提取数据矩阵
expr_matrix <- seurat[["RNA"]]$data
avg_expression <- rowMeans(expr_matrix[, cell_barcodes_seurat])

top100_genes <- names(sort(avg_expression, decreasing = TRUE))[1:100]
print(top100_genes[1:100])

# 精简版：只保留湿实验最需要的信息
A_priority_clones <- read.csv("~/2026sc/WetLab_A_最高优先（显著扩增）.csv")
A_priority_simplified <- A_priority_clones %>%
  select(
    Sample,           # 样本
    CDR3_comb,        # 克隆型ID
    Count,            # 克隆大小
    cell_barcodes,    # 细胞列表（可用于追溯）
    TRAV, TRBV,       # V基因
    CDR3a, CDR3b,     # CDR3序列
    Ja, Jb,           # J基因
    TRA_fwr1, TRA_cdr1, TRA_fwr2, TRA_cdr2, TRA_fwr3, TRA_fwr4,  # TRA 框架区
    TRB_fwr1, TRB_cdr1, TRB_fwr2, TRB_cdr2, TRB_fwr3, TRB_fwr4,  # TRB 框架区
    TRBV_Complexity,  # 双V基因标记
    WetLab_Ready      # 湿实验就绪状态
  ) %>%
  arrange(Sample, desc(Count))

write.csv(A_priority_simplified, "~/2026sc/A_最高优先_显著扩增克隆_精简版.csv", row.names = FALSE)

head(A_priority_simplified)

# 1. 加载完整的 Seurat 对象
#seurat <- readRDS("total_105cluster.rds")
# 读取之前保存的 clone_info_fix
clone_info_fix <- read.csv("~/2026sc/clone_info_final_fix.csv")


# 1. 提取各图层的表达矩阵（保持稀疏格式）
expr_matrices <- list()
expr_matrices[["data.1"]] <- GetAssayData(seurat, assay = "RNA", layer = "data.1")
expr_matrices[["data.2"]] <- GetAssayData(seurat, assay = "RNA", layer = "data.2")
expr_matrices[["data.3"]] <- GetAssayData(seurat, assay = "RNA", layer = "data.3")
expr_matrices[["data.4"]] <- GetAssayData(seurat, assay = "RNA", layer = "data.4")
cat("图层提取完成\n")

# 2. 获取所有 A_优先克隆的 barcode
clone_barcodes <- clone_info_fix %>%
  filter(Count >= 50) %>%
  group_by(Sample, CDR3_comb, Count) %>%
  summarise(
    cell_barcodes = paste(barcode_full, collapse = ";"),
    n_cells = n(),
    .groups = "drop"
  )

cat("A_优先克隆数:", nrow(clone_barcodes), "\n")

# 3. 获取所有 Seurat 细胞
all_cells <- rownames(seurat@meta.data)

# 4. 存储结果
successful_clones <- list()

# 5. 处理每个克隆
for(i in 1:nrow(clone_barcodes)) {
  clone_info <- clone_barcodes[i, ]
  clone_name <- clone_info$CDR3_comb
  sample_name <- clone_info$Sample
  clone_size <- clone_info$Count
  
  # 每10个克隆打印一次进度
  if(i %% 5 == 1) {
    cat("\n处理克隆", i, "/", nrow(clone_barcodes), ":", sample_name, "-", substr(clone_name, 1, 40), "...\n")
  }
  
  barcodes <- strsplit(clone_info$cell_barcodes, ";")[[1]]
  barcodes_seurat <- barcodes[barcodes %in% all_cells]
  
  if(length(barcodes_seurat) < 3) {
    next
  }
  
  # 确定该克隆属于哪个样本组
  group <- ""
  if(sample_name %in% c("G1_71", "G1_84", "G1_87")) {
    group <- "data.1"
  } else if(sample_name %in% c("G2_63", "G2_73", "G2_79")) {
    group <- "data.2"
  } else if(sample_name %in% c("G3_68", "G3_80", "G3_88")) {
    group <- "data.3"
  } else if(sample_name %in% c("G4_83", "G4_90", "G4_92")) {
    group <- "data.4"
  } else {
    next
  }
  
  # 从对应图层提取
  expr_matrix_layer <- expr_matrices[[group]]
  cells_in_layer <- colnames(expr_matrix_layer)
  valid_cols <- which(cells_in_layer %in% barcodes_seurat)
  
  if(length(valid_cols) == 0) {
    next
  }
  
  # 计算平均表达量（稀疏矩阵直接计算，不转换）
  sub_matrix <- expr_matrix_layer[, valid_cols, drop = FALSE]
  avg_expr <- Matrix::rowMeans(sub_matrix)
  
  # 取前 100 个基因
  top100 <- head(names(sort(avg_expr, decreasing = TRUE)), 100)
  
  successful_clones[[i]] <- data.frame(
    Sample = sample_name,
    Clone_ID = clone_name,
    Clone_Size = clone_size,
    Cells_Analyzed = length(valid_cols),
    Rank = 1:100,
    Gene = top100,
    Avg_Expression = as.numeric(avg_expr[top100]),
    stringsAsFactors = FALSE
  )
}

# 6. 合并导出
successful_clones <- successful_clones[!sapply(successful_clones, is.null)]
all_top100 <- bind_rows(successful_clones)

write.csv(all_top100, "~/2026sc/All_A_priority_clones_top100_genes_COMPLETE_ALL.csv", row.names = FALSE)

cat("\n========== 完成 ==========\n")
cat("成功处理克隆数:", length(unique(all_top100$Clone_ID)), "/", nrow(clone_barcodes), "\n")

# 查看结果
clone_summary <- all_top100 %>%
  distinct(Clone_ID, Sample, Clone_Size, Cells_Analyzed) %>%
  arrange(desc(Clone_Size))

print(clone_summary)
