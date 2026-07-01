library(stringr)
library(dplyr)
setwd("~/2026sc/scRepertoire")

# Group 1

df_clone <- read.csv("WetLab_G1_Ready.csv", stringsAsFactors = F)

df_contig <- read.csv("71_filtered_contig_annotations.csv", stringsAsFactors = F)

# baseR分割下划线取末尾
split_last <- function(s){
  vec <- strsplit(s, "_")[[1]]
  vec[length(vec)]
}
df_clone$match_barcode <- sapply(df_clone$cell_barcodes, split_last)

contig_sub <- df_contig[, c("barcode", "d_gene", "c_gene")]
df_merge <- merge(df_clone, contig_sub, by.x="match_barcode", by.y="barcode", all.x=T)
df_merge$match_barcode <- NULL

write.csv(df_merge, "71合并后TCR克隆表.csv", row.names=F)


df_contig <- read.csv("84_filtered_contig_annotations.csv", stringsAsFactors = F)

# baseR分割下划线取末尾
split_last <- function(s){
  vec <- strsplit(s, "_")[[1]]
  vec[length(vec)]
}
df_clone$match_barcode <- sapply(df_clone$cell_barcodes, split_last)

contig_sub <- df_contig[, c("barcode", "d_gene", "c_gene")]
df_merge <- merge(df_clone, contig_sub, by.x="match_barcode", by.y="barcode", all.x=T)
df_merge$match_barcode <- NULL

write.csv(df_merge, "84合并后TCR克隆表.csv", row.names=F)



df_contig <- read.csv("87_filtered_contig_annotations.csv", stringsAsFactors = F)

# baseR分割下划线取末尾
split_last <- function(s){
  vec <- strsplit(s, "_")[[1]]
  vec[length(vec)]
}
df_clone$match_barcode <- sapply(df_clone$cell_barcodes, split_last)

contig_sub <- df_contig[, c("barcode", "d_gene", "c_gene")]
df_merge <- merge(df_clone, contig_sub, by.x="match_barcode", by.y="barcode", all.x=T)
df_merge$match_barcode <- NULL

write.csv(df_merge, "87合并后TCR克隆表.csv", row.names=F)

#修改文件名称后继续做
df1 <- read.csv("71合并后TCR克隆表.csv", stringsAsFactors = F)
df2 <- read.csv("84合并后TCR克隆表.csv", stringsAsFactors = F)
df3 <- read.csv("87合并后TCR克隆表.csv", stringsAsFactors = F)

# 2. 分步合并，优先填充非NA的d_gene、c_gene
# 第一步：df1左连df2，用coalesce替换NA
merge12 <- left_join(df1,
                     df2[, c("Sample", "CDR3_comb", "d_gene", "c_gene")],
                     by = c("Sample", "CDR3_comb"),
                     suffix = c("", "_y2")) %>%
  mutate(
    d_gene = coalesce(d_gene, d_gene_y2),
    c_gene = coalesce(c_gene, c_gene_y2)
  ) %>%
  select(-ends_with("_y2")) # 删除临时后缀列

# 第二步：合并结果再连df3，继续补全NA
final_df <- left_join(merge12,
                      df3[, c("Sample", "CDR3_comb", "d_gene", "c_gene")],
                      by = c("Sample", "CDR3_comb"),
                      suffix = c("", "_y3")) %>%
  mutate(
    d_gene = coalesce(d_gene, d_gene_y3),
    c_gene = coalesce(c_gene, c_gene_y3)
  ) %>%
  select(-ends_with("_y3")) %>%
  distinct() # 删除完全重复行

# 3. 导出最终合并文件
write.csv(final_df, "G1TCR_合并去NA最终表.csv", row.names = FALSE)

# 校验关键列
head(final_df[, c("Sample","CDR3_comb","d_gene","c_gene")])


# Group 2

df_clone <- read.csv("WetLab_G2_Ready.csv", stringsAsFactors = F)

df_contig <- read.csv("63_filtered_contig_annotations.csv", stringsAsFactors = F)

# baseR分割下划线取末尾
split_last <- function(s){
  vec <- strsplit(s, "_")[[1]]
  vec[length(vec)]
}
df_clone$match_barcode <- sapply(df_clone$cell_barcodes, split_last)

contig_sub <- df_contig[, c("barcode", "d_gene", "c_gene")]
df_merge <- merge(df_clone, contig_sub, by.x="match_barcode", by.y="barcode", all.x=T)
df_merge$match_barcode <- NULL

write.csv(df_merge, "63合并后TCR克隆表.csv", row.names=F)


df_contig <- read.csv("73_filtered_contig_annotations.csv", stringsAsFactors = F)

# baseR分割下划线取末尾
split_last <- function(s){
  vec <- strsplit(s, "_")[[1]]
  vec[length(vec)]
}
df_clone$match_barcode <- sapply(df_clone$cell_barcodes, split_last)

contig_sub <- df_contig[, c("barcode", "d_gene", "c_gene")]
df_merge <- merge(df_clone, contig_sub, by.x="match_barcode", by.y="barcode", all.x=T)
df_merge$match_barcode <- NULL

write.csv(df_merge, "73合并后TCR克隆表.csv", row.names=F)



df_contig <- read.csv("79_filtered_contig_annotations.csv", stringsAsFactors = F)

# baseR分割下划线取末尾
split_last <- function(s){
  vec <- strsplit(s, "_")[[1]]
  vec[length(vec)]
}
df_clone$match_barcode <- sapply(df_clone$cell_barcodes, split_last)

contig_sub <- df_contig[, c("barcode", "d_gene", "c_gene")]
df_merge <- merge(df_clone, contig_sub, by.x="match_barcode", by.y="barcode", all.x=T)
df_merge$match_barcode <- NULL

write.csv(df_merge, "79合并后TCR克隆表.csv", row.names=F)

#修改文件名称后继续做
df1 <- read.csv("63合并后TCR克隆表.csv", stringsAsFactors = F)
df2 <- read.csv("73合并后TCR克隆表.csv", stringsAsFactors = F)
df3 <- read.csv("79合并后TCR克隆表.csv", stringsAsFactors = F)

# 2. 分步合并，优先填充非NA的d_gene、c_gene
# 第一步：df1左连df2，用coalesce替换NA
merge12 <- left_join(df1,
                     df2[, c("Sample", "CDR3_comb", "d_gene", "c_gene")],
                     by = c("Sample", "CDR3_comb"),
                     suffix = c("", "_y2")) %>%
  mutate(
    d_gene = coalesce(d_gene, d_gene_y2),
    c_gene = coalesce(c_gene, c_gene_y2)
  ) %>%
  select(-ends_with("_y2")) # 删除临时后缀列

# 第二步：合并结果再连df3，继续补全NA
final_df <- left_join(merge12,
                      df3[, c("Sample", "CDR3_comb", "d_gene", "c_gene")],
                      by = c("Sample", "CDR3_comb"),
                      suffix = c("", "_y3")) %>%
  mutate(
    d_gene = coalesce(d_gene, d_gene_y3),
    c_gene = coalesce(c_gene, c_gene_y3)
  ) %>%
  select(-ends_with("_y3")) %>%
  distinct() # 删除完全重复行

# 3. 导出最终合并文件
write.csv(final_df, "G2TCR_合并去NA最终表.csv", row.names = FALSE)

# 校验关键列
head(final_df[, c("Sample","CDR3_comb","d_gene","c_gene")])


# Group 3
df_clone <- read.csv("WetLab_G3_Ready.csv", stringsAsFactors = F)

df_contig <- read.csv("68_filtered_contig_annotations.csv", stringsAsFactors = F)

# baseR分割下划线取末尾
split_last <- function(s){
  vec <- strsplit(s, "_")[[1]]
  vec[length(vec)]
}
df_clone$match_barcode <- sapply(df_clone$cell_barcodes, split_last)

contig_sub <- df_contig[, c("barcode", "d_gene", "c_gene")]
df_merge <- merge(df_clone, contig_sub, by.x="match_barcode", by.y="barcode", all.x=T)
df_merge$match_barcode <- NULL

write.csv(df_merge, "68合并后TCR克隆表.csv", row.names=F)


df_contig <- read.csv("80_filtered_contig_annotations.csv", stringsAsFactors = F)

# baseR分割下划线取末尾
split_last <- function(s){
  vec <- strsplit(s, "_")[[1]]
  vec[length(vec)]
}
df_clone$match_barcode <- sapply(df_clone$cell_barcodes, split_last)

contig_sub <- df_contig[, c("barcode", "d_gene", "c_gene")]
df_merge <- merge(df_clone, contig_sub, by.x="match_barcode", by.y="barcode", all.x=T)
df_merge$match_barcode <- NULL

write.csv(df_merge, "80合并后TCR克隆表.csv", row.names=F)



df_contig <- read.csv("88_filtered_contig_annotations.csv", stringsAsFactors = F)

# baseR分割下划线取末尾
split_last <- function(s){
  vec <- strsplit(s, "_")[[1]]
  vec[length(vec)]
}
df_clone$match_barcode <- sapply(df_clone$cell_barcodes, split_last)

contig_sub <- df_contig[, c("barcode", "d_gene", "c_gene")]
df_merge <- merge(df_clone, contig_sub, by.x="match_barcode", by.y="barcode", all.x=T)
df_merge$match_barcode <- NULL

write.csv(df_merge, "88合并后TCR克隆表.csv", row.names=F)

#修改文件名称后继续做
df1 <- read.csv("68合并后TCR克隆表.csv", stringsAsFactors = F)
df2 <- read.csv("80合并后TCR克隆表.csv", stringsAsFactors = F)
df3 <- read.csv("88合并后TCR克隆表.csv", stringsAsFactors = F)

# 2. 分步合并，优先填充非NA的d_gene、c_gene
# 第一步：df1左连df2，用coalesce替换NA
merge12 <- left_join(df1,
                     df2[, c("Sample", "CDR3_comb", "d_gene", "c_gene")],
                     by = c("Sample", "CDR3_comb"),
                     suffix = c("", "_y2")) %>%
  mutate(
    d_gene = coalesce(d_gene, d_gene_y2),
    c_gene = coalesce(c_gene, c_gene_y2)
  ) %>%
  select(-ends_with("_y2")) # 删除临时后缀列

# 第二步：合并结果再连df3，继续补全NA
final_df <- left_join(merge12,
                      df3[, c("Sample", "CDR3_comb", "d_gene", "c_gene")],
                      by = c("Sample", "CDR3_comb"),
                      suffix = c("", "_y3")) %>%
  mutate(
    d_gene = coalesce(d_gene, d_gene_y3),
    c_gene = coalesce(c_gene, c_gene_y3)
  ) %>%
  select(-ends_with("_y3")) %>%
  distinct() # 删除完全重复行

# 3. 导出最终合并文件
write.csv(final_df, "G3TCR_合并去NA最终表.csv", row.names = FALSE)

# 校验关键列
head(final_df[, c("Sample","CDR3_comb","d_gene","c_gene")])


# Group 4

df_clone <- read.csv("WetLab_G4_Ready.csv", stringsAsFactors = F)

df_contig <- read.csv("83_filtered_contig_annotations.csv", stringsAsFactors = F)

# baseR分割下划线取末尾
split_last <- function(s){
  vec <- strsplit(s, "_")[[1]]
  vec[length(vec)]
}
df_clone$match_barcode <- sapply(df_clone$cell_barcodes, split_last)

contig_sub <- df_contig[, c("barcode", "d_gene", "c_gene")]
df_merge <- merge(df_clone, contig_sub, by.x="match_barcode", by.y="barcode", all.x=T)
df_merge$match_barcode <- NULL

write.csv(df_merge, "83合并后TCR克隆表.csv", row.names=F)


df_contig <- read.csv("90_filtered_contig_annotations.csv", stringsAsFactors = F)

# baseR分割下划线取末尾
split_last <- function(s){
  vec <- strsplit(s, "_")[[1]]
  vec[length(vec)]
}
df_clone$match_barcode <- sapply(df_clone$cell_barcodes, split_last)

contig_sub <- df_contig[, c("barcode", "d_gene", "c_gene")]
df_merge <- merge(df_clone, contig_sub, by.x="match_barcode", by.y="barcode", all.x=T)
df_merge$match_barcode <- NULL

write.csv(df_merge, "90合并后TCR克隆表.csv", row.names=F)



df_contig <- read.csv("92_filtered_contig_annotations.csv", stringsAsFactors = F)

# baseR分割下划线取末尾
split_last <- function(s){
  vec <- strsplit(s, "_")[[1]]
  vec[length(vec)]
}
df_clone$match_barcode <- sapply(df_clone$cell_barcodes, split_last)

contig_sub <- df_contig[, c("barcode", "d_gene", "c_gene")]
df_merge <- merge(df_clone, contig_sub, by.x="match_barcode", by.y="barcode", all.x=T)
df_merge$match_barcode <- NULL

write.csv(df_merge, "92合并后TCR克隆表.csv", row.names=F)

#修改文件名称后继续做
df1 <- read.csv("83合并后TCR克隆表.csv", stringsAsFactors = F)
df2 <- read.csv("90合并后TCR克隆表.csv", stringsAsFactors = F)
df3 <- read.csv("92合并后TCR克隆表.csv", stringsAsFactors = F)

# 2. 分步合并，优先填充非NA的d_gene、c_gene
# 第一步：df1左连df2，用coalesce替换NA
merge12 <- left_join(df1,
                     df2[, c("Sample", "CDR3_comb", "d_gene", "c_gene")],
                     by = c("Sample", "CDR3_comb"),
                     suffix = c("", "_y2")) %>%
  mutate(
    d_gene = coalesce(d_gene, d_gene_y2),
    c_gene = coalesce(c_gene, c_gene_y2)
  ) %>%
  select(-ends_with("_y2")) # 删除临时后缀列

# 第二步：合并结果再连df3，继续补全NA
final_df <- left_join(merge12,
                      df3[, c("Sample", "CDR3_comb", "d_gene", "c_gene")],
                      by = c("Sample", "CDR3_comb"),
                      suffix = c("", "_y3")) %>%
  mutate(
    d_gene = coalesce(d_gene, d_gene_y3),
    c_gene = coalesce(c_gene, c_gene_y3)
  ) %>%
  select(-ends_with("_y3")) %>%
  distinct() # 删除完全重复行

# 3. 导出最终合并文件
write.csv(final_df, "G4TCR_合并去NA最终表.csv", row.names = FALSE)

# 校验关键列
head(final_df[, c("Sample","CDR3_comb","d_gene","c_gene")])





