library(Seurat)
library(scRepertoire)
library(devtools)
#library(circlize)
library(scales)
library(dplyr)
library(tidyr)
#devtools::install_github("ncborcherding/scRepertoire@dev")
#manual: https://ncborcherding.github.io/vignettes/vignette.html
#can install from bioconductor （download *tar.gz file,and install from toold）,it also from github. 
#2.0 changes cloneType -> cloneSize compareClonaltypes -> clonalCompare highlightClonotypes -> highlightClones
# alluvialClonotypes -> alluvialClones
# split.by was combined into group.by
# expression2List was deprecated, but can be used with flag "force = TRUE"
# clonalCompare arg "numbers" is now "top.clones"
# StartracDiversity "sample" is not an arg any longer, replaced by "group.by"; "by" is deprecated but default is all 3 (expa, migr, tran)
# clonalNetwork "identity" is not an arg any longer, replaced by "group.by"

setwd("~/2026sc/scRepertoire")

G2_63 <- read.csv("~/2026sc/scRepertoire/63_filtered_contig_annotations.csv")
G3_68 <- read.csv("~/2026sc/scRepertoire/68_filtered_contig_annotations.csv")
G1_71 <- read.csv("~/2026sc/scRepertoire/71_filtered_contig_annotations.csv")
G2_73 <- read.csv("~/2026sc/scRepertoire/73_filtered_contig_annotations.csv")
G2_79 <- read.csv("~/2026sc/scRepertoire/79_filtered_contig_annotations.csv")
G3_80 <- read.csv("~/2026sc/scRepertoire/80_filtered_contig_annotations.csv")
G4_83 <- read.csv("~/2026sc/scRepertoire/83_filtered_contig_annotations.csv")
G1_84 <- read.csv("~/2026sc/scRepertoire/84_filtered_contig_annotations.csv")
G1_87 <- read.csv("~/2026sc/scRepertoire/87_filtered_contig_annotations.csv")
G3_88 <- read.csv("~/2026sc/scRepertoire/88_filtered_contig_annotations.csv")
G4_90 <- read.csv("~/2026sc/scRepertoire/90_filtered_contig_annotations.csv")
G4_92 <- read.csv("~/2026sc/scRepertoire/92_filtered_contig_annotations.csv")
contig_list <- list( G1_71, G1_84, G1_87, G2_63, G2_73, G2_79, G3_68, G3_80, G3_88, G4_83, G4_90, G4_92)

combined <- combineTCR(contig_list, 
                       samples = c("G1_71","G1_84", "G1_87", "G2_63",  "G2_73", "G2_79", "G3_68","G3_80",  "G3_88", "G4_83","G4_90", "G4_92"), 
)

clonalQuant(combined, cloneCall="strict", scale = T)


##1.合并12个原始contig，批量加Sample标签（和你samples顺序一一对应）
samp_order = c("G1_71","G1_84", "G1_87", "G2_63",  "G2_73", "G2_79", "G3_68","G3_80",  "G3_88", "G4_83","G4_90", "G4_92")

#逐个给数据框添加Sample列
raw_list = list(G2_63, G3_68, G1_71, G2_73, G2_79, G3_80,
                G4_83, G1_84, G1_87, G3_88, G4_90, G4_92)
for(i in seq_along(raw_list)){
  raw_list[[i]]$Sample = samp_order[i]
}
contig_raw = bind_rows(raw_list)

##2.构建克隆注释库：TRAV/TRBV+CDR3αβ+克隆Count（用来匹配vizGenes/clonalcomp）
clone_info = contig_raw %>%
  filter(productive == "true", chain %in% c("TRA","TRB")) %>%
  select(Sample, barcode, chain, v_gene, j_gene, cdr3) %>%
  pivot_wider(names_from = chain, values_from = c(v_gene,j_gene,cdr3)) %>%
  rename(TRAV = v_gene_TRA, TRBV = v_gene_TRB,
         CDR3a = cdr3_TRA, CDR3b = cdr3_TRB,
         Ja = j_gene_TRA, Jb = j_gene_TRB) %>%
  mutate(CDR3_comb = paste(CDR3a, CDR3b, sep = "_")) %>%
  filter(!is.na(TRAV) & !is.na(TRBV))

#统计每个CDR3克隆细胞数
clone_count = clone_info %>%
  group_by(Sample, CDR3_comb) %>%
  summarise(Count = n(), .groups = "drop")

clone_info_final = left_join(clone_info, clone_count, by = c("Sample", "CDR3_comb"))

clone_info_fix <- clone_info_final %>%
  mutate(across(c(TRAV,TRBV,Ja,Jb,CDR3a,CDR3b,CDR3_comb),
                ~sapply(.x, function(z) paste0(z,collapse = ";"))))

# 导出
write.csv(clone_info_fix,"~/2026sc/clone_info_final_fix.csv",row.names = F)

# 匹配键
clone_key <- clone_info_fix %>% select(Sample, CDR3_comb, TRAV, TRBV, CDR3a, CDR3b, Count, Ja, Jb)
write.csv(clone_key,"~/2026sc/clone_matching_key.csv",row.names=F)

viztable<-vizGenes(combined, 
                   x.axis = "TRAV", 
                   plot = "barplot", 
                   summary.fun = "proportion",
                   order = "variance",
                   exportTable = TRUE)
write.csv(viztable, file = "~/2026sc/viztable_prop.csv")

viztable<-vizGenes(combined, 
                   x.axis = "TRAV", 
                   plot = "barplot", 
                   summary.fun = "count",
                   order = "variance",
                   exportTable = TRUE)
write.csv(viztable, file = "~/2026sc/viztable_ct.csv")

clonalcomp<-clonalCompare(combined, top.clones = 100, samples = c("G1_71","G1_84", "G1_87", "G2_63",  "G2_73", "G2_79", "G3_68","G3_80",  "G3_88", "G4_83","G4_90", "G4_92"), cloneCall="aa", graph = "alluvial", exportTable = TRUE)
write.csv(clonalcomp, file = "~/2026sc/clonalcomp.csv")

clonalcomp<-clonalCompare(combined, top.clones = 100, samples = c("G1_71","G1_84", "G1_87", "G2_63",  "G2_73", "G2_79", "G3_68","G3_80",  "G3_88", "G4_83","G4_90", "G4_92"),  cloneCall="aa", graph = "alluvial", proportion = FALSE, exportTable = TRUE)
write.csv(clonalcomp, file = "~/2026sc/clonalcomp_ct.csv")

##=========1.TRAV proportion 表+注释=========
viztable_prop <- vizGenes(combined,
                          x.axis = "TRAV",
                          plot = "barplot",
                          summary.fun = "proportion",
                          order = "variance",
                          exportTable = TRUE)
# matrix转长格式df
viztable_prop <- as.data.frame(viztable_prop) %>%
  mutate(TRAV = rownames(.)) %>%
  pivot_longer(cols = -TRAV, names_to = "Sample", values_to = "proportion")
# 匹配V、CDR3、Count信息
viztable_prop_anno <- left_join(viztable_prop, clone_info_fix, by = c("TRAV","Sample"))
# 导出原始+注释两张
write.csv(viztable_prop, file = "~/2026sc/viztable_prop.csv", row.names=F)
write.csv(viztable_prop_anno, file = "~/2026sc/viztable_prop_anno.csv", row.names=F)

##=========2.TRAV count 表+注释=========
viztable_ct <- vizGenes(combined,
                        x.axis = "TRAV",
                        plot = "barplot",
                        summary.fun = "count",
                        order = "variance",
                        exportTable = TRUE)
viztable_ct <- as.data.frame(viztable_ct) %>%
  mutate(TRAV = rownames(.)) %>%
  pivot_longer(cols = -TRAV, names_to = "Sample", values_to = "count")
viztable_ct_anno <- left_join(viztable_ct, clone_info_fix, by = c("TRAV","Sample"))
write.csv(viztable_ct, file = "~/2026sc/viztable_ct.csv", row.names=F)
write.csv(viztable_ct_anno, file = "~/2026sc/viztable_ct_anno.csv", row.names=F)

samp_vec = c("G1_71","G1_84", "G1_87", "G2_63",  "G2_73", "G2_79", "G3_68","G3_80",  "G3_88", "G4_83","G4_90", "G4_92")

##=====proportion=====
clonalcomp <- clonalCompare(combined, top.clones = 100, samples = samp_vec, cloneCall="aa", graph = "alluvial", exportTable = TRUE) %>% as.data.frame()
clonalcomp <- rename(clonalcomp, CDR3_comb = clones)
clonalcomp_anno <- left_join(clonalcomp, clone_info_fix, by = c("CDR3_comb","Sample"))

# 替换write.csv，用write.table
write.table(clonalcomp, "~/2026sc/clonalcomp.csv", sep = ",", row.names = F, col.names = T)
write.table(clonalcomp_anno, "~/2026sc/clonalcomp_anno.csv", sep = ",", row.names = F, col.names = T)

##=====count=====
clonalcomp_ct <- clonalCompare(combined, top.clones = 100, samples = samp_vec, cloneCall="aa", graph = "alluvial", proportion = FALSE, exportTable = TRUE) %>% as.data.frame()
clonalcomp_ct <- rename(clonalcomp_ct, CDR3_comb = clones)
clonalcomp_ct_anno <- left_join(clonalcomp_ct, clone_info_fix, by = c("CDR3_comb","Sample"))

write.table(clonalcomp_ct, "~/2026sc/clonalcomp_ct.csv", sep = ",", row.names = F, col.names = T)
write.table(clonalcomp_ct_anno, "~/2026sc/clonalcomp_ct_anno.csv", sep = ",", row.names = F, col.names = T)

#TRBV proportion
viz_trbv_prop <- vizGenes(combined,x.axis="TRBV",summary.fun="proportion",exportTable=T)
viz_trbv_prop <- as.data.frame(viz_trbv_prop) %>%
  mutate(TRBV=rownames(.)) %>% pivot_longer(-TRBV,names_to="Sample",values_to="proportion")
viz_trbv_anno <- left_join(viz_trbv_prop,clone_info_fix,by=c("TRBV","Sample"))

write.table(viz_trbv_prop,"~/2026sc/vizTRBV_prop.csv",sep=",",row.names=F,col.names=T)
write.table(viz_trbv_anno,"~/2026sc/vizTRBV_prop_anno.csv",sep=",",row.names=F,col.names=T)

#TRBV count
viz_trbv_ct <- vizGenes(combined,x.axis="TRBV",summary.fun="count",exportTable=T)
viz_trbv_ct <- as.data.frame(viz_trbv_ct) %>%
  mutate(TRBV=rownames(.)) %>% pivot_longer(-TRBV,names_to="Sample",values_to="count")
viz_trbv_ct_anno <- left_join(viz_trbv_ct,clone_info_fix,by=c("TRBV","Sample"))

write.table(viz_trbv_ct,"~/2026sc/vizTRBV_ct.csv",sep=",",row.names=F,col.names=T)
write.table(viz_trbv_ct_anno,"~/2026sc/vizTRBV_ct_anno.csv",sep=",",row.names=F,col.names=T)


#提取seurat元数据
meta_raw <- seurat@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  select(orig.ident, barcode, celltype) %>%
  rename(Sample = orig.ident)

#1 把meta_raw的barcode去掉前缀，只保留-1结尾部分
meta_raw$barcode = sub("^.*_", "", meta_raw$barcode)

#2 再合并#Sample+barcode双主键合并TCR信息
final_all <- left_join(meta_raw, clone_info_fix, by = c("Sample","barcode"))
# 查看匹配数量
table(!is.na(final_all$TRAV))

final_all <- final_all %>%
  mutate(Clone_Grade = case_when(
    Count >= 22 ~ "A_优先湿实验",
    Count >= 1 ~ "B_备选",
    TRUE         ~ "C_留存暂不做"
  ))

#全量汇总表
write.table(final_all,"~/2026sc/TCR_All_Final.csv",sep=",",row.names=F,col.names=T)
#仅A级候选（交付实验室）
A_TCR <- subset(final_all,Clone_Grade=="A_优先湿实验")
write.table(A_TCR,"~/2026sc/A级_TCR候选序列.csv",sep=",",row.names=F,col.names=T)
# 备选B级单独存表
B_TCR <- subset(final_all,Clone_Grade=="B_备选")
write.table(B_TCR,"~/2026sc/B级_TCR备选.csv",sep=",",row.names=F,col.names=T)

A_unique <- A_TCR %>%
  distinct(TRAV,TRBV,CDR3a,CDR3b,.keep_all = TRUE)

write.table(A_unique,"~/2026sc/A级去重_TCR构建清单.csv",sep=",",row.names=F)


#1读取seurat，提取meta信息
obj <- readRDS("only_lymphoid_from_total.rds")
unique(obj$orig.ident)
cells_in_G4 <- colnames(obj)[obj$orig.ident == "G4_90"]
# 模糊匹配条码
grep("AAACGGCCACTGATTA", cells_in_G4, value = T)
#[1] "G4_90_AAACGGCCACTGATTA-1"

cell_real = "G4_92_ACTGAAGCAGTAAGAT-1"
cell_find = obj[, cell_real]


# log标准化数据 layer=data
expr <- GetAssayData(cell_find, assay = "RNA", layer = "data")

# 排序取TOP100
top100 <- sort(expr[,1], decreasing = TRUE)[1:100]
df_top100 <- data.frame(Gene=names(top100), LogExpr=top100)

# 输出前20预览
head(df_top100,100)

