setwd("~/2026sc/important_data/rds")

library(Seurat)
library(BiocParallel)
library(BiocNeighbors)
library(DEsingle)
library(calibrate) # for textxy function

Dendritic<-readRDS("dc_annotated01.rds")
h<-Dendritic@meta.data
write.table(h,"h")
class(Dendritic)
Dendritic[["RNA"]] <- JoinLayers(Dendritic[["RNA"]])
Dendritic.sec<-as.SingleCellExperiment(Dendritic)
Dendritic$group <- gsub("_.*", "", colnames(Dendritic))

# 查看分组
table(Dendritic$group)

# 只保留 G3 和 G4
Dendritic_sub <- subset(Dendritic, group %in% c("G3", "G4"))

counts <- GetAssayData(Dendritic_sub, assay = "RNA", layer = "counts")
group <- factor(Dendritic_sub$group)

# 验证长度一致
cat("counts 列数:", ncol(counts), "\n")
cat("group 长度:", length(group), "\n")

results<-DEsingle(counts=counts,group=group)
results.classified <- DEtype(results = results, threshold = 0.05)
write.table(results.classified,"results.classified_Dendritic")

library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
res <- read.csv("volcano.txt", header=TRUE,sep="\t")
head(res)
res$p_adj <- ifelse(res$p == 0, 1e-18, res$p)
with(res, plot(log2fc, -log10(p_adj), 
               pch=20, 
               main="Volcano plot", 
               xlim=c(-4,6.5), 
               col="grey"))

with(subset(res, p<.05 & abs(log2fc)>1), points(log2fc, -log10(p), pch=20, col="blue"))
with(subset(res, p<.05 & abs(log2fc)>1), 
     points(log2fc, -log10(p_adj), pch=20, col="blue"))

with(subset(res, p<.05 & log2fc>1), points(log2fc, -log10(p), pch=20, col="red"))
with(subset(res, p<.05 & log2fc>1), 
     points(log2fc, -log10(p_adj), pch=20, col="red"))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)

target_genes <- c(# === 抗原加工/溶酶体 ===
  "Cxcl1", "Cxcl14", "Il10", "Ccl8",  "Pf4",
  "C4b", "Cd300ld3")
target <- subset(res, gene %in% target_genes)
with(target, textxy(log2fc, -log10(p), labs = gene, cex = 1.0, col = "black",font=2))
with(target, textxy(log2fc, -log10(p_adj), labs = gene, cex = 1.0, col = "black",font=2))

legend("topright",   # 位置：topright, topleft, bottomright, bottomleft
       legend = c("Up (G4 high)", "Down (G3 high)", "Not Significant"),
       col = c("red", "blue", "grey"),
       pch = 20,
       cex = 0.8,
       bty = "n")   # 去掉图例边框



target_genes <- c(
  # === 趋化因子/免疫 ===
  "Cxcl1", "Cxcl14", "Il10", "Ccl6", "Ccl8", "Ccl9", "Pf4",
  "C4b", "Il1rapl1", "Ms4a7", "Cd300ld3",
  # === DC 受体/抗原呈递 ===
  "Trem2", "Clec10a", "Clec4b1", "Mrc1", "Cd36", "Stab1", 
  "Cadm1", "Mgl2", "Folr2", "Dab2",
  # === 抗原加工/溶酶体 ===
  "Ctsb", "Ctsd", "Ctsl", "Ctso", "Lgmn", "Acp5",
  "Manba", "Tpcn2", "Steap3",
  # === DNA 损伤/细胞周期 ===
  "H2ax", "E2f1", "Bub1", "Plk1", "Mki67", "Ccne1", "Ccnb2",
  "Aurkb", "Mad2l1", "Cenpe", "Mcm2", "Pbk",
  # === 补体/免疫调节 ===
  "C1qa", "C1qb", "C1qc", "Apoe", "Gdf15",
  # === 泛素连接酶/免疫调控 ===
  "Rnf128", "Rnf121",
  # === DC 功能/代谢 ===
  "Cd63", "Gpr34", "Arg1", "Gas6", "Pdgfa", "Pdgfc",
  "Alox5", "Impdh1", "Timp2", "Pdpn", "Itgb5",
  # === 转录/表观遗传 ===
  "Dnmt3a", "Kat2a", "Kmt2d", "Bach2", "Nr1h3",
  # === 自噬/应激 ===
  "Atg14", "Sesn1", "Hvcn1", "Mcoln1"
)
