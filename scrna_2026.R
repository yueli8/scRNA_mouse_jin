library(Seurat)
library(ggplot2)
library(cowplot)
library(scater)
library(scran)
library(BiocParallel)
library(BiocNeighbors)
library(data.table)
library(dplyr)
library(Matrix)
library(clustree)

setwd("~/2026sc/rds")

P63 <- readRDS("63.rds")
P68 <- readRDS("68.rds")
P71 <- readRDS("71.rds")
P73 <- readRDS("73.rds")
P79 <- readRDS("79.rds")
P80 <- readRDS("80.rds")
P83 <- readRDS("83.rds")
P84 <- readRDS("84.rds")
P87 <- readRDS("87.rds")
P88 <- readRDS("88.rds")
P90 <- readRDS("90.rds")
P92 <- readRDS("92.rds")

# 只设置 tech / celltype，不做 RenameCells
P63@meta.data$tech <- "G2_63"; P63@meta.data$celltype <- "G2"
P68@meta.data$tech <- "G3_68"; P68@meta.data$celltype <- "G3"
P71@meta.data$tech <- "G1_71"; P71@meta.data$celltype <- "G1"
P73@meta.data$tech <- "G2_73"; P73@meta.data$celltype <- "G2"
P79@meta.data$tech <- "G2_79"; P79@meta.data$celltype <- "G2"
P80@meta.data$tech <- "G3_80"; P80@meta.data$celltype <- "G3"
P83@meta.data$tech <- "G4_83"; P83@meta.data$celltype <- "G4"
P84@meta.data$tech <- "G1_84"; P84@meta.data$celltype <- "G1"
P87@meta.data$tech <- "G1_87"; P87@meta.data$celltype <- "G1"
P88@meta.data$tech <- "G3_88"; P88@meta.data$celltype <- "G3"
P90@meta.data$tech <- "G4_90"; P90@meta.data$celltype <- "G4"
P92@meta.data$tech <- "G4_92"; P92@meta.data$celltype <- "G4"

obj_list <- list(
  P63 = P63, P68 = P68, P71 = P71, P73 = P73,
  P79 = P79, P80 = P80, P83 = P83, P84 = P84,
  P87 = P87, P88 = P88, P90 = P90, P92 = P92
)

# 合并：只在这里加一次前缀
merged <- merge(
  x = obj_list[[1]],
  y = obj_list[-1],
  add.cell.ids = names(obj_list)
)

# 修复 orig.ident：从细胞名前缀提取样本名
merged$orig.ident <- sub("_.*", "", colnames(merged))

# 检查
merged
table(merged$orig.ident)
head(colnames(merged))   # 形如 P63_AAACCAAAGACCCCAG-1

saveRDS(merged, file = "merged.rds")

hms <- readRDS("merged.rds")  
#before integrate
hms[["percent.mt"]] <- PercentageFeatureSet(hms, pattern = "^Mt-")
VlnPlot(hms, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
pancreas <- NormalizeData(object = hms, normalization.method = "LogNormalize", scale.factor = 1e4)
pancreas <- FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
pancreas <- ScaleData(pancreas, verbose = FALSE)
pancreas <- RunPCA(pancreas, npcs = 30, verbose = FALSE)
pancreas <- RunUMAP(pancreas, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) + 
  NoLegend()
plot_grid(p1,p2)
saveRDS(pancreas, file = "pancreas.rds")
pancreas <- readRDS("pancreas.rds")  
#integrate
pancreas.list <- SplitObject(pancreas, split.by = "tech")

#saveRDS(pancreas.list, file = "pancreas_list.rds") 
#pancreas.list <- readRDS("pancreas_list.rds") 
for (i in 1: length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000, 
                                             verbose = FALSE)
}

reference.list <- pancreas.list[c("G2_63","G3_68","G1_71","G2_73",
                                  "G2_79","G3_80", "G4_83", "G1_84","G1_87","G3_88", "G4_90", "G4_92" )]

# 选参考样本（用 reference.list 里的真实名字）
ref_names <- c("G2_79", "G1_71", "G2_73")
ref_idx <- match(ref_names, names(reference.list))
ref_idx   # 输出 5 3 4

pancreas.anchors <- FindIntegrationAnchors(
  object.list = reference.list,
  dims = 1:30,
  reduction = "rpca",
  reference = ref_idx
)

# 接下去
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)










pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30,reduction = "rpca")

# 1. 保存 anchors（最关键，重跑最耗时）
saveRDS(pancreas.anchors, file = "pancreas_anchors.rpca.rds")

# 2. 保存 reference.list（整合时还需要它）
saveRDS(reference.list, file = "reference_list.rds")













pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
DefaultAssay(pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype")
plot_grid(p1,p2)
saveRDS(pancreas.integrated, file = "hms_after_integrated.rds")
saveRDS(pancreas.integrated, file = "P57586061_after_integrated.rds")

hms_individual_integrated<-readRDS(file="hms_after_integrated.rds")
p1 <- DimPlot(hms_individual_integrated, reduction = "umap", group.by = "celltype")
p1


hms_neighbor<- FindNeighbors(hms_individual_integrated, dims = 1:30)
obj <- FindClusters(hms_neighbor, resolution = seq(0.5,1.2,by=0.1))
#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
clustree(obj)
hms_cluster <- FindClusters( hms_neighbor, resolution = 1.2)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:30)
DimPlot(hms_cluster, reduction = "umap")
saveRDS(hms_cluster, file = "hms_cluster_test_1.2.rds")

scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.05, 
                                logfc.threshold = 0.05
)
write.table(scRNA.markers,file="cellMarkers.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) 
write.csv(file="top20_cell_markers.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(file="top20_marker_genes_1.2.csv",top20_table,row.names=F)

#细胞及细胞中基因与RNA数量
slotNames(hms_cluster)
#assay
hms_cluster@assays
dim(hms_cluster@meta.data)
View(hms_cluster@meta.data)

hms_cluster<-readRDS(file="hms_cluster_test_1.2.rds")
DimPlot(hms_cluster, reduction = "umap")
