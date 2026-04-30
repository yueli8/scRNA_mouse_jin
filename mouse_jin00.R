library(Seurat)
library(ggplot2)
library(cowplot)
library(BiocParallel)
library(BiocNeighbors)
library(data.table)
library(dplyr)
library(Matrix)
library(readr)

setwd("~/2026sc")
total<-readRDS(file="mts_cluster_test_5.rds")
lymphoid<-readRDS(file="lymphoid.rds")
hms_neighbor<- FindNeighbors(lymphoid, dims = 1:30)
obj <- FindClusters(hms_neighbor, resolution = seq(0.5,1.2,by=0.1))
hms_cluster <- FindClusters( hms_neighbor, resolution = 10)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:30)
DimPlot(hms_cluster, reduction = "umap",label = TRUE)
saveRDS(hms_cluster, file = "T_test_10.rds")

scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.00000001, 
                                logfc.threshold = 0.00000001)
write.table(scRNA.markers,file="TcellMarkers.txt",sep="\t",row.names=F,quote=F)

scRNA.markers <- read.table("TcellMarkers.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 查看数据
head(data)

#挑选每个细胞亚群中特意高表达的20个基因
top100 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) 
write.csv(file="top100_cell_markers.csv",top100)
top100_table=unstack(top100, gene ~ cluster)
names(top100_table)=gsub("X","cluster",names(top100_table))
write.csv(file="Top100_marker_genes_10.csv",top100_table,row.names=F)

hms_cluster<-readRDS(file="T_test_10.rds")
DimPlot(hms_cluster, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster, reduction = "umap", label = FALSE, pt.size = 0.5) 
new.cluster.ids <- c("t","t","t","t","stromal","t","t","t","t","t","t",
                     "t","t","b","t","myeloid","t","t","t","myeloid","t",
                     "t","t","t","myeloid","t","t","t","stromal","t","t",
                     "myeloid","t","nk","t","t","t","t","nk","ilc","t",
                     "t","t","nk","t","t","nk","t","myeloid","t","t",
                     "t","t","t","myeloid","t","t","t","b","t","nk",
                     "t","t","t","b","t","b","t","t","t","t",
                     "b","t","t","myeloid","myeloid","myeloid","t","stromal","t","b",
                     "t","t","t","t","myeloid","myeloid","myeloid","t","b","myeloid",
                     "t","t","myeloid","t","t","b","t","stromal","t","t",
                     "t","stromal","t","t","stromal","t","t","myeloid","myeloid","t",
                     "t","myeloid","t","t","b","myeloid") 

names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(hms_cluster_id, file = "hms_cluster_id_only_lymphoid.rds")
Only_lymphoid<-subset(hms_cluster_id, idents=c("t", "b","nk","ilc"))
DimPlot(Only_lymphoid, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(Only_lymphoid, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(Only_lymphoid, file = "Only_lymphoid_cluster_id_test.rds")

table(Idents(Only_lymphoid), Only_lymphoid$orig.ident)
t<-subset(hms_cluster_id, idents=c('t'))
DimPlot(t, reduction = "umap")
saveRDS(t, file="t.rds")

