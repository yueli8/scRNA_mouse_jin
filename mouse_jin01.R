library(Seurat)
library(ggplot2)
library(cowplot)
library(BiocParallel)
library(BiocNeighbors)
library(data.table)
library(dplyr)
library(Matrix)

setwd("~/2026sc")
lymphoid<-readRDS(file="Lymphoid.rds")
hms_neighbor<- FindNeighbors(lymphoid, dims = 1:30)
obj <- FindClusters(hms_neighbor, resolution = seq(0.5,1.2,by=0.1))
hms_cluster <- FindClusters( hms_neighbor, resolution = 10)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:30)
DimPlot(hms_cluster, reduction = "umap",label = TRUE)
saveRDS(hms_cluster, file = "T_test_10.rds")
hms_cluster<-readRDS(file="T_test_10.rds")
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
write.csv(file="Ttop100_marker_genes_10.csv",top100_table,row.names=F)




