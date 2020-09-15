require(Seurat)
require(dplyr)
require(Matrix)
require(RColorBrewer)
require(pdist)
require(plyr)
require(reshape)
require(ggplot2)
require(reshape2)
library(openxlsx)
require(gplots)
require(corrplot)
library(viridis)
require(reshape2)
require(gridExtra)
library(cowplot)
library(reticulate)
library(umap)

load('TLS_120h.Robj')
Idents(TLS) <- TLS@meta.data$TLS_cluster


# Subclustering
obj <- SubsetData(TLS, ident.use =c('Seurat_6'))

s.genes <- cc.genes$s.genes
s.genes <- paste(toupper(substr(s.genes, 1, 1)), substr(tolower(s.genes), 2, nchar(s.genes)), sep= " ")
g2m.genes <- cc.genes$g2m.genes
g2m.genes <- paste(toupper(substr(g2m.genes, 1, 1)), substr(tolower(g2m.genes), 2, nchar(g2m.genes)), sep="")
obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
obj <- ScaleData(obj, vars.to.regress = c('S.Score', 'G2M.Score'), features = rownames(obj))
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:10)
obj <- FindClusters(obj)
obj <- RunUMAP(object = obj, dims = 1:30, n.neighbors=10)
obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
obj <- FindNeighbors(obj, dims = 1:10)
obj <- FindClusters(obj, resolution = 0.4)

obj@meta.data$orig.ident <- factor(obj@meta.data$orig.ident, levels=c("TLS_120h","Gastruloid","TLSCL"))
plots <- DimPlot(obj, group.by = c("experiment", "seurat_clusters"), combine = FALSE)
plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3))))
CombinePlots(plots)
DimPlot(obj,c group.by = "TLS_cluster", split.by = "orig.ident", ncol = 3, label=T)
DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", split.by = "orig.ident", ncol = 3, label=T)
neuralsub.markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


# Subclustering Somitic Cells
obj <- subset(TLS, cells=rownames(subset(TLS@meta.data, orig.ident == c('TLS_120h', 'Gastruloid','TLSCL'))))
Idents(TLS) <- TLS@meta.data$TLS_cluster

obj <- SubsetData(TLS, ident.use =c('Seurat_1','Seurat_8','Seurat_3','Seurat_0'))
s.genes <- cc.genes$s.genes
s.genes <- paste(toupper(substr(s.genes, 1, 1)), substr(tolower(s.genes), 2, nchar(s.genes)), sep= " ")
g2m.genes <- cc.genes$g2m.genes
g2m.genes <- paste(toupper(substr(g2m.genes, 1, 1)), substr(tolower(g2m.genes), 2, nchar(g2m.genes)), sep="")
obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
obj <- ScaleData(obj, vars.to.regress = c('S.Score', 'G2M.Score'), features = rownames(obj))
obj <- FindVariableFeatures(obj, verbose = FALSE)
somitesubs <- read.csv('/project/Organoid/Organoids/markergenessomitesubs.csv',header=FALSE)
somitesubs <- somitesubs[,1]
somitesubs
obj <- RunPCA(obj, features = VariableFeatures(object = obj))
#obj <- RunPCA(obj, verbose = FALSE, features = somitesubs)
obj <- FindNeighbors(obj, dims = 1:10)
obj <- FindClusters(obj)
obj <- RunUMAP(object = obj, dims = 1:30, n.neighbors=10)
obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
obj <- FindNeighbors(obj, dims = 1:10)
obj <- FindClusters(obj, resolution = 0.8)

obj@meta.data$orig.ident <- factor(obj@meta.data$orig.ident, levels=c("TLS_120h","Gastruloid","TLSCL"))
plots <- DimPlot(obj, group.by = c("experiment", "seurat_clusters"), combine = FALSE)
plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3))))
CombinePlots(plots)
DimPlot(obj,c group.by = "TLS_cluster", split.by = "orig.ident", ncol = 3, label=T)
DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", split.by = "orig.ident", ncol = 3, label=T)
