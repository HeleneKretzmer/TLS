require(Seurat)
require(ComplexHeatmap)
require(cluster)
require(circlize)

# load data
load("TLS.Robj")
load("Mouse_ref.Robj")

# Subcluster in vivo mouse PSM (30) and Somites states
mouse <- subset(Mouse, cells=rownames(subset(Mouse@meta.data, cluster %in% c("30","18"))))
mouse <- RunPCA(mouse)
mouse <- FindNeighbors(mouse, dims = 1:10)
mouse <- FindClusters(mouse, resolution=0.14)
mouse <- RunUMAP(object = mouse, dims = 1:10)

mouse.markers <- FindAllMarkers(mouse, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- mouse.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(mouse, features = top10$gene) + NoLegend()

# split Somites into aPSM (18.1) and Somites (18)
aPSM_cells <- rownames(subset(mouse@meta.data, seurat_clusters==3))
mouse@meta.data$cluster[rownames(mouse@meta.data) %in% aPSM_cells] <- "18.1"


# remove in vivo mouse clusters that don't matter
mouse <- subset(mouse, cells=rownames(subset(mouse@meta.data, !cluster %in% c('0','25','28','14','29','40'))))
mouse <- subset(mouse, cells=rownames(subset(mouse@meta.data, !cluster %in% c('8','19','1','11','10','24','35'))))
mouse <- subset(mouse, cells=rownames(subset(mouse@meta.data, !cluster %in% c('20','34'))))
mouse <- subset(mouse, cells=rownames(subset(mouse@meta.data, stage %in% c('WT_75','WT_80','WT_85'))))


# Normalize data and find variable features
Mouse <- SCTransform(Mouse, verbose = FALSE)
TLS <- SCTransform(TLS, verbose = FALSE)

# Cell type classification using an integrated reference
anchors <- FindTransferAnchors(reference = Mouse, query = TLS, dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = Mouse$cluster, dims = 1:30)
TLS <- AddMetaData(TLS, metadata = predictions)


# prediction and assignment scores
pred.scores <- ddply(TLS@meta.data[,c(5,17:44)], .(seurat_clusters), colwise(mean))
rownames(pred.scores) <- pred.scores$seurat_clusters
colnames(pred.scores) <- gsub("prediction.score.","WT_",colnames(pred.scores))

assigned.states <- data.frame(dcast(data.frame(table(TLS@meta.data[,c(5,16)])), seurat_clusters~predicted.id))
rownames(assigned.states) <- assigned.states$seurat_clusters
colnames(assigned.states) <- gsub("X","WT_",colnames(assigned.states))
assigned.states[,-1] <- 100*assigned.states[,-1]/rowSums(assigned.states[,-1])

Heatmap(pred.scores[,-1], cluster_rows = diana(pred.scores[,-1]), cluster_columns = agnes(t(pred.scores[,-1])), col=colorRamp2(c(0,0.7), c("white", "orange")))
Heatmap(assigned.states[,-1], cluster_rows = diana(assigned.states[,-1]), cluster_columns = agnes(t(assigned.states[,-1])), col=colorRamp2(c(0,100), c("white", "grey20")))


# merged object for TLS/Mouse marker identification
anchors <- FindIntegrationAnchors(object.list = c(TLS, Mouse), dims = 1:20)
Mouse_TLS_combined <- IntegrateData(anchorset = anchors, dims = 1:20)

CellType <- data.frame(CellType=Mouse_TLS_combined@meta.data$seurat_clusters)
CellType$CellType[is.na(CellType$CellType)] <- Mouse_TLS_combined@meta.data$predicted.id[is.na(CellType$CellType)]
CellType$study <- gsub('embryos','WT',gsub('_.*','',Mouse_TLS_combined@meta.data$orig.ident))
CellType$ID <- paste(CellType$CellType, CellType$study)
rownames(CellType) <- rownames(Mouse_TLS_combined@meta.data)
Mouse_TLS_combined <- AddMetaData(Mouse_TLS_combined, metadata = CellType)


# conserved and differential markers per cell state, e.g. Seurat_0
Idents(Mouse_TLS_combined) <- Mouse_TLS_combined@meta.data$CellType
conserved.markers.Seurat_13 <- FindConservedMarkers(Mouse_TLS_combined, ident.1 = "Seurat_0", grouping.var = "study")
Idents(Mouse_TLS_combined) <- Mouse_TLS_combined@meta.data$study
differential.markers <- FindMarkers(subset(Mouse_TLS_combined, cells=rownames(subset(Mouse_TLS_combined@meta.data, CellType %in% c('Seurat_0')))), ident.1 = "TLS", ident.2 = "WT", verbose = FALSE)
