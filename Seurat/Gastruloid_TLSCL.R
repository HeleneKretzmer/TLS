require(Seurat)

# Load data
data_10X <- Read10X(data.dir = 'Gastruloid/outs/filtered_feature_bc_matrix')
Gastruloid <- CreateSeuratObject(counts=data_10X, project='Gastruloid', min.cells=3, min.features=200)

data_10X <- Read10X(data.dir = 'TLSCL/outs/filtered_feature_bc_matrix')
TLSCL <- CreateSeuratObject(counts=data_10X, project='TLSCL', min.cells=3, min.features=200)


# QC, normalize and scale
Gastruloid[["percent.mt"]] <- PercentageFeatureSet(Gastruloid, pattern = "^mt-")
TLSCL[["percent.mt"]] <- PercentageFeatureSet(TLSCL, pattern = "^mt-")
Gastruloid <- subset(Gastruloid, subset = nFeature_RNA > 3000 & nCount_RNA > 10000 & nCount_RNA < 80000 & percent.mt < 5)
TLSCL <- subset(TLSCL, subset = nFeature_RNA > 3000 & nCount_RNA > 10000 & nCount_RNA < 80000 & percent.mt < 5)

mGast.list <- list(Gastruloid, TLSCL)
mGast.list <- lapply(X = mGast.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
    x <- ScaleData(x, verbose = FALSE)
    x <- RunPCA(x, verbose = FALSE)
    x <- RunUMAP(x, dims = 1:30, n.neighbors=10)
    x <- FindNeighbors(x, dims = 1:20)
    x <- FindClusters(x, resolution = 0.5)

})
Gastruloid <- mGast.list[[1]]
TLSCL <- mGast.list[[2]]

# Cell cycle
s.genes <- cc.genes$s.genes
s.genes <- paste(toupper(substr(s.genes, 1, 1)), substr(tolower(s.genes), 2, nchar(s.genes)), sep="")
g2m.genes <- cc.genes$g2m.genes
g2m.genes <- paste(toupper(substr(g2m.genes, 1, 1)), substr(tolower(g2m.genes), 2, nchar(g2m.genes)), sep="")

Gastruloid <- CellCycleScoring(Gastruloid, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
TLSCL <- CellCycleScoring(TLSCL, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Gastruloid <- ScaleData(Gastruloid, vars.to.regress = c('S.Score', 'G2M.Score'), features = rownames(Gastruloid))
Gastruloid <- FindVariableFeatures(Gastruloid, verbose = FALSE)
Gastruloid <- RunPCA(Gastruloid, verbose = FALSE)
Gastruloid <- FindNeighbors(Gastruloid, dims = 1:10)
Gastruloid <- FindClusters(Gastruloid)
Gastruloid <- RunUMAP(object = Gastruloid, dims = 1:10)
Gastruloid <- CellCycleScoring(Gastruloid, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

TLSCL <- ScaleData(TLSCL, vars.to.regress = c('S.Score', 'G2M.Score'), features = rownames(TLSCL))
TLSCL <- FindVariableFeatures(TLSCL, verbose = FALSE)
TLSCL <- RunPCA(TLSCL, verbose = FALSE)
TLSCL <- FindNeighbors(TLSCL, dims = 1:10)
TLSCL <- FindClusters(TLSCL)
TLSCL <- RunUMAP(object = TLSCL, dims = 1:10)
TLSCL <- CellCycleScoring(TLSCL, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)



################################################################################
# Integration with TLS 120h
load("TLS.Robj")
TLS <- subset(TLS, cells=rownames(subset(TLS@meta.data, orig.ident == 'mGast_120h')))

# Normalize data and find variable features
Gastruloid <- SCTransform(Gastruloid, verbose = FALSE)
TLSCL <- SCTransform(TLSCL, verbose = FALSE)
TLS <- SCTransform(TLS, verbose = FALSE)

# Cell type classification using an integrated reference: TLS to Gastruloid
Idents(Gastruloid) <- Gastruloid@meta.data$seurat_clusters
Idents(TLS) <- TLS@meta.data$seurat_clusters

anchors <- FindTransferAnchors(reference = Gastruloid, query = TLS, dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = Gastruloid$seurat_clusters, dims = 1:30)
colnames(predictions) <- paste0('Gastruloid_', colnames(predictions))
TLS <- AddMetaData(TLS, metadata = predictions)


# Cell type classification using an integrated reference
anchors <- FindTransferAnchors(reference = TLS, query = Gastruloid, dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = TLS$seurat_clusters, dims = 1:30)
Gastruloid <- AddMetaData(Gastruloid, metadata = predictions)

anchors <- FindTransferAnchors(reference = TLS, query = TLSCL, dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = TLS$seurat_clusters, dims = 1:30)
TLSCL <- AddMetaData(TLSCL, metadata = predictions)


# Integration
TLS <- AddMetaData(TLS, metadata = data.frame(experiment = TLS@meta.data$orig.ident, row.names=rownames(TLS@meta.data)))
Gastruloid <- AddMetaData(Gastruloid, metadata = data.frame(experiment = Gastruloid@meta.data$orig.ident, row.names=rownames(Gastruloid@meta.data)))
TLSCL <- AddMetaData(TLSCL, metadata = data.frame(experiment = TLSCL@meta.data$orig.ident, row.names=rownames(TLSCL@meta.data)))

TLS <- AddMetaData(TLS, metadata = data.frame(TLS_cluster = TLS@meta.data$seurat_clusters, row.names=rownames(TLS@meta.data)))
Gastruloid <- AddMetaData(Gastruloid, metadata = data.frame(TLS_cluster = Gastruloid@meta.data$predicted.id, row.names=rownames(Gastruloid@meta.data)))
TLSCL <- AddMetaData(TLSCL, metadata = data.frame(TLS_cluster = TLSCL@meta.data$predicted.id, row.names=rownames(TLSCL@meta.data)))


# Normalization
mGast.list <- list(TLS, Gastruloid, TLSCL)
mGast.list <- lapply(X = mGast.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

# Integration
anchors <- FindIntegrationAnchors(object.list = mGast.list, reference = c(1), reduction = "rpca", dims = 1:30)
TLS_Gast_TLSCL <- IntegrateData(anchorset = anchors, dims = 1:30)

# Processing
TLS_Gast_TLSCL <- ScaleData(TLS_Gast_TLSCL, verbose = FALSE)
TLS_Gast_TLSCL <- RunPCA(TLS_Gast_TLSCL, verbose = FALSE)
TLS_Gast_TLSCL <- RunUMAP(TLS_Gast_TLSCL, dims = 1:30, n.neighbors=10)
TLS_Gast_TLSCL <- FindNeighbors(TLS_Gast_TLSCL, dims = 1:20)
TLS_Gast_TLSCL <- FindClusters(TLS_Gast_TLSCL, resolution = 0.5)

# Cell cycle
s.genes <- cc.genes$s.genes
s.genes <- paste(toupper(substr(s.genes, 1, 1)), substr(tolower(s.genes), 2, nchar(s.genes)), sep="")
g2m.genes <- cc.genes$g2m.genes
g2m.genes <- paste(toupper(substr(g2m.genes, 1, 1)), substr(tolower(g2m.genes), 2, nchar(g2m.genes)), sep="")

TLS_Gast_TLSCL <- CellCycleScoring(TLS_Gast_TLSCL, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
TLS_Gast_TLSCL <- ScaleData(TLS_Gast_TLSCL, vars.to.regress = c('S.Score', 'G2M.Score'), features = rownames(TLS_Gast_TLSCL))
TLS_Gast_TLSCL <- FindVariableFeatures(TLS_Gast_TLSCL, verbose = FALSE)
TLS_Gast_TLSCL <- RunPCA(TLS_Gast_TLSCL, verbose = FALSE)
TLS_Gast_TLSCL <- FindNeighbors(TLS_Gast_TLSCL, dims = 1:10)
TLS_Gast_TLSCL <- FindClusters(TLS_Gast_TLSCL)
TLS_Gast_TLSCL <- RunUMAP(object = TLS_Gast_TLSCL, dims = 1:10)
TLS_Gast_TLSCL <- CellCycleScoring(TLS_Gast_TLSCL, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


# Marker gene matrix
Idents(TLS_Gast_TLSCL) <- TLS_Gast_TLSCL@meta.data$seurat_clusters
TLS_Gast_TLSCL.markers <- FindAllMarkers(TLS_Gast_TLSCL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- TLS_Gast_TLSCL.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top50 <- TLS_Gast_TLSCL.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
