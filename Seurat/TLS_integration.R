require(Seurat)

# Load data
data_10X <- Read10X(data.dir = "TLS_96h/outs/filtered_feature_bc_matrix")
TLS_96h <- CreateSeuratObject(counts=data_10X, project="TLS_96h", min.cells=3, min.features=200)

data_10X <- Read10X(data.dir = "TLS_108h/outs/filtered_feature_bc_matrix")
TLS_108h <- CreateSeuratObject(counts=data_10X, project="TLS_108h", min.cells=3, min.features=200)

data_10X <- Read10X(data.dir = "TLS_120h/outs/filtered_feature_bc_matrix/")
TLS_120h <- CreateSeuratObject(counts=data_10X, project="TLS_120h", min.cells=3, min.features=200)

# Normalize and scale
TLS.list <- list(TLS_120h, TLS_108h, TLS_96h)

TLS.list <- lapply(X = TLS.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

# Integrate different time points
anchors <- FindIntegrationAnchors(object.list = TLS.list, reference = c(1, 2, 3), reduction = "rpca", dims = 1:30)
TLS.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
TLS.integrated <- ScaleData(TLS.integrated, verbose = FALSE)
TLS.integrated <- RunPCA(TLS.integrated, verbose = FALSE)
TLS.integrated <- RunUMAP(TLS.integrated, dims = 1:30, n.neighbors=10)

TLS <- TLS.integrated
TLS@meta.data$orig.ident <- factor(TLS@meta.data$orig.ident, levels=c('TLS_96h','TLS_108h','TLS_120h'))


# Regress out cell cycle genes, run dimension reductions and find clusters
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

TLS <- CellCycleScoring(TLS, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
TLS <- RunPCA(TLS, features = c(s.genes, g2m.genes))
TLS <- ScaleData(TLS, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(TLS))
TLS <- FindNeighbors(TLS, dims = 1:10)
TLS <- FindClusters(TLS, resolution = 0.5)
TLS <- RunUMAP(object = TLS, dims = 1:10)
TLS <- RunPCA(TLS, features = c(s.genes, g2m.genes))


# Filter for cells that survived scanpy based QC and remove cluster 10 and 11
BCs <- read.table("good_cells.csv")
TLS <- subset(TLS, cells=BCs$V1)
TLS <- subset(TLS, cells=rownames(subset(TLS@meta.data, !seurat_clusters %in% c("Seurat_10","Seurat_11"))))
