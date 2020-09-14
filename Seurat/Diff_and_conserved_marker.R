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

#load Robj

load("invivoWT_TLS.Robj")

#set default assay to RNA
DefaultAssay(combined) <- "RNA"

#remove XY genes
#this is done to prevent that XY linked genes will hijack the DEG analysis
XY_genes <- read.table('X_Y_genes.txt')
XY_genes <- XY_genes$V1

combined@assays$RNA@counts <- combined@assays$RNA@counts[!rownames(combined@assays$RNA@counts) %in% XY_genes, ]
combined@assays$RNA@data   <- combined@assays$RNA@data[!rownames(combined@assays$RNA@data) %in% XY_genes, ]

#Subset in vivo celltypes

combined2 <- subset(combined, cells=rownames(subset(combined@meta.data, Mouse_CellType %in% c("18","18.1","7","26","31","39","30","27","2"))))

##Define idents

Idents(combined2) <- combined2@meta.data$Mouse_CellType

combined2

##Select for genes expressed above threshold X
combined3 <- combined2
DefaultAssay(combined3) <- "RNA"

Idents(combined3) <- combined2@meta.data$Mouse_CellType

combined4 <- combined2
DefaultAssay(combined4) <- "RNA"

Idents(combined4) <- combined4@meta.data$study

#calculate average:
subset18 <- subset(combined3, idents = "18")
Idents(subset18) <- c("TLS","WT")
avg.subset18 <- AverageExpression(subset18, verbose = FALSE)$RNA

#call max
avg.subset18[, "max"] <- apply(avg.subset18, 1, max)
head(avg.subset18)

#subset based on max
df18 <- avg.subset18
df18 <- subset(df18, df18$max > 1)
genes.to.use <- row.names(df18)

#calculate differential markers all expressed

differential.markers.18 <- FindMarkers(subset(combined4, cells=rownames(subset(combined4@meta.data, Mouse_CellType %in% c('18')))),
ident.1 = "TLS", ident.2 = "WT", verbose = FALSE, assay.type="RNA",logfc.threshold=1.0, min.diff.pct=0.25, pseudocount.use=0.1)
nrow(differential.markers.18)

#calculate total number of expressed genes in cluster

subset18 <- subset(combined3, idents = "18")
Idents(subset18) <- c("TLS","WT")
avg.subset18 <- AverageExpression(subset18, verbose = FALSE)$RNA
#call max
avg.subset18[, "max"] <- apply(avg.subset18, 1, max)
head(avg.subset18)
df18 <- avg.subset18
df18 <- subset(df18, df18$max > 0)
nrow(df18)

#calculate differential markers top 20 in vivo markers

genes.to.use18 <- read.table('/project/Organoid/markers18.txt', row.names=1)
genes.to.use18 <- genes.to.use18[!rownames(genes.to.use18) %in% XY_genes,]
differential.markers.18 <- FindMarkers(subset(combined4, cells=rownames(subset(combined4@meta.data, Mouse_CellType %in% c('18')))),
features = rownames(genes.to.use18[1:20,]), ident.1 = "TLS", ident.2 = "WT", verbose = FALSE, assay.type="RNA",logfc.threshold=1.0, pseudocount.use=0.1, min.diff.pct=0.25)
nrow(genes.to.use18)
nrow(differential.markers.18)

#calculate conserved markers
Conserved.18 <- FindConservedMarkers(combined3, ident.1="18", grouping.var="study", min.diff.pct=0.25)
