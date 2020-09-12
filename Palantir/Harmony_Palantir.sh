```R preparation
library(Matrix)
require(Seurat)

# Load organoid data
load('TLS.Robj')

TLS_120h <- subset(TLS, cells=rownames(subset(TLS@meta.data, orig.ident == 'TLS_120h')))
TLS_108h <- subset(TLS, cells=rownames(subset(TLS@meta.data, orig.ident == 'TLS_108h')))
TLS_96h <- subset(TLS, cells=rownames(subset(TLS@meta.data, orig.ident == 'TLS_96h')))

write.csv(file='TLS_cluster.csv', data.frame(BC=rownames(TLS@meta.data), Cluster=TLS@meta.data$seurat_clusters), quote=F)


# TLS 120h
matrix_dir = "TLS_120h/outs/filtered_feature_bc_matrix/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V2

TLS_120h_counts <- data.frame(t(mat))

rownames(TLS_120h_counts) <- gsub('-','_',rownames(TLS_120h_counts))
TLS_120h_counts <- TLS_120h_counts[which(rownames(TLS_120h_counts) %in% colnames(TLS_120h@assays$RNA@counts)), which(colnames(TLS_120h_counts) %in% rownames(TLS_120h@assays$RNA@counts))]

write.csv(file='TLS_120h_counts.csv', TLS_120h_counts, quote=F)


# TLS 108h
matrix_dir = "TLS_108h/outs/filtered_feature_bc_matrix/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V2

TLS_108h_counts <- data.frame(t(mat))
rownames(TLS_108h_counts) <- gsub('-1','_2',rownames(TLS_108h_counts))
TLS_108h_counts <- TLS_108h_counts[which(rownames(TLS_108h_counts) %in% colnames(TLS_108h@assays$RNA@counts)), which(colnames(TLS_108h_counts) %in% rownames(TLS_108h@assays$RNA@counts))]

write.csv(file='TLS_108h_counts.csv', TLS_108h_counts, quote=F)


# TLS 96h
matrix_dir = "TLS_96h/outs/filtered_feature_bc_matrix/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V2

TLS_96h_counts <- data.frame(t(mat))
rownames(TLS_96h_counts) <- gsub('-1','_3',rownames(TLS_96h_counts))
TLS_96h_counts <- TLS_96h_counts[which(rownames(TLS_96h_counts) %in% colnames(TLS_96h@assays$RNA@counts)), which(colnames(TLS_96h_counts) %in% rownames(TLS_96h@assays$RNA@counts))]

write.csv(file='TLS_96h_counts.csv', TLS_96h_counts, quote=F)
```



```ipython
import harmony
import palantir

# Plotting and miscellaneous imports
import os
import pandas as pd
from pandas.core.index import RangeIndex
import matplotlib
import matplotlib.pyplot as plt

# Initialize random seed
import random
random.seed(101)


# load data
csv_files = ['TLS_96h_counts.csv',
            'TLS_108h_counts.csv',
            '120h_counts.csv']
sample_names = ['TLS_96h', 'TLS_108h', 'TLS_120h']

counts = harmony.utils.load_from_csvs(csv_files, sample_names)


# Normalization, var gene selection and log transformation
norm_df = harmony.utils.normalize_counts(counts)
hvg_genes = harmony.utils.hvg_genes(norm_df)
data_df = harmony.utils.log_transform(norm_df.loc[:,hvg_genes])

# Harmony augmented affinity matrix
tp = pd.Series(index=data_df.index)
for t in ['TLS_96h', 'TLS_108h', 'TLS_120h']:
    cells = data_df.index[data_df.index.str.contains(t)]
    tp[cells] = t

timepoint_connections = pd.DataFrame(columns=[0, 1])
index = 0
timepoint_connections.loc[index, :] = ['TLS_96h', 'TLS_108h']; index += 1
timepoint_connections.loc[index, :] = ['TLS_108h', 'TLS_120h']; index += 1
timepoint_connections

aug_aff, aff = harmony.core.augmented_affinity_matrix(data_df, tp, timepoint_connections)

# Visualization using force directed layouts
layout = harmony.plot.force_directed_layout(aug_aff, data_df.index)
cluster = pd.read_csv('TLS_cluster.csv')
cluster.set_index('BC', inplace=True)
tmp = [c.split("h_")[1] for c in tp.index]
seurat_cluster = [cluster['Cluster']][0][tmp]
seurat_cluster.index = tp.index

# Palantir trajectory detection
dm_res = palantir.utils.run_diffusion_maps(aug_aff)
ms_data = palantir.utils.determine_multiscale_space(dm_res)
ms_data.index = data_df.index


## Run Palantir NMP 96h start cell
start_cell = 'TLS_96h_GACCCAGAGTTGGCGA_3'

r_res = palantir.core.run_palantir(ms_data, start_cell, num_waypoints=2500)
pr_res.branch_probs.columns = ['Neural', 'Somite']
palantir.plot.plot_palantir_results(pr_res, layout)


# MAGIC imputation
imp_df = palantir.utils.run_magic_imputation(data_df, dm_res)

## Visualizing gene expression
palantir.plot.plot_gene_expression(imp_df, layout, ['Sox2', 'T'])


# Gene expression trends
genes = ['Sox2', 'T',]
gene_trends = palantir.presults.compute_gene_trends(pr_res, imp_df.loc[:, genes])
palantir.plot.plot_gene_trends(gene_trends)

genes = ['Sox2', 'T', 'Sox1', 'Tbx6', 'Rspo3', 'Wnt3a', 'Cyp26a1', 'Hoxc10', 'Hes7', 'Pax6', 'Fgf8', 'Cdx2']
gene_trends = palantir.presults.compute_gene_trends( pr_res, imp_df.loc[:, genes])
palantir.plot.plot_gene_trends(gene_trends)


# Branch probabilities per cluster and time point
branch_probs = pd.concat([pd.DataFrame(pr_res2.branch_probs), pd.DataFrame(pr_res2.entropy, columns=['diff_pot']), pd.DataFrame(seurat_cluster), pd.DataFrame(tp, columns=['TP']), pd.DataFrame(layout['x']), pd.DataFrame(layout['y']), pd.DataFrame(counts[['Sox2','T']])], axis=1, sort=False)

# Differentiation potential per cell, cluster and time point
diff_potential = pd.concat([pd.DataFrame(pr_res2.branch_probs), pd.DataFrame(pr_res2.entropy, columns=['diff_pot']), pd.DataFrame(seurat_cluster), pd.DataFrame(tp, columns=['TP']), pd.DataFrame(layout['x']), pd.DataFrame(layout['y'])], axis=1, sort=False)
```
