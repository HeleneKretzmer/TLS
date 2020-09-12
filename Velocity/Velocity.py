#!/usr/bin/env python
# coding: utf-8

# # Load ibraries
import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb
from gprofiler import gprofiler
import rpy2.rinterface_lib.callbacks
import logging
from rpy2.robjects import pandas2ri
import anndata2ri

rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)

pandas2ri.activate()
anndata2ri.activate()
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')

plt.rcParams['figure.figsize']=(8,8) #rescale figures
sc.settings.verbosity = 3
sc.logging.print_versions()

get_ipython().run_cell_magic('R', '', '# Load all the R libraries we will be using in the notebook\nlibrary(scran)\nlibrary(RColorBrewer)\nlibrary(slingshot)\nlibrary(monocle)\nlibrary(gam)\nlibrary(clusterExperiment)\nlibrary(ggplot2)\nlibrary(plyr)\nlibrary(MAST)')


# Load data
data_file = 'TLS_120h/outs/filtered_feature_bc_matrix/matrix.mtx.gz'
barcode_file = 'TLS_120h/outs/filtered_feature_bc_matrix/barcodes.tsv.gz'
gene_file = 'TLS_120h/outs/filtered_feature_bc_matrix/features.tsv.gz'

adata = sc.read(data_file, cache=True)
adata = adata.transpose()
adata.X = adata.X.toarray()
barcodes = pd.read_csv(barcode_file, header=None, sep='\t')
genes = pd.read_csv(gene_file, header=None, sep='\t')
genes = genes.drop(2, axis=1)

barcodes.rename(columns={0:'barcode'}, inplace=True)
barcodes.set_index('barcode', inplace=True)
adata.obs = barcodes
adata.obs['sample'] = 'TLS_120h'
adata.obs['region'] = 'TLS'
adata.obs['donor'] = 'TLS_120h'

genes.rename(columns={0:'gene_id', 1:'gene_symbol'}, inplace=True)
genes.set_index('gene_symbol', inplace=True)
adata.var = genes

adata_120h = adata


# Load data
data_file = 'TLS_96h/outs/filtered_feature_bc_matrix/matrix.mtx.gz'
barcode_file = 'TLS_96h/outs/filtered_feature_bc_matrix/barcodes.tsv.gz'
gene_file = 'TLS_96h/outs/filtered_feature_bc_matrix/features.tsv.gz'

adata = sc.read(data_file, cache=True)
adata = adata.transpose()
adata.X = adata.X.toarray()
barcodes = pd.read_csv(barcode_file, header=None, sep='\t')
genes = pd.read_csv(gene_file, header=None, sep='\t')
genes = genes.drop(2, axis=1)

barcodes.rename(columns={0:'barcode'}, inplace=True)
barcodes.set_index('barcode', inplace=True)
adata.obs = barcodes
adata.obs['sample'] = 'TLS_96h'
adata.obs['region'] = 'TLS'
adata.obs['donor'] = 'TLS_96h'

genes.rename(columns={0:'gene_id', 1:'gene_symbol'}, inplace=True)
genes.set_index('gene_symbol', inplace=True)
adata.var = genes

adata_96h = adata


# Load data
data_file = 'TLS_108h/outs/filtered_feature_bc_matrix/matrix.mtx.gz'
barcode_file = 'TLS_108h/outs/filtered_feature_bc_matrix/barcodes.tsv.gz'
gene_file = 'TLS_108h/outs/filtered_feature_bc_matrix/features.tsv.gz'

#Load data
adata = sc.read(data_file, cache=True)
adata = adata.transpose()
adata.X = adata.X.toarray()
barcodes = pd.read_csv(barcode_file, header=None, sep='\t')
genes = pd.read_csv(gene_file, header=None, sep='\t')
genes = genes.drop(2, axis=1)

#Annotate data
barcodes.rename(columns={0:'barcode'}, inplace=True)
barcodes.set_index('barcode', inplace=True)
adata.obs = barcodes
adata.obs['sample'] = 'TLS_108h'
adata.obs['region'] = 'TLS'
adata.obs['donor'] = 'TLS_108h'

genes.rename(columns={0:'gene_id', 1:'gene_symbol'}, inplace=True)
genes.set_index('gene_symbol', inplace=True)
adata.var = genes

adata_108h = adata


# Concatenate to main adata object
adata = adata_120h.concatenate(adata_108h, batch_key='sample')
adata.var['gene_id'] = adata.var['gene_id-1']
adata.var.drop(columns = ['gene_id-1', 'gene_id-0'], inplace=True)
adata.obs_names = [c.split("-")[0] for c in adata.obs_names]
adata.obs_names_make_unique(join='_')
adata.var_names_make_unique()
adata.obs_names_make_unique()

adata = adata.concatenate(adata_96h, batch_key='sample')
adata.var['gene_id'] = adata.var['gene_id-1']
adata.var.drop(columns = ['gene_id-1', 'gene_id-0'], inplace=True)
adata.obs_names = [c.split("-")[0] for c in adata.obs_names]
adata.obs_names_make_unique(join='_')
adata.var_names_make_unique()
adata.obs_names_make_unique()

adata.obs['sample'] = adata.obs['donor']
adata.var_names_make_unique()
adata.obs_names_make_unique()


# Add UMAP from Seurat clustering
meta = pd.read_csv('UMAP.tsv', sep='\t')
adata = adata[[i for i,x in enumerate(adata.obs.index) if x in np.array([meta['BC']])],]
meta.set_index('BC', inplace=True)
adata.obs['UMAP1'] = [meta['UMAP_1']][0][adata.obs.index]
adata.obs['UMAP2'] = [meta['UMAP_2']][0][adata.obs.index]


# Pre-processing and visualization
adata.obs['n_counts'] = adata.X.sum(1)
adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
adata.obs['n_genes'] = (adata.X > 0).sum(1)

mt_gene_mask = [gene.startswith('mt-') for gene in adata.var_names]
adata.obs['mt_frac'] = adata.X[:, mt_gene_mask].sum(1)/adata.obs['n_counts']

# Filter cells according to identified QC thresholds:
sc.pp.filter_cells(adata, min_counts = 10000)
sc.pp.filter_cells(adata, max_counts = 40000)
adata = adata[adata.obs['mt_frac'] < 0.1]
sc.pp.filter_cells(adata, min_genes = 3000)
adata.var_names_make_unique()
sc.pp.filter_genes(adata, min_cells=20)


# Normalization
adata_pp = adata.copy()
sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
sc.pp.log1p(adata_pp)
sc.pp.pca(adata_pp, n_comps=15)
sc.pp.neighbors(adata_pp)
sc.tl.louvain(adata_pp, key_added='groups', resolution=1)

input_groups = adata_pp.obs['groups']
data_mat = adata.X.T
get_ipython().run_cell_magic('R', '-i data_mat -i input_groups -o size_factors', '\nsize_factors = computeSumFactors(data_mat, clusters=input_groups, min.mean=0.1)')
del adata_pp

adata.obs['size_factors'] = size_factors
adata.layers["counts"] = adata.X.copy()

adata.X /= adata.obs['size_factors'].values[:,None]
sc.pp.log1p(adata)

adata.raw = adata


# Highly Variable Genes
sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=4000)

# Visualization
sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack')
sc.pp.neighbors(adata)

sc.tl.tsne(adata)
sc.tl.umap(adata)
sc.tl.diffmap(adata)
sc.tl.draw_graph(adata)

adata.obsm['X_umap'] = np.array(adata.obs[['UMAP1','UMAP2']])


# Cell cycle scoring
cc_genes = pd.read_table('Macosko_cell_cycle_genes.txt', delimiter='\t')
s_genes = cc_genes['S'].dropna()
g2m_genes = cc_genes['G2.M'].dropna()

s_genes_mm = [gene.lower().capitalize() for gene in s_genes]
g2m_genes_mm = [gene.lower().capitalize() for gene in g2m_genes]

s_genes_mm_ens = adata.var_names[np.in1d(adata.var_names, s_genes_mm)]
g2m_genes_mm_ens = adata.var_names[np.in1d(adata.var_names, g2m_genes_mm)]

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes_mm_ens, g2m_genes=g2m_genes_mm_ens)

sc.pl.umap(adata, color=['S_score', 'G2M_score'], use_raw=False)
sc.pl.umap(adata, color='phase', use_raw=False)
sc.pp.regress_out(adata, ['S_score', 'G2M_score'])

sc.pp.scale(adata)
sc.pp.pca(adata, n_comps=50, svd_solver='arpack')
sc.pp.neighbors(adata)

sc.tl.tsne(adata)
sc.tl.umap(adata)
sc.tl.diffmap(adata)
sc.tl.draw_graph(adata)

adata.obsm['X_umap'] = np.array(adata.obs[['UMAP1','UMAP2']])


# ## 3.1 Clustering
# Perform clustering - using highly variable genes
sc.tl.louvain(adata, resolution=1, key_added='louvain_r1')
adata.obs['louvain_r1'].value_counts()

# use seurat clustering
adata.obs['louvain']  = adata.obs['louvain_r1']

meta = pd.read_csv('TLS_cluster.tsv', sep='\t')
meta.set_index('BC', inplace=True)
adata.obs['louvain'] = [meta['seurat_clusters']][0][adata.obs.index]




# velocity
## Read and merge velocity
adata_loom_120 = scv.read("TLS_120h/velocyto/TLS_120h.loom", sparse=True, cache=True)
adata_loom_120.var_names_make_unique()
adata_loom_108 = scv.read("TLS_108h/velocyto/TLS_108h.loom", sparse=True, cache=True)
adata_loom_108.var_names_make_unique()
adata_loom_96 = scv.read("TLS_96h/velocyto/TLS_96h.loom", sparse=True, cache=True)
adata_loom_96.var_names_make_unique()

loom = adata_loom_120.concatenate(adata_loom_108)
loom.var_names_make_unique()
loom.obs_names_make_unique()
loom = loom.concatenate(adata_loom_96)
loom.var_names_make_unique()
loom.obs_names_make_unique()

## merge loom file into an already existing AnnData object
adata = scv.utils.merge(adata, loom)
adata.var_names_make_unique()
adata.obs_names_make_unique()

scv.utils.show_proportions(adata)

## Preprocess the data
scv.pp.filter_and_normalize(adata, min_counts=20, min_counts_u=10, n_top_genes=3000)
scv.pp.moments(adata, n_pcs=50, n_neighbors=30)

## Compute velocity and velocity graph
scv.tl.velocity(adata)

## Compute the (cosine) correlation of potential cell transitions with the velocity vector in high dimensional space
scv.tl.velocity_graph(adata)

## Project the velocity graph onto an embedding
scv.tl.velocity_embedding(adata, basis='umap')

## Plot results
scv.pl.velocity_embedding_stream(adata, legend_loc='on data', alpha=.05, color='louvain')
scv.pl.velocity_embedding(adata, basis='umap', dpi=600, color='louvain')
scv.pl.velocity_embedding_grid(adata, color='Tmsb10', layer=['velocity', 'spliced'], colorbar=True)
