##load scanpy
import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb
from gprofiler import gprofiler
import scvelo as scv
scv.logging.print_version()
import rpy2.rinterface_lib.callbacks
import logging
from rpy2.robjects import pandas2ri
import anndata2ri
# Automatically convert rpy2 outputs to pandas dataframes
pandas2ri.activate()
anndata2ri.activate()
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')
plt.rcParams['figure.figsize']=(8,8) #rescale figures
plt.rcParams.update({'font.size': 40})#rescale fontsize
sc.settings.verbosity = 3
#sc.set_figure_params(dpi=200, dpi_save=300)
sc.logging.print_versions()
get_ipython().run_cell_magic('R', '', '# Load all the R libraries we will be using in the notebook\nlibrary(scran)\nlibrary(RColorBrewer)\nlibrary(slingshot)\nlibrary(monocle)\nlibrary(gam)\nlibrary(clusterExperiment)\nlibrary(ggplot2)\nlibrary(plyr)\nlibrary(MAST)')
# scVelo
scv.settings.set_figure_params('scvelo')
# colors
colors2 = plt.cm.Reds(np.linspace(0, 1, 128))
colors3 = plt.cm.Greys_r(np.linspace(0.7,0.8,20))
colorsComb = np.vstack([colors3, colors2])
mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)

##load data
adata_all = sc.read(“TLS.h5ad", sparse=True)

adata_120h = adata_all[adata_all.obs['donor'].isin(['TLS_120h'])]
adata_108h = adata_all[adata_all.obs['donor'].isin(['TLS_108h'])]
adata_96h = adata_all[adata_all.obs['donor'].isin(['TLS_96h'])]

remove = ['Seurat_10', 'Seurat_11']
adata_all_final = adata_all[~adata_all.obs['louvain'].isin(remove)]
sc.pl.umap(adata_all_final, color=['louvain'], palette=sc.pl.palettes.default_64)

remove = ['Seurat_10', 'Seurat_11']
adata_96h_final = adata_96h[~adata_96h.obs['louvain'].isin(remove)]
sc.pl.umap(adata_96h_final, color=['louvain'], palette=sc.pl.palettes.default_64)

remove = ['Seurat_11']
adata_108h_final = adata_108h[~adata_108h.obs['louvain'].isin(remove)]
sc.pl.umap(adata_108h_final, color=['louvain'], palette=sc.pl.palettes.default_64)


remove = ['Seurat_10', 'Seurat_11']
adata_120h_final = adata_120h[~adata_120h.obs['louvain'].isin(remove)]
sc.pl.umap(adata_120h_final, color=['louvain'], palette=sc.pl.palettes.default_64)

# Pseudotime reconstruction
## PAGA P1 = adata_120h
scv.tl.paga(P1, groups='louvain', use_rna_velocity=False)
# Reconstructing gene changes along PAGA paths for a given set of genes
##120HOURS woSubclustering and with merging somites
##mergesomites
sc.tl.louvain(adata_120h, restrict_to=('louvain', ["Seurat_1","Seurat_3", "Seurat_8", "Seurat_0"]), resolution=0, key_added='louvain2')
sc.pl.umap(adata_120h, color=['louvain2'], palette=sc.pl.palettes.default_64,legend_fontsize=20)

## PAGA P1 = adata_120h
scv.tl.paga(P1, groups='louvain2', use_rna_velocity=False)
# Reconstructing gene changes along PAGA paths for a given set of genes
## Somitogenesis Trajectory P1.uns['iroot'] = np.flatnonzero(P1.obs['louvain2']  == "Seurat_4")[0] sc.tl.dpt(P1) sc.pl.diffmap(P1, components='1,2', color='dpt_pseudotime')  c.pl.umap(P1, color=['dpt_pseudotime'], legend_loc='on data') 
paths = [('p1', ["Seurat_4", "Seurat_7", "Seurat_2", "Seurat_5","Seurat_1-Seurat_3-Seurat_8-Seurat_0,0"])]
P1.obs['distance'] = P1.obs['dpt_pseudotime'] P1.raw = sc.AnnData(sc.pp.scale(P1.X, copy=True), var=P1.raw.var)  import matplotlib.pyplot as pl

gene_names = ["Hoxc10","Nkx1-2", "Sox2", "Cyp26a1","Fgf8", "T","Tbx6", "Hes7","Rspo3","Dll1","Mesp2","Ripply2","Cer1","Aldh1a2","Meox1", "Tcf15","Pax3","Pax1"]
##Genes To Plot Somitic Trajectory
_, axs = pl.subplots(ncols=2, figsize=(6, 2.5), gridspec_kw={'wspace': 0.05, 'left': 0.12}) pl.subplots_adjust(left=0.05, right=0.98, top=0.82, bottom=0.2) for ipath, (descr, path) in enumerate(paths):   _, axs = pl.subplots(ncols=2, figsize=(6, 2.5), gridspec_kw={'wspace': 0.05, 'left': 0.12})  p.subplots_adjust(left=0.05, right=0.98, top=0.82, bottom=0.2)  fo ipath, (descr, path) in enumerate(paths):  _, data = sc.pl.paga_path(         P1, path, gene_names,         show_node_names=False,         ax=axs[ipath],         ytick_fontsize=12,         left_margin=0.5,          _avg=50,          anotations=['distance'],          shw_yticks=True if ipath==0 else False,          sho_colorbar=True, 
normalize_to_zero_one=True,
         color_map='magma',         color_maps_annotations={'distance': 'viridis'},         title='{} path'.format(descr),         return_data=True,         show=False)  p.show()

## Neural Trajectory
paths = [('p1', ["Seurat_4", "Seurat_9", "Seurat_6"])] P1.obs['distance'] = P1.obs['dpt_pseudotime'] P1.raw = sc.AnnData(sc.pp.scale(P1.X, copy=True), var=P1.raw.var)
gene_names = ["T", "Sox2","Pax6", "Dbx1","Dbx2"]
##Genes To Plot Neural Trajectory
_, axs = pl.subplots(ncols=2, figsize=(6, 2.5), gridspec_kw={'wspace': 0.05, 'left': 0.12}) pl.subplots_adjust(left=0.05, right=0.98, top=0.82, bottom=0.2) for ipath, (descr, path) in enumerate(paths):   _, axs = pl.subplots(ncols=2, figsize=(6, 2.5), gridspec_kw={'wspace': 0.05, 'left': 0.12})  p.subplots_adjust(left=0.05, right=0.98, top=0.82, bottom=0.2)  fo ipath, (descr, path) in enumerate(paths):  _, data = sc.pl.paga_path(         P1, path, gene_names,         show_node_names=False,         ax=axs[ipath],         ytick_fontsize=12,         left_margin=0.5,          _avg=50,          anotations=['distance'],          shw_yticks=True if ipath==0 else False,          sho_colorbar=True, 
normalize_to_zero_one=True,
         color_map='magma',         color_maps_annotations={'distance': 'viridis'},         title='{} path'.format(descr),         return_data=True,         show=False)  p.show()

##Define Hox Genes
hoxa = ['Hoxa1','Hoxa2','Hoxa3','Hoxa4','Hoxa5','Hoxa6','Hoxa7','Hoxa9','Hoxa10']
hoxb = ['Hoxb1','Hoxb2','Hoxb3','Hoxb4','Hoxb5','Hoxb6','Hoxb7','Hoxb8','Hoxb9','Hoxb13']
hoxc = ['Hoxc4','Hoxc5','Hoxc6','Hoxc8','Hoxc9','Hoxc10']
hoxd = ['Hoxd1','Hoxd3','Hoxd4','Hoxd8','Hoxd9','Hoxd10','Hoxd13']

##Plot Hox Genes pseudotime for NMPs, meso and neuro progenitors
include = ["Seurat_9", "Seurat_4", "Seurat_7"]
adata_prog = adata_all_final[adata_all_final.obs['louvain'].isin(include)]

P1 = adata_prog
scv.tl.paga(P1, groups='donor', use_rna_velocity=False)
 P1.uns['iroot'] = np.flatnonzero(P1.obs['donor']  == 'TLS_96h' )[0]
sc.tl.dpt(P1) sc.pl.diffmap(P1, components='1,2', color='dpt_pseudotime') sc.pl.umap(P1, color=['dpt_pseudotime'], legend_loc='on data') 
paths = [('Time', ['TLS_96h', 'TLS_108h', 'TLS_120h'])]
P1.obs['distance'] = P1.obs['dpt_pseudotime'] P1.raw = sc.AnnData(sc.pp.scale(P1.X, copy=True), var=P1.raw.var)

## example HoxC
_, axs = pl.subplots(ncols=2, figsize=(6, 2.5), gridspec_kw={'wspace': 0.05, 'left': 0.12}) pl.subplots_adjust(left=0.05, right=0.98, top=0.82, bottom=0.2) for ipath, (descr, path) in enumerate(paths):   _, axs = pl.subplots(ncols=2, figsize=(6, 2.5), gridspec_kw={'wspace': 0.05, 'left': 0.12})  p.subplots_adjust(left=0.05, right=0.98, top=0.82, bottom=0.2)  fo ipath, (descr, path) in enumerate(paths):  _, data = sc.pl.paga_path(         P1, path, hoxc,         show_node_names=False,         ax=axs[ipath],         ytick_fontsize=12,         left_margin=0.5,          _avg=50,          anotations=['distance'],          shw_yticks=True if ipath==0 else False,          sho_colorbar=True, 
normalize_to_zero_one=True,
         color_map='magma',         color_maps_annotations={'distance': 'viridis'},         title='{} path'.format(descr),         return_data=True,         show=False)  p.show()

##Subclustering Neural Tube Cluster
neuro2 = ['Seurat_6']

adata_neuro2 = adata_all_final[adata_all_final.obs['louvain'].isin(neuro2)]

sc.pl.umap(adata_neuro2, color=['louvain'], palette=sc.pl.palettes.default_64)

adata_neuro2 = adata_neuro2[adata_neuro2.obs['donor'].isin(['TLS_120h'])]
sc.tl.louvain(adata_neuro2, restrict_to=('louvain', ['Seurat_6']), resolution=0.65, key_added='louvain_sub')

sc.pl.umap(adata_neuro2, color=['louvain_sub'], palette=sc.pl.palettes.default_64, legend_fontsize=20)

g2 = ["T","Nkx1-2","Igfbp5","Apela","Cdkn1c","Greb1l","Hoxc9","Irx3","Msx1","Msx3","Pax3","Id1","Id2","Id3","Wnt4","Prdm13","Zic1","Zic5","Fzd10","Pax7","Pax6","Prdm8","Dbx1","Dbx2","Pou3f2","Mir124-2hg","Sp8","Prdm12","Hes3","Nkx6-2","Ptch1","Olig2","Nkx6-1","Neurog2","Bmp7"]

sc.pl.matrixplot(adata_neuro2, var_names=g2, groupby='louvain_sub', standard_scale='var', cmap='Reds')

##Subclustering Somite Cluster 120h
somite = ['Seurat_0','Seurat_1','Seurat_8','Seurat_3']
adata_somite = adata_all_final[adata_all_final.obs['louvain'].isin(somite)]
sc.pl.umap(adata_somite, color=['louvain'], palette=sc.pl.palettes.default_64)

sc.tl.louvain(adata_somite, restrict_to=('louvain', ['Seurat_0','Seurat_1','Seurat_8','Seurat_3']), resolution=0.3, key_added='louvain_sub')
sc.pl.umap(adata_somite, color=['louvain_sub'], palette=sc.pl.palettes.default_64, legend_fontsize=20)

sc.tl.umap(adata_somite)
sc.pl.umap(adata_somite, color='louvain_sub')

adata_somite_120 = adata_somite[adata_somite.obs['donor'].isin(['TLS_120h'])]
sc.tl.louvain(adata_somite_120, restrict_to=('louvain', ['Seurat_0','Seurat_1','Seurat_8','Seurat_3']), resolution=0.4, key_added='louvain_sub')
sc.pl.umap(adata_somite_120, color=['louvain_sub'], palette=sc.pl.palettes.default_64, legend_fontsize=20)
sc.tl.umap(adata_somite_120)
sc.pl.umap(adata_somite_120, color='louvain_sub')

##DEGs per subcluster
sc.tl.rank_genes_groups(adata_somite_120, groupby='louvain_sub', key_added='test',  rankby_abs=False, use_raw=True)
sc.pl.rank_genes_groups_matrixplot(adata_somite_120, key='test', dendrogram=False, standard_scale='var', cmap='Reds')

##Subclustering Somite Cluster 96h
adata_somite_96 = adata_somite[adata_somite.obs['donor'].isin(['TLS_96h'])]
sc.tl.louvain(adata_somite_96, restrict_to=('louvain', ['Seurat_0','Seurat_1','Seurat_8','Seurat_3']), resolution=0.3, key_added='louvain_sub')
sc.pl.umap(adata_somite_96, color=['louvain_sub'], palette=sc.pl.palettes.default_64, legend_fontsize=20)
sc.tl.umap(adata_somite_96)
sc.pl.umap(adata_somite_96, color='louvain_sub')
remove96 = ['Seurat_0-Seurat_1-Seurat_8-Seurat_3,2','Seurat_0-Seurat_1-Seurat_8-Seurat_3,3']
adata_somite_96f = adata_somite_96 [~adata_somite_96.obs['louvain_sub'].isin(remove96)]
sc.pl.umap(adata_somite_96f, color=['louvain_sub'], palette=sc.pl.palettes.default_64)
