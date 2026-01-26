import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.pyplot import rc_context

out_dir = './NM_mice/out'

adata = sc.read_h5ad(os.path.join(out_dir, 'integrated_T_NM_mice.h5ad'))
adata

marker_genes = [
    'Cd3d', 'Cd4', 'Cd8a', 
    'Ccr7', 'Lef1', 'Sell', 'Ccl5', 'Cxcr3', #naive 
    'Prf1', 'Gzmm', 'Gzmk', 'Gzmb', #effectors
    'Trdc', 'Rorc', #GD
    'Bcl6', 'Il21', 'Cxcr5', #Tfh
    'Gata3', 'Cxcr6', 'Il4', 'Pparg',  #Th2
    'Foxp3', 'Il2ra', 'Ctla4', #Treg
    'Tcf7', 'Bcl2', 'Sell', 'S1pr1', 'Ccr7', #nTregs
    'Tnfrsf4', 'Tnfrsf9', 'Cd83', 'Pdcd1', 'Ikzf2', #rTregs
    'Lgals1' ,'Itgae', 'S100a4', #NLT-like
    'Top2a', 'Mki67', #cycling
]

with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(
        adata, 
        color="leiden_scVI_0.6", 
        legend_loc='on data'
        )
sc.pl.dotplot(
  adata,
  marker_genes,
  groupby='leiden_scVI_0.6',
  dendrogram=False,
  cmap='Reds',
  swap_axes=False,
  use_raw=True,
  standard_scale='var'
 )

def subcluster_and_umap(adata, cluster_column, cluster_names, subcluster_res=0.1, pt_size=0.5):
  subcluster_column = f"{cluster_column}_{'_'.join(cluster_names)}"
  sc.tl.leiden(
    adata,
    restrict_to=(cluster_column, cluster_names),
    key_added=subcluster_column,
    resolution=subcluster_res
  )
  sc.pl.umap(
    adata,
    color=subcluster_column,
    frameon=False,
    legend_loc='on data',
    size=pt_size,
    palette='tab20'
  )

with rc_context({'figure.figsize': (4, 4)}):
  subcluster_and_umap(adata, 'leiden_scVI_0.6', ['2'], 0.2, pt_size=6)
with rc_context({'figure.figsize': (4, 4)}):
  subcluster_and_umap(adata, 'leiden_scVI_0.6', ['3'], 0.2, pt_size=6)

adata.obs['cell_type_V2'] = 'T other'
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.6'] == '0', 'CD4N', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.6'] == '1', 'CD8N', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.6_2'] == '2,0', 'CD4EM', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.6_2'] == '2,1', 'Tfh', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.6_3'] == '3,0', 'nTregs', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.6_3'] == '3,1', 'eTregs', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.6_3'] == '3,2', 'NLT-like Tregs', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.6'] == '4', 'CD8EM Gzmk+', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.6'] == '5', 'CD8EM', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.6'] == '6', 'T proliferating', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.6'] == '7', 'Th2', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.6'] == '9', 'γδ T', adata.obs['cell_type_V2'])

with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(
        adata, 
        color="cell_type_V2", 
        legend_loc='on data'
        )
sc.pl.dotplot(
  adata,
  marker_genes,
  groupby='cell_type_V2',
  dendrogram=False,
  cmap='Reds',
  swap_axes=False,
  use_raw=True,
  standard_scale='var'
 )

adata.obs.to_csv(os.path.join(out_dir, 'annotated_T_metadata_NM_mice.tsv'), sep='\t')
adata.write(os.path.join(out_dir, 'T_NM_mice.h5ad'))
