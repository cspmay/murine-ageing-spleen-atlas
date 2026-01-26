import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.pyplot import rc_context

out_dir = './atlas/out/'

adata = sc.read_h5ad(os.path.join(out_dir, 'integrated_T_cleaned_datasets.h5ad'))

T_marker_genes = [
    'Cd3d', 'Cd4', 'Cd8a', 
    'Ccr7', 'Lef1', 'Sell', 'Ccl5', 'Cxcr3', #naive 
    'Prf1', 'Gzmm', 'Gzmk', 'Gzmb', #effectors
    'Trdc', 'Rorc', #GD
    'Bcl6', 'Il21', 'Cxcr5', #Tfh
    'Gata3', 'Cxcr6', 'Il4', 'Pparg',  #Th2
    'Foxp3', 'Il2ra', 'Ctla4', #Treg
    'Tcf7', 'Bcl2', 'Sell', 'S1pr1', 'Ccr7', #cTregs
    'Tnfrsf4', 'Tnfrsf9', 'Cd83', 'Pdcd1', 'Ikzf2', #eTregs
    'Lgals1' ,'Itgae', 'S100a4', #NLT-like
    'Top2a', 'Mki67', #cycling
]

with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(
        adata, 
        color="leiden_scVI_0.4", 
        legend_loc='on data',
        )
sc.pl.dotplot(
  adata,
  marker_genes,
  groupby='leiden_scVI_0.4',
  dendrogram=False,
  cmap='Reds',
  swap_axes=False,
  use_raw=True,
  standard_scale='var',
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
  subcluster_and_umap(adata, 'leiden_scVI_0.4', ['2'], 0.3, pt_size=6)
with rc_context({'figure.figsize': (4, 4)}):
  subcluster_and_umap(adata, 'leiden_scVI_0.4', ['4'], 0.2, pt_size=6)

adata.obs['cell_type_V2'] = 'T other'
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.4'] == '0', 'CD8N', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.4'] == '1', 'CD4N', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.4_2'] == '2,0', 'CD4EM', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.4_2'] == '2,1', 'CD4EM', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.4_2'] == '2,2', 'Tfh', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.4_2'] == '2,3', 'CD4EM', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.4_2'] == '2,4', 'CD4EM', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.4_2'] == '2,5', 'CD4EM', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.4'] == '3', 'CD8EM', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.4_4'] == '4,0', 'eTregs', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.4_4'] == '4,1', 'nTregs', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.4_4'] == '4,2', 'NLT-like Tregs', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.4'] == '5', 'CD8EM Gzmk+', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.4'] == '6', 'Th2', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.4'] == '8', 'T proliferating', adata.obs['cell_type_V2'])

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

adata.obs.to_csv(os.path.join(out_dir, 'annotated_T_metadata.tsv'), sep='\t')
adata.write(os.path.join(out_dir, 'T_cell_atlas.h5ad'))
