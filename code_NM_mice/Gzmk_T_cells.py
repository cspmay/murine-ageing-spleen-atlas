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

adata = sc.read_h5ad(os.path.join(out_dir, 'integrated_Gzmk_NM_mice.h5ad')
adata

marker_genes = [
    'Gzmk', 'Pdcd1', 'Tox', 'Itga4', 'Ccl5', 'Ikzf3', 'Runx2', 'Chst11', 'Myo1f', 'Fyn', 'Atxn1', 'Zeb1',
    'Zeb2', 'Hip1', 'Chsy1', 'Sytl3', 'Dock5', 'Il18rap', 'Rora', 'Itgal'
]

with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(
        adata, 
        color="leiden_scVI_0.3", 
        legend_loc='on data',
        legend_fontsize=8, 
        )
sc.pl.dotplot(
  adata,
  marker_genes,
  groupby='leiden_scVI_0.3',
  dendrogram=False,
  cmap='Reds',
  swap_axes=False,
  use_raw=True,
  standard_scale='var',
 )

adata.obs['cell_type_V2'] = 'T other'
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.3'] == '0', 'CD8EM Gzmk-lo', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.3'] == '1', 'CD8EM Gzmk-hi', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.3'] == '2', 'CD8EM Gzmk-hi', adata.obs['cell_type_V2'])

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

adata.obs.to_csv(os.path.join(out_dir, 'annotated_Gzmk_metadata_NM_mice.tsv'), sep='\t')
adata.write(os.path.join(out_dir, 'Gzmk_NM_mice.h5ad'))
