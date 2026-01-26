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
from adjustText import adjust_text

out_dir = './atlas/out/'

adata = sc.read_h5ad(os.path.join(out_dir, 'integrated_B_datasets.h5ad'))

B_marker_genes = [
	'Cd79a', 'Cd19',
    'Apoe', 'B2m', 'Plac8', #Apoe
    'Zeb2', 'Itgax', 'Tbx21', #ABCs 
    'Fcer2a', 'Cxcr4', #Fo
    'Cd83', 'Hspd1', 'Myc', #GC
    'Ccdc28b', 'Ptpn22', 'Zbtb20', 'Fcrl5', 'Il10', 'Ctla4', 'Sox5', 'Myc', 'Atf3', #Bregs
    'Cd9', 'Mzb1', 'Cr2', 'Cd1d1', #MZ
    'Pax5', 'Bcl2', 'Cd80', #Mem
    'Cd93', 'Cd24a', 'Vpreb3', #Tr
    'Top2a', 'Mki67', 'Pcna', #cycling
    ]

with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(
        adata, 
        color="leiden_scVI_0.5", 
        legend_loc='on data'
        )
sc.pl.dotplot(
  adata,
  B_marker_genes,
  groupby='leiden_scVI_0.5',
  dendrogram=False,
  cmap='Reds',
  swap_axes=False,
  use_raw=True,
  standard_scale='var',
 )

adata.obs['cell_type_V2'] = 'B other'
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5'] == '0', 'FoB', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5'] == '1', 'Memory B', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5'] == '2', 'FoB', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5'] == '3', 'MZB', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5'] == '4', 'GC B', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5'] == '5', 'TrB', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5'] == '6', 'Apoe-hi B', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5'] == '7', 'B proliferating', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5'] == '10', 'FoB', adata.obs['cell_type_V2'])

with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(
        adata, 
        color="cell_type_V2", 
        legend_loc='on data'
        )
sc.pl.dotplot(
  adata,
  B_marker_genes,
  groupby='cell_type_V2',
  dendrogram=False,
  cmap='Reds',
  swap_axes=False,
  use_raw=True,
  standard_scale='var'
 )

adata.obs.to_csv(os.path.join(out_dir, 'annotated_B_metadata.tsv'), sep='\t')
adata.write(os.path.join(out_dir, 'B_cell_atlas.h5ad'))
