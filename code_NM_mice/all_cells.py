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

out_dir = './NM_mice/out/'

adata = sc.read_h5ad(os.path.join(out_dir, 'integrated_all_NM_mice.h5ad'))
adata

marker_genes = [
	'Ptprc',
    'Cd79a', 'Cd19', 'Ighd', #B cells
    'Cd200r3', 'Mcpt8', #Basophils    
    'Csf1r', 'Ms4a6c', #Mono
    'S100a4', 'Ly6c2', #CM
    'Cd300e', 'Cd36', 'Fcgr4', #NCM
    'Gata2', 'Sox4', 'Igfbp7', 'Plvap', #Endothelial
    'Csf3r', 'S100a8', 'S100a9', #Mature Neutrophils
    'Camp', 'Chil3', #Immature Neutrophils
    'C1qb', 'C1qc', 'C1qa', #Macrophages
    'Mertk', 'Mrc1', 'Marco', 'Cxcl9', 'Cxcl10', #Exposed Macrophages
    'Nkg7', 'Klrd1', 'Ccl5', #NK cells
    'Jchain', 'Mzb1', 'Sdc1', #PCs
    'Cd3e', 'Cd4', 'Cd8a', #T cells
    'Kit', 'Naaa', 'Plbd1', #cDCs
    'Xcr1', 'Irf8', #cDC1
    'Dtx1', 'Tbx21', 'Ccr6', #cDC2
    'Ccr7', 'Nudt17', #migDC
    'Siglech', 'Klk1', 'Ctsl', #pDCs
    'Hbb-bt', 'Hbb-bs', 'Hba-a1', #RBCs 
    'Top2a', 'Mki67' #Cycling
]

with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(
        adata, 
        color="leiden_scVI_0.5", 
        legend_loc='on data'
        )
sc.pl.dotplot(
  adata,
  marker_genes,
  groupby='leiden_scVI_0.5',
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
  subcluster_and_umap(adata, 'leiden_scVI_0.5', ['6'], 0.1, pt_size=6)
with rc_context({'figure.figsize': (4, 4)}):
  subcluster_and_umap(adata, 'leiden_scVI_0.5', ['8'], 0.1, pt_size=6)
with rc_context({'figure.figsize': (4, 4)}):
  subcluster_and_umap(adata, 'leiden_scVI_0.5', ['12'], 0.1, pt_size=6)
with rc_context({'figure.figsize': (4, 4)}):
  subcluster_and_umap(adata, 'leiden_scVI_0.5', ['17'], 0.1, pt_size=6)

adata.obs['cell_type_V2'] = 'Other'
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5'] == '0', 'B cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5'] == '1', 'B cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5'] == '2', 'T cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5'] == '3', 'Doublets', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5'] == '4', 'T cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5'] == '5', 'B cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5_6'] == '6,0', 'Classical Monocytes', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5_6'] == '6,1', 'Non-Classical Monocytes', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5_6'] == '6,2', 'Classical Monocytes', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5'] == '7', 'B cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5_8'] == '8,0', 'Mature Neutrophils', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5_8'] == '8,1', 'Immature Neutrophils', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5_8'] == '8,2', 'Immature Neutrophils', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5_8'] == '8,3', 'Immature Neutrophils', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5'] == '9', 'RBCs proliferating', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5'] == '10', 'T proliferating', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5'] == '11', 'B proliferating', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5_12'] == '12,0', 'cDC2', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5_12'] == '12,1', 'pDCs', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5_12'] == '12,2', 'cDC1', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5_12'] == '12,3', 'migDCs', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5'] == '13', 'B cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5'] == '14', 'NK cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5'] == '15', 'Macrophages', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5'] == '16', 'Plasma cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5_17'] == '17,0', 'Proliferating Monocytes', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5_17'] == '17,1', 'Endothelial cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.5_17'] == '17,2', 'Basophils', adata.obs['cell_type_V2'])

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

adata.obs.to_csv(os.path.join(out_dir, 'annotated_metadata_NM_mice.tsv'), sep='\t')
adata.write(os.path.join(out_dir, 'all_NM_mice.h5ad'))
