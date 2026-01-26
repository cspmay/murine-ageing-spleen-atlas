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

out_dir = './atlas/out/'

adata = sc.read_h5ad(os.path.join(out_dir, 'integrated_datasets.h5ad'))
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
        color="leiden_scVI_0.7", 
        legend_loc='on data'
        )
sc.pl.dotplot(
  adata,
  marker_genes,
  groupby='leiden_scVI_0.7',
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
  subcluster_and_umap(adata, 'leiden_scVI_0.7', ['13'], 0.2, pt_size=6)

with rc_context({'figure.figsize': (4, 4)}):
  subcluster_and_umap(adata, 'leiden_scVI_0.7', ['15'], 0.1, pt_size=6)

with rc_context({'figure.figsize': (4, 4)}):
  subcluster_and_umap(adata, 'leiden_scVI_0.7', ['23'], 0.2, pt_size=6)

with rc_context({'figure.figsize': (4, 4)}):
  subcluster_and_umap(adata, 'leiden_scVI_0.7', ['27'], 0.1, pt_size=6)

adata.obs['cell_type_V2'] = 'Other'
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7'] == '0', 'B cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7'] == '1', 'B cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7'] == '2', 'B cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7'] == '3', 'T cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7'] == '4', 'T cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7'] == '5', 'T cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7'] == '6', 'B cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7'] == '7', 'B cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7'] == '8', 'B cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7'] == '9', 'T cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7'] == '10', 'B cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7'] == '11', 'Doublets', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7'] == '12', 'B cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7_13'] == '13,0', 'Classical Monocytes', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7_13'] == '13,1', 'Non-Classical Monocytes', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7_13'] == '13,2', 'Non-Classical Monocytes', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7_13'] == '13,3', 'Resident Monocytes', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7_13'] == '13,4', 'Classical Monocytes', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7'] == '14', 'NK cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7_15'] == '15,0', 'cDC2', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7_15'] == '15,1', 'cDC1', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7_15'] == '15,2', 'migDCs', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7_15'] == '15,3', 'cDC1', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7'] == '16', 'T cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7'] == '17', 'Immature Neutrophils', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7'] == '18', 'Proliferating cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7'] == '19', 'Mature Neutrophils', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7'] == '20', 'B proliferating', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7'] == '21', 'RBCs', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7'] == '22', 'Plasma cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7_23'] == '23,0', 'Macrophages', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7_23'] == '23,1', 'Macrophages', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7_23'] == '23,2', 'Exposed Macrophages', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7_23'] == '23,3', 'Exposed Macrophages', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7'] == '24', 'T cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7'] == '25', 'B cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7'] == '26', 'pDCs', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7_27'] == '27,0', 'Basophils', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7_27'] == '27,1', 'Endothelial cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7_27'] == '27,2', 'Endothelial cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7_27'] == '27,3', 'Endothelial cells', adata.obs['cell_type_V2'])
adata.obs['cell_type_V2'] = np.where(adata.obs['leiden_scVI_0.7'] == '28', 'RBCs', adata.obs['cell_type_V2'])

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

adata.obs.to_csv(os.path.join(out_dir, 'annotated_metadata.tsv'), sep='\t')
adata.write(os.path.join(out_dir, 'Splenocyte_atlas.h5ad'))
