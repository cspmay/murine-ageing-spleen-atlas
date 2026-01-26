import os, sys
import pandas as pd
import anndata as ad
import scanpy as sc

sys.path.append('../utils_py')
import genes
import misc_utils

out_dir = './NM_mice/out'

metadata_V2 = pd.read_csv(
  os.path.join(out_dir, 'annotated_metadata_NM_mice.tsv'),
  sep='\t',
  index_col=0
)
metadata_V2.head(3)

T_barcodes = metadata_V2.index[metadata_V2['cell_type_V2'].isin(['T cells'])].tolist()
len(T_barcodes)

B_barcodes = metadata_V2.index[metadata_V2['cell_type_V2'].isin(['B cells'])].tolist()
len(B_barcodes)


atlas_dir = './atlas/out'
adata = sc.read_h5ad(os.path.join(atlas_dir, 'combined_atlas.h5ad'))
adata

adata_NM_mice = adata[adata.obs['dataset'] == 'NM_mice']
adata_NM_mice

genes_to_remove = genes.IG_genes + genes.TCR_genes

adata_NM_mice = misc_utils.remove_genes(adata_NM_mice, genes_to_remove)
adata_NM_mice
adata_NM_mice.write(os.path.join(out_dir, 'prepared_NM_mice.h5ad'))


T_barcodes_to_include = adata_NM_mice.obs_names[adata_NM_mice.obs_names.isin(T_barcodes)]
T_adata_NM_mice = adata_NM_mice[T_barcodes_to_include]
T_adata_NM_mice

T_adata_NM_mice = misc_utils.remove_genes(T_adata_NM_mice, genes_to_remove)
T_adata_NM_mice
T_adata_NM_mice.write(os.path.join(out_dir, 'prepared_T_NM_mice.h5ad'))


Gzmk_barcodes_to_include = adata_NM_mice.obs_names[adata_NM_mice.obs_names.isin(Gzmk_barcodes)]
Gzmk_adata_NM_mice = adata_NM_mice[Gzmk_barcodes_to_include]
Gzmk_adata_NM_mice

Gzmk_adata_NM_mice = misc_utils.remove_genes(Gzmk_adata_NM_mice, genes_to_remove)
Gzmk_adata_NM_mice
Gzmk_adata_NM_mice.write(os.path.join(out_dir, 'prepared_Gzmk_NM_mice.h5ad'))


B_barcodes_to_include = adata_NM_mice.obs_names[adata_NM_mice.obs_names.isin(B_barcodes)]
B_adata_NM_mice = adata_NM_mice[B_barcodes_to_include]
B_adata_NM_mice

B_adata_NM_mice = misc_utils.remove_genes(B_adata_NM_mice, genes_to_remove)
B_adata_NM_mice
B_adata_NM_mice.write(os.path.join(out_dir, 'prepared_B_NM_mice.h5ad'))
