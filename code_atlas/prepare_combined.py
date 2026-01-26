import os
import sys
import scanpy as sc

sys.path.append('../utils_py')
import genes
import misc_utils

out_dir = './out'
print(out_dir)

adata = sc.read_h5ad(os.path.join(out_dir, 'combined_atlas.h5ad'))
adata

genes_to_remove = genes.IG_genes + genes.TCR_genes
adata = misc_utils.remove_genes(adata, genes_to_remove)
adata

adata.write(os.path.join(out_dir, 'prepared_all_datasets.h5ad'))


T_adata = sc.read_h5ad(os.path.join(out_dir, 'combined_T.h5ad'))
T_adata

genes_to_remove = genes.IG_genes + genes.TCR_genes
T_adata = misc_utils.remove_genes(T_adata, genes_to_remove)
T_adata

T_adata.write(os.path.join(out_dir, 'prepared_T_datasets.h5ad'))
