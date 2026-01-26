import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

data_dir = './data/NM_mice/'
print(data_dir)

out_dir = './out'
os.makedirs(out_dir, exist_ok=True)
print(out_dir)

adata = sc.read_mtx(data_dir + 'DGE.mtx')
gene_data = pd.read_csv(data_dir + 'all_genes.csv')
cell_meta = pd.read_csv(data_dir + 'cell_metadata.csv')

gene_data = gene_data[gene_data.gene_name.notnull()]
notNa = gene_data.index
notNa = notNa.to_list()

adata = adata[:,notNa]
adata.var = gene_data
adata.var.set_index('gene_name', inplace=True)
adata.var.index.name = None
adata.var_names_make_unique()

adata.obs = cell_meta
adata.obs.set_index('bc_wells', inplace=True)
adata.obs.index.name = None
adata.obs_names_make_unique()
adata

adata = adata[~adata.obs['sample'].str.startswith('tenxP')]
adata.obs['sample_name'] = adata.obs['sample']
adata

adata.obs['age_group'] = np.where(
	adata.obs['sample_name'].str.startswith('CVS'),
	'old',
	'young'
)
adata.obs['condition'] = np.where(
	adata.obs['sample_name'].str.startswith('YD'),
	'control',
	np.where(
		adata.obs['sample_name'].str.startswith('CVS_C'),
		'control',
		'stress'
	)
)

adata.write(os.path.join(out_dir, 'combined.h5ad'))
