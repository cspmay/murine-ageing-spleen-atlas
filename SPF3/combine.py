import os
import scanpy as sc
import anndata as ad

data_dir = './data/SPF3/'
print(data_dir)

out_dir = './out'
os.makedirs(out_dir, exist_ok=True)
print(out_dir)

adata_old = sc.read_10x_mtx(os.path.join(data_dir, 'old'))
adata_old.obs['sample_name'] = 'MO'
adata_old.obs['condition'] = 'control'
adata_old

adata_young = sc.read_10x_mtx(os.path.join(data_dir, 'young'))
adata_young.obs['sample_name'] = 'MY'
adata_young.obs['condition'] = 'control'
adata_young

combined_adata = ad.concat(
	{'old':adata_old, 'young':adata_young},
  label='age_group',
  index_unique='_'
)

combined_adata.write(os.path.join(out_dir, 'combined.h5ad'))
