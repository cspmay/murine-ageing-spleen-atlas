import os
import scanpy as sc
import anndata as ad

data_dir = './data/SPF5/'
print(data_dir)

out_dir = './out'
os.makedirs(out_dir, exist_ok=True)
print(out_dir)

adata_young_saline = sc.read_10x_mtx(os.path.join(data_dir, 'young_saline', 'filtered_feature_bc_matrix'))
adata_young_saline.obs['condition'] = 'control'
adata_young_saline.obs['age_group'] = 'young'
adata_young_saline

adata_old_saline = sc.read_10x_mtx(os.path.join(data_dir, 'aged_saline', 'filtered_feature_bc_matrix'))
adata_old_saline.obs['condition'] = 'control'
adata_old_saline.obs['age_group'] = 'old'
adata_old_saline

adata_old_PF4 = sc.read_10x_mtx(os.path.join(data_dir, 'aged_PF4', 'filtered_feature_bc_matrix'))
adata_old_PF4.obs['condition'] = 'PF4'
adata_old_PF4.obs['age_group'] = 'old'
adata_old_PF4

combined_adata = ad.concat(
	{'SYS':adata_young_saline, 'SAS':adata_old_saline, 'SAPF4': adata_old_PF4},
  label='sample_name',
  index_unique='_'
)

combined_adata.write(os.path.join(out_dir, 'combined.h5ad'))
