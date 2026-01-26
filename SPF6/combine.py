import os
import scanpy as sc
import anndata as ad

data_dir = './data/SPF6/'
print(data_dir)

out_dir = './out'
os.makedirs(out_dir, exist_ok=True)
print(out_dir)

adatas = {}

for sample_name in sorted(os.listdir(data_dir)):
	age_group = 'young' if sample_name.startswith('young') else 'old'
	print(f"Sample: {sample_name}, age: {age_group}")
	adata = sc.read_10x_mtx(os.path.join(data_dir, sample_name, 'filtered_feature_bc_matrix'))
	adata.var_names_make_unique()
	adata.obs['condition'] = 'control'
	adata.obs['age_group'] = age_group
	adatas[sample_name] = adata

combined_adata = ad.concat(
	adatas,
  label='sample_name',
  index_unique='_'
)

combined_adata.write(os.path.join(out_dir, 'combined.h5ad'))
