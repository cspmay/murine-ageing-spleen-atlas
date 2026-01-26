import os
import scanpy as sc
import anndata as ad

data_dir = './data/SPF4/'
print(data_dir)

out_dir = './out'
os.makedirs(out_dir, exist_ok=True)
print(out_dir)

adatas = {}

for file_name in sorted(os.listdir(data_dir)):
	condition = file_name.split('_')[0]
	age_group = file_name.split('_')[1]
	age_group = 'young' if age_group == 'Y' else 'old'
	print(f"Sample: {file_name}, age: {age_group}, condition: {condition}")
	adata = sc.read_10x_mtx(os.path.join(data_dir, file_name, 'filtered_feature_bc_matrix'))
	adata.var_names_make_unique()
	adata.obs['age_group'] = age_group
	adata.obs['condition'] = condition
	adatas[file_name] = adata

combined_adata = ad.concat(
	adatas,
  label='sample_name',
  index_unique='_'
)

combined_adata.write(os.path.join(out_dir, 'combined.h5ad'))
