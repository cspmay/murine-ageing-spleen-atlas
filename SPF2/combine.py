import os
import scanpy as sc
import anndata as ad

data_dir = './data/SPF2/'
print(data_dir)

out_dir = './out'
os.makedirs(out_dir, exist_ok=True)
print(out_dir)

adatas = {}

for file_name in sorted(os.listdir(data_dir)):
	sample_name = file_name[:4]
	age_group = 'young' if sample_name.startswith('Y') else 'old'
	print(f"Sample: {sample_name}, age: {age_group}")
	adata = sc.read_10x_h5(os.path.join(data_dir, file_name))
	adata.var_names_make_unique()
	adata.obs['age_group'] = age_group
	adata.obs['condition'] = 'control'
	adatas[sample_name] = adata

combined_adata = ad.concat(
	adatas,
  label='sample_name',
  index_unique='_'
)

combined_adata.write(os.path.join(out_dir, 'combined.h5ad'))
