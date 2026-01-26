import os
import sys
import anndata as ad
import scanpy as sc

out_dir = './out'
os.makedirs(out_dir, exist_ok=True)
print(out_dir)

def load_filtered(out_dir):
  return sc.read_h5ad(os.path.join(out_dir, 'filtered.h5ad'))

adatas = {
    'SPF1': load_filtered('../SPF1/out'),
    'SPF2': load_filtered('../SPF2/out'),
    'SPF3': load_filtered('../SPF3/out'),	
    'SPF4': load_filtered('../SPF4/out'),
    'SPF5': load_filtered('../SPF5/out'),
    'NM mice': load_filtered('../NM_mice/out'),
}

combined_adata = ad.concat(adatas, label='dataset')
combined_adata

combined_adata.obs.groupby(['dataset', 'age_group', 'condition']).size().reset_index(name='counts')

combined_adata.write(os.path.join(out_dir, 'combined_atlas.h5ad'))
