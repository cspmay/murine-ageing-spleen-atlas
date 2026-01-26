import os
import sys
import scvi
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context

sys.path.append('../utils_py')
import misc_utils

scvi.settings.num_threads = 8
scvi.settings.seed = 123

out_dir = './out'

model_path = os.path.join(out_dir, 'scvi_model_NM_mice_5000')
# model_path = os.path.join(out_dir, 'scvi_model_NM_mice_T_cleaned_3000')
# model_path = os.path.join(out_dir, 'scvi_model_NM_mice_Gzmk_cleaned_2000')
# model_path = os.path.join(out_dir, 'scvi_model_NM_mice_B_3000')
print(model_path)

adata = sc.read_h5ad(os.path.join(out_dir, 'prepared_NM_mice.h5ad'))
# adata = sc.read_h5ad(os.path.join(out_dir, 'prepared_T_NM_mice.h5ad'))
# adata = sc.read_h5ad(os.path.join(out_dir, 'prepared_B_NM_mice.h5ad'))
adata

# adata = misc_utils.remove_barcodes(
  # adata,
  # os.path.join(out_dir, 'contaminating_barcodes.txt')
# )

# T_annotations = pd.read_csv(os.path.join(out_dir, 'annotated_T_metadata_NM_mice.tsv'), sep='\t', index_col=0)
# T_annotations.head()
# adata.obs['cell_type_V2'] = T_annotations['cell_type_V2']
# adata = adata[adata.obs['cell_type_V2'] == 'CD8EM Gzmk+']

adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(
  adata,
  n_top_genes=5000,
  # n_top_genes=3000,
  # n_top_genes=2000,
  subset=True,
  layer="counts",
  flavor="seurat_v3"
)

scvi.model.SCVI.setup_anndata(
  adata,
  layer='counts',
  categorical_covariate_keys=['sample_name'],
  continuous_covariate_keys=['pct_counts_mt']
)

model = scvi.model.SCVI(adata)
model

n_epochs = 100

get_ipython().system(' rm -rf $model_path')
model.train(max_epochs=n_epochs, early_stopping=True)
model.save(model_path)


model = scvi.model.SCVI.load(model_path, adata)
latent = model.get_latent_representation()
adata.obsm['X_scVI'] = latent


sc.pp.neighbors(adata, use_rep='X_scVI')
sc.tl.umap(adata)

with rc_context({'figure.figsize': (6, 4)}):
  elbo_train_set = model.history["elbo_train"]["elbo_train"]
  elbo_val_set = model.history["elbo_validation"]["elbo_validation"]
  x = np.linspace(0, n_epochs, (len(elbo_train_set)))
  plt.plot(x, elbo_train_set, label="train")
  plt.plot(x, elbo_val_set, label="val")
  plt.show()

for current_resolution in np.round(np.arange(0.2, 0.8, 0.1), 1):
# for current_resolution in np.round(np.arange(0.3, 0.6, 0.1), 1):
  clustering_key = 'leiden_scVI_%0.1f' % current_resolution
  sc.tl.leiden(
    adata,
    key_added=clustering_key,
    resolution=current_resolution
  )

adata.write(os.path.join(out_dir, 'integrated_NM_mice.h5ad'))
# adata.write(os.path.join(out_dir, 'integrated_T_NM_mice.h5ad'))
# adata.write(os.path.join(out_dir, 'integrated_B_NM_mice.h5ad'))
# adata.write(os.path.join(out_dir, 'integrated_Gzmk_NM_mice.h5ad'))
