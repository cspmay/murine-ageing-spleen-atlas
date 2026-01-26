import os
import sys
import pandas as pd
import scanpy as sc
import scipy

out_dir = '/gpfs/helios/home/chinnasu/parse_scrnaseq/NM_mice/out'

all_annotated_out_dir = os.path.join(out_dir, 'all_annotated')
# all_annotated_out_dir = os.path.join(out_dir, 'T_annotated')
# all_annotated_out_dir = os.path.join(out_dir, 'B_annotated')
# all_annotated_out_dir = os.path.join(out_dir, 'Gzmk_annotated')
get_ipython().system(' mkdir -p $all_annotated_out_dir')
print(all_annotated_out_dir)

adata = sc.read_h5ad(os.path.join(out_dir, 'prepared_all_NM_mice.h5ad'))
# adata = sc.read_h5ad(os.path.join(out_dir, 'prepared_T_NM_mice.h5ad'))
# adata = sc.read_h5ad(os.path.join(out_dir, 'prepared_B_NM_mice.h5ad'))
# adata = sc.read_h5ad(os.path.join(out_dir, 'prepared_Gzmk_NM_mice.h5ad'))
adata

scipy.io.mmwrite(
  os.path.join(all_annotated_out_dir, 'matrix.mtx'),
  adata.X
)
with open(os.path.join(all_annotated_out_dir, 'features.txt'), 'w') as f:
  for feature in adata.var.index.to_list():
    f.write(f"{feature}\n")
with open(os.path.join(all_annotated_out_dir, 'barcodes.txt'), 'w') as f:
  for barcode in adata.obs.index.to_list():
    f.write(f"{barcode}\n")
adata.obs.to_csv(os.path.join(all_annotated_out_dir, 'metadata.tsv'), sep='\t')

print(f"✓ Matrix shape: {adata.X.shape} (genes x cells)") 
print(f"✓ Features: {len(adata.var.index)} genes")
print(f"✓ Barcodes: {len(adata.obs.index)} cells")
