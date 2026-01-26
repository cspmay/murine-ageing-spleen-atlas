import os
import sys
import pandas as pd
import anndata as ad
import scanpy as sc

sys.path.append('../utils_py')
import misc_utils

out_dir = './out'
print(out_dir)

T_out_dir = os.path.join(out_dir, 'T_cells')
os.makedirs(T_out_dir, exist_ok=True)
print(T_out_dir)

metadata_V2 = pd.read_csv(
  os.path.join(out_dir, 'annotated_metadata.tsv'),
  sep='\t',
  index_col=0
)
metadata_V2.head(3)

T_barcodes = metadata_V2.index[metadata_V2['leiden_scVI_0.7'].isin([3, 4, 5, 9, 16, 24])].tolist()
len(T_barcodes)

adatas = {
  'Kimmel': misc_utils.load_filtered('../Kimmel/out', T_barcodes),
	'Mogilenko': misc_utils.load_filtered('../Mogilenko/out', T_barcodes),
	'parse_WT': misc_utils.load_filtered('../parse_WT/out', T_barcodes),
	'TMS': misc_utils.load_filtered('../TMS/out', T_barcodes),
	'Schroer': misc_utils.load_filtered('../Schroer/out', T_barcodes),
	'Ma': misc_utils.load_filtered('../Ma/out', T_barcodes),
	'Tyrrell': misc_utils.load_filtered('../Tyrrell/out')
}

combined_adata = ad.concat(adatas, label='dataset')
combined_adata

cells_df = combined_adata.obs.groupby(['dataset', 'age_group', 'condition']).size().reset_index(name='counts')
cells_df[cells_df['counts'] != 0]

combined_adata.write(os.path.join(T_out_dir, 'combined_T.h5ad'))
