import os
import sys
import re
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
sys.path.append('../utils_py')
import misc_utils

out_dir = './out'

DE_out_dir = os.path.join(out_dir, 'scanpy_DE')
os.makedirs(DE_out_dir, exist_ok=True)
DE_out_dir

adata = sc.read_h5ad(os.path.join(out_dir, 'prepared_all_datasets.h5ad'))
adata

T_adata = sc.read_h5ad(os.path.join(out_dir, 'combined_T_datasets.h5ad'))
T_adata

T_adata = misc_utils.remove_barcodes(
  T_adata,
  os.path.join(out_dir, 'T_contaminating_barcodes.txt')
)
T_adata

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.normalize_total(T_adata, target_sum=1e4)
sc.pp.log1p(T_adata)

all_metadata = pd.read_csv(
  os.path.join(out_dir, 'combined_metadata_cleaned.tsv'),
  sep='\t',
  low_memory=False,
  index_col=0
)
all_metadata.head(3)
T_metadata = pd.read_csv(
  os.path.join(out_dir, 'annotated_T_metadata.tsv'),
  sep='\t',
  low_memory=False,
  index_col=0
)
T_metadata.head(3)

adata.obs['cell_type_V2'] = all_metadata['cell_type_V2']
adata.obs.head(3)
T_adata.obs['cell_type_V2'] = T_metadata['cell_type_V2']
T_adata.obs.head(3)

def dataset_DE(adata, dataset, out_dir, min_counts=10):
	dataset_out_dir = os.path.join(out_dir, dataset)
	if not os.path.exists(dataset_out_dir):
		os.mkdir(dataset_out_dir)
	print(dataset_out_dir)

	dataset_adata = adata[adata.obs['dataset'] == dataset]

	condition = 'Iso' if dataset == 'SPF4' else 'control'
	dataset_adata = dataset_adata[dataset_adata.obs['condition'] == condition]

	print(f"Total cells: {len(dataset_adata)}")
	group_counts = dataset_adata.obs.groupby(['age_group', 'cell_type_V2']).size().reset_index(name='counts')

	cell_type_counts = pd.pivot(group_counts, index='cell_type_V2', columns='age_group', values='counts').reset_index()
	cell_type_counts['do_DE'] = cell_type_counts.apply(
    lambda row: row['old'] > min_counts and row['young'] > min_counts,
    axis=1
	)

	cell_types_for_DE = cell_type_counts[cell_type_counts['do_DE']]['cell_type_V2'].tolist()

	for cell_type in cell_types_for_DE:
		if cell_type in ['Other']:
			continue
		cell_type_adata = dataset_adata[dataset_adata.obs['cell_type_V2'] == cell_type].copy()
		print(f"+++ {cell_type}: {len(cell_type_adata)}")
		sc.tl.rank_genes_groups(
	    cell_type_adata,
	    groupby='age_group',
	    method='wilcoxon',
	    reference='young',
	    pts=True
	  )
		de_genes = sc.get.rank_genes_groups_df(cell_type_adata, group='old')
		de_genes = de_genes[de_genes['pvals_adj'] < 0.05]
		de_genes = de_genes[abs(de_genes['logfoldchanges']) > 0.25]
		de_genes = de_genes[de_genes['pct_nz_group'] > 0.1]
		de_genes.rename(
			columns={'names':'name', 'pvals':'pval', 'pvals_adj':'adj_pval', 'logfoldchanges':'lfc'},
			inplace=True
		)
		tsv_out_path = os.path.join(
			dataset_out_dir,
			f"{cell_type}_sig.tsv"
		)
		de_genes.to_csv(tsv_out_path, sep='\t', index=False)

dataset_DE(adata, 'SPF1', DE_out_dir)
dataset_DE(adata, 'SPF2', DE_out_dir)
dataset_DE(adata, 'SPF3', DE_out_dir)
dataset_DE(adata, 'SPF4', DE_out_dir)
dataset_DE(adata, 'SPF5', DE_out_dir)
dataset_DE(T_adata, 'SPF6', DE_out_dir)
dataset_DE(adata, 'NM mice', DE_out_dir)
