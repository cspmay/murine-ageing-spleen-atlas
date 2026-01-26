import os
import scanpy as sc
import anndata as ad
import pandas as pd

def remove_genes(adata, gene_list):
  result = adata[:, ~adata.var_names.isin(list(gene_list))].copy() 
  return result

def remove_barcodes(adata, barcode_path):
  with open(barcode_path) as f:
    barcode_list = [line.rstrip() for line in f]
    print(len(barcode_list))
  result = adata[~adata.obs.index.isin(barcode_list)].copy() 
  return result

def load_filtered(out_dir, barcodes=[]):
	adata = sc.read_h5ad(os.path.join(out_dir, 'filtered.h5ad'))
	if barcodes:
		barcodes_to_include = adata.obs_names[adata.obs_names.isin(barcodes)]
		adata = adata[barcodes_to_include]
	return adata.copy()
