import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.pyplot import rc_context

out_dir = './atlas/out/'
print(out_dir)

plot_dir = './atlas/out/figures/all'
os.makedirs(plot_dir, exist_ok=True)
plot_dir

adata = sc.read_h5ad(os.path.join(out_dir, 'not_integrated_datasets.h5ad'))
adata

datasets_palette = {
	'SPF1': '#00A087',
    'SPF2': '#E64B35',
	'SPF3': '#4DBBD5',
	'SPF4': '#8491B4',
    'SPF5': '#3C5488',
	'NM mice': '#F39B7F',	
}

with rc_context({'figure.figsize': (3, 3)}):
    sc.pl.umap(adata, color='dataset', palette=datasets_palette, frameon=False, title='Before integration')
with rc_context({'figure.figsize': (3, 3)}):
    sc.pl.umap(adata, color='dataset', groups=['SPF1'], palette=datasets_palette, frameon=False, title='SPF1', legend_loc=None)
with rc_context({'figure.figsize': (3, 3)}):
	sc.pl.umap(adata, color='dataset', groups=['SPF2'], palette=datasets_palette, frameon=False, title='SPF2', legend_loc=None)
with rc_context({'figure.figsize': (3, 3)}):
	sc.pl.umap(adata, color='dataset', groups=['SPF3'], palette=datasets_palette, frameon=False, title='SPF3', legend_loc=None)
with rc_context({'figure.figsize': (3, 3)}):
	sc.pl.umap(adata, color='dataset', groups=['SPF4'], palette=datasets_palette, frameon=False, title='SPF4', legend_loc=None)
with rc_context({'figure.figsize': (3, 3)}):
	sc.pl.umap(adata, color='dataset', groups=['SPF5'], palette=datasets_palette, frameon=False, title='SPF5', legend_loc=None)
with rc_context({'figure.figsize': (3, 3)}):
	sc.pl.umap(adata, color='dataset', groups=['NM mice'], palette=datasets_palette, frameon=False, title='NM mice', legend_loc=None)
