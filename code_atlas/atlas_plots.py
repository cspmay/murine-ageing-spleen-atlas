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
from adjustText import adjust_text

out_dir = './atlas/out/'
print(out_dir)

plot_dir = './atlas/out/figures'
os.makedirs(plot_dir, exist_ok=True)
plot_dir

adata = sc.read_h5ad(os.path.join(out_dir, 'Splenocyte_atlas.h5ad'))
# adata = sc.read_h5ad(os.path.join(out_dir, 'B_cell_atlas.h5ad'))
# adata = sc.read_h5ad(os.path.join(out_dir, 'T_cell_atlas.h5ad'))
adata

### DATASET UMAPs AFTER INTEGRATION ###
datasets_palette = {
	'SPF1': '#00A087',
    'SPF2': '#E64B35',
	'SPF3': '#4DBBD5',
	'SPF4': '#8491B4',
    'SPF5': '#3C5488',
	'NM mice': '#F39B7F',	
}

with rc_context({'figure.figsize': (3, 3)}):
    sc.pl.umap(adata, color='dataset', palette=datasets_palette, frameon=False, title='After integration')
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

### ATLAS UMAPs ### 
all_celltype_palette = {
	'B cells': '#3182bd',
    'B proliferating': '#17becf', 
    'Doublets': '#c6dbef',
    'T cells': '#fc9272',
	'NK cells': '#a4ac73', 
	'Macrophages': '#2ca02c',
    'Exposed Macrophages': '#98df8a',
	'Classical Monocytes': '#fd8d3c',	
    'Non-Classical Monocytes': '#fcbba1',	
    'Resident Monocytes': '#8c510a',	
    'Immature Neutrophils': '#c5b0d5',
    'Mature Neutrophils': '#9467bd', 
    'cDC1': '#a50f15',	
    'cDC2': '#ce6dbd',
    'migDCs': '#c49c94',
    'pDCs': '#7b4173', 
    'Plasma cells': '#54278f', 
    'Basophils': '#636363', 
    'Endothelial cells': '#a63603',
    'RBCs': '#de2d26',
    'Proliferating cells': '#c7c7c7',
    'Other': '#8c564b',
}
B_celltype_palette = {
	'FoB': '#3182bd',
    'Memory B': '#fc9272',
    'Apoe-hi B': '#98df8a',
    'MZB': '#ce6dbd',	
    'TrB': '#fdd0a2', 
    'GC B': '#17becf', 
    'B other': '#c7c7c7',
    'B proliferating': '#c49c94'
}
T_celltype_palette = {
	'CD8N': '#3182bd',
    'CD4N': '#fc9272',
	'CD8EM': '#fdd0a2', 
    'CD4EM': '#98df8a',
    'CD8EM Gzmk+': '#ce6dbd',	
    'nTregs': '#17becf',
    'eTregs': '#c6dbef', 
    'NLT-like Tregs': '#a63603',
    'Tfh': '#9467bd', 
    'Th2': '#636363', 
    'T other': '#c7c7c7',
    'T proliferating': '#c49c94'
}

sc.pl.umap(
    adata, 
    color="cell_type_V2", 
    palette=all_celltype_palette,
    # palette=B_celltype_palette,
    # palette=T_celltype_palette,
    frameon=False,
    title="Splenocyte atlas",
    # title="B cell atlas",
    # title="T cell atlas",
    ax=ax,
    show=False
)
save_path = os.path.join(plot_dir, 'Splenocyte_atlas.png')
# save_path = os.path.join(plot_dir, 'B_cell_atlas.png')
# save_path = os.path.join(plot_dir, 'T_cell_atlas.png')
plt.savefig(save_path, dpi=600, bbox_inches='tight')

### CELL TYPE PROPORTIONS ###
def calculate_proportions(adata, group_col, cell_type_col):
    df = adata.obs[[group_col, cell_type_col]].copy()
    counts = df.groupby([group_col, cell_type_col]).size().reset_index(name='count')
    totals = df.groupby(group_col).size().reset_index(name='total')
    merged = counts.merge(totals, on=group_col)
    merged['proportion'] = merged['count'] / merged['total']
    return merged
plot_data = adata[adata.obs['condition'].isin(['control', 'Iso'])].copy()
plot_data.obs['dataset_label'] = plot_data.obs['dataset'].astype(str)
plot_data.obs['age_group_clean'] = plot_data.obs['age_group'].str.capitalize()
plot_data.obs['group_label'] = (plot_data.obs['age_group_clean'] + '_' + 
                                 plot_data.obs['dataset_label'])
exclude_cell_types = ['Doublets', 'RBCs', 'Other']  
# exclude_cell_types = ['T cells', 'B cells', 'Doublets', 'RBCs', 'Other']  
# exclude_cell_types = ['B other', 'B proliferating']  
# exclude_cell_types = ['T other', 'T proliferating'] 
plot_data_filtered = plot_data[~plot_data.obs['cell_type_V2'].isin(exclude_cell_types)].copy()
prop_data_all = calculate_proportions(plot_data_filtered, 'group_label', 'cell_type_V2')
pivot_all = prop_data_all.pivot(index='cell_type_V2', columns='group_label', values='proportion').fillna(0)
young_cols = [col for col in pivot_all.columns if col.startswith('Young_')]
old_cols = [col for col in pivot_all.columns if col.startswith('Old_')]
pivot_all = pivot_all[young_cols + old_cols]
fig, ax = plt.subplots(figsize=(5, 4.75))
pivot_all.T.plot(kind='bar', stacked=True, ax=ax, 
                 color=[all_celltype_palette[ct] for ct in pivot_all.index], 
                 # color=[B_celltype_palette[ct] for ct in pivot_all.index],
                 # color=[T_celltype_palette[ct] for ct in pivot_all.index],
                 width=0.8, legend=False)
n_datasets = len(young_cols)
labels = [col.replace('Young_', '') for col in young_cols] + [col.replace('Old_', '') for col in old_cols]
ax.set_xticklabels(labels, rotation=90, ha='center', fontsize=10)

ax.set_ylabel('Proportion')
ax.set_frame_on(False)
ax.grid(False)
plt.tight_layout()
save_path = os.path.join(plot_dir, 'all_cell_type_proportions.png')
# save_path = os.path.join(plot_dir, 'B_cell_type_proportions.png')
# save_path = os.path.join(plot_dir, 'T_cell_type_proportions.png')
plt.savefig(save_path, dpi=600, bbox_inches='tight')
