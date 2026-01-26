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

out_dir = './NM_mice/out/'
plot_dir = './out/figures'
os.makedirs(plot_dir, exist_ok=True)
plot_dir

adata = sc.read_h5ad(os.path.join(out_dir, 'all_NM_mice.h5ad'))
# adata = sc.read_h5ad(os.path.join(out_dir, 'T_NM_mice.h5ad'))
# adata = sc.read_h5ad(os.path.join(out_dir, 'B_NM_mice.h5ad'))
# adata = sc.read_h5ad(os.path.join(out_dir, 'Gzmk_NM_mice.h5ad'))
adata

all_celltype_palette = {
	'B cells': '#3182bd',
    'B proliferating': '#17becf', 
    'Doublets': '#c6dbef',
    'T cells': '#fc9272',
	'NK cells': '#a4ac73', 
	'Macrophages': '#2ca02c',
	'Classical Monocytes': '#f16913',	
    'Non-Classical Monocytes': '#fcbba1',	
    'Proliferating Monocytes': '#8c510a',	
    'Immature Neutrophils': '#c5b0d5',
    'Mature Neutrophils': '#9467bd', 
    'cDC1': '#a50f15',	
    'cDC2': '#ce6dbd',
    'migDCs': '#c49c94',
    'pDCs': '#7b4173', 
    'Plasma cells': '#54278f',
    'Basophils': '#636363', 
    'Endothelial cells': '#98df8a',
    'RBCs proliferating': '#de2d26',
    'T proliferating': '#c7c7c7'
}
T_celltype_palette = {
	'CD8N': '#3182bd',
    'CD4N': '#fc9272',
	'CD8EM': '#fdd0a2', 
    'CD4EM': '#98df8a',
    'CD8EM Gzmk+': '#ce6dbd',	
    'nTregs': '#17becf',
    'eTregs': '#c6dbef', 
    'NLT-like Tregs': '#54278f',
    'Tfh': '#9467bd', 
    'Th2': '#636363', 
    'γδ T': '#de2d26',
    'T other': '#c7c7c7',
    'T proliferating': '#c49c94'
}
B_celltype_palette = {
	'FoB': '#3182bd',
    'Memory B': '#fc9272',
    'Apoe-hi B': '#fdd0a2',
    'MZB': '#ce6dbd',	
    'TrB': '#98df8a', 
    'GC B': '#17becf',
    'ABCs': '#de2d26',
    'Bregs': '#a4ac73',
    'B other': '#c7c7c7',
    'B proliferating': '#c49c94'
}
Gzmk_celltype_palette = {
	'CD8EM Gzmk-lo': '#3182bd',
    'CD8EM Gzmk-hi': '#c6dbef', 
}

### UMAP PLOT ###
sc.pl.umap(
    adata, 
    color="cell_type_V2", 
    palette=all_celltype_palette,
    # palette=T_celltype_palette,
    # palette=B_celltype_palette,
    # palette=Gzmk_celltype_palette,
    frameon=False,
    ax=ax,
    show=False
)
plt.tight_layout()
save_path = os.path.join(plot_dir, 'all_cells.png')
# save_path = os.path.join(plot_dir, 'T_cells.png')
# save_path = os.path.join(plot_dir, 'B_cells.png')
# save_path = os.path.join(plot_dir, 'Gzmk_cells.png')
plt.savefig(save_path, dpi=600, bbox_inches='tight')


### MATRIX PLOT ###
all_marker_genes = [
	'Ptprc',
    'Cd79a', 'Cd19', 'Ighd', #B cells
    'Cd200r3', 'Mcpt8', #Basophils    
    'Csf1r', 'Ms4a6c', #Mono
    'S100a4', 'Ly6c2', #CM
    'Cd300e', 'Cd36', 'Fcgr4', #NCM
    'Gata2', 'Sox4', 'Igfbp7', 'Plvap', #Endothelial
    'Csf3r', 'S100a8', 'S100a9', #Mature Neutrophils
    'Camp', 'Chil3', #Immature Neutrophils
    'C1qb', 'C1qc', 'C1qa', #Macrophages
    'Nkg7', 'Klrd1', 'Ccl5', #NK cells
    'Jchain', 'Mzb1', 'Sdc1', #PCs
    'Cd3e', 'Cd4', 'Cd8a', #T cells
    'Kit', 'Naaa', 'Plbd1', #cDCs
    'Xcr1', 'Irf8', #cDC1
    'Dtx1', 'Tbx21', 'Ccr6', #cDC2
    'Ccr7', 'Nudt17', #migDC
    'Siglech', 'Klk1', 'Ctsl', #pDCs
    'Hbb-bt', 'Hbb-bs', 'Hba-a1', #RBCs 
    'Top2a', 'Mki67' #Cycling
]
T_marker_genes = [
    'Cd3d', 'Cd4', 'Cd8a', 
    'Lef1', 'Sell', 'Ccr7', #naive 
    'Prf1', 'Cxcr3', 'Gzmm', 'Gzmk', 'Ccl5', #effectors
    'Bcl6', 'Il21', 'Pdcd1', 'Cxcr5', #Tfh
    'Gata3', 'Cxcr6', 'Il4', #Th2
    'Foxp3', 'Il2ra', 'Ctla4', #nTreg
    'Itgae', 'S100a4', #NLT-like
    'Tnfrsf4', 'Tnfrsf9', 'Cd83', #eTregs
    'Trdc', 'Rorc', #GD
    'Top2a', 'Mki67', #cycling
]
B_marker_genes = [
	'Cd79a', 'Cd19',
    'Zeb2', 'Itgax', 'Tbx21', #ABCs 
    'Apoe', 'B2m', 'Tcf4', #Apoe-hi
    'Fcrl5', 'Zbtb20', 'Il10', 'Ctla4', #Bregs
    'Fcer2a', 'Cxcr4', #Fo
    'Cd83', 'Hspd1', 'Myc', #GC
    'Cr2', 'Mzb1', 'Cd9', 'Cd1d1', #MZ
    'Bcl2', 'Cd80', 'Zbtb32', #Mem
    'Cd93', 'Cd24a', 'Vpreb3', #Tr
    'Top2a', 'Mki67', #cycling
]
sc.pl.matrixplot(adata, 
                 all_marker_genes,
                 # T_marker_genes,
                 # B_marker_genes,
                 groupby='cell_type_V2',
                 cmap='viridis',
                 standard_scale='var',
                 use_raw=True,
                 show=False)
save_path = os.path.join(plot_dir, 'all_cells_geneplot.png')
# save_path = os.path.join(plot_dir, 'T_cells_geneplot.png')
# save_path = os.path.join(plot_dir, 'B_cells_geneplot.png')
plt.savefig(save_path, dpi=600, bbox_inches='tight')


### PROPORTION PLOT ###
def calculate_proportions(adata, group_col, cell_type_col):
    df = adata.obs[[group_col, cell_type_col]].copy()
    counts = df.groupby([group_col, cell_type_col]).size().reset_index(name='count')
    totals = df.groupby(group_col).size().reset_index(name='total')

    merged = counts.merge(totals, on=group_col)
    merged['proportion'] = merged['count'] / merged['total']
    return merged

adata.obs['combined_group'] = adata.obs['age_group'].astype(str).str.capitalize() + ' ' + adata.obs['condition'].astype(str)

fig, ax = plt.subplots(figsize=(3, 5))

plot_data = adata[adata.obs['combined_group'].isin(['Young control', 'Old control', 'Old stress'])].copy()
plot_data.obs['combined_group'] = plot_data.obs['combined_group'].astype('category')
plot_data.obs['combined_group'] = plot_data.obs['combined_group'].cat.reorder_categories(
    ['Young control', 'Old control', 'Old stress']
)
exclude_cell_types = ['Doublets', 'RBCs proliferating']   
# exclude_cell_types = ['T cells', 'B cells', 'T proliferating', 'B proliferating', 'Doublets', 'RBCs proliferating']  
# exclude_cell_types = ['T proliferating', 'T other']  
# exclude_cell_types = ['B proliferating', 'B other']  
plot_data_filtered = plot_data[~plot_data.obs['cell_type_V2'].isin(exclude_cell_types)].copy()
prop_data_all = calculate_proportions(plot_data_filtered, 'combined_group', 'cell_type_V2')
pivot_all = prop_data_all.pivot(index='cell_type_V2', columns='combined_group', values='proportion').fillna(0)
pivot_all.T.plot(kind='bar', stacked=True, ax=ax, 
                 color=[all_celltype_palette[ct] for ct in pivot_all.index], 
                 # color=[T_celltype_palette[ct] for ct in pivot_all.index], 
                 # color=[B_celltype_palette[ct] for ct in pivot_all.index], 
                 # color=[Gzmk_celltype_palette[ct] for ct in pivot_all.index], 
                 width=0.7, legend=False)
ax.set_xlabel(' ', fontsize=12)
ax.set_ylabel('Proportion', fontsize=12)
ax.set_xticklabels(['Young', 'Old', 'Old\nstressed'], rotation=0, ha='center')
ax.set_frame_on(False)
ax.grid(False)
plt.tight_layout()
save_path = os.path.join(plot_dir, 'all_cell_type_proportions.png')
# save_path = os.path.join(plot_dir, 'all_cell_type_proportions_woTB.png')
# save_path = os.path.join(plot_dir, 'T_cell_type_proportions.png')
# save_path = os.path.join(plot_dir, 'B_cell_type_proportions.png')
# save_path = os.path.join(plot_dir, 'Gzmk_cell_type_proportions.png')
plt.savefig(save_path, dpi=600, bbox_inches='tight')
plt.show()


### UMAP DENSITY PLOTS ###
label_mapping = {
    'Young control': 'Young',
    'Old control': 'Old',
    'Old stress': 'Old stressed'
}
plot_data.obs['combined_group'] = plot_data.obs['combined_group'].replace(label_mapping)
sc.tl.embedding_density(plot_data, basis='umap', groupby='combined_group')
density_plot = sc.pl.embedding_density(plot_data, basis='umap', groupby='combined_group', 
                                       ncols=3, 
                                       show=False, 
                                       return_fig=True)
save_path = os.path.join(plot_dir, 'all_cells_umap_density.png')
# save_path = os.path.join(plot_dir, 'T_cells_umap_density.png')
# save_path = os.path.join(plot_dir, 'B_cells_umap_density.png')
# save_path = os.path.join(plot_dir, 'Gzmk_cells_umap_density.png')
plt.savefig(save_path, dpi=600, bbox_inches='tight')
plt.show()


### BOX PLOTS ###
condition_1_filter = (adata.obs['age_group'] == 'young') & (adata.obs['condition'] == 'control')
condition_2_filter = (adata.obs['age_group'] == 'old') & (adata.obs['condition'] == 'control')
condition_1_label = 'Young'
condition_2_label = 'Old'
output_subdir = 'age'

# condition_1_filter = (adata.obs['age_group'] == 'old') & (adata.obs['condition'] == 'control')
# condition_2_filter = (adata.obs['age_group'] == 'old') & (adata.obs['condition'] == 'stress')
# condition_1_label = 'Old'
# condition_2_label = 'Old stressed'
# output_subdir = 'stress'

all_cell_types = adata.obs['cell_type_V2'].unique()
def calculate_proportions_per_sample(adata, cell_type_column='cell_type_V2', sample_column='sample_name'):
    samples = adata.obs[sample_column].unique()
    cell_types = adata.obs[cell_type_column].unique()
    results = pd.DataFrame(index=samples, columns=cell_types)
    for sample in samples:
        sample_mask = adata.obs[sample_column] == sample
        total_cells = np.sum(sample_mask)
        if total_cells > 0:
            for cell_type in cell_types:
                cell_type_mask = adata.obs[cell_type_column] == cell_type
                cell_count = np.sum(sample_mask & cell_type_mask)
                results.loc[sample, cell_type] = (cell_count / total_cells) * 100
    return results.fillna(0)

def plot_single_cell_type(cond1_data, cond2_data, cell_type, plot_dir, fdr_value, 
                          cond1_label, cond2_label):
    cond1_values = cond1_data[cell_type].values
    cond2_values = cond2_data[cell_type].values
    fig, ax = plt.subplots(figsize=(3, 5))
    plot_data = pd.DataFrame({
        'Group': [cond1_label] * len(cond1_values) + [cond2_label] * len(cond2_values),
        'Percentage': np.concatenate([cond1_values, cond2_values])
    })
    color_palette = ['#8ECAE6', '#FF9B85']  # age
    # color_palette = ['#FF9B85', '#90D5A5']  # stress 
    sns.boxplot(x='Group', y='Percentage', data=plot_data, width=0.5, 
                palette=color_palette, showfliers=False, ax=ax)
    sns.stripplot(x='Group', y='Percentage', data=plot_data, jitter=0.2, 
                  color='black', size=5, alpha=0.7, ax=ax)

    ax.set_title(f'{cell_type}', fontsize=12, fontweight='bold')
    ax.set_ylabel('% of cells', fontsize=11)
    ax.set_xlabel('')
    y_max = max(np.max(cond1_values), np.max(cond2_values))
    ax.set_ylim(0, y_max * 1.25)
    if not np.isnan(fdr_value):
        y_pos = y_max * 1.1
        ax.plot([0, 1], [y_pos, y_pos], 'k-', linewidth=1.5)
        ax.plot([0, 0], [y_pos, y_pos * 0.98], 'k-', linewidth=1.5)
        ax.plot([1, 1], [y_pos, y_pos * 0.98], 'k-', linewidth=1.5)
        if fdr_value < 0.0001:
            p_text = f'FDR = {fdr_value:.2e}'
        elif fdr_value < 0.001:
            p_text = f'FDR = {fdr_value:.4f}'
        else:
            p_text = f'FDR = {fdr_value:.3f}'
        ax.text(0.5, y_pos, p_text, ha='center', va='bottom', 
                fontsize=10, fontweight='bold')
    sns.despine(ax=ax)        
    plt.tight_layout()
    file_name = cell_type.replace(' ', '_').replace('/', '_').lower()
    save_path = os.path.join(plot_dir, f'{file_name}.png')
    plt.savefig(save_path, dpi=600, bbox_inches='tight')
    plt.close()
individual_plots_dir = os.path.join(plot_dir, output_subdir)
os.makedirs(individual_plots_dir, exist_ok=True)

condition_1_samples = adata.obs.loc[condition_1_filter, 'sample_name'].unique()
condition_1_adata = adata[adata.obs['sample_name'].isin(condition_1_samples)]
condition_1_proportions = calculate_proportions_per_sample(condition_1_adata)

condition_2_samples = adata.obs.loc[condition_2_filter, 'sample_name'].unique()
condition_2_adata = adata[adata.obs['sample_name'].isin(condition_2_samples)]
condition_2_proportions = calculate_proportions_per_sample(condition_2_adata)

all_stats = []
for cell_type in all_cell_types:
    cond1_values = condition_1_proportions[cell_type].values
    cond2_values = condition_2_proportions[cell_type].values
    cond1_mean = np.mean(cond1_values)
    cond1_sem = stats.sem(cond1_values) if len(cond1_values) > 1 else 0
    cond2_mean = np.mean(cond2_values)
    cond2_sem = stats.sem(cond2_values) if len(cond2_values) > 1 else 0
    if len(cond1_values) > 1 and len(cond2_values) > 1:
        t_stat, p_val = stats.ttest_ind(cond1_values, cond2_values, equal_var=False)
    else:
        p_val = np.nan
    all_stats.append({
        'cell_type': cell_type,
        f'{condition_1_label}_mean': cond1_mean,
        f'{condition_1_label}_sem': cond1_sem,
        f'{condition_2_label}_mean': cond2_mean,
        f'{condition_2_label}_sem': cond2_sem,
        'p_value': p_val
    })
stats_df = pd.DataFrame(all_stats)
stats_df.set_index('cell_type', inplace=True)
valid_p = ~np.isnan(stats_df['p_value'])
fdr_results = np.full(len(stats_df), np.nan)
if valid_p.sum() > 0:
    _, fdr_corrected, _, _ = multipletests(stats_df.loc[valid_p, 'p_value'], 
                                            method='fdr_bh')
    fdr_results[valid_p] = fdr_corrected
stats_df['fdr'] = fdr_results
stats_df['significant'] = stats_df['fdr'] < 0.05
stats_df.sort_values('fdr', inplace=True)
for cell_type in stats_df.index:
    fdr_value = stats_df.loc[cell_type, 'fdr']
    plot_single_cell_type(condition_1_proportions, condition_2_proportions, cell_type, 
                          individual_plots_dir, fdr_value, condition_1_label, condition_2_label)
print(f"Created {len(all_cell_types)} individual cell type plots")
