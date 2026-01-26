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
from scipy import stats
from statsmodels.stats.multitest import multipletests

out_dir = './atlas/out/'
adata = sc.read_h5ad(os.path.join(out_dir, 'Splenocyte_atlas.h5ad'))
# adata = sc.read_h5ad(os.path.join(out_dir, 'B_cell_atlas.h5ad'))
# adata = sc.read_h5ad(os.path.join(out_dir, 'T_cell_atlas.h5ad'))
adata

meta_plot_dir = './out/figures/meta_analysis'
os.makedirs(meta_plot_dir, exist_ok=True)
meta_plot_dir

control_conditions = ['control', 'Iso']  
adata_filtered = adata[adata.obs['condition'].isin(control_conditions)].copy()


def calculate_proportions_per_study(adata, study_col='dataset', age_col='age_group', 
                                     celltype_col='cell_type_V2', sample_col='sample_name',
                                     condition_col='condition'):
    results = []
    for study in sorted(adata.obs[study_col].unique()):
        study_data = adata[adata.obs[study_col] == study]
        conditions_in_dataset = sorted(study_data.obs[condition_col].unique())
        young_samples = study_data[study_data.obs[age_col] == 'young'].obs[sample_col].unique()
        old_samples = study_data[study_data.obs[age_col] == 'old'].obs[sample_col].unique()
        
        has_replicates = (len(young_samples) > 1) and (len(old_samples) > 1)
        for celltype in study_data.obs[celltype_col].unique():
            young_props = []
            old_props = []
            for sample in young_samples:
                mask = (study_data.obs[sample_col] == sample)
                total = mask.sum()
                celltype_count = (mask & (study_data.obs[celltype_col] == celltype)).sum()
                young_props.append(celltype_count / total if total > 0 else 0)
            for sample in old_samples:
                mask = (study_data.obs[sample_col] == sample)
                total = mask.sum()
                celltype_count = (mask & (study_data.obs[celltype_col] == celltype)).sum()
                old_props.append(celltype_count / total if total > 0 else 0)
            mean_young = np.mean(young_props)
            mean_old = np.mean(old_props)
            logfc = np.log2((mean_old + 1e-6) / (mean_young + 1e-6))
            if has_replicates:
                t_stat, p_val = stats.ttest_ind(old_props, young_props)
                se = np.sqrt(np.var(young_props)/len(young_props) + 
                            np.var(old_props)/len(old_props))
                if mean_young > 1e-10:
                    se_logfc = se / (mean_young * np.log(2))
                else:
                    se_logfc = np.nan
                    p_val = np.nan
            else:
                se_logfc = np.nan
                p_val = np.nan

            results.append({
                'dataset': study,
                'cell_type': celltype,
                'logFC': logfc,
                'se': se_logfc,
                'n_young': len(young_samples),
                'n_old': len(old_samples),
                'has_replicates': has_replicates,
                'p_value': p_val,
                'mean_young_prop': mean_young,
                'mean_old_prop': mean_old,
                'conditions': ', '.join(conditions_in_dataset)
            })
    return pd.DataFrame(results)


def meta_analysis_celltype(df_celltype):
    with_reps = df_celltype[df_celltype['has_replicates'] == True].copy()
    without_reps = df_celltype[df_celltype['has_replicates'] == False].copy()
    if len(with_reps) == 0:
        return {
            'cell_type': df_celltype['cell_type'].iloc[0],
            'meta_logFC': np.nan,
            'meta_se': np.nan,
            'meta_pval': np.nan,
            'I2': np.nan,
            'tau2': np.nan,
            'Q': np.nan,
            'Q_pval': np.nan,
            'n_studies_meta': 0,
            'n_studies_descriptive': len(without_reps)
        }
    with_reps = with_reps[with_reps['se'].notna()].copy()
    if len(with_reps) == 0:
        return {
            'cell_type': df_celltype['cell_type'].iloc[0],
            'meta_logFC': np.nan,
            'meta_se': np.nan,
            'meta_pval': np.nan,
            'I2': np.nan,
            'tau2': np.nan,
            'Q': np.nan,
            'Q_pval': np.nan,
            'n_studies_meta': 0,
            'n_studies_descriptive': len(without_reps)
        }
    weights = 1 / (with_reps['se'] ** 2)
    fe_estimate = np.sum(weights * with_reps['logFC']) / np.sum(weights)
    Q = np.sum(weights * (with_reps['logFC'] - fe_estimate) ** 2)
    df_Q = len(with_reps) - 1
    Q_pval = 1 - stats.chi2.cdf(Q, df_Q) if df_Q > 0 else np.nan
    C = np.sum(weights) - np.sum(weights**2) / np.sum(weights)
    tau2 = max(0, (Q - df_Q) / C) if df_Q > 0 and C > 0 else 0
    re_weights = 1 / (with_reps['se'] ** 2 + tau2)
    re_estimate = np.sum(re_weights * with_reps['logFC']) / np.sum(re_weights)
    re_se = np.sqrt(1 / np.sum(re_weights))
    z = re_estimate / re_se
    p_val = 2 * (1 - stats.norm.cdf(abs(z)))
    I2 = max(0, 100 * (Q - df_Q) / Q) if Q > 0 else 0
    return {
        'cell_type': df_celltype['cell_type'].iloc[0],
        'meta_logFC': re_estimate,
        'meta_se': re_se,
        'meta_pval': p_val,
        'meta_ci_lower': re_estimate - 1.96 * re_se,
        'meta_ci_upper': re_estimate + 1.96 * re_se,
        'I2': I2,
        'tau2': tau2,
        'Q': Q,
        'Q_pval': Q_pval,
        'n_studies_meta': len(with_reps),
        'n_studies_descriptive': len(without_reps),
        'studies_with_reps': ', '.join(sorted(with_reps['dataset'].tolist())),
        'studies_without_reps': ', '.join(sorted(without_reps['dataset'].tolist())) if len(without_reps) > 0 else 'None'
    }
prop_df = calculate_proportions_per_study(adata_filtered, 
                                          study_col='dataset',
                                          age_col='age_group', 
                                          celltype_col='cell_type_V2',
                                          sample_col='sample_name',
                                          condition_col='condition')
rep_summary = prop_df.groupby('dataset', observed=True).agg({
    'has_replicates': 'first',
    'n_young': 'first',
    'n_old': 'first'
}).reset_index()
with_reps = rep_summary[rep_summary['has_replicates'] == True]
without_reps = rep_summary[rep_summary['has_replicates'] == False]
meta_results = []
for celltype in prop_df['cell_type'].unique():
    ct_data = prop_df[prop_df['cell_type'] == celltype]
    meta_results.append(meta_analysis_celltype(ct_data))
meta_df = pd.DataFrame(meta_results)
valid_pvals = meta_df['meta_pval'].notna()
if valid_pvals.sum() > 0:
    meta_df.loc[valid_pvals, 'meta_fdr'] = stats.false_discovery_control(
        meta_df.loc[valid_pvals, 'meta_pval']
    )
else:
    meta_df['meta_fdr'] = np.nan


def create_forest_plot(prop_df, meta_df, celltype, save_path=None):
    ct_data = prop_df[prop_df['cell_type'] == celltype].copy()
    ct_data = ct_data.sort_values('dataset').reset_index(drop=True)
    meta_row = meta_df[meta_df['cell_type'] == celltype].iloc[0]
    fig, ax = plt.subplots(figsize=(4, max(2, len(ct_data) * 0.3 + 1.5)))
    y_pos = range(len(ct_data))
    for i, (_, row) in enumerate(ct_data.iterrows()):
        color = 'black' if row['has_replicates'] else 'gray'
        marker = 'D' if row['has_replicates'] else 'o'
        alpha = 1.0 if row['has_replicates'] else 0.5
        ax.plot(row['logFC'], i, marker=marker, color=color, markersize=6, alpha=alpha)
        if row['has_replicates'] and not np.isnan(row['se']):
            ci_lower = row['logFC'] - 1.96 * row['se']
            ci_upper = row['logFC'] + 1.96 * row['se']
            ax.plot([ci_lower, ci_upper], [i, i], color=color, linewidth=2, alpha=alpha)
    study_labels = ct_data['dataset'].tolist()
    if not np.isnan(meta_row['meta_logFC']):
        ax.axhline(len(ct_data) - 0.5, color='black', linestyle='--', linewidth=1)

        meta_y = len(ct_data)
        ax.plot(meta_row['meta_logFC'], meta_y, marker='D', color='red', 
                markersize=8, label='Meta-analysis')
        ax.plot([meta_row['meta_ci_lower'], meta_row['meta_ci_upper']], 
                [meta_y, meta_y], color='red', linewidth=2)
        study_labels.append('META-ANALYSIS')
    ax.axvline(0, color='black', linestyle='-', linewidth=0.5)
    ax.set_yticks(range(len(study_labels)))
    ax.set_yticklabels(study_labels)
    ax.set_xlabel('Log2 Fold Change (Old vs Young)', fontsize=10)
    fdr = meta_row["meta_fdr"]
    if pd.notna(fdr):
        if fdr < 0.0001:
            fdr_str = f'{fdr:.2e}'
        else:
            fdr_str = f'{fdr:.4f}'
    else:
        fdr_str = 'N/A'
    ax.set_title(f'{celltype}\nI² = {meta_row["I2"]:.1f}%, FDR = {fdr_str}', 
                 fontsize=10)
    ax.grid(axis='x', alpha=0.3)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=600, bbox_inches='tight')
    plt.show()
    plt.close()

sig_results = meta_df[meta_df['meta_fdr'] < 0.05].sort_values('meta_pval')
high_het = meta_df[meta_df['I2'] > 75].sort_values('I2', ascending=False)
prop_df = prop_df.sort_values(['cell_type', 'dataset']).reset_index(drop=True)
meta_df.to_csv(os.path.join(meta_plot_dir, 'meta_analysis_results_control_iso.csv'), index=False)
prop_df.to_csv(os.path.join(meta_plot_dir, 'dataset_level_results_control_iso.csv'), index=False)
if len(sig_results) > 0:
    for celltype in sig_results['cell_type']:
        print(f"  Creating plot for {celltype}...")
        save_path = os.path.join(meta_plot_dir, f'forest_plot_{celltype.replace("/", "_").replace(" ", "_")}.png')
        create_forest_plot(prop_df, meta_df, celltype, save_path=save_path)
    non_sig = meta_df[(meta_df['meta_fdr'] >= 0.05) & (meta_df['n_studies_meta'] > 0)]
    if len(non_sig) > 0:
        for celltype in non_sig['cell_type']:
            print(f"  Creating plot for {celltype}...")
            save_path = os.path.join(meta_plot_dir, f'forest_plot_{celltype.replace("/", "_").replace(" ", "_")}.png')
            create_forest_plot(prop_df, meta_df, celltype, save_path=save_path)
else:
    plottable = meta_df[meta_df['n_studies_meta'] > 0]
    for celltype in plottable['cell_type']:
        print(f"  Creating plot for {celltype}...")
        save_path = os.path.join(meta_plot_dir, f'forest_plot_{celltype.replace("/", "_").replace(" ", "_")}.png')
        create_forest_plot(prop_df, meta_df, celltype, save_path=save_path)
