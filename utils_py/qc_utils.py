import pandas as pd
import numpy as np
import scanpy as sc

def get_qc_df(adata, qc_cutoffs):
  qc_df = pd.merge(
    adata.obs,
    qc_cutoffs,
    left_on='sample_name',
    right_on='sample_name',
    how='left',
    sort=False
  )
  qc_df['barcode'] = adata.obs.index
  qc_df['qc'] = qc_df.apply(filter_counts, axis=1)
  print(qc_df['qc'].value_counts())  
  return qc_df

def format_qc_df(qc_df):
	result =  qc_df.groupby(['sample_name', 'qc']).size().reset_index(name='count')
	result = result.pivot(index='sample_name', columns='qc', values='count').reset_index()
	result['total'] = result['fail'] + result['pass']
	result['failed_pct'] = 100 * (result['fail'] / result['total'])
	return result

def qc_plot(adata, qc_column, qc_cutoffs, ylabel):
  ax = sc.pl.violin(
    adata,
    [qc_column],
    groupby='sample_name',
    rotation=90,
    stripplot=False,
		ylabel=ylabel,
    show=False,
		order=qc_cutoffs['sample_name'].tolist()
  )
  n_samples = len(np.unique(adata.obs['sample_name']))
  x_step = 1 / n_samples
  x_min = 0
  x_max = x_step
  ax.grid(False)
  for i in range(0, n_samples):
    y_min = qc_cutoffs.iloc[i][qc_column + '_min']
    ax.axhline(y=y_min, xmin=x_min, xmax=x_max, color='red') 
    y_max = qc_cutoffs.iloc[i][qc_column + '_max']
    ax.axhline(y=y_max, xmin=x_min, xmax=x_max, color='red')
    ax.grid(False)
    x_min = x_min + x_step
    x_max = x_max + x_step
