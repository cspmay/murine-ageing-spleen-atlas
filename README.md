# Murine ageing spleen atlas

We profiled splenocytes from young, old, and chronically stressed old mice with natural microbiota using single-cell RNA sequencing. A comprehensive atlas was built by integrating our dataset with six publicly available datasets. 

## Data availability 
NM mice dataset analysed in this study is available in GEO under accession number GSE253406. Additional publicly available datasets were obtained from GEO (GSE132901, GSE155006, GSE179095, and GSE210719), AWS S3 bucket (https://registry.opendata.aws/tabula-muris-senis/), and GSA (CRA004660).

## Directory details
The dataset-specific folders NM_mice and SPF1-6 contain scripts to import and combine samples into a single Anndata file (combine.py) and to perform quality control (qc.py).

For atlas, the code_atlas folder contains scripts to combine the datasets into single Anndata object (combine_datasets.py), clean the Anndata object (prepare_combined.py), run scVI integration (integrate.py), manual annotations (Splenocyte_atlas_annotation.py), visualisation (atlas_plots_no_integration.py, atlas_plots.py, atlas_meta_analysis.py, combine_annotation.r, metadata_plots.r), and differential gene expression analysis (scanpy_DE.py, shared_DEGs.r).    
Scripts for re-analysis of T cells (T_combine_datasets.py, T_integrate.py, T_cell_atlas_annotation.py) and B cells (B_integrate.py, B_cell_atlas_annotation.py) have also been included. 

For NM mice, the code_NM_mice folder contains scripts to clean the Anndata object (prepare_combined.py), run scVI integration (integrate.py), manual annotations (all_cells.py, T_cells.py, Gzmk_T_cells.py, B_cells.py), visualisation (NM_plots.py), and differential gene expression analysis (DEG_export_prepare.py, DEG_pseudobulk_prepare_age.r, DEG_pseudobulk_prepare_stress.r, DESeq2_age.r, DESeq2_stress.r, DEG_venn_diagram.r).    

Shared Python utility functions have been included in the utils_py folder.

## Acknowledgements
This study was supported by the European Union's Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement No: 955321 and Estonian Research Council grants, PRG2011 and TARISTU24-TK22.
