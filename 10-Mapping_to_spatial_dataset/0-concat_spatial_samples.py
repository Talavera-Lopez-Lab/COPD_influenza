from pathlib import Path
import scanpy as sc
import squidpy as sq
import pandas as pd

local_folder_data = Path.cwd() / '.data'
metadata = pd.read_csv(local_folder_data / 'hs_visium_metadata.tsv', sep='\t', index_col='sample_id')

adatas = []
for folder in (local_folder_data / 'raw').iterdir():
    adata = sq.read.visium(folder)
    adata.obs['sample_id'] = folder.name
    adatas.append(adata)
print('samples Loaded')
for adata in adatas:
    adata.var_names_make_unique()
adata = sc.concat(adatas, uns_merge='unique')
adata.obs = adata.obs.reset_index().merge(metadata, on='sample_id', how='left').set_index('index')
adata.obs_names_make_unique()
sc.pp.filter_cells(adata, min_counts=200)
sc.pp.filter_genes(adata, min_cells=3)

filename = 'Franzen_L_2024_human.h5ad'
adata.write_h5ad(local_folder_data / filename)
print(f'{filename} written')