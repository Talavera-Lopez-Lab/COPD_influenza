from pathlib import Path
import scanpy as sc
import squidpy as sq

local_folder_data = Path.cwd() / '.data'

adatas = [sq.read.visium(folder) for folder in (local_folder_data / 'raw').iterdir()]
print('samples Loaded')
for adata in adatas:
    adata.var_names_make_unique()
adata = sc.concat(adatas)
adata.obs_names_make_unique()
del adatas

filename = 'Franzen_L_2024_humand.h5ad'
adata.write_h5ad(local_folder_data / filename)
print(f'{filename} written')
