import liana as li
import scanpy as sc
from collections import defaultdict
import os
import pandas as pd

data_dir = os.path.join(os.getcwd(), '.data')
os.makedirs(data_dir, exist_ok=True)
adata_path = os.path.join(data_dir, 'Marburg_cell_states_locked_ctl240709.raw.h5ad')
print('load data')
adata_all = sc.read_h5ad(adata_path) 

sample_key = 'batch'
condition_key = 'group'
groupby = 'cell_compartment'

sc.pp.normalize_total(adata_all, target_sum = 1e6, exclude_highly_expressed = True)
sc.pp.log1p(adata_all)

context_dict = adata_all.obs[[sample_key, condition_key]].drop_duplicates()
context_dict = dict(zip(context_dict[sample_key], context_dict[condition_key]))
context_dict = defaultdict(lambda: 'Unknown', context_dict)
context_dict

print('run liana rank aggregate by sample')
li.mt.rank_aggregate.by_sample(
    adata_all,
    groupby=groupby,
    sample_key=sample_key, # sample key by which we which to loop
    use_raw=False,
    seed=1789,
    min_cells=5,
    verbose=True, # use 'full' to show all verbose information
    #n_perms=100, # reduce permutations for speed
    return_all_lrs=True, # return all LR values
    )
print('successfulyl ran liana aggregate rank by sample')

liana_res = adata_all.uns['liana_res']
filename = 'liana_res_compartment.csv'
liana_res.to_csv(os.path.join(data_dir, filename))
print(f'written {filename} to {data_dir}')