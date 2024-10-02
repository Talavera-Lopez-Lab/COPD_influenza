'''
Using squidpy-env environment 
'''
import liana as li
import scanpy as sc
import cell2cell as c2c
from collections import defaultdict
import os
from tqdm import tqdm


data_dir = os.path.join(os.getcwd(), '.data')
os.makedirs(data_dir, exist_ok=True)
adata_path = os.path.join(data_dir, 'Marburg_cell_states_locked_ctl240709.raw.h5ad')
print('load data')
adata_all = sc.read_h5ad(adata_path) 

conditions = list(adata_all.obs['group'].unique())
for condition in tqdm(conditions):
    adata = adata_all[adata_all.obs['group'] == condition]
    sample_key = 'batch'
    condition_key = 'group'
    groupby = 'cell_compartment'
    
    sc.pp.normalize_total(adata, target_sum = 1e6, exclude_highly_expressed = True)
    sc.pp.log1p(adata)
    
    context_dict = adata.obs[[sample_key, condition_key]].drop_duplicates()
    context_dict = dict(zip(context_dict[sample_key], context_dict[condition_key]))
    context_dict = defaultdict(lambda: 'Unknown', context_dict)
    context_dict
    
    print('run liana rank aggregate by sample')
    li.mt.rank_aggregate.by_sample(
        adata,
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
    
    file_name = 'anndata_liana.h5ad'
    file_path = os.path.join(data_dir, file_name)
    # Write
    #print(f'writing {file_name} to {file_path}')
    #adata.write_h5ad(file_path)
    #print(f'succesfully written {file_name} to {file_path}')
    
    print('loading initial tensor')
    tensor = li.multi.to_tensor_c2c(adata,
                                    sample_key=sample_key,
                                    score_key='magnitude_rank', # can be any score from liana
                                    how='outer_cells' # how to join the samples
                                    )
    print('initial tensor created')
    
    file_name = 'initial_tensor_compartment.pkl'
    file_path = os.path.join(data_dir, file_name)
    
    # Save
    #print(f'writing {file_name} to {file_path}')
    #c2c.io.export_variable_with_pickle(tensor, file_path)
    #print(f'written {file_name} to {file_path}')
    
    
    print('create tensor metadata')
    context_dict = adata.obs[[sample_key, condition_key]].drop_duplicates()
    context_dict = dict(zip(context_dict[sample_key], context_dict[condition_key]))
    context_dict = defaultdict(lambda: 'Unknown', context_dict)
    
    tensor_meta = c2c.tensor.generate_tensor_metadata(interaction_tensor=tensor,
                                                      metadata_dicts=[context_dict, None, None, None],
                                                      fill_with_order_elements=True
                                                      )
    print('running tensor cell2cell pipeline')
    tensor = c2c.analysis.run_tensor_cell2cell_pipeline(tensor,
                                                        tensor_meta,
                                                        copy_tensor=True, # Whether to output a new tensor or modifying the original
                                                        rank=None, # Number of factors to perform the factorization. If None, it is automatically determined by an elbow analysis. Here, it was precomuputed.
                                                        tf_optimization='regular', # To define how robust we want the analysis to be.
                                                        random_state=1789, # Random seed for reproducibility
                                                        backend='pytorch', # This enables a banckend that supports using a GPU.
                                                        device='cuda:1', # Device to use. If using GPU and PyTorch, use 'cuda'. For CPU use 'cpu'
                                                        elbow_metric='error', # Metric to use in the elbow analysis.
                                                        smooth_elbow=False, # Whether smoothing the metric of the elbow analysis.
                                                        upper_rank=25, # Max number of factors to try in the elbow analysis
                                                        tf_init='random', # Initialization method of the tensor factorization
                                                        tf_svd='numpy_svd', # Type of SVD to use if the initialization is 'svd'
                                                        cmaps=None, # Color palettes to use in color each of the dimensions. Must be a list of palettes.
                                                        sample_col='Element', # Columns containing the elements in the tensor metadata
                                                        group_col='Category', # Columns containing the major groups in the tensor metadata
                                                        fig_fontsize=14,
                                                        output_fig=True, # Whether to output the figures. If False, figures won't be saved a files if a folder was passed in output_folder.
                                                        )
    print('succesfully created final tensor')
    file_name = f'final_tensor_compartment_{condition}.pkl'
    file_path = os.path.join(data_dir, file_name)
    
    print(f'writing {file_name} to {file_path}')
    # Save
    c2c.io.export_variable_with_pickle(tensor, file_path)
    print(f'succesfully written {file_name} to {file_path}')