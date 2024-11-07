"""
This script is identical to the identically named notebook
"""

from pathlib import Path
import scanpy as sc
import cell2location
import numpy as np
from datetime import datetime

global_repo_data = Path.cwd() / '..' / '.data'
local_folder_data = Path.cwd() / '.data'
figures_dir = Path.cwd() / 'figures'
figures_dir.mkdir(exist_ok=True)
results_folder = local_folder_data / 'results' / 'COPD_to_Spatial' 
ref_run_name = results_folder / 'reference_signatures'
run_name = results_folder / 'cell2location_map'

# Set this bool if you want to redo the reference setup
rerun_reference: bool= False
if rerun_reference:

    cell_label = "cell_compartment"
    adata_copd = sc.read_h5ad(global_repo_data / 'Marburg_cell_states_locked_ctl240709.raw.h5ad')
    adata_copd = adata_copd[adata_copd.obs['group'] == 'healthy_ctrl'].copy()
    adata_copd.obs["celltype"] = adata_copd.obs[cell_label]

    adata_habermann = sc.read_h5ad(local_folder_data / 'GSE135893_ILD_annotated_fullsize.h5ad')
    adata_habermann = adata_habermann[adata_habermann.obs['Diagnosis'] == "Control"].copy()
    adata_habermann.obs["batch"] = adata_habermann.obs["Sample_Name"]

    epithelial_cell_types = [
        'Basal', 'Ciliated', 'Differentiating Ciliated', 'SCGB3A2+',
        'SCGB3A2+ SCGB1A1+',  'MUC5AC+ High', 'MUC5B+', 'AT1',
        'AT2', 'Proliferating Epithelial Cells', 'Transitional AT2', 'KRT5-/KRT17+',
    ]
    adata_habermann = adata_habermann[~adata_habermann.obs["celltype"].isin(epithelial_cell_types)]

    adata_reference = sc.concat([adata_copd, adata_habermann])

    selected = cell2location.utils.filtering.filter_genes(adata_reference, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
    adata_reference = adata_reference[:, selected].copy()

    cell2location.models.RegressionModel.setup_anndata(
        adata=adata_reference,
        batch_key="batch",
        labels_key="celltype",
    )
    mod = cell2location.models.RegressionModel(
        adata_reference,
    )
    mod.train(max_epochs=250)

    adata_reference = mod.export_posterior(
        adata_reference, sample_kwargs={'num_samples': 1000, 'batch_size': 2500}
    )
    mod.save(ref_run_name, overwrite=True)
    adata_file = ref_run_name / f'{datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}_reference.h5ad'
    adata_reference.write(adata_file)

adata_file = max(ref_run_name.glob("*.h5ad"), key=lambda f: f.stem.split('_reference')[0])
adata_reference = sc.read_h5ad(adata_file)
mod = cell2location.models.RegressionModel.load(ref_run_name, adata_reference)

if 'means_per_cluster_mu_fg' in adata_reference.varm.keys():
    inf_aver = adata_reference.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_reference.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_reference.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_reference.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_reference.uns['mod']['factor_names']

adata_spatial = sc.read_h5ad(local_folder_data / 'Franzen_L_2024_human.h5ad')
filtered_columns = list(adata_spatial.obs.columns[~adata_spatial.obs.columns.str.startswith("c2l")])
adata_spatial.obs = adata_spatial.obs[filtered_columns]
adata_spatial.obs['batch'] = adata_spatial.obs['sample_id']
samples = [sample for sample in list(adata_spatial.obs['batch'].unique()) if isinstance(sample, str)]
for sample in samples:
    adata_spatial_sample = adata_spatial[adata_spatial.obs['batch'] == sample]
    intersect = np.intersect1d(adata_spatial_sample.var_names, inf_aver.index)
    adata_spatial_sample = adata_spatial_sample[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()
    cell2location.models.Cell2location.setup_anndata(adata=adata_spatial_sample)
    mod = cell2location.models.Cell2location(
        adata_spatial_sample, cell_state_df=inf_aver,
        N_cells_per_location=30,
        detection_alpha=20
    )
    mod.train(max_epochs=30000,
          batch_size=None,
          train_size=1,
    )
    adata_spatial_sample = mod.export_posterior(
        adata_spatial_sample, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs}
    )
    mod.save(run_name, overwrite=True)
    adata_file = run_name / f"{sample}.h5ad"
    adata_spatial_sample.write(adata_file)
