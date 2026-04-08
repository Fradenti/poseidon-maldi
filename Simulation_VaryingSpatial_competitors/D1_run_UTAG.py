import numpy as np
import pandas as pd

from utag import utag
import scanpy as sc
import squidpy as sq
import anndata
import matplotlib.pyplot as plt

i_sim = 1
i_str = 1

strength = [0.5,1,2,10] # 0.5 edit

# import sys
# i_sim = int(sys.argv[1])
# print(f"Running simulation with i_sim = {i_sim}")

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--i_sim", type=int, required=True)
parser.add_argument("--i_str", type=int, required=True)
args = parser.parse_args()
i_sim = args.i_sim
i_str = args.i_str
print(f"Running with i_sim={i_sim}, i_str={i_str}")

stre = strength[i_str-1]

full_data = pd.read_csv(f"data/simulated_dataC/simulated_data_str{stre}_{i_sim}.csv")
metadata_cols = ['cell_id', 'x', 'y', 'true_cluster', 'sample_id'] 
obs = full_data[metadata_cols].copy()
expression_cols = [col for col in full_data.columns if col not in metadata_cols]
X = full_data[expression_cols].values  # convert to numpy array
var = pd.DataFrame(index=expression_cols)
var['feature'] = expression_cols 

minimal_adata = anndata.AnnData(
    X = np.array(X).astype(np.float64),
    obs = obs,
    var = var
)

minimal_adata.obsm['spatial'] = np.array(minimal_adata.obs[['x', 'y']])
del minimal_adata.obs['y']
del minimal_adata.obs['x']


# Run UTAG on provided data
utag_single_results = utag(
    minimal_adata,
    slide_key="sample_id",
    max_dist=1.5, # equivalent to QUEEN adjacency matrix
    normalization_mode='l1_norm',
    apply_clustering=True,
    clustering_method = 'leiden', 
    resolutions = [0.2, 0.3, 0.4, 0.5] # [0.3,1.0,2.0]
)

true_clusters = minimal_adata.obs['true_cluster']
utag_labels_02 = utag_single_results.obs['UTAG Label_leiden_0.2']
utag_labels_03 = utag_single_results.obs['UTAG Label_leiden_0.3']
utag_labels_1 = utag_single_results.obs['UTAG Label_leiden_0.4']
utag_labels_2 = utag_single_results.obs['UTAG Label_leiden_0.5']
cell_ids = utag_single_results.obs['cell_id']
sample_ids = utag_single_results.obs['sample_id']
spatial_coords = utag_single_results.obsm['spatial']
x_coords = spatial_coords[:, 0]
y_coords = spatial_coords[:, 1]

results_df = pd.DataFrame({
    'cell_id': cell_ids,
    'sample_id': sample_ids,
    'true_cluster': true_clusters,
    'utag_cluster_02': utag_labels_02,
    'utag_cluster_03': utag_labels_03,
    'utag_cluster_04': utag_labels_1,
    'utag_cluster_05': utag_labels_2,
    'x_coord': x_coords,
    'y_coord': y_coords
})

# Export to CSV
results_df.to_csv(f'results/RData/UTAG/utag_results_str{stre}_{i_sim}.csv', index=False)
