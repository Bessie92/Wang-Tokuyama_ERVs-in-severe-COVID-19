# Import packages 
import warnings
warnings.filterwarnings("ignore")

import logging
import random

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
#import pertpy
import scanpy as sc
import seaborn as sns

sc.settings.verbosity = 0

# Load in adata
adata = sc.read_h5ad("../AnnData/merged_adata_filtered_with_cell_types_and_clustering.h5ad")

# We need raw counts
adata.X = adata.layers["raw"].copy()

# set categorical metadata to be indeed categorical to create pseudobulks
adata.obs["condition"] = adata.obs["condition"].astype("category")
adata.obs["sample"] = adata.obs["sample"].astype("category")
adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")


def aggregate_and_filter(
    adata,
    cell_identity,
    donor_key="sample",
    condition_key="label", # in our case 'condition'
    cell_identity_key="cell_type",
    obs_to_keep=None,  # which additional metadata to keep in the pseudobulk adata
    replicates_per_patient=1, # how many pseudoreplicates to make for each patient
):
    # 1.subset adata to the given cell identity
    if obs_to_keep is None:
        obs_to_keep = []
    adata_cell_pop = adata[adata.obs[cell_identity_key] == cell_identity].copy()
    
    # 2.check which donors to keep according to the number of cells specified with NUM_OF_CELL_PER_DONOR
    size_by_donor = adata_cell_pop.obs.groupby([donor_key]).size()
    donors_to_drop = [
        donor
        for donor in size_by_donor.index
        if size_by_donor[donor] <= NUM_OF_CELL_PER_DONOR
    ]
    if len(donors_to_drop) > 0:
        print("Dropping the following samples:")
        print(donors_to_drop)
    
    # 3.prepare empty DataFrame to hold the aggregated counts
    df = pd.DataFrame(columns=[*adata_cell_pop.var_names, *obs_to_keep])
    
    # 4.loop through donors/samples
    adata_cell_pop.obs[donor_key] = adata_cell_pop.obs[donor_key].astype("category")
    for i, donor in enumerate(donors := adata_cell_pop.obs[donor_key].cat.categories):
        print(f"\tProcessing donor {i+1} out of {len(donors)}...", end="\r")
        if donor not in donors_to_drop:
            #### *** Randomly shuffles and splits the donor's cells into replicates_per_patient (e.g. 2) equally sized subgroups.
            adata_donor = adata_cell_pop[adata_cell_pop.obs[donor_key] == donor]
            indices = list(adata_donor.obs_names)
            random.shuffle(indices)
            indices = np.array_split(np.array(indices), replicates_per_patient)
            for i, rep_idx in enumerate(indices):
                adata_replicate = adata_donor[rep_idx]
                # specify how to aggregate: sum gene expression for each gene for each donor and also keep the condition information
                agg_dict = {gene: "sum" for gene in adata_replicate.var_names}
                for obs in obs_to_keep:
                    agg_dict[obs] = "first"
                # create a df with all genes, donor and condition info
                df_donor = pd.DataFrame(adata_replicate.X.A)
                df_donor.index = adata_replicate.obs_names
                df_donor.columns = adata_replicate.var_names
                df_donor = df_donor.join(adata_replicate.obs[obs_to_keep])
                # aggregate
                df_donor = df_donor.groupby(donor_key).agg(agg_dict)
                df_donor[donor_key] = donor
                df.loc[f"donor_{donor}_{i}"] = df_donor.loc[donor]
    print("\n")
    # create AnnData object from the df
    adata_cell_pop = sc.AnnData(
        df[adata_cell_pop.var_names], obs=df.drop(columns=adata_cell_pop.var_names)
    )
    return adata_cell_pop


NUM_OF_CELL_PER_DONOR = 30
obs_to_keep = ["condition", "sample","cell_type"]  # keep this metadata in pseudobulked object

# process first cell type separately
cell_type = adata.obs["cell_type"].cat.categories[0]
print(f'Processing {cell_type} (1 out of {len(adata.obs["cell_type"].cat.categories)})...')
adata_pb = aggregate_and_filter(
    adata,
    cell_identity=cell_type,
    donor_key="sample",
    condition_key="condition",
    cell_identity_key="cell_type",
    obs_to_keep=obs_to_keep,
    replicates_per_patient=2
)

# process remaining cell types and concatenate
for i, cell_type in enumerate(adata.obs["cell_type"].cat.categories[1:]):
    print(f'Processing {cell_type} ({i+2} out of {len(adata.obs["cell_type"].cat.categories)})...')
    adata_cell_type = aggregate_and_filter(
        adata,
        cell_identity=cell_type,
        donor_key="sample",
        condition_key="condition",
        cell_identity_key="cell_type",
        obs_to_keep=obs_to_keep,
        replicates_per_patient=2
    )
    adata_pb = adata_pb.concatenate(adata_cell_type)

# save the raw counts in the 'counts' layer, then normalize the counts and calculate the PCA coordinates for the normalized pseudobulk counts.
adata_pb.layers['counts'] = adata_pb.X.copy()
sc.pp.normalize_total(adata_pb, target_sum=1e6)
sc.pp.log1p(adata_pb)
sc.pp.pca(adata_pb)

# look at created pseudo-replicates on a PCA plot and color by all the available metadata to see if there are any confounding factors that we might want to include in the design matrix
# also add a lib_size and log_lib_size columns to check if there is a correlation between library size and PC component
sc.pl.pca(adata_pb, color=adata_pb.obs, ncols=1, size=300)

# edgeR takes raw counts as input, so we put counts back into the .X field before we proceed
adata_pb.X = adata_pb.layers['counts'].copy()

# Subset to one cell type
adata_neu = adata_pb[adata_pb.obs["cell_type"] == "Neutrophils"].copy()
print(adata_neu)

# Import proviral ERV names 
ERVs_list = pd.read_csv("proviral_ERV_list.csv")
common_ERVs = list(set(ERVs_list.erv_id) & set(adata.var_names))
ERV_adata = adata[:, common_ERVs] 
ERV_adata

# 2. extract replicate info
adata_neu.obs["replicate"] = (
    adata_neu.obs_names
    .str.split("_").str[-1]
    .str.split("-").str[0]
    .astype(str))

adata_neu_clean = adata_neu.copy()

if sparse.issparse(adata_neu_clean.X):
    adata_neu_clean.X = adata_neu_clean.X.tocsr().astype(np.float32)
else:
    adata_neu_clean.X = np.asarray(adata_neu_clean.X).astype(np.float32)


for col in adata_neu_clean.obs.columns:
    if pd.api.types.is_categorical_dtype(adata_neu_clean.obs[col]):
        adata_neu_clean.obs[col] = adata_neu_clean.obs[col].astype(str)


for slot in ["obsm", "obsp", "varm", "layers"]:
    getattr(adata_neu_clean, slot).clear()


adata_neu_clean.write("~/SRP279746_sc/Pseudobulk_EdgeR/pseudobulk_neutrophils.h5ad")