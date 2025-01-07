# %% Define paths ----

# ==> To be set by user <==

# Path to working directory
# work_dir_path = "path/to/working/directory/"

# Path to output directory
# out_dir_path = "path/to/output/directory/" 

# Path to directory containing animal-specific count matrices
# => Available at https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE129788&format=file
# cnt_dir_path = "path/to/GSE129788_RAW/"

# Path to tissue metadata file (csv)
# => Produced by scripts/Ximerakis_2019_pmid_31551601/01_make_sample_metadata.R
# tissue_metadata_path = "path/to/tissue_metadata.csv"

# Path to sample metadata file (csv)
# => Produced by scripts/Ximerakis_2019_pmid_31551601/01_make_sample_metadata.R
# sample_metadata_path = "path/to/sample_metadata.csv"

# Our paths (for internal use only)
# work_dir_path = "/scratch/ben/rnaseq/"
# out_dir_path = "data/public/Ximerakis_2019_pmid_31551601/output/build_adata/"
# cnt_dir_path = "seq_data/public/Ximerakis_2019_pmid_31551601/supplementary_files/GSE129788_RAW/"
# tissue_metadata_path = "data/public/Ximerakis_2019_pmid_31551601/output/make_sample_metadata/tissue_metadata.csv"
# sample_metadata_path = "data/public/Ximerakis_2019_pmid_31551601/output/make_sample_metadata/sample_metadata.csv"

# %% Set up ----

# Import required libraries
import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

# Set working directory
os.chdir(work_dir_path)

# %% Load metadata ----

# Load tissue metadata
tissue_metadata = pd.read_csv(tissue_metadata_path)

# Load sample metadata
sample_metadata = pd.read_csv(sample_metadata_path)

# %% Build adata ----

# Build adata from tissue-specific count matrices
adatas = {}

for i in np.arange(len(tissue_metadata)):

  adata = sc.read_text(cnt_dir_path+tissue_metadata["sample_id"][i]+"_"+tissue_metadata["animal_id"][i]+"_10X.txt.gz").transpose()
  obs_metadata = tissue_metadata.iloc[[i],:].loc[:, ['run_id', 'sample_id', 'animal_id_description', 'animal_id', 'age', 'age_status']].copy()

  for j in obs_metadata.columns:
    adata.obs[j] = obs_metadata[j].iloc[0]

  adatas[tissue_metadata["sample_id"][i]] = adata

adata = sc.concat(adatas = adatas)

# Convert counts to sparse matrix 
adata.X = sparse.csr_matrix(adata.X.copy()) 

# Set loge CP10K+1 layer
adata.layers["loge_cp10k1p"] = adata.X.copy()

# %% Add sample metadata ----

# Filter/reorder observations to match sample metadata
adata = adata.copy()[sample_metadata["animal_num_barcode"]]

# Add sample metadata to adata.obs
adata.obs = pd.merge(left=adata.obs, 
                     how="inner",
                     right=sample_metadata.set_index("animal_num_barcode", drop=False).drop(columns="age_status"),
                     left_index=True, right_index=True)

## Cell type description
# OPC: oligodendrocyte precursor cells 
# OLG: oligodendrocytes 
# OEG: olfactory ensheathing glia 
# NSC: neural stem cells 
# ARP: astrocyte-restricted precursors
# ASC: astrocytes 
# NRP: neuronal-restricted precursors 
# ImmN: immature neurons
# mNEUR: mature neurons 
# NendC: neuroendocrine cells
# EPC: ependymocytes 
# HypEPC: hypendymal cells 
# TNC: tanycytes 
# CPC: choroid plexus epithelial cells
# EC: endothelial cells 
# PC: pericytes 
# VSMC: vascular smooth muscle cells 
# Hb-VC: hemoglobin-expressing vascular cells 
# VLMC: vascular and leptomeningeal cells 
# ABC: arachnoid barrier cells
# MG: microglia 
# MNC: monocytes 
# MAC: macrophages 
# DC: dendritic cells
# NEUT: neutrophils 

## Lineage (cell class) description
# OLG_Lin: oligodendrocyte lineage, 
# ASC_Lin: astrocyte lineage and stem cells, 
# NEURON_Lin: neuronal lineage, 
# EPC_Lin: ependymal cells, 
# VASC_Lin: vasculature cells 
# IMMUNE_Lin: immune cells

# Set cell class as category 
adata.obs["cell_class"] = adata.obs["cell_class"].astype("category")
adata.obs["cell_class"] = adata.obs["cell_class"].cat.reorder_categories(
  ['NEURON_Lin', 'ASC_Lin', 'EPC_Lin', 'OLG_Lin', 'IMMUNE_Lin', 'VASC_Lin']
  )

# Set age_status as category and provide colors for plotting 
age_status_colors = {"young": "#5d76cb", "old": "#ca3767"}
adata.obs["age_status"] = adata.obs["age_status"].astype("category")
adata.obs["age_status"] = adata.obs["age_status"].cat.reorder_categories(age_status_colors.keys())
adata.uns["age_status_colors"] = list(age_status_colors.values())

# Write adata to disk
adata.copy().write_h5ad(filename=out_dir_path+"raw_adata.h5ad") 

# %%
