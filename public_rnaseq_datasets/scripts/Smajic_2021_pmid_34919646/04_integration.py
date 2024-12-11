# %% Define paths ----

# ==> To be set by user <==

# Path to working directory
# work_dir_path <- "path/to/working/directory/"

# Path to output directory
# out_dir_path <- "path/to/output/directory/" 

# Path to raw adata file (h5ad)
# => Produced by scripts/Smajic_2021_pmid_34919646/03_build_adata.py
# adata_path <- "path/to/raw_adata.h5ad"

# Our paths
work_dir_path = "/scratch/ben/rnaseq/"
out_dir_path = "data/public/Smajic_2021_pmid_34919646/output/integration/"
adata_path = "data/public/Smajic_2021_pmid_34919646/output/build_adata/raw_adata.h5ad"

# %% Set up ----

# Import required libraries
import os
import numpy as np
import scanpy as sc
import scvi
import torch

# Set working directory
os.chdir(work_dir_path)

# Define path to output files
sc.settings.figdir = out_dir_path
out_dir_path = str(sc.settings.figdir)+"/"

torch.set_float32_matmul_precision('high')

# %% Define analysis parameters ----

# Define batch variable
batch_key = "donor_id"

# Define layer to use for plotting
layer = "log2_cpm1p"

# %% Load adata ----

# Load adata
adata = sc.read_h5ad(filename=adata_path)

# %% scVI integration and clustering ----

# Find HVG on each batch separately
sc.pp.highly_variable_genes(adata=adata, 
                            flavor="cell_ranger", 
                            layer=layer,
                            n_top_genes=3000,
                            batch_key=batch_key)

# Prepare adata for scVI integration
# => Create a copy of adata containing only the HVG
adata_hvg = adata[:, adata.var["highly_variable"]].copy()

# => Set up the HVG-only adata
scvi.model.SCVI.setup_anndata(adata=adata_hvg, 
                              layer="raw_count", 
                              batch_key=batch_key)

scvi.model.SCVI.setup_anndata(adata=adata_hvg, 
                              layer="raw_count", 
                              batch_key=batch_key,
                              continuous_covariate_keys=["total_counts", "pct_counts_ribo"])

# Build scVI model
model_scvi = scvi.model.SCVI(adata=adata_hvg)
model_scvi.view_anndata_setup()

# Train scVI model
model_scvi.train()

# Write ref_model_scvi to disk
model_scvi.save(dir_path=out_dir_path+"scvi_model", overwrite=True)

# Load ref_model_scvi from disk
model_scvi = scvi.model.SCVI.load(dir_path=out_dir_path+"scvi_model", adata=adata_hvg)

# Pass the scVI-normalised counts to the HVG-only adata
adata_hvg.layers["layers_scvi"] = model_scvi.get_normalized_expression(library_size=1e4)

# Pass the scVI-corrected embedding to the HVG-only adata
adata_hvg.obsm["X_scvi"] = model_scvi.get_latent_representation()
adata_hvg.obsm["X_scvi_mde"] = scvi.model.utils.mde(adata_hvg.obsm["X_scvi"])

# Transfer the scVI-corrected embedding to the full adata
for i in ["X_scvi", "X_scvi_mde"]:
  adata.obsm[i] = adata_hvg.obsm[i]

# Compute neighbors using the scVI-corrected embedding
sc.pp.neighbors(adata, use_rep="X_scvi", key_added="neighbors_scvi")

# Compute UMAP
sc.tl.umap(adata, neighbors_key="neighbors_scvi")

# Find clusters
leiden_key = "leiden scVI "
for i in np.round(np.arange(0.1, 1, 0.1), 1):
  sc.tl.leiden(adata, 
               neighbors_key="neighbors_scvi",
               key_added=leiden_key+str(i),
               resolution=i)

# Write scVI adata to disk
adata.copy().write_h5ad(filename=out_dir_path+"scvi_adata.h5ad") 

# %% Plot cell clustering results ----

# Expression level of marker genes on UMAP
sc.pl.umap(adata, 
           neighbors_key="neighbors_scvi",
           color=["cell_type", "donor_id", "MAP2", "TH", "CADPS2",
                  "SLC17A6", "GAD2", "GRIK1", "AQP4", "FOXJ1",
                  "VCAN", "MOBP", "CD74", "CLDN5", "PDGFRB"],
           title=["cell type", "donor ID", "MAP2", "TH", "CADPS2",
                  "SLC17A6", "GAD2", "GRIK1", "AQP4", "FOXJ1",
                  "VCAN", "MOBP", "CD74", "CLDN5", "PDGFRB"],
           layer=layer, sort_order=False, frameon=False, edges=False, cmap="Reds", ncols=5,
           gene_symbols="gene_name",  vmin=0, vmax="p99", wspace=0.3,
           legend_fontsize="medium", save="_"+layer+"_marker_genes.pdf")

# %%
