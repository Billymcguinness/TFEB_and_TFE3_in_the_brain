# %% Define paths ----

# ==> To be set by user <==

# Path to working directory
# work_dir_path <- "path/to/working/directory/"

# Path to output directory
# out_dir_path <- "path/to/output/directory/" 

# Path to raw adata file (h5ad)
# => Produced by scripts/Ximerakis_2019_pmid_31551601/02_build_adata.py
# adata_path <- "path/to/raw_adata.h5ad"

# Our paths (for internal use only)
# work_dir_path = "/scratch/ben/rnaseq/"
# out_dir_path = "data/public/Ximerakis_2019_pmid_31551601/output/integration/"
# adata_path = "data/public/Ximerakis_2019_pmid_31551601/output/build_adata/raw_adata.h5ad"

# %% Set up ----

# Import required libraries
import os
import scanpy as sc

# Set working directory
os.chdir(work_dir_path)

# Define path to output files
sc.settings.figdir = out_dir_path
out_dir_path = str(sc.settings.figdir)+"/"

# %% Define analysis parameters ----

# Define batch variable
batch_key = "animal_id"

# %% Load adata ----

# Load adata
adata = sc.read_h5ad(filename=adata_path)

# %% Cluster Harmony-integrated data ----

# Find HVG on each batch separately
sc.pp.highly_variable_genes(adata=adata, flavor="seurat", layer="loge_cp10k1p", batch_key=batch_key)

# Perform Z-score transformation
sc.pp.scale(adata)

# Compute PCA
sc.pp.pca(data=adata, use_highly_variable=True)

# Perform Harmony integration
sc.external.pp.harmony_integrate(adata, key=batch_key)
 
# Compute neighbors using the Harmony-corrected embedding
sc.pp.neighbors(adata, use_rep="X_pca_harmony")

# Compute UMAP
sc.tl.umap(adata)

# Write harmony adata to disk
adata.copy().write_h5ad(filename=out_dir_path+"harmony_adata.h5ad") 

# %% Plot cell clustering results ----

# Expression level of marker genes on UMAP
sc.pl.umap(adata, 
           color=["animal_id", "cell_class", "cell_type",
                  "Snap25", "Cldn10", "Dynlrb2", "Sox10", "C1qc", "Esam"],
           title=["animal ID", "lineage", "cell type",
                  "Snap25", "Cldn10", "Dynlrb2", "Sox10", "C1qc", "Esam"],
           layer="loge_cp10k1p", sort_order=False, frameon=False, edges=False, cmap="Reds", ncols=3,
           gene_symbols="gene_name",  vmin=0, vmax="p99", wspace=0.2,
           legend_fontsize="medium", save="_loge_cp10k1p_marker_genes.pdf")

sc.pl.umap(adata, 
           color=["cell_class", "cell_type", "age_status",
                  "Snap25", "Cldn10", "Dynlrb2", "Sox10", "C1qc", "Esam"],
           title=["lineage", "cell type", "age status",
                  "Snap25", "Cldn10", "Dynlrb2", "Sox10", "C1qc", "Esam"],
           layer="loge_cp10k1p", sort_order=False, frameon=False, edges=False, cmap="Reds", ncols=3,
           gene_symbols="gene_name",  vmin=0, vmax="p99", wspace=0.2,
           legend_fontsize="medium", save="_loge_cp10k1p_marker_genes_age_res.pdf")

# %%
