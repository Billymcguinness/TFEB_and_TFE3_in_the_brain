# %% Define paths ----

# ==> To be set by user <==

# Path to working directory
# work_dir_path = "path/to/working/directory/"

# Path to output directory
# out_dir_path = "path/to/output/directory/" 

# Path to count matrix file (h5)
# => Produced by scripts/Smajic_2021_pmid_34919646/02_prepare_count_matrix.R
# cnt_path = "path/to/raw_count.h5"

# Path to sample metadata file (csv)
# => Produced by scripts/Smajic_2021_pmid_34919646/01_make_sample_metadata.R
# sample_metadata_path = "path/to/sample_metadata.csv"

# Path to gene metadata file (csv)
# => Produced by scripts/feature_metadata/01_make_feature_metadata.R
# gene_metadata_path = "path/to/gene_metadata.csv"

# Our paths (for internal use only)
# work_dir_path = "/scratch/ben/rnaseq/"
# out_dir_path = "data/public/Smajic_2021_pmid_34919646/output/build_adata/"
# cnt_path = "data/public/Smajic_2021_pmid_34919646/output/prepare_count_matrix/raw_count.h5"
# sample_metadata_path = "data/public/Smajic_2021_pmid_34919646/output/make_sample_metadata/sample_metadata.csv" 
# gene_metadata_path = "ref_data/feature_metadata/gene_metadata.csv"

# %% Set up ----

# Import required libraries
import os
import numpy as np
import pandas as pd
import scanpy as sc

# Set working directory
os.chdir(work_dir_path)

# %% Build adata ----

# Construct adata from raw count matrix
adata = sc.read_10x_h5(filename=cnt_path)

# Set raw count layer
adata.layers["raw_count"] = adata.X.copy()

# %% Add CPM data ----

# Compute CPM
adata.layers["cpm"] = sc.pp.normalize_total(adata=adata, target_sum=1e6, copy=True).X

# Compute log2 CPM+1
adata.layers["log2_cpm1p"] = sc.pp.log1p(adata.layers["cpm"].copy(), base=2)

# %% Add sample/gene metadata ----

# Load sample metadata
sample_metadata = pd.read_csv(sample_metadata_path)
sample_metadata = sample_metadata.set_index(keys="barcode", drop=False)

# Reorder rows to match adata.obs_names
sample_metadata = sample_metadata.loc[adata.obs_names]

# Add sample metadata to adata
adata.obs = sample_metadata

# Set cell type as category and reorder levels
adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")
adata.obs["cell_type"] = adata.obs["cell_type"].cat.reorder_categories(
['DaNs', 'CADPS2+ neurons', 'Excitatory', 
 'Inhibitory', 'GABA', 'Astrocytes',
 'Ependymal', 'OPCs', 'Oligodendrocytes', 
 'Microglia','Endothelial cells', 'Pericytes']
 )

# Define broad cell type and set as category
adata.obs["broad_cell_type"] = np.where(
  adata.obs["cell_type"].isin(['DaNs', 'CADPS2+ neurons', 'Excitatory', 'Inhibitory', 'GABA']), 
  "Neurons", 
  adata.obs["cell_type"]
  )

adata.obs["broad_cell_type"] = adata.obs["broad_cell_type"].astype("category")
adata.obs["broad_cell_type"] = adata.obs["broad_cell_type"].cat.reorder_categories(
  ['Neurons', 'Astrocytes',
 'Ependymal', 'OPCs', 'Oligodendrocytes', 
 'Microglia','Endothelial cells', 'Pericytes']
  )

# Load gene metadata
gene_metadata = pd.read_csv(gene_metadata_path)
gene_metadata = gene_metadata.set_index(keys="ensembl_gene_id", drop=False)

# Update gene metadata with ensembl gene IDs present only in adata
gene_metadata = pd.DataFrame(
  data={"ensembl_gene_id": adata.var_names}
  ).merge(
    right=gene_metadata.reset_index(drop=True),
    how="left",
    on="ensembl_gene_id"
    ).set_index(
      keys="ensembl_gene_id",
      drop=False
      )

# Annotate mitochondrial and ribosomal genes
gene_metadata["mt"] = gene_metadata["gene_name"].str.startswith("MT-", na=False)
gene_metadata["ribo"] = gene_metadata["gene_name"].str.startswith((("RPS", "RPL")), na=False)

# Annotate genes located on chromosomes X, Y and M
gene_metadata["chr_y"] = gene_metadata["chr_name"] == "chrY"
gene_metadata["chr_x"] = gene_metadata["chr_name"] == "chrX"
gene_metadata["chr_m"] = gene_metadata["chr_name"] == "chrM"

# Add gene metadata to adata
adata.var = gene_metadata
adata.var_names = adata.var["ensembl_gene_id"]
for i in [x for x in adata.var.columns if "marker" in x]:
  adata.var[i] = adata.var[i].astype("category")

# Compute QC metrics 
sc.pp.calculate_qc_metrics(adata=adata,
                           expr_type="counts",
                           var_type="genes",
                           qc_vars=["mt", "ribo", "chr_x", "chr_y", "chr_m"],
                           percent_top=None,
                           inplace=True,
                           log1p=False)

# Add XIST counts to observations metadata 
adata.obs["xist_count"] = adata.copy().X[:,adata.var["gene_name"].str.match("XIST", na=False)].toarray()

# Write raw adata to disk
adata.copy().write_h5ad(filename=out_dir_path+"raw_adata.h5ad") 

# %%
