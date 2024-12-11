# %% Define paths ----

# ==> To be set by user <==

# Path to working directory
# work_dir_path <- "path/to/working/directory/"

# Path to output directory
# out_dir_path <- "path/to/output/directory/" 

# Path to scVI adata file (h5ad)
# => Produced by scripts/Smajic_2021_pmid_34919646/04_integration.py
# adata_path <- "path/to/scvi_adata.h5ad"

# Our paths
work_dir_path = "/scratch/ben/rnaseq/"
out_dir_path = "data/billy/Smajic_2021_pmid_34919646/output/sc_umap_dotplot_tfes_oct24/"
adata_path = "data/public/Smajic_2021_pmid_34919646/output/integration/scvi_adata.h5ad"

# %% Set up ----

# Import required libraries
import os
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt

# Set working directory
os.chdir(work_dir_path)

# Define path to output files
sc.settings.figdir = out_dir_path
out_dir_path = str(sc.settings.figdir)+"/"

# %% Define analysis parameters ----

# Define layer to use for plotting
layer = "log2_cpm1p"
# layer = "cpm"

# %% Load scVI adata ----

# Load scVI adata
adata = sc.read_h5ad(filename=adata_path)

# Filter out CADPS2+ neurons
adata = adata[adata.obs["cell_type"] != 'CADPS2+ neurons']

# Preserve identical non-neuronal cell type colors across cell type and broad cell type categories 
cell_type_colors = ['#1f77b4', '#279e68', '#d62728',
                    '#aa40fc', '#8c564b', '#e377c2', '#b5bd61',
                    '#17becf', '#aec7e8', '#ffbb78', '#98df8a']

adata.uns["cell_type_colors"] = cell_type_colors
adata.uns["broad_cell_type_colors"] = [adata.uns["cell_type_colors"][0]]+adata.uns["cell_type_colors"][4:]

# %% Plot cell clustering results ----

# Expression level of TFE3 and TFEB on UMAP
fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(40, 6))
for i, j, k in zip(np.arange(4), 
                   ["cell_type", "broad_cell_type", "TFE3", "TFEB"],
                   ["cell type", "cell type - low resolution", "TFE3", "TFEB"]):
    sc.pl.umap(adata,  
               neighbors_key="neighbors_scvi",
               color=j,
               layer=layer, sort_order=True, frameon=False, edges=False, cmap="Reds", ncols=4,
               gene_symbols="gene_name",  vmin=0, vmax="p99", show=False, ax=axs[i])
    if k == "TFE3":
       axs[i].set_title(label=k, fontdict={"fontsize": 16, "weight": "bold", "color":"#0077b6"})
    elif k == "TFEB":
       axs[i].set_title(label=k, fontdict={"fontsize": 16, "weight": "bold", "color":"#7f7f7f"})
    else:
       axs[i].set_title(label=k, fontdict={"fontsize": 16, "weight": "bold"})
       axs[i].legend(fontsize=14, bbox_to_anchor=(0.85, 1), loc="upper left", frameon=False)
fig.suptitle("Smajić et al., 2022", fontsize=24, fontdict={"style": "italic"}, x = 0.2, y = 1.05)
plt.savefig(out_dir_path+"umap_"+layer+"_TFE3_TFEB_low_res_included.pdf", bbox_inches="tight", pad_inches=0)
 
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(30, 6))
for i, j, k in zip(np.arange(3), 
                   ["cell_type", "TFE3", "TFEB"],
                   ["cell type", "TFE3", "TFEB"]):
    sc.pl.umap(adata,  
               neighbors_key="neighbors_scvi",
               color=j,
               layer=layer, sort_order=True, frameon=False, edges=False, cmap="Reds", ncols=4,
               gene_symbols="gene_name",  vmin=0, vmax="p99", show=False, ax=axs[i])
    if k == "TFE3":
       axs[i].set_title(label=k, fontdict={"fontsize": 16, "weight": "bold", "color":"#0077b6"})
    elif k == "TFEB":
       axs[i].set_title(label=k, fontdict={"fontsize": 16, "weight": "bold", "color":"#7f7f7f"})
    else:
       axs[i].set_title(label=k, fontdict={"fontsize": 16, "weight": "bold"})
       axs[i].legend(fontsize=14, bbox_to_anchor=(0.85, 1), loc="upper left", frameon=False)
fig.suptitle("Smajić et al., 2022", fontsize=24, fontdict={"style": "italic"}, x = 0.2, y = 1.05)
plt.savefig(out_dir_path+"umap_"+layer+"_TFE3_TFEB.pdf", bbox_inches="tight", pad_inches=0)

# Expression level of TFE3, TFEB and MITF on UMAP
fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(40, 6))
for i, j, k in zip(np.arange(4), 
                   ["cell_type", "TFE3", "TFEB", "MITF"],
                   ["cell type", "TFE3", "TFEB", "MITF"]):
    sc.pl.umap(adata,  
               neighbors_key="neighbors_scvi",
               color=j,
               layer=layer, sort_order=True, frameon=False, edges=False, cmap="Reds", ncols=4,
               gene_symbols="gene_name",  vmin=0, vmax="p99", show=False, ax=axs[i])
    if k == "TFE3":
       axs[i].set_title(label=k, fontdict={"fontsize": 16, "weight": "bold", "color":"#0077b6"})
    elif k == "TFEB":
       axs[i].set_title(label=k, fontdict={"fontsize": 16, "weight": "bold", "color":"#7f7f7f"})
    elif k == "MITF":
       axs[i].set_title(label=k, fontdict={"fontsize": 16, "weight": "bold", "color":"#ff6d00"})
    else:
       axs[i].set_title(label=k, fontdict={"fontsize": 16, "weight": "bold"})
       axs[i].legend(fontsize=14, bbox_to_anchor=(0.85, 1), loc="upper left", frameon=False)
fig.suptitle("Smajić et al., 2022", fontsize=24, fontdict={"style": "italic"}, x = 0.2, y = 1.05)
plt.savefig(out_dir_path+"umap_"+layer+"_TFE3_TFEB_MITF.pdf", bbox_inches="tight", pad_inches=0)

# Expression levels of TFEs on dotplot (cell type resolution)
if layer == "log2_cpm1p":
  colorbar_title = "log\u2082(CPM+1)"
elif layer == "cpm":
  colorbar_title = "CPM"

# => TFE3 and TFEB only
dp = sc.pl.dotplot(
  adata=adata,
  var_names=["TFE3", "TFEB"],
  groupby="cell_type",
  gene_symbols="gene_name",
  layer=layer,
  mean_only_expressed=True,
  swap_axes=True,
  colorbar_title=colorbar_title,
  size_title="fraction of cells (%)",
  return_fig=True,
  figsize=(6,2.5)
)
dp.add_totals(color="#5d76cb")
dp.style(dot_edge_color='black',
         dot_edge_lw=0.5, 
         cmap="magma", 
         grid=True)
ax = dp.get_axes()["mainplot_ax"]
plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
for label in ax.get_yticklabels():
    label.set_weight("bold")
    if label.get_text() == "TFE3":
       label.set_color("#0077b6")
    else:
       label.set_color("#7f7f7f")   
ax.set_title(label="Smajić et al., 2022", fontdict={"fontsize": 10, "style": "italic"}, x = 0.1, y = 1.25)
plt.savefig(out_dir_path+"dotplot_"+layer+"_TFE3_TFEB.pdf", bbox_inches="tight", pad_inches=0)
 
# => TFE3, TFEB and MITF
dp = sc.pl.dotplot(
  adata=adata,
  var_names=["TFE3", "TFEB", "MITF"],
  groupby="cell_type",
  gene_symbols="gene_name",
  layer=layer,
  mean_only_expressed=True,
  swap_axes=True,
  colorbar_title=colorbar_title,
  size_title="fraction of cells (%)",
  return_fig=True,
  figsize=(6,2.5)
)
dp.add_totals(color="#5d76cb")
dp.style(dot_edge_color='black',
         dot_edge_lw=0.5, 
         cmap="magma", 
         grid=True)
ax = dp.get_axes()["mainplot_ax"]
plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
for label in ax.get_yticklabels():
    label.set_weight("bold")
    if label.get_text() == "TFE3":
       label.set_color("#0077b6")
    elif label.get_text() == "TFEB":
       label.set_color("#7f7f7f")  
    else:
       label.set_color("#ff6d00")   
ax.set_title(label="Smajić et al., 2022", fontdict={"fontsize": 10, "style": "italic"}, x = 0.1, y = 1.25)
plt.savefig(out_dir_path+"dotplot_"+layer+"_TFE3_TFEB_MITF.pdf", bbox_inches="tight", pad_inches=0)
# %%
