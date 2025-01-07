# %% Define paths ----

# ==> To be set by user <==

# Path to working directory
# work_dir_path = "path/to/working/directory/"

# Path to output directory
# out_dir_path = "path/to/output/directory/" 

# Path to Harmony adata file (h5ad)
# => Produced by scripts/Ximerakis_2019_pmid_31551601/03_integration.py
# adata_path = "path/to/harmony_adata.h5ad"

# Our paths (for internal use only)
# work_dir_path = "/scratch/ben/rnaseq/"
# out_dir_path = "data/billy/Ximerakis_2019_pmid_31551601/output/sc_umap_dotplot_tfes_oct24/"
# adata_path = "data/public/Ximerakis_2019_pmid_31551601/output/integration/harmony_adata.h5ad"

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

# %% Load Harmony adata ----

# Load Harmony adata
adata = sc.read_h5ad(filename=adata_path)

# %% Plot cell clustering results ----

# Expression level of Tfe3 and Tfeb on UMAP
fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(40, 6))
for i, j, k in zip(np.arange(4), 
                   ["cell_type", "cell_class", "Tfe3", "Tfeb"],
                   ["cell type", "lineage", "Tfe3", "Tfeb"]):
    sc.pl.umap(adata,  
               neighbors_key="neighbors_scvi",
               color=j,
               layer="loge_cp10k1p", sort_order=True, frameon=False, edges=False, cmap="Reds", ncols=4,
               gene_symbols="gene_name",  vmin=0, vmax="p99", show=False, ax=axs[i])
    if k == "Tfe3":
      axs[i].set_title(label=k, fontdict={"fontsize": 16, "weight": "bold", "color":"#0077b6"})
    elif k == "Tfeb":
      axs[i].set_title(label=k, fontdict={"fontsize": 16, "weight": "bold", "color":"#7f7f7f"})
    else:
      axs[i].set_title(label=k, fontdict={"fontsize": 16, "weight": "bold"})
      if k == "cell type":
        axs[i].legend(fontsize=14, bbox_to_anchor=(0.90, 1), loc="upper left", frameon=False, ncol=2, columnspacing=0.1)
      else:
        axs[i].legend(fontsize=14, bbox_to_anchor=(0.90, 1), loc="upper left", frameon=False)
fig.suptitle("Ximerakis et al., 2019", fontsize=24, fontdict={"style": "italic"}, x = 0.2, y = 1.05)
plt.savefig(out_dir_path+"umap_loge_cp10k1p_Tfe3_Tfeb.pdf", bbox_inches="tight", pad_inches=0)

# Expression level of Tfe3, Tfeb and Mitf on UMAP
fig, axs = plt.subplots(nrows=1, ncols=5, figsize=(50, 6))
for i, j, k in zip(np.arange(5), 
                   ["cell_type", "cell_class", "Tfe3", "Tfeb", "Mitf"],
                   ["cell type", "lineage", "Tfe3", "Tfeb", "Mitf"]):
    sc.pl.umap(adata,  
               neighbors_key="neighbors_scvi",
               color=j,
               layer="loge_cp10k1p", sort_order=True, frameon=False, edges=False, cmap="Reds", ncols=4,
               gene_symbols="gene_name",  vmin=0, vmax="p99", show=False, ax=axs[i])
    if k == "Tfe3":
      axs[i].set_title(label=k, fontdict={"fontsize": 16, "weight": "bold", "color":"#0077b6"})
    elif k == "Tfeb":
      axs[i].set_title(label=k, fontdict={"fontsize": 16, "weight": "bold", "color":"#7f7f7f"})
    elif k == "Mitf":
       axs[i].set_title(label=k, fontdict={"fontsize": 16, "weight": "bold", "color":"#ff6d00"})
    else:
      axs[i].set_title(label=k, fontdict={"fontsize": 16, "weight": "bold"})
      if k == "cell type":
        axs[i].legend(fontsize=14, bbox_to_anchor=(0.90, 1), loc="upper left", frameon=False, ncol=2, columnspacing=0.1)
      else:
        axs[i].legend(fontsize=14, bbox_to_anchor=(0.90, 1), loc="upper left", frameon=False)
fig.suptitle("Ximerakis et al., 2019", fontsize=24, fontdict={"style": "italic"}, x = 0.2, y = 1.05)
plt.savefig(out_dir_path+"umap_loge_cp10k1p_Tfe3_Tfeb_Mitf.pdf", bbox_inches="tight", pad_inches=0)

# Expression levels of TFEs on dotplot (lineage resolution)
# => Tfe3 and Tfeb only
dp = sc.pl.dotplot(
  adata=adata,
  var_names=["Tfe3", "Tfeb"],
  groupby="cell_class",
  layer="loge_cp10k1p",
  mean_only_expressed=True,
  swap_axes=True,
  colorbar_title="log\u2091(CP10K+1)",
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
    if label.get_text() == "Tfe3":
       label.set_color("#0077b6")
    else:
       label.set_color("#7f7f7f")   
ax.set_title(label="Ximerakis et al., 2019", fontdict={"fontsize": 10, "style": "italic"}, x = 0.1, y = 1.25)
plt.savefig(out_dir_path+"dotplot_loge_cp10k1p_Tfe3_Tfeb.pdf", bbox_inches="tight", pad_inches=0)

# => Tfe3, Tfeb and Mitf
dp = sc.pl.dotplot(
  adata=adata,
  var_names=["Tfe3", "Tfeb", "Mitf"],
  groupby="cell_class",
  layer="loge_cp10k1p",
  mean_only_expressed=True,
  swap_axes=True,
  colorbar_title="log\u2091(CP10K+1)",
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
    if label.get_text() == "Tfe3":
       label.set_color("#0077b6")
    elif label.get_text() == "Tfeb":
       label.set_color("#7f7f7f")  
    else:
       label.set_color("#ff6d00")      
ax.set_title(label="Ximerakis et al., 2019", fontdict={"fontsize": 10, "style": "italic"}, x = 0.1, y = 1.25)
plt.savefig(out_dir_path+"dotplot_loge_cp10k1p_Tfe3_Tfeb_Mitf.pdf", bbox_inches="tight", pad_inches=0)

# %%
