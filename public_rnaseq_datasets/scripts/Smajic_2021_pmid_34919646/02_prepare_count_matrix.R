# Define paths ------------------------------------------------------------

# ==> To be set by user <==

# Path to working directory
# work_dir_path <- "path/to/working/directory/"

# Path to output directory
# out_dir_path <- "path/to/output/directory/" 

# Path to genes ID file (tsv)
# => Available at https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157783/suppl/GSE157783%5FIPDCO%5Fhg%5Fmidbrain%5Fgenes.tar.gz 
# gene_id_path <- "path/to/IPDCO_hg_midbrain_genes.tsv"

# Path to count matrix file (tsv)
# => Available at https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157783/suppl/GSE157783%5FIPDCO%5Fhg%5Fmidbrain%5FUMI.tar.gz
# cnt_path <- "path/to/IPDCO_hg_midbrain_UMI.tsv"

# Our paths
work_dir_path <- "/scratch/ben/rnaseq/"
out_dir_path <- "data/public/Smajic_2021_pmid_34919646/output/prepare_count_matrix/"
gene_id_path <- "seq_data/public/Smajic_2021_pmid_34919646/supplementary_files/IPDCO_hg_midbrain_genes.tsv"
cnt_path <- "seq_data/public/Smajic_2021_pmid_34919646/supplementary_files/IPDCO_hg_midbrain_UMI.tsv"

# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)
library(data.table)
library(Matrix)
library(DropletUtils)

# Set working directory
setwd(work_dir_path)

# Construct count matrix --------------------------------------------------

# Load raw UMI count matrix
raw_count <- as.matrix(fread(cnt_path))

# Set row names to Ensembl gene IDs 
rownames(raw_count) <- read_tsv(gene_id_path) %>% 
  arrange(row) %>% 
  pull(gene)

# Convert to sparse matrix
raw_count <- Matrix(raw_count, sparse = TRUE)

# Write count matrix to file
write10xCounts(path = paste0(out_dir_path, "raw_count.h5"),
               x = raw_count,
               barcodes = colnames(raw_count),
               gene.id = rownames(raw_count),
               overwrite = T)

