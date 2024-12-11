# Define paths ------------------------------------------------------------

# ==> To be set by user <==

# Path to working directory
# work_dir_path <- "path/to/working/directory/"

# Path to output directory
# out_dir_path <- "path/to/output/directory/" 

# Path to sample metadata file (rds)
# => Produced by scripts/foundinpd/02_make_sample_metadata.R
# sample_metadata_path <- "path/to/filtered_sample_metadata.rds"

# Path to transcript metadata file (csv)
# => Produced by scripts/foundinpd/01_make_feature_metadata.R
# transcript_metadata_path <- "path/to/transcript_metadata.csv"

# Our paths (for internal use only)
# work_dir_path <- "/scratch/ben/rnaseq/"
# out_dir_path <- "data/public/foundinpd/output/summarise_to_gene_level/" 
# sample_metadata_path <- "data/public/foundinpd/output/make_sample_metadata/filtered_sample_metadata.rds"
# transcript_metadata_path <- "data/public/foundinpd/output/feature_metadata/transcript_metadata.csv"

# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)
library(tximport)
library(data.table)

# Set working directory
setwd(work_dir_path)

# Import metadata --------------------------------------------------

# Import filtered sample metadata
sample_metadata <- read_rds(sample_metadata_path) 

# Import transcript metadata
transcript_metadata <- read_csv(transcript_metadata_path)

# Summarise transcript data to gene level ---------------------------------

# Construct txi object
txi <- tximport(files = sample_metadata$files,
                type = "salmon",
                txOut = FALSE,
                countsFromAbundance = "no",
                tx2gene = transcript_metadata)

write_rds(x = txi,
          file = paste0(out_dir_path, "txi.rds"))

# Extract estimated TPMs --------------------------------------------------

txi <- read_rds(file = paste0(out_dir_path, "txi.rds"))

tpm <- as_tibble(txi$abundance, rownames = "ensembl_gene_id_version")

fwrite(x = tpm, file = paste0(out_dir_path, "tpm.csv"))

# Extract estimated raw read counts ---------------------------------------

raw_count <- as_tibble(txi$counts, rownames = "ensembl_gene_id_version")

fwrite(x = raw_count, file = paste0(out_dir_path, "raw_count.csv"))

# Extract transcript length -----------------------------------------------

transcript_length <- as_tibble(txi$length, rownames = "ensembl_gene_id_version")

fwrite(x = transcript_length,
       file = paste0(out_dir_path, "transcript_length.csv"))



