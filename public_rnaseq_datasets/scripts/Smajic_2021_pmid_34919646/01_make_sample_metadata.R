# Define paths ------------------------------------------------------------

# ==> To be set by user <==

# Path to working directory
# work_dir_path <- "path/to/working/directory/"

# Path to output directory
# out_dir_path <- "path/to/output/directory/" 

# Path to donor metadata file (csv)
# => Provided at data/Smajic_2021_pmid_34919646/donor_metadata.csv
# donor_metadata_path <- "path/to/donor_metadata.csv"

# Path to barcode metadata file (tsv)
# => Available at https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157783/suppl/GSE157783%5FIPDCO%5Fhg%5Fmidbrain%5Fcell.tar.gz
# barcode_metadata_path <- "path/to/salmon/IPDCO_hg_midbrain_cell.tsv"

# Our paths (for internal use only)
# work_dir_path <- "/scratch/ben/rnaseq/"
# out_dir_path <- "data/public/Smajic_2021_pmid_34919646/output/make_sample_metadata/" 
# donor_metadata_path <- "seq_data/public/Smajic_2021_pmid_34919646/supplementary_files/donor_metadata.csv"
# barcode_metadata_path <- "seq_data/public/Smajic_2021_pmid_34919646/supplementary_files/IPDCO_hg_midbrain_cell.tsv"

# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)

# Set working directory
setwd(work_dir_path)

# Construct sample metadata -----------------------------------------------

# Load donor metadata
donor_metadata <- read_csv(donor_metadata_path)

# Construct sample metadata
sample_metadata <- read_tsv(barcode_metadata_path) %>% 
  select(barcode, 
         cell_type = cell_ontology,
         donor_id = patient) %>% 
  mutate(donor_id = str_replace(donor_id, "C", "CTRL")) %>% 
  left_join(donor_metadata, by = join_by(donor_id))

write_csv(sample_metadata, paste0(out_dir_path, "sample_metadata.csv"))


