# Define paths ------------------------------------------------------------

# ==> To be set by user <==

# Path to working directory
# work_dir_path <- "path/to/working/directory/"

# Path to output directory
# out_dir_path <- "path/to/output/directory/" 

# Path to SraRunTable file (txt)
# => Provided at data/Ximerakis_2019_pmid_31551601/SraRunTable.txt
# sra_run_table_path <- "path/to/SraRunTable.txt"

# Path to animal metadata file (csv)
# => Provided at data/Ximerakis_2019_pmid_31551601/sample_id_animal_id.csv
# animal_metadata_path <- "path/to/sample_id_animal_id.csv"

# Path to barcode metadata file
# => Available at https://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129788/suppl/GSE129788%5FSupplementary%5Fmeta%5Fdata%5FCell%5FTypes%5FEtc.txt.gz
# barcode_metadata_path <- "path/to/GSE129788_Supplementary_meta_data_Cell_Types_Etc.txt.gz"

# Our paths (for internal use only)
# work_dir_path <- "/scratch/ben/rnaseq/"
# out_dir_path <- "data/public/Ximerakis_2019_pmid_31551601/output/make_sample_metadata/"
# sra_run_table_path <- "seq_data/public/Ximerakis_2019_pmid_31551601/SraRunTable.txt"
# animal_metadata_path <- "seq_data/public/Ximerakis_2019_pmid_31551601/supplementary_files/sample_id_animal_id.csv"
# barcode_metadata_path <- "seq_data/public/Ximerakis_2019_pmid_31551601/supplementary_files/GSE129788_Supplementary_meta_data_Cell_Types_Etc.txt.gz"

# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)

# Set working directory
setwd(work_dir_path)

# Construct sample metadata -----------------------------------------------

sample_metadata <- read_delim(barcode_metadata_path,
                              delim = "\t", show_col_types = F) %>%
  filter(NAME != "TYPE") %>% 
  mutate(animal_num_barcode = str_remove(NAME, "^Aging_mouse_brain_portal_data_"),
         barcode = str_remove(NAME, "^Aging_mouse_brain_portal_data_\\d{1,}_")) %>% 
  select(animal_num_barcode, barcode, 
         nGene, nUMI, 
         cell_type = cluster, cell_class = cell_classes, 
         age_status = animal_type, cell_type_age_status = cell_type_age)

write_csv(sample_metadata, paste0(out_dir_path, "sample_metadata.csv"))

# Construct tissue metadata -----------------------------------------------

# Load list of sample IDs - donor IDs
sample_id_donor_id <- read_csv(animal_metadata_path) %>% 
  mutate(animal_id = str_remove(animal_id_description, ":.*$"))

tissue_metadata <- read_delim(file = sra_run_table_path,
                              delim = ",") %>%
  select(run_id = Run,
         sample_id = `Sample Name`,
         age = AGE,
         LibraryLayout,
         AvgSpotLen,
         Instrument) %>%
  mutate(age_status = ifelse(str_detect(age, "21-22"), "old", "young"), .after = age) %>% 
  left_join(sample_id_donor_id, by = join_by(sample_id))

write_csv(tissue_metadata, paste0(out_dir_path, "tissue_metadata.csv"))
