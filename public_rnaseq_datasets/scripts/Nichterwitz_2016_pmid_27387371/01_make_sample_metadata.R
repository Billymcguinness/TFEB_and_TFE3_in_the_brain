# Define paths ------------------------------------------------------------

# ==> To be set by user <==

# Path to working directory
# work_dir_path <- "path/to/working/directory/"

# Path to output directory
# out_dir_path <- "path/to/output/directory/" 

# Path to sample ID description file (csv)
# => Provided at data/Nichterwitz_2016_pmid_27387371/sample_id_description.csv
# sample_id_description_path <- "path/to/sample_id_description.csv"

# Path to SraRunTable file (txt)
# => Provided at data/Nichterwitz_2016_pmid_27387371/SraRunTable.txt
# sra_run_table_path <- "path/to/SraRunTable.txt"

# Path to directory containing Salmon quant files
# => Available upon request or see methods to produce from FASTQ files
# salmon_quant_dir_path <- "path/to/salmon/quant/directory/"

# Our paths (for internal use only)
# work_dir_path <- "/scratch/ben/rnaseq/"
# out_dir_path <- "data/public/Nichterwitz_2016_pmid_27387371/output/make_sample_metadata/" 
# sample_id_description_path <- "seq_data/public/Nichterwitz_2016_pmid_27387371/supplementary_files/sample_id_description.csv"
# sra_run_table_path <- "seq_data/public/Nichterwitz_2016_pmid_27387371/SraRunTable.txt"
# salmon_quant_dir_path <- "seq_data/public/Nichterwitz_2016_pmid_27387371/03_salmon_selective_alignment_trim_galore/"

# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)

# Set working directory
setwd(work_dir_path)

# Load list of sample ID - sample description -----------------------------

sample_id_sample_description <- read_csv(sample_id_description_path)

# Construct sample metadata -----------------------------------------------

sample_metadata <- read_delim(file = sra_run_table_path,
                              delim = ",") %>% 
  filter(Organism == "Homo sapiens") %>%
  select(run_id = Run,
         sample_id = `Sample Name`,
         cell_number = Cell_number,
         source_name,
         tissue,
         LibraryLayout,
         AvgSpotLen,
         Instrument) %>%
  mutate(batch = "Nichterwitz",
         species = "human", 
         disease_status = "CTRL",
         cell_number = as_factor(cell_number),
         files = setNames(
           object = map2_chr(.x = species,
                             .y = run_id,
                             .f = ~ list.files(path = paste0(salmon_quant_dir_path, .y),
                                               pattern = "quant.sf",
                                               full.names = TRUE)), 
           nm = sample_id)
  ) %>% 
  left_join(sample_id_sample_description, by = join_by(sample_id)) %>% 
  select(files, batch, run_id:species, disease_status, cell_number:tissue, sample_description, LibraryLayout:Instrument) 

write_rds(sample_metadata, paste0(out_dir_path, "sample_metadata.rds"))
