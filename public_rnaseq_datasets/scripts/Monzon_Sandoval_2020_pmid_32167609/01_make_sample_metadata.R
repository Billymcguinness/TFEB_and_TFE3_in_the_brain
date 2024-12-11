# Define paths ------------------------------------------------------------

# ==> To be set by user <==

# Path to working directory
# work_dir_path <- "path/to/working/directory/"

# Path to output directory
# out_dir_path <- "path/to/output/directory/" 

# Path to SraRunTable file (txt)
# => Provided at data/Monzon_Sandoval_2020_pmid_32167609/SraRunTable.txt
# sra_run_table_path <- "path/to/SraRunTable.txt"

# Path to directory containing Salmon quant files
# => Available upon request or see methods to produce from FASTQ files
# salmon_quant_dir_path <- "path/to/salmon/quant/directory/"

# Our paths
work_dir_path <- "/scratch/ben/rnaseq/"
out_dir_path <- "data/public/Monzon_Sandoval_2020_pmid_32167609/output/make_sample_metadata/" 
sra_run_table_path <- "seq_data/public/Monzon_Sandoval_2020_pmid_32167609/SraRunTable.txt"
salmon_quant_dir_path <- "seq_data/public/Monzon_Sandoval_2020_pmid_32167609/03_salmon_selective_alignment_trim_galore/"

# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)

# Set working directory
setwd(work_dir_path)

# Construct sample metadata -----------------------------------------------

sample_metadata <- read_delim(file = sra_run_table_path,
                              delim = ",") %>%
  select(run_id = Run,
         sample_id = `Sample Name`,
         donor_id = Origin,
         donor_sex = sex,
         donor_age = AGE,
         tier,
         rin = RIN,
         LibraryLayout,
         AvgSpotLen,
         Instrument) %>%
  mutate(batch = "Monzon-Sandoval",
         species = "human",
         disease_status = "CTRL",
         donor_sex = ifelse(donor_sex == "female", "F", "M"),
         tissue = "SNpc",
         precise_tissue = paste0(tier, " tier ", tissue),
         tier = fct_relevel(tier, c("ventral", "dorsal")),
         sample_name = paste0(donor_id, " ", tier, " ", tissue, " DaNs"),
         files = setNames(
           object = map_chr(.x = run_id,
                            .f = ~ list.files(path = paste0(salmon_quant_dir_path, .x),
                                              pattern = "quant.sf",
                                              full.names = TRUE)), 
           nm = sample_id)
  ) %>% 
  select(files, batch, run_id:donor_id, sample_name, 
         species, disease_status, donor_sex, donor_age, 
         tissue, tier, precise_tissue, rin, LibraryLayout:Instrument) 

write_rds(sample_metadata, paste0(out_dir_path, "sample_metadata.rds"))
