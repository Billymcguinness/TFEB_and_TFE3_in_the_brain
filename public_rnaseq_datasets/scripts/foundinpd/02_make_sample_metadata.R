# Define paths ------------------------------------------------------------

# ==> To be set by user <==

# Path to working directory
# work_dir_path <- "path/to/working/directory/"

# Path to output directory
# out_dir_path <- "path/to/output/directory/" 

# Path to participant metadata file
# => Available from the Parkinson’s Progression Markers Initiative and provided at data/foundinpd/Participant_Status_17Sep2024.csv
# participant_metadata_path <- "path/to/Participant_Status_17Sep2024.csv"

# Path to iPSC metadata file
# => Available from the Parkinson’s Progression Markers Initiative and provided at data/foundinpd/iPSC_Catalog_Metadata_17Sep2024.csv
# ipsc_metadata_path <- "path/to/iPSC_Catalog_Metadata_17Sep2024.csv"

# Path to directory containing Salmon quant files
# => Available from the Parkinson’s Progression Markers Initiative
# salmon_quant_dir_path <- "path/to/salmon/quant/directory/"

# Our paths
work_dir_path <- "/scratch/ben/rnaseq/"
out_dir_path <- "data/public/foundinpd/output/make_sample_metadata/"
participant_metadata_path <- "data/public/foundinpd/input/subject_characteristics/Participant_Status_17Sep2024.csv"
ipsc_metadata_path <- "data/public/foundinpd/input/biosample_inventory/iPSC_Catalog_Metadata_17Sep2024.csv"
salmon_quant_dir_path <- "seq_data/public/foundinpd/processed/RNAB/salmon_quant/"

# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)

# Set working directory
setwd(work_dir_path)

# Construct participant metadata ------------------------------------------

participant_metadata <- read_csv(participant_metadata_path) %>% 
  select(patno = PATNO,
         disease_status = COHORT_DEFINITION,
         starts_with("ENRL")) %>% 
  mutate(patno = as.character(patno),
         disease_status = case_when(disease_status == "Healthy Control" ~ "CTRL",
                                    disease_status == "Parkinson's Disease" ~ "PD",
                                    disease_status %in% c("Prodromal", "SWEDD") ~ disease_status,
                                    TRUE ~ NA_character_) %>% 
           fct_relevel(c("CTRL", "Prodromal", "PD", "SWEDD"))) %>% 
  pivot_longer(cols = starts_with("ENRL"), names_to = "risk_name", values_to = "risk_value") %>% 
  mutate(risk_factor = map2_chr(.x = risk_value,
                                .y = risk_name,
                                .f = ~ ifelse(.x == 1, 
                                              str_remove(.y, "ENRL"), 
                                              NA_character_))) %>% 
  select(-c(risk_name, risk_value)) %>% 
  nest(data = risk_factor) %>% 
  mutate(risk_factor = map_chr(.x = data,
                               .f = ~ paste0(na.omit(.x$risk_factor), collapse = " / ")) %>% 
           {ifelse(. == "", NA_character_, .)}) %>% 
  select(-data)

# Construct iPSC metadata -------------------------------------------------

ipsc_metdata <- read_csv(ipsc_metadata_path) %>% 
  select(patno = PATNO,
         disease_status = COHORT,
         gender = GENDER,
         age_at_diagnosis = AGE_AT_DIAGNOSIS,
         age_at_sample_collection = AGE_AT_SAMPLE_COLLECTION,
         karyotype = KARYOTYPE) %>% 
  mutate(patno = as.character(patno),
         disease_status = case_when(disease_status == "Healthy Control" ~ "CTRL",
                                    disease_status == "Parkinson's Disease" ~ "PD",
                                    disease_status %in% c("Prodromal", "SWEDD") ~ disease_status,
                                    TRUE ~ NA_character_) %>% 
           fct_relevel(c("CTRL", "Prodromal", "PD", "SWEDD"))) %>% 
  unique() %>% 
  filter(karyotype != "Abnormal")

# Construct sample metadata -----------------------------------------------

sample_metadata <- tibble(sample_spe_dir = list.dirs(path = salmon_quant_dir_path, 
                                                     full.names = F, recursive = F)) %>% 
  mutate(files = setNames(
    object = map_chr(.x = sample_spe_dir,
                     .f = ~ list.files(path = paste0(salmon_quant_dir_path, .x),
                                       pattern = "quant.sf",
                                       full.names = TRUE)), 
    nm = sample_spe_dir),
    split_name = str_split(string = sample_spe_dir, pattern = "_"),
    patno = map_chr(.x = split_name, .f = ~ .x[[2]]) %>% 
      str_remove("^PPMI") %>% str_remove("B.*$"),
    cdi_number = map_chr(.x = split_name, .f = ~ .x[[3]]),
    patno_cdi_number = paste0(patno, "_", cdi_number),
    batch = map_chr(.x = split_name, .f = ~ .x[[2]]) %>% 
      str_extract("B.*$") %>% {ifelse(is.na(.), "B1", .)},
    div = map_chr(.x = split_name, .f = ~ .x[[4]]) %>% 
      str_replace("^[a-z]+", "DIV"),
    version = map_chr(.x = split_name, .f = ~ .x[[5]])) %>% 
  select(-split_name) %>% 
  left_join(ipsc_metdata, by = join_by(patno)) %>% 
  left_join(participant_metadata, by = join_by(patno, disease_status)) 

write_rds(sample_metadata, paste0(out_dir_path, "sample_metadata.rds"))

# Subset sample metadata 
# => For each combination of patno and DIV, keep only the sample-specific directories:
# - matching the "da" differentiation protocol
# - the sample with smallest version index among samples with smallest batch index
# => Keep only the samples with disease status CTRL, Prodromal or PD
filtered_files <- sample_metadata %>% 
  filter(str_detect(files, "_da\\d{1,2}_v\\d{1,}/quant.sf$")) %>%
  select(patno, batch, div, version, files) %>% 
  unique() %>% 
  arrange(patno, batch, div, version) %>% 
  nest(data = c(batch, version, files)) %>% 
  mutate(data = map(data, ~ .x[1,])) %>% 
  unnest(data) %>% 
  pull(files)

filtered_sample_metadata <- sample_metadata %>% 
  filter(files %in% filtered_files,
         disease_status %in% c("CTRL", "Prodromal", "PD"))

write_rds(filtered_sample_metadata, paste0(out_dir_path, "filtered_sample_metadata.rds"))