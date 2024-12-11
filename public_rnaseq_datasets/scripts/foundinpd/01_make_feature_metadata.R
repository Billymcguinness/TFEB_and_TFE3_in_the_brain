# Define paths ------------------------------------------------------------

# ==> To be set by user <==

# Path to working directory
# work_dir_path <- "path/to/working/directory/"

# Path to output directory
# out_dir_path <- "path/to/output/directory/" 

# Path to GENCODE annotation file (gtf)
# => Can be produced as described at https://github.com/FOUNDINPD/annotation-RNA
# gencode_annotation_path <- "path/to/gencode_annotation.gtf"

# Our paths (for internal use only)
# work_dir_path <- "/scratch/ben/rnaseq/"
# out_dir_path <- "data/public/foundinpd/output/feature_metadata/"
# gencode_annotation_path <- "data/public/foundinpd/input/generate_annotation/gencode_v29.lncipedia_v5_2_hc.annotation.gtf"

# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)

# Set working directory
setwd(work_dir_path)

# Load GENCODE annotation -------------------------------------------------

gencode_annotation <- as_tibble(rtracklayer::import(gencode_annotation_path))

write_rds(gencode_annotation, paste0(out_dir_path, "gencode_annotation.rds"))

# Construct gene metadata from GENCODE annotation -------------------------

gene_metadata <- unique(gencode_annotation[!is.na(gencode_annotation$transcript_id), 
                                           c("seqnames", "gene_id", "gene_type", "gene_name")]) %>% 
  rename(chr_name = seqnames,
         ensembl_gene_id_version = gene_id) %>% 
  mutate(ensembl_gene_id = str_remove(ensembl_gene_id_version, "\\..*$")) %>% 
  select(chr_name, ensembl_gene_id_version, ensembl_gene_id, gene_type, gene_name)

write_csv(gene_metadata, paste0(out_dir_path, "gene_metadata.csv"))

# Construct transcript metadata from GENCODE annotation -------------------

# Construct transcript metadata
transcript_metadata <- unique(gencode_annotation[!is.na(gencode_annotation$transcript_id), 
                                                 c("transcript_id", "gene_id")]) %>% 
  rename(ensembl_transcript_id_version = transcript_id,
         ensembl_gene_id_version = gene_id) %>% 
  arrange(ensembl_gene_id_version) %>% 
  group_by(ensembl_gene_id_version) %>% 
  mutate(n_transcripts_for_gene = length(ensembl_transcript_id_version)) %>% 
  ungroup()

write_csv(transcript_metadata, paste0(out_dir_path, "transcript_metadata.csv"))


