# Define paths ------------------------------------------------------------

# ==> To be set by user <==

# Path to working directory
# work_dir_path <- "path/to/working/directory/"

# Path to output directory
# out_dir_path <- "path/to/output/directory/" 

# Path to sample metadata file (rds)
# => Produced by scripts/Nichterwitz_2016_pmid_27387371/01_make_sample_metadata.R
# sample_metadata_path <- "path/to/sample_metadata.rds"

# Path to gene metadata file (csv)
# => Produced by scripts/feature_metadata/01_make_feature_metadata.R
# gene_metadata_path <- "path/to/gene_metadata.csv"

# Path to TPM file (csv)
# => Produced by scripts/Nichterwitz_2016_pmid_27387371/02_summarise_to_gene_level.R
# tpm_path <- "path/to/tpm.csv"

# Our paths (for internal use only)
# work_dir_path <- "/scratch/ben/rnaseq/"
# out_dir_path <- "data/billy/Nichterwitz_2016_pmid_27387371/output/cnt_tfes/" 
# sample_metadata_path <- "data/public/Nichterwitz_2016_pmid_27387371/output/make_sample_metadata/sample_metadata.rds"
# gene_metadata_path <- "ref_data/feature_metadata/gene_metadata.csv"
# tpm_path <- "data/public/Nichterwitz_2016_pmid_27387371/output/summarise_to_gene_level/tpm.csv"

# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)
library(data.table)
library(ggpubr)

# Set working directory
setwd(work_dir_path)

# Import metadata --------------------------------------------------

# Construct human dopaminergic neuron-specific sample metadata
sample_metadata <- read_rds(sample_metadata_path) %>% 
  filter(!tissue %in% c("spinal motor", "motor neurons")) %>% 
  mutate(donor_id = str_extract(sample_description, "^Case\\d") %>% 
           str_replace("Case", "HC") %>% 
           str_replace("5", "1") %>% 
           str_replace("6", "2") %>% 
           str_replace("7", "3"),
         region = case_when(str_detect(sample_description, "_SNc_") ~ "SNpc",
                            str_detect(sample_description, "_VTA_") ~ "VTA",
                            TRUE ~ NA_character_) %>% 
           fct_relevel(c("SNpc", "VTA")),
         lcm_seq_sample_id = case_when(cell_number == 120 ~ str_extract(sample_description, "LCM\\d+$"),
                                       cell_number == 75 ~ str_extract(sample_description, "^Case.*\\dmin"),
                                       TRUE ~ NA_character_),
         batch = paste0(batch, " et al., 2016"),
         donor_id_cell_number = paste0(donor_id, " - ", cell_number, "c"))

# Import gene metadata 
gene_metadata <- read_csv(gene_metadata_path)

# Load TPMs ---------------------------------------------------------------

tpm <- fread(tpm_path) %>% 
  select(ensembl_gene_id_version, all_of(sample_metadata$sample_id))

# Construct count data ----------------------------------------------------

count_data <- tpm %>% 
  pivot_longer(cols = -ensembl_gene_id_version, names_to = "sample_id", values_to = "tpm") %>% 
  mutate(log2_tpm1p = log(x = tpm+1L, base = 2)) %>% 
  left_join(sample_metadata, by = join_by(sample_id)) %>% 
  left_join(gene_metadata, by = join_by(ensembl_gene_id_version)) 

# Plot TFEs counts --------------------------------------------------------

# Perform Wilcoxon signed-rank test to compare the expression level of TFE3 and TFEB
wilcox_df <- count_data %>%
  filter(gene_name %in% c("TFE3", "TFEB")) %>%
  select(gene_name, donor_id_cell_number, region, lcm_seq_sample_id, log2_tpm1p) %>%
  pivot_wider(names_from = "gene_name", values_from = "log2_tpm1p")

wilcox_res <- wilcox.test(x = wilcox_df$TFEB,
                          y = wilcox_df$TFE3,
                          alternative = "two.sided",
                          paired = TRUE)

# Display log2 TPM+1 for TFE3 and TFEB in LCM dopaminergic neurons (SNpc / VTA)
ggplot(data = count_data %>% 
         filter(gene_name %in% c("TFE3", "TFEB")) %>% 
         mutate(gene_name = fct_relevel(gene_name, c("TFE3", "TFEB"))),
       mapping = aes(x = gene_name, y = log2_tpm1p)) +
  geom_boxplot(outliers = F, fill = "#e3d5ca", alpha = 0.6) +
  geom_line(mapping = aes(group = lcm_seq_sample_id),
            position = position_dodge(width = 0.4),
            color = "grey", linewidth = 0.5, alpha = 0.5) +
  geom_point(mapping = aes(fill = donor_id_cell_number, shape = region, group = lcm_seq_sample_id), 
             position = position_dodge(width = 0.4),
             size = 4, alpha = 0.6) +
  geom_label(data = tibble(gene_name = "TFEB",
                           log2_tpm1p = 1.1*max(count_data[count_data$gene_name == "TFEB", "log2_tpm1p"][[1]],
                                                na.rm = T),
                           label = paste0("p = ", formatC(wilcox_res$p.value, format = "e", digits = 1))),
             mapping = aes(label = label),
             size = 5) +
  facet_wrap(~batch) +
  scale_fill_manual(values = c("#bbdefb", "#0d47a1", "#e30b5d", "#758bfd", "#fb6107","#74c69d")) +
  scale_shape_manual(values = c(22, 21)) +
  theme_pubr(legend = "right") +
  theme(strip.text = element_text(size = 16, face="italic"),
        axis.text.x = element_text(size = 16, face="bold", color = c("#0077b6", "#7f7f7f")),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)) +
  labs(x = NULL, 
       y = paste0("log", utf8::utf8_print("\u2082"), "(TPM+1)"), 
       fill = "donor ID - n cells", 
       shape = "region") +
  guides(fill = guide_legend(override.aes = list(shape = 22, size = 5), order = 1L),
         shape = guide_legend(override.aes = list(size = 5), order = 2L)) 

ggsave(plot = last_plot(),
       filename = paste0(out_dir_path, "TFE3_TFEB_log2_tpm1p.pdf"),
       width = 10,
       height = 8,
       device = cairo_pdf)
