# Define paths ------------------------------------------------------------

# ==> To be set by user <==

# Path to working directory
# work_dir_path <- "path/to/working/directory/"

# Path to output directory
# out_dir_path <- "path/to/output/directory/" 

# Path to sample metadata file (rds)
# => Produced by scripts/foundinpd/02_make_sample_metadata.R
# sample_metadata_path <- "path/to/filtered_sample_metadata.rds"

# Path to gene metadata file (csv)
# => Produced by scripts/foundinpd/01_make_feature_metadata.R
# gene_metadata_path <- "path/to/gene_metadata.csv"

# Path to TPM file (csv)
# => Produced by scripts/foundinpd/03_summarise_to_gene_level.R
# tpm_path <- "path/to/tpm.csv"

# Our paths (for internal use only)
# work_dir_path <- "/scratch/ben/rnaseq/"
# out_dir_path <- "data/billy/foundinpd/output/cnt_tfes/" 
# sample_metadata_path <- "data/public/foundinpd/output/make_sample_metadata/filtered_sample_metadata.rds"
# gene_metadata_path <- "data/public/foundinpd/output/feature_metadata/gene_metadata.csv"
# tpm_path <- "data/public/foundinpd/output/summarise_to_gene_level/tpm.csv"

# Set up ------------------------------------------------------------------

# Import required libraries
library(tidyverse)
library(data.table)
library(ggpubr)

# Set working directory
setwd(work_dir_path)

# Import metadata --------------------------------------------------

# Import filtered sample metadata
sample_metadata <- read_rds(sample_metadata_path)

# Import gene metadata 
gene_metadata <- read_csv(gene_metadata_path)

# Load TPMs ---------------------------------------------------------------

tpm <- fread(tpm_path)

# Construct count data ----------------------------------------------------

count_data <- tpm %>%
  select(ensembl_gene_id_version, all_of(sample_metadata$sample_spe_dir)) %>% 
  filter(ensembl_gene_id_version %in% gene_metadata[gene_metadata$gene_name %in% c("TFE3", "TFEB"), 
                                                    "ensembl_gene_id_version"][[1]]) %>% 
  pivot_longer(cols = -ensembl_gene_id_version, names_to = "sample_spe_dir", values_to = "tpm") %>% 
  left_join(sample_metadata, by = join_by(sample_spe_dir)) %>% 
  left_join(gene_metadata, by = join_by(ensembl_gene_id_version)) %>% 
  mutate(log2_tpm1p = log(x = tpm+1L, base = 2),
         div_gene_name = paste0(div, "\n", gene_name) %>% 
           fct_relevel(c("DIV0\nTFE3", "DIV0\nTFEB", "DIV25\nTFE3", "DIV25\nTFEB", "DIV65\nTFE3", "DIV65\nTFEB")))

# Plot TFEs counts --------------------------------------------------------

# Perform Wilcoxon signed-rank test to compare the expression level of TFE3 and TFEB at each time point
wilcox_dfs <- count_data %>%
  select(gene_name, patno, div, log2_tpm1p) %>%
  pivot_wider(names_from = "gene_name", values_from = "log2_tpm1p") %>% 
  split(f = .$div)

wilcox_res <- map(.x = wilcox_dfs,
                  .f = function(x) { 
                    
                    wilcox.test(x = x$TFEB,
                                y = x$TFE3,
                                alternative = "two.sided",
                                paired = TRUE)
                  })

# Extract p values and correct for multiple testing
ori_pvals <- map_dbl(.x = wilcox_res, .f = ~ .x$p.value)
adj_pvals <- p.adjust(p = ori_pvals, method = "bonferroni")

# Display log2 TPM+1 for TFE3 and TFEB in DIV0, DIV25 and DIV65 iPSC-DaNs  
# => Show original p values on plot
ggplot(data = count_data %>% 
         filter(disease_status %in% c("CTRL", "PD")),
       mapping = aes(x = gene_name, 
                     y = log2_tpm1p)) +
  geom_boxplot(outliers = F, fill = "#e3d5ca", alpha = 0.6) +
  geom_line(mapping = aes(group = interaction(patno, div)),
            position = position_dodge(width = 0.6),
            color = "grey", linewidth = 0.5, alpha = 0.5) +
  geom_point(mapping = aes(fill = disease_status, 
                           shape = disease_status, 
                           group = interaction(patno, div)), 
             position = position_dodge(width = 0.6),
             size = 4, alpha = 0.6) +
  geom_label(data = tibble(gene_name = rep(x = "TFEB", times = 3),
                           div = names(ori_pvals),
                           log2_tpm1p = map_dbl(.x = names(ori_pvals),
                                                .f = ~ 1.2*max(count_data[count_data$gene_name == "TFEB" & count_data$div == .x, 
                                                                          "log2_tpm1p"][[1]],
                                                               na.rm = T)),
                           label = map_chr(.x = ori_pvals,
                                           .f = ~ paste0("p = ", formatC(.x, format = "e", digits = 1)))),
             mapping = aes(label = label),
             size = 5) +
  facet_wrap(~div) +
  scale_fill_manual(values = c("#5d76cb", "#ca3767")) +
  scale_shape_manual(values = c(22, 21)) +
  theme_pubr(legend = "right") +
  theme(strip.text = element_text(size = 16),
        axis.text.x = element_text(size = 16, face="bold", color = c("#0077b6", "#7f7f7f")),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)) +
  labs(x = NULL, 
       y = paste0("log", utf8::utf8_print("\u2082"), "(TPM+1)"), 
       fill = NULL,
       shape = NULL)

ggsave(plot = last_plot(),
       filename = paste0(out_dir_path, "TFE3_TFEB_log2_tpm1p_oripval.pdf"),
       width = 10,
       height = 8,
       device = cairo_pdf)

# => Show Bonferroni-adjusted p values on plot
ggplot(data = count_data %>% 
         filter(disease_status %in% c("CTRL", "PD")),
       mapping = aes(x = gene_name, 
                     y = log2_tpm1p)) +
  geom_boxplot(outliers = F, fill = "#e3d5ca", alpha = 0.6) +
  geom_line(mapping = aes(group = interaction(patno, div)),
            position = position_dodge(width = 0.6),
            color = "grey", linewidth = 0.5, alpha = 0.5) +
  geom_point(mapping = aes(fill = disease_status, 
                           shape = disease_status, 
                           group = interaction(patno, div)), 
             position = position_dodge(width = 0.6),
             size = 4, alpha = 0.6) +
  geom_label(data = tibble(gene_name = rep(x = "TFEB", times = 3),
                           div = names(adj_pvals),
                           log2_tpm1p = map_dbl(.x = names(adj_pvals),
                                                .f = ~ 1.2*max(count_data[count_data$gene_name == "TFEB" & count_data$div == .x, 
                                                                          "log2_tpm1p"][[1]],
                                                               na.rm = T)),
                           label = map_chr(.x = adj_pvals,
                                           .f = ~ paste0("adj. p = ", formatC(.x, format = "e", digits = 1)))),
             mapping = aes(label = label),
             size = 4.6) +
  facet_wrap(~div) +
  scale_fill_manual(values = c("#5d76cb", "#ca3767")) +
  scale_shape_manual(values = c(22, 21)) +
  theme_pubr(legend = "right") +
  theme(strip.text = element_text(size = 16),
        axis.text.x = element_text(size = 16, face="bold", color = c("#0077b6", "#7f7f7f")),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)) +
  labs(x = NULL, 
       y = paste0("log", utf8::utf8_print("\u2082"), "(TPM+1)"), 
       fill = NULL,
       shape = NULL)

ggsave(plot = last_plot(),
       filename = paste0(out_dir_path, "TFE3_TFEB_log2_tpm1p_adjpval.pdf"),
       width = 10,
       height = 8,
       device = cairo_pdf)




