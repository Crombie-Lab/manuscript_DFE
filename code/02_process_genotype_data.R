library(tidyverse)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load parental genotypes
parent_genos <- data.table::fread("data/raw/parent_genos.csv")

# line order
line_order <- data.table::fread("data/raw/line_order.csv")

# get the unimputed genotypes for all lines
full_unimp <- data.table::fread("data/raw/M_RILS_unimp_s1.csv") %>% # this file is in alphanumeric order not genomic order, need to rearrange
  dplyr::select(full_id, all_of(line_order$line)) %>%
  dplyr::slice(-(1:3)) %>%
  tidyr::separate(full_id, into = c("chrom", "pos", "var_type"), remove = F) %>%
  dplyr::mutate(pos = as.numeric(pos)) %>%
  dplyr::arrange(chrom, pos) %>%
  dplyr::mutate(N2_G0_ANC = 0) %>%
  dplyr::left_join(., parent_genos, by = c("full_id" = "loci")) %>%
  dplyr::select(chrom, pos, var_type, full_var_id = full_id, 
                N2_G0_ANC, MA530, MA563, everything())

# get the set of 517 lines that were genotyped
genotyped_lines <- full_unimp %>%
  dplyr::select(dplyr::where(~!all(is.na(.x)))) %>%
  dplyr::select(any_of(line_order$line))
# get the names of lines
genotyped_lines2 <- tibble::tibble(line = names(genotyped_lines))

# get the imputation probabilities redacted lines (517 RI(AI)Ls)
p_matrix <- data.table::fread("data/processed/P-Matrix.csv")

# transpose p matrix
p_matrix2 <- as_tibble(t(p_matrix)) %>%
  janitor::row_to_names(row_number = 1)

# bind to var ids
full_p <- full_unimp %>%
  dplyr::select(1:4) %>%
  dplyr::bind_cols(p_matrix2)
  
# export the processed data
rio::export(full_unimp, file = "data/processed/02_DFE_unimputed_genotypes.csv")
rio::export(full_p, file = "data/processed/03_DFE_prob_genotypes.csv")
rio::export(genotyped_lines2, file = "data/processed/genotyped_lines.csv")

