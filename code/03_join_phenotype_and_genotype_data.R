library(tidyverse)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

#=====================================================#
# Part 1: load data
#=====================================================#
# load phenotype data
pheno <- data.table::fread("data/processed/01_DFE_phenotypes.csv")

# load unimputed genotype data and transpose to join
geno_unimp <-  data.table::fread("data/processed/02_DFE_unimputed_genotypes.csv") %>% # read data
  dplyr::select(-(1:3)) %>% # drop unwanted vars
  data.table::transpose(., keep.names = "temp") %>% # transpose the data with temp full_id
  janitor::row_to_names(row_number = 1) %>% # push first row to header
  dplyr::select(full_id = full_var_id, everything()) # rename to full_id

# load fully imputed genotype data and transpose to join
geno_imp_p <-  data.table::fread("data/processed/03_DFE_prob_genotypes.csv") %>% # read data
  dplyr::select(-(1:3)) %>% # drop unwanted vars
  data.table::transpose(., keep.names = "temp") %>% # transpose the data with temp full_id
  janitor::row_to_names(row_number = 1) %>% # push first row to header
  dplyr::select(full_id = full_var_id, everything()) # rename to full_id

#=====================================================#
# Part 2: Join data
#=====================================================#
# join the unimptuted data
join_unimp <- pheno %>%
  dplyr::mutate(join = case_when(stringr::str_detect(full_id, pattern = "N2_G0_ANC|MA530|MA563") ~ stringr::str_replace(full_id, pattern = "(.[^_]+$)", replacement = ""),
                                 TRUE ~ full_id)) %>% # fix full_id to join genotypes
  dplyr::left_join(geno_unimp, by = c("join" = "full_id")) %>%
  dplyr::select(-join)

# join the imputed data with p matrix
join_imp_p <- pheno %>%
  dplyr::mutate(join = case_when(stringr::str_detect(full_id, pattern = "N2_G0_ANC|MA530|MA563") ~ stringr::str_replace(full_id, pattern = "(.[^_]+$)", replacement = ""),
                                 TRUE ~ full_id)) %>% # fix full_id to join genotypes
  dplyr::left_join(geno_imp_p, by = c("join" = "full_id")) %>%
  dplyr::select(-join)

#=====================================================#
# Part 3: Export joined data
#=====================================================#
# export data with all genotypes as .csv files for manuscript
rio::export(join_unimp, file = "data/processed/04_DFE_joined_unimputed_geno_pheno.csv")
rio::export(join_imp_p, file = "data/processed/05_DFE_joined_imputed_prob_geno_pheno.csv")

# #=====================================================#
# # Part 4: Explore low representation loci - lots of missing genotypes
# #=====================================================#
# # get distribution of loci genotype representation in RILs
# lr_loci <- join_imp %>%
#   dplyr::filter(grepl(full_id, pattern = "RIL_|RIAIL_")) %>%
#   dplyr::filter(rowSums(is.na(.)) != 169) %>% # remove lines with all missing genotypes
#   dplyr::select(-p_focal, -p_comp) %>%
#   dplyr::distinct(full_id, .keep_all = T) %>%
#   dplyr::select(-1:-4) %>%
#   summarise(across(everything(), ~ sum(!is.na(.)))) %>%
#   tidyr::pivot_longer(cols = everything())
# 
# # plot it
# lr_loci_hist <- ggplot(lr_loci) +
#   aes(x = value) +
#   geom_histogram(bins = 60) +
#   xlim(0, 520) +
#   theme_bw() +
#   geom_vline(xintercept = 120, color = "red", linetype = 2) +
#   annotate("text", x = 90, y = 10, label = "n = 120", color = "red") +
#   labs(x = "n lines with genotypes at locus", y = "n loci")
# 
# # save it
# cowplot::ggsave2(lr_loci_hist, filename = "plots/loci_genotype_representation_histogram.png", width = 5, height = 5)
# 
# # Get loci with too few genotyped lines n < 120
# bad_loci <- lr_loci %>%
#   dplyr::filter(value <120) %>%
#   dplyr::pull(name)
# # consider removing these loci from analyses?
