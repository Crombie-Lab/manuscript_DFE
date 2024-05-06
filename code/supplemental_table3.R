library(tidyverse)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# read in supp table 3
st3 <- readxl::read_xlsx("data/processed/Supplemental Table 3_list of mutations & effects.xlsx")

# get effects
post.all <- data.table::fread("data/processed/Table_posterior_mutation_effects_RIL+RIAIL.csv")
mar.all <- data.table::fread("data/processed/Table_marg_effects_RIL+RIAIL.csv") %>%
  dplyr::rename_with(~paste0(.,"_mar"))

# get variant impacts - hmm, hard to match here without more info
vi1 <- data.table::fread("data/raw/RI(AI)Ls_BS_MasterFile_VEP.txt")
vi2 <- data.table::fread("data/raw/RI(AI)Ls_Indel_MasterFile_VEP.txt")
vi3 <- data.table::fread("data/raw/RI(AI)Ls_BS_severe_VEP.txt") %>%
  tidyr::separate(Location, into = c("Chromosome", "Position", "end"), remove = T)
vi4 <- data.table::fread("data/raw/RI(AI)Ls_Indel_severe_VEP.txt") %>%
  tidyr::separate(Location, into = c("Chromosome", "Position", "end"), remove = T)

# get all
vibind <- bind_rows(vi3, vi4) %>%
  dplyr::mutate(Position = as.double(Position)) %>%
  dplyr::select(-end)

# add vars to table 3
out <- st3 %>%
  dplyr::select(-`Effect (marginal)`, -`Effect (Bayesian posterior mean)`) %>%
  dplyr::bind_cols(`Effect (marginal)` = mar.all$mean_mar, `Effect (Bayesian posterior mean)` = post.all$mean) %>%
  dplyr::left_join(., vibind)

# save it
rio::export(out, file = "data/processed/Supplemental_table_3.csv")
