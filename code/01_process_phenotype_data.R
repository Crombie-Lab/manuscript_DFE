library(tidyverse)

# set the working dir relative to repo
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# Load competitive fitness data
cf <- as.data.frame(read.csv(file = "data/raw/20170707_RILsCompetitiveFitnessData.csv", header = T))

# Load lookup table for setting strain names and reps
lt <- as.data.frame((read.csv(file = "data/raw/20180223_RandomizationRILslookuptabel.csv", header = T)))

# Process complete dataframe 
df <- left_join(cf, lt, by = c("block", "rand_id")) %>% # join lookup and competitive fitness data
  dplyr::select(-X, -X.1, -X.2, -X.3, -line_id.x) %>% # remove join artifacts
  dplyr::rename(line_id = line_id.y) %>% 
  tidyr::separate(line_id, into = c("type", "ps", "id"), sep = "_") %>% # parse line ids
  dplyr::mutate(id = ifelse(ps != "Ps", ps, id)) %>%
  dplyr::mutate(type = ifelse(type == "A.10", "RIAIL", type),
                type = ifelse(type == "N2 GO ANC", "N2_G0_ANC", type)) %>%
  tidyr::separate(id, into = c("line_id", "rep")) %>%
  dplyr::select(-ps) %>%
  tidyr::unite(full_id, c("type", "line_id"), remove = FALSE) %>%
  dplyr::mutate(full_id = ifelse(full_id == "NA_NA", NA, full_id)) %>%
  dplyr::mutate(p_focal = ((L1_MA + L1_plus_MA) / (L1_VP604 + L1_plus_VP604 + L1_MA + L1_plus_MA))) %>% # calculate proportion focal 
  dplyr::mutate(p_comp = ((L1_VP604 + L1_plus_VP604) / (L1_VP604 + L1_plus_VP604 + L1_MA + L1_plus_MA))) %>% # calculate proportion competitor NEED TO UPDATE TO USE BOTH GATES!
  dplyr::mutate(ci = p_focal/p_comp) %>% # calculate ci
  dplyr::mutate(ln_ci = logb(ci, base = exp(1))) %>% #calculate natural log of ci
  dplyr::mutate(type = factor(type, levels = c("N2_G0_ANC", "MA530", "MA563", "RIAIL", "RIL"))) %>% # set contrasts for line type
  dplyr::mutate(rep = as.double(rep),
                full_id_rep = paste0(full_id, "_", rep))

# filter data
df_proc <- df %>%
  dplyr::mutate(remove1 = case_when(full_id == "RIL_327" & rep == 5 ~ "remove",
                                   TRUE ~ "keep")) %>% # remove duplicate RIL_327 rep 5
  dplyr::mutate(remove2 = ifelse(p_focal < 0.95 & p_focal > 0.05, "keep", "remove")) %>% # filter outliers at p_focal > 0.95 and < 0.05
  dplyr::filter(!is.na(type), is.na(censor) & remove1 == "keep" & remove2 == "keep") %>% # apply filters
  dplyr::filter(L1_MA >=100 | L1_VP604 >=100) %>%
  dplyr::select(-remove1, -remove2)

# Save data frame and .csv of processed data
df_proc_export <- df_proc %>%
  dplyr::arrange(type, as.numeric(line_id), block, rep) %>%
  dplyr::select(full_id, block, rep, p_focal, p_comp, ln_ci) 

rio::export(df_proc_export, file = "data/processed/01_DFE_phenotypes.csv")
