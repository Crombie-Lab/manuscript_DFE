library(tidyverse)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# get the joined phenotype and imputed genotype data - P-matrix
d1 <- data.table::fread("data/processed/05_DFE_joined_imputed_prob_geno_pheno.csv")

# get the list of genotyped lines
genotyped_lines <- data.table::fread("data/processed/genotyped_lines.csv")

#------------------------------------------------------------------------------#
# Step 1: Bootstrap across all lines and use weighted average from P-matrix
#------------------------------------------------------------------------------#

# hold the output
out <- NULL

# set the seed
set.seed(99)
for(i in 1:1000) {
# resample lines from the genotyped lines
line.samp <- sample(genotyped_lines$line, replace = TRUE)

# build full dataset with these lines
dat.samp.list <- NULL
for(j in 1:length(line.samp)){
  # get the line name
  line <- line.samp[j]
  # get rep data
  dat.samp <- d1 %>%
    dplyr::filter(full_id == line)
  # add to list
  dat.samp.list[[j]] <-  dat.samp
}

# bind sampled line data together and assign new id
dat.samp.df <- data.table::rbindlist(dat.samp.list, idcol = "line.samp.id") %>%
  dplyr::mutate(full_id2 = paste(full_id, line.samp.id, sep = "_"), .after = full_id)

# filter the observations to RILs with phenotypes and non-missing genotypes
d2 <- dat.samp.df %>%
  dplyr::mutate(type = stringr::str_extract(full_id, pattern = "^[^_]+")) %>% # make type from full_id
  dplyr::group_by(full_id2, block) %>%
  dplyr::mutate(line_mean_fitness = mean(ln_ci)) %>% # calculate the mean fitness within all lines ln_ci
  dplyr::select(-rep) %>% # drop rep because no longer relevant
  dplyr::ungroup() %>%
  dplyr::distinct(full_id2, block, .keep_all = T) %>% # just keep line means
  dplyr::mutate(reg.line_mean_fitness = residuals(lm(line_mean_fitness ~ as.factor(block)))) %>%
  dplyr::mutate(overall_mean_ln_ci = mean(line_mean_fitness),
                reg.overall_mean_ln_ci = mean(reg.line_mean_fitness)) %>% # get the mean of all lines
  dplyr::ungroup() %>%
  dplyr::select(full_id2, type, full_id, block:ln_ci, line_mean_fitness, reg.line_mean_fitness, overall_mean_ln_ci, reg.overall_mean_ln_ci, everything(), -line.samp.id)

# reshape to long format
d3 <- d2 %>%
  tidyr::gather(variant, genotype, !(full_id2:reg.overall_mean_ln_ci)) #%>% # reshape the variant data to long format
  #dplyr::mutate(genotype = dplyr::case_when(genotype > 0.9 ~ 1,
  #                                          genotype < 0.1 ~ 0,
  #                                          TRUE ~ genotype)) # set high confidence genotypes

# calculate variant effects
mar_effects <- d3 %>%
  dplyr::mutate(reg.line_mean_fitness_mean_scaled = reg.line_mean_fitness - reg.overall_mean_ln_ci) %>% # calculate the mean scaled fitness for each line
  dplyr::mutate(r_genotype = 1-genotype,
                mut = genotype * reg.line_mean_fitness_mean_scaled,
                wt = r_genotype * reg.line_mean_fitness_mean_scaled) %>%
  dplyr::group_by(variant) %>%
  dplyr::mutate(reg.mar_eff = mean(mut) - mean(wt)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(variant, .keep_all = T) %>%
  dplyr::select(variant, reg.mar_eff) %>%
  dplyr::mutate(boot = glue::glue("{i}"))

  out[[i]] <-  mar_effects
  message(glue::glue("{i}"))
}

# bind these up
mar_eff_boots <- data.table::rbindlist(out)

# export these data.
rio::export(mar_eff_boots, file = 'data/processed/marginal_effects_boots.csv')

#------------------------------------------------------------------------------#
# Step 2: Bootstrap across RIAILs and use weighted average from P-matrix
#------------------------------------------------------------------------------#
# get the RIAILs
genotyped_riails <- genotyped_lines %>%
  dplyr::filter(grepl(line, pattern = "RIAIL_"))

# hold the output
out.riails <- NULL

# set the seed
set.seed(138)
for(i in 1:1000) {
  # resample lines from the genotyped lines
  line.samp <- sample(genotyped_riails$line, replace = TRUE)
  
  # build full dataset with these lines
  dat.samp.list <- NULL
  for(j in 1:length(line.samp)){
    # get the line name
    line <- line.samp[j]
    # get rep data
    dat.samp <- d1 %>%
      dplyr::filter(full_id == line)
    # add to list
    dat.samp.list[[j]] <-  dat.samp
  }
  
  # bind sampled line data together and assign new id
  dat.samp.df <- data.table::rbindlist(dat.samp.list, idcol = "line.samp.id") %>%
    dplyr::mutate(full_id2 = paste(full_id, line.samp.id, sep = "_"), .after = full_id)
  
  # filter the observations to RILs with phenotypes and non-missing genotypes
  d2 <- dat.samp.df %>%
    dplyr::mutate(type = stringr::str_extract(full_id, pattern = "^[^_]+")) %>% # make type from full_id
    dplyr::group_by(full_id2, block) %>%
    dplyr::mutate(line_mean_fitness = mean(ln_ci)) %>% # calculate the mean fitness within all lines ln_ci
    dplyr::select(-rep) %>% # drop rep because no longer relevant
    dplyr::ungroup() %>%
    dplyr::distinct(full_id2, block, .keep_all = T) %>% # just keep line means
    dplyr::mutate(reg.line_mean_fitness = residuals(lm(line_mean_fitness ~ as.factor(block)))) %>%
    dplyr::mutate(overall_mean_ln_ci = mean(line_mean_fitness),
                  reg.overall_mean_ln_ci = mean(reg.line_mean_fitness)) %>% # get the mean of all lines
    dplyr::ungroup() %>%
    dplyr::select(full_id2, type, full_id, block:ln_ci, line_mean_fitness, reg.line_mean_fitness, overall_mean_ln_ci, reg.overall_mean_ln_ci, everything(), -line.samp.id)
  
  # reshape to long format
  d3 <- d2 %>%
    tidyr::gather(variant, genotype, !(full_id2:reg.overall_mean_ln_ci)) #%>% # reshape the variant data to long format
  #dplyr::mutate(genotype = dplyr::case_when(genotype > 0.9 ~ 1,
  #                                          genotype < 0.1 ~ 0,
  #                                          TRUE ~ genotype)) # set high confidence genotypes
  
  # calculate variant effects
  mar_effects <- d3 %>%
    dplyr::mutate(reg.line_mean_fitness_mean_scaled = reg.line_mean_fitness - reg.overall_mean_ln_ci) %>% # calculate the mean scaled fitness for each line
    dplyr::mutate(r_genotype = 1-genotype,
                  mut = genotype * reg.line_mean_fitness_mean_scaled,
                  wt = r_genotype * reg.line_mean_fitness_mean_scaled) %>%
    dplyr::group_by(variant) %>%
    dplyr::mutate(reg.mar_eff = mean(mut) - mean(wt)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(variant, .keep_all = T) %>%
    dplyr::select(variant, reg.mar_eff) %>%
    dplyr::mutate(boot = glue::glue("{i}"))
  
  out.riails[[i]] <-  mar_effects
  message(glue::glue("{i}"))
}

# bind these up
mar_eff_boots_riails <- data.table::rbindlist(out.riails)

# export these data.
rio::export(mar_eff_boots_riails, file = 'data/processed/marginal_effects_boots_riails.csv')

#------------------------------------------------------------------------------#
# Step 3: Bootstrap across RILs and use weighted average from P-matrix
#------------------------------------------------------------------------------#
# get the RIAILs
genotyped_rils <- genotyped_lines %>%
  dplyr::filter(grepl(line, pattern = "RIL_"))

# hold the output
out.rils <- NULL

# set the seed
set.seed(37)
for(i in 1:1000) {
  # resample lines from the genotyped lines
  line.samp <- sample(genotyped_rils$line, replace = TRUE)
  
  # build full dataset with these lines
  dat.samp.list <- NULL
  for(j in 1:length(line.samp)){
    # get the line name
    line <- line.samp[j]
    # get rep data
    dat.samp <- d1 %>%
      dplyr::filter(full_id == line)
    # add to list
    dat.samp.list[[j]] <-  dat.samp
  }
  
  # bind sampled line data together and assign new id
  dat.samp.df <- data.table::rbindlist(dat.samp.list, idcol = "line.samp.id") %>%
    dplyr::mutate(full_id2 = paste(full_id, line.samp.id, sep = "_"), .after = full_id)
  
  # filter the observations to RILs with phenotypes and non-missing genotypes
  d2 <- dat.samp.df %>%
    dplyr::mutate(type = stringr::str_extract(full_id, pattern = "^[^_]+")) %>% # make type from full_id
    dplyr::group_by(full_id2, block) %>%
    dplyr::mutate(line_mean_fitness = mean(ln_ci)) %>% # calculate the mean fitness within all lines ln_ci
    dplyr::select(-rep) %>% # drop rep because no longer relevant
    dplyr::ungroup() %>%
    dplyr::distinct(full_id2, block, .keep_all = T) %>% # just keep line means
    dplyr::mutate(reg.line_mean_fitness = residuals(lm(line_mean_fitness ~ as.factor(block)))) %>%
    dplyr::mutate(overall_mean_ln_ci = mean(line_mean_fitness),
                  reg.overall_mean_ln_ci = mean(reg.line_mean_fitness)) %>% # get the mean of all lines
    dplyr::ungroup() %>%
    dplyr::select(full_id2, type, full_id, block:ln_ci, line_mean_fitness, reg.line_mean_fitness, overall_mean_ln_ci, reg.overall_mean_ln_ci, everything(), -line.samp.id)
  
  # reshape to long format
  d3 <- d2 %>%
    tidyr::gather(variant, genotype, !(full_id2:reg.overall_mean_ln_ci)) #%>% # reshape the variant data to long format
  #dplyr::mutate(genotype = dplyr::case_when(genotype > 0.9 ~ 1,
  #                                          genotype < 0.1 ~ 0,
  #                                          TRUE ~ genotype)) # set high confidence genotypes
  
  # calculate variant effects
  mar_effects <- d3 %>%
    dplyr::mutate(reg.line_mean_fitness_mean_scaled = reg.line_mean_fitness - reg.overall_mean_ln_ci) %>% # calculate the mean scaled fitness for each line
    dplyr::mutate(r_genotype = 1-genotype,
                  mut = genotype * reg.line_mean_fitness_mean_scaled,
                  wt = r_genotype * reg.line_mean_fitness_mean_scaled) %>%
    dplyr::group_by(variant) %>%
    dplyr::mutate(reg.mar_eff = mean(mut) - mean(wt)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(variant, .keep_all = T) %>%
    dplyr::select(variant, reg.mar_eff) %>%
    dplyr::mutate(boot = glue::glue("{i}"))
  
  out.rils[[i]] <-  mar_effects
  message(glue::glue("{i}"))
}

# bind these up
mar_eff_boots_rils <- data.table::rbindlist(out.rils)

# export these data.
rio::export(mar_eff_boots_rils, file = 'data/processed/marginal_effects_boots_rils.csv')

#------------------------------------------------------------------------------#
# Step 2: Bootstrap across RIAILs and use weighted average from P-matrix
#------------------------------------------------------------------------------#
# get the RIAILs
genotyped_riails <- genotyped_lines %>%
  dplyr::filter(grepl(line, pattern = "RIAIL_"))

# hold the output
out.riails <- NULL

# set the seed
set.seed(138)
for(i in 1:1000) {
  # resample lines from the genotyped lines
  line.samp <- sample(genotyped_riails$line, replace = TRUE)
  
  # build full dataset with these lines
  dat.samp.list <- NULL
  for(j in 1:length(line.samp)){
    # get the line name
    line <- line.samp[j]
    # get rep data
    dat.samp <- d1 %>%
      dplyr::filter(full_id == line)
    # add to list
    dat.samp.list[[j]] <-  dat.samp
  }
  
  # bind sampled line data together and assign new id
  dat.samp.df <- data.table::rbindlist(dat.samp.list, idcol = "line.samp.id") %>%
    dplyr::mutate(full_id2 = paste(full_id, line.samp.id, sep = "_"), .after = full_id)
  
  # filter the observations to RILs with phenotypes and non-missing genotypes
  d2 <- dat.samp.df %>%
    dplyr::mutate(type = stringr::str_extract(full_id, pattern = "^[^_]+")) %>% # make type from full_id
    dplyr::group_by(full_id2, block) %>%
    dplyr::mutate(line_mean_fitness = mean(ln_ci)) %>% # calculate the mean fitness within all lines ln_ci
    dplyr::select(-rep) %>% # drop rep because no longer relevant
    dplyr::ungroup() %>%
    dplyr::distinct(full_id2, block, .keep_all = T) %>% # just keep line means
    dplyr::mutate(reg.line_mean_fitness = residuals(lm(line_mean_fitness ~ as.factor(block)))) %>%
    dplyr::mutate(overall_mean_ln_ci = mean(line_mean_fitness),
                  reg.overall_mean_ln_ci = mean(reg.line_mean_fitness)) %>% # get the mean of all lines
    dplyr::ungroup() %>%
    dplyr::select(full_id2, type, full_id, block:ln_ci, line_mean_fitness, reg.line_mean_fitness, overall_mean_ln_ci, reg.overall_mean_ln_ci, everything(), -line.samp.id)
  
  # reshape to long format
  d3 <- d2 %>%
    tidyr::gather(variant, genotype, !(full_id2:reg.overall_mean_ln_ci)) #%>% # reshape the variant data to long format
  #dplyr::mutate(genotype = dplyr::case_when(genotype > 0.9 ~ 1,
  #                                          genotype < 0.1 ~ 0,
  #                                          TRUE ~ genotype)) # set high confidence genotypes
  
  # calculate variant effects
  mar_effects <- d3 %>%
    dplyr::mutate(reg.line_mean_fitness_mean_scaled = reg.line_mean_fitness - reg.overall_mean_ln_ci) %>% # calculate the mean scaled fitness for each line
    dplyr::mutate(r_genotype = 1-genotype,
                  mut = genotype * reg.line_mean_fitness_mean_scaled,
                  wt = r_genotype * reg.line_mean_fitness_mean_scaled) %>%
    dplyr::group_by(variant) %>%
    dplyr::mutate(reg.mar_eff = mean(mut) - mean(wt)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(variant, .keep_all = T) %>%
    dplyr::select(variant, reg.mar_eff) %>%
    dplyr::mutate(boot = glue::glue("{i}"))
  
  out.riails[[i]] <-  mar_effects
  message(glue::glue("{i}"))
}

# bind these up
mar_eff_boots_riails <- data.table::rbindlist(out.riails)

# export these data.
rio::export(mar_eff_boots_riails, file = 'data/processed/marginal_effects_boots_riails.csv')

#------------------------------------------------------------------------------#
# Step 4: Bootstrap across RILs and use weighted average from P-matrix NO ASSAY REG
#------------------------------------------------------------------------------#
# get the RIAILs
genotyped_rils <- genotyped_lines %>%
  dplyr::filter(grepl(line, pattern = "RIL_"))

# hold the output
out.rils2 <- NULL

# set the seed
set.seed(37)
for(i in 1:1000) {
  # resample lines from the genotyped lines
  line.samp <- sample(genotyped_rils$line, replace = TRUE)
  
  # build full dataset with these lines
  dat.samp.list <- NULL
  for(j in 1:length(line.samp)){
    # get the line name
    line <- line.samp[j]
    # get rep data
    dat.samp <- d1 %>%
      dplyr::filter(full_id == line)
    # add to list
    dat.samp.list[[j]] <-  dat.samp
  }
  
  # bind sampled line data together and assign new id
  dat.samp.df <- data.table::rbindlist(dat.samp.list, idcol = "line.samp.id") %>%
    dplyr::mutate(full_id2 = paste(full_id, line.samp.id, sep = "_"), .after = full_id)
  
  # filter the observations to RILs with phenotypes and non-missing genotypes
  d2 <- dat.samp.df %>%
    dplyr::mutate(type = stringr::str_extract(full_id, pattern = "^[^_]+")) %>% # make type from full_id
    dplyr::group_by(full_id2, block) %>%
    dplyr::mutate(line_mean_fitness = mean(ln_ci)) %>% # calculate the mean fitness within all lines ln_ci
    dplyr::select(-rep) %>% # drop rep because no longer relevant
    dplyr::ungroup() %>%
    dplyr::distinct(full_id2, block, .keep_all = T) %>% # just keep line means
    #dplyr::mutate(reg.line_mean_fitness = residuals(lm(line_mean_fitness ~ as.factor(block)))) %>%
    dplyr::mutate(overall_mean_ln_ci = mean(line_mean_fitness)) %>% # get the mean of all lines
    dplyr::ungroup() %>%
    dplyr::select(full_id2, type, full_id, block:ln_ci, line_mean_fitness, overall_mean_ln_ci, everything(), -line.samp.id)
  
  # reshape to long format
  d3 <- d2 %>%
    tidyr::gather(variant, genotype, !(full_id2:overall_mean_ln_ci)) #%>% # reshape the variant data to long format
  #dplyr::mutate(genotype = dplyr::case_when(genotype > 0.9 ~ 1,
  #                                          genotype < 0.1 ~ 0,
  #                                          TRUE ~ genotype)) # set high confidence genotypes
  
  # calculate variant effects
  mar_effects <- d3 %>%
    dplyr::mutate(line_mean_fitness_mean_scaled = line_mean_fitness - overall_mean_ln_ci) %>% # calculate the mean scaled fitness for each line
    dplyr::mutate(r_genotype = 1-genotype,
                  mut = genotype * line_mean_fitness_mean_scaled,
                  wt = r_genotype * line_mean_fitness_mean_scaled) %>%
    dplyr::group_by(variant) %>%
    dplyr::mutate(mar_eff = mean(mut) - mean(wt)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(variant, .keep_all = T) %>%
    dplyr::select(variant, mar_eff) %>%
    dplyr::mutate(boot = glue::glue("{i}"))
  
  out.rils2[[i]] <-  mar_effects
  message(glue::glue("{i}"))
}

# bind these up
mar_eff_boots_rils2 <- data.table::rbindlist(out.rils2)

# export these data.
rio::export(mar_eff_boots_rils2, file = 'data/processed/marginal_effects_boots_rils_NOREG.csv')
