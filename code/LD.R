library(genetics) 
library(tidyverse)
#library(LDheatmap)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# # setup data for LD
# gm <- data.table::fread("data/processed/04_DFE_imputed_genotypes_10kbMax.csv") %>%
#   dplyr::mutate(snp_id = paste0(chrom, "_", pos)) %>%
#   dplyr::select(snp_id, everything(), -chrom:-MA563) %>%
#   tidyr::pivot_longer(cols = -snp_id) %>%
#   dplyr::group_by(name) %>%
#   dplyr::mutate(n.na = sum(is.na(value))) %>%
#   dplyr::filter(n.na<146) %>% # is this right???????????????????????
#   dplyr::ungroup() %>%
#   dplyr::select(-n.na) %>%
#   tidyr::pivot_wider()
# 
# # make a list of genetics genotype objects from genotype matrix
# go_list <- list()
# for (i in 1:nrow(gm)) {
#   go_list[[i]] <- genetics::genotype(as.character(gsub(1, "T/T",
#                                                        gsub(0, "A/A", gm[i, 4:ncol(gm)])))) # missing lines! ugh should be 2:ncol(gm)
# }

# assign genotypes based on parents for unimputed data
# setup data for LD
gm <- data.table::fread("data/processed/02_DFE_unimputed_genotypes.csv") %>%
  dplyr::mutate(snp_id = paste0(chrom, "_", pos)) %>%
  dplyr::select(snp_id, everything(), -chrom:-N2_G0_ANC) %>%
  tidyr::pivot_longer(cols = -snp_id:-MA563, names_to = "strain", values_to = "genotype") %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(n.na = sum(is.na(genotype))) %>%
  dplyr::filter(n.na<146) %>% # is this the right na filter?
  dplyr::ungroup() %>%
  dplyr::select(-n.na) %>%
  dplyr::mutate(genotype = case_when(genotype == MA563 ~ "L/L",
                                      genotype == MA530 ~ "H/H",
                                      TRUE ~ NA_character_)) %>%
  dplyr::select(-MA530, -MA563) %>%
  tidyr::pivot_wider(names_from = "strain", values_from = "genotype")

go_list <- NULL
for (i in 1:nrow(gm)) {
  go_list[[i]] <- genetics::genotype(as.character(gm[i, 2:ncol(gm)]))
}

# make a dataframe from list
go_df <- data.frame(go_list)

# add informative column names to dataframe
colnames(go_df) <- (gm$snp_id)

# get mean LD in D
ldcalc_D <- t(genetics::LD(go_df)[[2]])
mean_LD_D <- mean(ldcalc_D, na.rm = T) # 0.0089
sd_LD_D <- sd(ldcalc_D, na.rm = T) # 0.0516

# write pairwise LD table
write.table(ldcalc_D, file = "data/processed/FULL_LD_D_calc_unimputed.tsv", quote=F, row.names = T, col.names = NA, sep="\t")

#=========================================================
# Calculate for RILs and RIAILs separately
#========================================================#
gm_RILs <- gm %>%
  tidyr::pivot_longer(col = 2:ncol(gm), names_to = "strain") %>%
  dplyr::filter(stringr::str_detect(strain, pattern = "RIL")) %>%
  tidyr::pivot_wider(names_from = strain)

gm_RIAILs <- gm %>%
  tidyr::pivot_longer(col = 2:ncol(gm), names_to = "strain") %>%
  dplyr::filter(stringr::str_detect(strain, pattern = "RIAIL")) %>%
  tidyr::pivot_wider(names_from = strain)

# 1) RILs
# make a list of genetics genotype objects from genotype matrix
go_list_RILs <- list()
for (i in 1:nrow(gm_RILs)) {
  go_list_RILs[[i]] <- genetics::genotype(as.character(gm[i, 2:ncol(gm_RILs)]))
}

# make a dataframe from list
go_df_RILs <- data.frame(go_list_RILs)

# add informative column names to dataframe
colnames(go_df_RILs) <- (gm_RILs$snp_id)

# perform ld calculation for all pairwise combinations and output as a matrix 
ldcalc_RILs <- t(genetics::LD(go_df_RILs)[[2]])

# write pairwise LD table
write.table(ldcalc_RILs, file = "data/processed/RILs_LD_D_calc_unimputed.tsv", quote=F, row.names = T, col.names = NA, sep="\t")

# 2) RIAILs
# make a list of genetics genotype objects from genotype matrix
go_list_RIAILs <- list()
for (i in 1:nrow(gm_RIAILs)) {
  go_list_RIAILs[[i]] <- genetics::genotype(as.character(gm[i, 2:ncol(gm_RIAILs)]))
}

# make a dataframe from list
go_df_RIAILs <- data.frame(go_list_RIAILs)

# add informative column names to dataframe
colnames(go_df_RIAILs) <- (gm_RIAILs$snp_id)

# perform ld calculation for all pairwise combinations and output as a matrix 
ldcalc_RIAILs <- t(genetics::LD(go_df_RIAILs)[[2]])

# write pairwise LD table
write.table(ldcalc_RIAILs, file = "data/processed/RIAILs_LD_D_calc_unimputed.tsv", quote=F, row.names = T, col.names = NA, sep="\t")

#===========================================================================#
# Calculate avg pairwise LD within chormosome for RILs and RIAILs seperately
#===========================================================================#

#=========================================================
# Calculate overall mean LD, sd of LD, and SEM LD
#========================================================#
mean_ld = mean(as.matrix(ldcalc_D), na.rm = T)
sd_ld = sd(as.matrix(ldcalc_D), na.rm = T)
sem_ld =  sd(as.vector(t(ldcalc_D[!is.na(ldcalc_D)]))) / sqrt(length(as.vector(t(ldcalc_D[!is.na(ldcalc_D)]))))

# Full RILs and RIAILs
I <- as.data.frame(ldcalc_D) %>%
  dplyr::select(starts_with("I_")) %>% # filter columns to chrom 
  tibble::rownames_to_column(., var = "POS") %>% # get the row names to a column so we can filter on it
  dplyr::filter(stringr::str_detect(POS, pattern = "^I_")) %>% #filter to chrom
  dplyr::select(-POS)

I_mean_ld <- mean(as.matrix(I), na.rm = T)
I_sd_ld <- sd(as.matrix(I), na.rm = T)
I_sem_ld <- sd(as.vector(t(I[!is.na(I)]))) / sqrt(length(as.vector(t(I[!is.na(I)]))))

II <- as.data.frame(ldcalc_D) %>%
  dplyr::select(starts_with("II_")) %>% # filter columns to chrom 
  tibble::rownames_to_column(., var = "POS") %>% # get the row names to a column so we can filter on it
  dplyr::filter(stringr::str_detect(POS, pattern = "^II_")) %>% #filter to chrom
  dplyr::select(-POS)

II_mean_ld <- mean(as.matrix(II), na.rm = T)
II_sd_ld <- sd(as.matrix(II), na.rm = T)
II_sem_ld <- sd(as.vector(t(II[!is.na(II)]))) / sqrt(length(as.vector(t(II[!is.na(II)]))))

III <- as.data.frame(ldcalc_D) %>%
  dplyr::select(starts_with("III_")) %>% # filter columns to chrom 
  tibble::rownames_to_column(., var = "POS") %>% # get the row names to a column so we can filter on it
  dplyr::filter(stringr::str_detect(POS, pattern = "^III_")) %>% #filter to chrom
  dplyr::select(-POS)

III_mean_ld <- mean(as.matrix(III), na.rm = T)
III_sd_ld <- sd(as.matrix(III), na.rm = T)
III_sem_ld <- sd(as.vector(t(III[!is.na(III)]))) / sqrt(length(as.vector(t(III[!is.na(III)]))))

IV <- as.data.frame(ldcalc_D) %>%
  dplyr::select(starts_with("IV_")) %>% # filter columns to chrom 
  tibble::rownames_to_column(., var = "POS") %>% # get the row names to a column so we can filter on it
  dplyr::filter(stringr::str_detect(POS, pattern = "^IV_")) %>% #filter to chrom
  dplyr::select(-POS)

IV_mean_ld <- mean(as.matrix(IV), na.rm = T)
IV_sd_ld <- sd(as.matrix(IV), na.rm = T)
IV_sem_ld <- sd(as.vector(t(IV[!is.na(IV)]))) / sqrt(length(as.vector(t(IV[!is.na(IV)]))))

V <- as.data.frame(ldcalc_D) %>%
  dplyr::select(starts_with("V_")) %>% # filter columns to chrom 
  tibble::rownames_to_column(., var = "POS") %>% # get the row names to a column so we can filter on it
  dplyr::filter(stringr::str_detect(POS, pattern = "^V_")) %>% #filter to chrom
  dplyr::select(-POS)

V_mean_ld <- mean(as.matrix(V), na.rm = T)
V_sd_ld <- sd(as.matrix(V), na.rm = T)
V_sem_ld <- sd(as.vector(t(V[!is.na(V)]))) / sqrt(length(as.vector(t(V[!is.na(V)]))))

X <- as.data.frame(ldcalc_D) %>%
  dplyr::select(starts_with("X_")) %>% # filter columns to chrom 
  tibble::rownames_to_column(., var = "POS") %>% # get the row names to a column so we can filter on it
  dplyr::filter(stringr::str_detect(POS, pattern = "^X_")) %>% #filter to chrom
  dplyr::select(-POS)

X_mean_ld <- mean(as.matrix(X), na.rm = T)
X_sd_ld <- sd(as.matrix(X), na.rm = T)
X_sem_ld <- sd(as.vector(t(X[!is.na(X)]))) / sqrt(length(as.vector(t(X[!is.na(X)]))))

intra_ld <-tibble::tibble(I_mean_ld = I_mean_ld,
                          II_mean_ld = II_mean_ld,
                          III_mean_ld = III_mean_ld,
                          IV_mean_ld = IV_mean_ld,
                          V_mean_ld = V_mean_ld,
                          X_mean_ld = X_mean_ld) %>%
  tidyr::pivot_longer(cols = everything(),  values_to = c("mean_intrachrom_pairwise_ld")) %>%
  tidyr::separate(name, into = c("chr", "type")) 

intra_ld_sem <- tibble::tibble(I_sem_ld = I_sem_ld,
                               II_sem_ld = II_sem_ld,
                               III_sem_ld = III_sem_ld,
                               IV_sem_ld = IV_sem_ld,
                               V_sem_ld = V_sem_ld,
                               X_sem_ld = X_sem_ld) %>%
  tidyr::pivot_longer(cols = everything(),  values_to = c("sem_intrachrom_pairwise_ld")) %>%
  tidyr::separate(name, into = c("chr", "type"))

intra_ld2 <- intra_ld %>%
  cbind(intra_ld_sem$sem_intrachrom_pairwise_ld) %>%
  dplyr::rename(sem_intrachrom_pairwise_ld = `intra_ld_sem$sem_intrachrom_pairwise_ld`) %>%
  dplyr::mutate(type = "all")

#===================================================#
# RILs: manual calc all b/c out of time to loop. ugh!
#===================================================#
I_ril <- as.data.frame(ldcalc_RILs) %>%
  dplyr::select(starts_with("I_")) %>% # filter columns to chrom 
  tibble::rownames_to_column(., var = "POS") %>% # get the row names to a column so we can filter on it
  dplyr::filter(stringr::str_detect(POS, pattern = "^I_")) %>% #filter to chrom
  dplyr::select(-POS)
I_ril_mean_ld <- mean(as.matrix(I_ril), na.rm = T)
I_ril_sem_ld <- sd(as.vector(t(I_ril[!is.na(I_ril)]))) / sqrt(length(as.vector(t(I_ril[!is.na(I_ril)]))))

II_ril <- as.data.frame(ldcalc_RILs) %>%
  dplyr::select(starts_with("II_")) %>% # filter columns to chrom 
  tibble::rownames_to_column(., var = "POS") %>% # get the row names to a column so we can filter on it
  dplyr::filter(stringr::str_detect(POS, pattern = "^II_")) %>% #filter to chrom
  dplyr::select(-POS)
II_ril_mean_ld <- mean(as.matrix(II_ril), na.rm = T)
II_ril_sem_ld <- sd(as.vector(t(II_ril[!is.na(II_ril)]))) / sqrt(length(as.vector(t(II_ril[!is.na(II_ril)]))))

III_ril <- as.data.frame(ldcalc_RILs) %>%
  dplyr::select(starts_with("III_")) %>% # filter columns to chrom 
  tibble::rownames_to_column(., var = "POS") %>% # get the row names to a column so we can filter on it
  dplyr::filter(stringr::str_detect(POS, pattern = "^III_")) %>% #filter to chrom
  dplyr::select(-POS)
III_ril_mean_ld <- mean(as.matrix(III_ril), na.rm = T)
III_ril_sem_ld <- sd(as.vector(t(III_ril[!is.na(III_ril)]))) / sqrt(length(as.vector(t(III_ril[!is.na(III_ril)]))))

IV_ril <- as.data.frame(ldcalc_RILs) %>%
  dplyr::select(starts_with("IV_")) %>% # filter columns to chrom 
  tibble::rownames_to_column(., var = "POS") %>% # get the row names to a column so we can filter on it
  dplyr::filter(stringr::str_detect(POS, pattern = "^IV_")) %>% #filter to chrom
  dplyr::select(-POS)
IV_ril_mean_ld <- mean(as.matrix(IV_ril), na.rm = T)
IV_ril_sem_ld <- sd(as.vector(t(IV_ril[!is.na(IV_ril)]))) / sqrt(length(as.vector(t(IV_ril[!is.na(IV_ril)]))))

V_ril <- as.data.frame(ldcalc_RILs) %>%
  dplyr::select(starts_with("V_")) %>% # filter columns to chrom 
  tibble::rownames_to_column(., var = "POS") %>% # get the row names to a column so we can filter on it
  dplyr::filter(stringr::str_detect(POS, pattern = "^V_")) %>% #filter to chrom
  dplyr::select(-POS)
V_ril_mean_ld <- mean(as.matrix(V_ril), na.rm = T)
V_ril_sem_ld <- sd(as.vector(t(V_ril[!is.na(V_ril)]))) / sqrt(length(as.vector(t(V_ril[!is.na(V_ril)]))))

X_ril <- as.data.frame(ldcalc_RILs) %>%
  dplyr::select(starts_with("X_")) %>% # filter columns to chrom 
  tibble::rownames_to_column(., var = "POS") %>% # get the row names to a column so we can filter on it
  dplyr::filter(stringr::str_detect(POS, pattern = "^X_")) %>% #filter to chrom
  dplyr::select(-POS)
X_ril_mean_ld <- mean(as.matrix(X_ril), na.rm = T)
X_ril_sem_ld <- sd(as.vector(t(X_ril[!is.na(X_ril)]))) / sqrt(length(as.vector(t(X_ril[!is.na(X_ril)]))))

ril_intra_ld <-tibble::tibble(I_ril_mean_ld = I_ril_mean_ld,
                              II_ril_mean_ld = II_ril_mean_ld,
                              III_ril_mean_ld = III_ril_mean_ld,
                              IV_ril_mean_ld = IV_ril_mean_ld,
                              V_ril_mean_ld = V_ril_mean_ld,
                              X_ril_mean_ld = X_ril_mean_ld) %>%
  tidyr::pivot_longer(cols = everything(),  values_to = c("mean_intrachrom_pairwise_ld")) %>%
  tidyr::separate(name, into = c("chr", "type"))

ril_intra_ld_sem <-tibble::tibble(I_ril_sem_ld = I_ril_sem_ld,
                                  II_ril_sem_ld = II_ril_sem_ld,
                                  III_ril_sem_ld = III_ril_sem_ld,
                                  IV_ril_sem_ld = IV_ril_sem_ld,
                                  V_ril_sem_ld = V_ril_sem_ld,
                                  X_ril_sem_ld = X_ril_sem_ld) %>%
  tidyr::pivot_longer(cols = everything(),  values_to = c("sem_intrachrom_pairwise_ld")) %>%
  tidyr::separate(name, into = c("chr", "type"))

ril_intra_ld2 <- ril_intra_ld %>%
  cbind(ril_intra_ld_sem$sem_intrachrom_pairwise_ld)%>%
  dplyr::rename(sem_intrachrom_pairwise_ld = `ril_intra_ld_sem$sem_intrachrom_pairwise_ld`)

#=======================================#
# RIAILs: now for the others
#======++++=============================#
I_riail <- as.data.frame(ldcalc_RIAILs) %>%
  dplyr::select(starts_with("I_")) %>% # filter columns to chrom 
  tibble::rownames_to_column(., var = "POS") %>% # get the row names to a column so we can filter on it
  dplyr::filter(stringr::str_detect(POS, pattern = "^I_")) %>% #filter to chrom
  dplyr::select(-POS)
I_riail_mean_ld <- mean(as.matrix(I_riail), na.rm = T)
I_riail_sem_ld <- sd(as.vector(t(I_riail[!is.na(I_riail)]))) / sqrt(length(as.vector(t(I_riail[!is.na(I_riail)]))))

II_riail <- as.data.frame(ldcalc_RIAILs) %>%
  dplyr::select(starts_with("II_")) %>% # filter columns to chrom 
  tibble::rownames_to_column(., var = "POS") %>% # get the row names to a column so we can filter on it
  dplyr::filter(stringr::str_detect(POS, pattern = "^II_")) %>% #filter to chrom
  dplyr::select(-POS)
II_riail_mean_ld <- mean(as.matrix(II_riail), na.rm = T)
II_riail_sem_ld <- sd(as.vector(t(II_riail[!is.na(II_riail)]))) / sqrt(length(as.vector(t(II_riail[!is.na(II_riail)]))))

III_riail <- as.data.frame(ldcalc_RIAILs) %>%
  dplyr::select(starts_with("III_")) %>% # filter columns to chrom 
  tibble::rownames_to_column(., var = "POS") %>% # get the row names to a column so we can filter on it
  dplyr::filter(stringr::str_detect(POS, pattern = "^III_")) %>% #filter to chrom
  dplyr::select(-POS)
III_riail_mean_ld <- mean(as.matrix(III_riail), na.rm = T)
III_riail_sem_ld <- sd(as.vector(t(III_riail[!is.na(III_riail)]))) / sqrt(length(as.vector(t(III_riail[!is.na(III_riail)]))))

IV_riail <- as.data.frame(ldcalc_RIAILs) %>%
  dplyr::select(starts_with("IV_")) %>% # filter columns to chrom 
  tibble::rownames_to_column(., var = "POS") %>% # get the row names to a column so we can filter on it
  dplyr::filter(stringr::str_detect(POS, pattern = "^IV_")) %>% #filter to chrom
  dplyr::select(-POS)
IV_riail_mean_ld <- mean(as.matrix(IV_riail), na.rm = T)
IV_riail_sem_ld <- sd(as.vector(t(IV_riail[!is.na(IV_riail)]))) / sqrt(length(as.vector(t(IV_riail[!is.na(IV_riail)]))))

V_riail <- as.data.frame(ldcalc_RIAILs) %>%
  dplyr::select(starts_with("V_")) %>% # filter columns to chrom 
  tibble::rownames_to_column(., var = "POS") %>% # get the row names to a column so we can filter on it
  dplyr::filter(stringr::str_detect(POS, pattern = "^V_")) %>% #filter to chrom
  dplyr::select(-POS)
V_riail_mean_ld <- mean(as.matrix(V_riail), na.rm = T)
V_riail_sem_ld <- sd(as.vector(t(V_riail[!is.na(V_riail)]))) / sqrt(length(as.vector(t(V_riail[!is.na(V_riail)]))))

X_riail <- as.data.frame(ldcalc_RIAILs) %>%
  dplyr::select(starts_with("X_")) %>% # filter columns to chrom 
  tibble::rownames_to_column(., var = "POS") %>% # get the row names to a column so we can filter on it
  dplyr::filter(stringr::str_detect(POS, pattern = "^X_")) %>% #filter to chrom
  dplyr::select(-POS)
X_riail_mean_ld <- mean(as.matrix(X_riail), na.rm = T)
X_riail_sem_ld <- sd(as.vector(t(X_riail[!is.na(X_riail)]))) / sqrt(length(as.vector(t(X_riail[!is.na(X_riail)]))))

# make a single dataframe
riails_intra_ld <-tibble::tibble(I_riail_mean_ld = I_riail_mean_ld,
                                 II_riail_mean_ld = II_riail_mean_ld,
                                 III_riail_mean_ld = III_riail_mean_ld,
                                 IV_riail_mean_ld = IV_riail_mean_ld,
                                 V_riail_mean_ld = V_riail_mean_ld,
                                 X_riail_mean_ld = X_riail_mean_ld) %>%
  tidyr::pivot_longer(cols = everything(),  values_to = c("mean_intrachrom_pairwise_ld")) %>%
  tidyr::separate(name, into = c("chr", "type"))

riails_intra_ld_sem <-tibble::tibble(I_riail_sem_ld = I_riail_sem_ld,
                                     II_riail_sem_ld = II_riail_sem_ld,
                                     III_riail_sem_ld = III_riail_sem_ld,
                                     IV_riail_sem_ld = IV_riail_sem_ld,
                                     V_riail_sem_ld = V_riail_sem_ld,
                                     X_riail_sem_ld = X_riail_sem_ld) %>%
  tidyr::pivot_longer(cols = everything(),  values_to = c("sem_intrachrom_pairwise_ld")) %>%
  tidyr::separate(name, into = c("chr", "type"))

riails_intra_ld2 <- riails_intra_ld %>%
  cbind(riails_intra_ld_sem$sem_intrachrom_pairwise_ld) %>%
  dplyr::rename(sem_intrachrom_pairwise_ld = `riails_intra_ld_sem$sem_intrachrom_pairwise_ld`)

# output the intra chromosomal ld data for rils and riails
rio::export(intra_ld2, file = "data/processed/FULL_intra_chrom_ld_D_unimputed.tsv")
rio::export(ril_intra_ld2, file = "data/processed/ril_intra_chrom_ld_D_unimputed.tsv")
rio::export(riails_intra_ld2, file = "data/processed/riail_intra_chrom_ld_D_unimputed.tsv")

#======================================================#
# Plot intrachromosomal LD for each
#======================================================#
intra_ld_full2 <- intra_ld2 %>%
  dplyr::full_join(ril_intra_ld2) %>%
  dplyr::full_join(riails_intra_ld2)

intra_ld_plot <- ggplot(intra_ld_full2) +
  aes(x = factor(chr, levels = c("I", "II", "III", "IV", "V", "X")), y = mean_intrachrom_pairwise_ld, fill = type) +
  geom_bar(stat = "identity", width = 0.5, position = "dodge") +
  scale_fill_manual(labels = c("all strains", "RIAILs", "RILs"), values = c("grey40", "#00B9F1", "#00A875")) +
  geom_errorbar(aes(ymin=mean_intrachrom_pairwise_ld -sem_intrachrom_pairwise_ld, ymax=mean_intrachrom_pairwise_ld + sem_intrachrom_pairwise_ld), width=.25, position=position_dodge(.5), size = 0.25) +
  #geom_col(fill = "dark blue", gr) +
  #geom_hline(yintercept = unique(inter_chr_ld$exp_ld), linetype = "dashed", color = "light blue") +
  theme_bw() +
  labs(x = "Chromosome", y = 'Mean intra-chromosomal LD (D)', fill = "") +
  theme(legend.position = "none")
intra_ld_plot

#cowplot::ggsave2(intra_ld_plot, filename = "figures/mean_intra_ld_plot_unimputed.png", width = 7.5, height = 3.5)

#=====================================================#
# Now calculate mean interchromosomal LDs
#=====================================================#
# mean interchrom ld function. Pairs is set of patterns to match chroms, see make pairs below function.
mean_ic_ld <- function(pairs, data) {
  # make a list
  lds <- NULL
  # loop through the pairs
  for(i in 1:nrow(pairs)){
    l <- pairs[[i,1]]
    r <- pairs[[i,2]]
    
    # shape data 
    d <- as.data.frame(data) %>%
      dplyr::select_if(grepl(l, names(.))) %>% # filter columns to chrom 
      tibble::rownames_to_column(., var = "POS") %>% # get the row names to a column so we can filter on it
      dplyr::filter(stringr::str_detect(POS, pattern = r)) %>% #filter to chrom
      dplyr::select(-POS)
    
    # calculate mean
    mean_d <- tibble::tibble(chr1 = stringr::str_replace_all(l, pattern = "\\^|_", replacement = ""),
                             chr2 = stringr::str_replace_all(r, pattern = "\\^|_", replacement = ""),
                             mean_ic_ld = mean(as.matrix(d), na.rm = T),
                             sd_ic_ld = sd(as.matrix(d), na.rm = T),
                             sem_ic_ld = sd_ic_ld / sqrt(length(as.vector(t(d[!is.na(d)])))))
    
    # Add to list
    lds[[i]] <- mean_d
  }
  # return the list as dataframe
  out <- data.table::rbindlist(lds)
  return(out)
}

# Make pairs
pairs <- as.data.frame(t(combn(c("^I_", "^II_", "^III_", "^IV_", "^V_", "^X_"),2)))

# run function for all rils and riails
ic_ld_D <- mean_ic_ld(pairs = pairs, data = ldcalc_D) %>%
  dplyr::mutate(type = "all")

# run function for rils
ril_ic_ld_D <- mean_ic_ld(pairs = pairs, data = ldcalc_RILs) %>%
  dplyr::mutate(type = "ril")

# run function for riails
riail_ic_ld_D <- mean_ic_ld(pairs = pairs, data = ldcalc_RIAILs) %>%
  dplyr::mutate(type = "riail")

# output these tables
rio::export(ic_ld_D, file = "data/processed/full_mean_interchrom_ld_D_unimputed.tsv")
rio::export(ril_ic_ld_D, file = "data/processed/ril_mean_interchrom_ld_D_unimputed.tsv")
rio::export(riail_ic_ld_D, file = "data/processed/riail_mean_interchrom_ld_D_unimputed.tsv")
#======================================================#
# Plot interchromosomal LD for each
#======================================================#
pair_order <- ic_ld_D %>%
  dplyr::mutate(chrom_pair = paste0(chr1, "-", chr2)) %>%
  dplyr::pull(chrom_pair)

ic_ld_full2 <- ic_ld_D %>%
  dplyr::full_join(ril_ic_ld_D) %>%
  dplyr::full_join(riail_ic_ld_D) %>%
  dplyr::mutate(chrom_pair = paste0(chr1, "-", chr2))

ic_ld_plot <- ggplot(ic_ld_full2) +
  aes(x = factor(chrom_pair, levels = c(pair_order)), y = mean_ic_ld, fill = type) +
  geom_bar(stat = "identity", width = 0.5, position = "dodge") +
  scale_fill_manual(labels = c("all strains", "RIAILs", "RILs"), values = c("grey40", "#00B9F1", "#00A875")) +
  #geom_col(fill = "dark blue", gr) +
  geom_errorbar(aes(ymin= mean_ic_ld - sem_ic_ld, ymax= mean_ic_ld + sem_ic_ld), width=.25, position=position_dodge(.5), size = 0.25) +
  #geom_hline(yintercept = unique(inter_chr_ld$exp_ld), linetype = "dashed", color = "light blue") +
  theme_bw() +
  labs(x = "Chromosome comparison", y = 'Mean inter-chromosomal LD (D)', fill = "") # expression('Mean inter-chromosomal\nlinkage disequillibrium r'^'2')
ic_ld_plot

#cowplot::ggsave2(ic_ld_plot, filename = "plots/mean_ic_ld_plot_unimputed.png", width = 7.5, height = 3.5)

sup_fig2 <- cowplot::plot_grid(intra_ld_plot, ic_ld_plot, ncol = 2, rel_widths = c(.8, 2), labels = c("a", "b"), align = "vh", axis = "bl")
cowplot::ggsave2(sup_fig2, filename = "figures/Supplemental_figure3.png", width = 10, height = 4)

# #==============================================================================#
# # Plot heatmaps for intra chromosomal LD in r, R^2, and D with ggplot
# #==============================================================================#
# # get parents genotypes
# parent_genos <- data.table::fread("data/raw/parent_genos.csv")
# pg.2 <- parent_genos %>%
#   dplyr::mutate(loci = stringr::str_replace(loci, pattern = regex("_[^_]*$"), replacement = ""),
#                 label = ifelse(MA530 == 1, "MA530", "MA563"))
# 
# # load the mutational effects bootstrap function
# source("code/functions.R")
# 
# # get ALL LD Measures for all Strains
# ld.all <- t(genetics::LD(go_df))
# # Export all LD measures
# rio::export(ld.all[[2]], file = "data/processed/LD_matrices/full_LD_unimputed_D.csv", row.names = T)
# rio::export(ld.all[[3]], file = "data/processed/LD_matrices/full_LD_unimputed_Dprime.csv", row.names = T)
# rio::export(ld.all[[4]], file = "data/processed/LD_matrices/full_LD_unimputed_r.csv", row.names = T)
# rio::export(ld.all[[5]], file = "data/processed/LD_matrices/full_LD_unimputed_r2.csv", row.names = T)
# rio::export(ld.all[[6]], file = "data/processed/LD_matrices/full_LD_unimputed_nlines.csv", row.names = T)
# rio::export(ld.all[[7]], file = "data/processed/LD_matrices/full_LD_unimputed_Chi-sq.csv", row.names = T)
# rio::export(ld.all[[8]], file = "data/processed/LD_matrices/full_LD_unimputed_p-value.csv", row.names = T)
# 
# # get it in D
# ld.all.D <- t(genetics::LD(go_df)[[2]])
# ld.RILs.D <- t(genetics::LD(go_df_RILs)[[2]])
# ld.RIAILs.D <- t(genetics::LD(go_df_RIAILs)[[2]])
# 
# # get the LD in r
# ld.all.r <- t(genetics::LD(go_df)[[4]])
# ld.RILs.r <- t(genetics::LD(go_df_RILs)[[4]])
# ld.RIAILs.r <- t(genetics::LD(go_df_RIAILs)[[4]])
# 
# # get the LD in R^2
# ld.all.r2 <- t(genetics::LD(go_df)[[5]])
# ld.RILs.r2 <- t(genetics::LD(go_df_RILs)[[5]])
# ld.RIAILs.r2 <- t(genetics::LD(go_df_RIAILs)[[5]])
# 
# # plot the intra chrom LD in D
# p1.D <- LDchroms(data = ld.all.D, parents = pg.2, type = 'D')
# p2.D <- LDchroms(data = ld.RILs.D, parents = pg.2, type = 'D')
# p3.D <- LDchroms(data = ld.RIAILs.D, parents = pg.2, type = 'D')
# 
# # plot the intra chrom LD in r
# p1 <- LDchroms(data = ld.all.r, parents = pg.2, type = 'r')
# p2 <- LDchroms(data = ld.RILs.r, parents = pg.2, type = 'r')
# p3 <- LDchroms(data = ld.RIAILs.r, parents = pg.2, type = 'r')
# 
# # plot the intra chrom LD in r2
# p1r2 <- LDchroms(data = ld.all.r2, parents = pg.2, type = 'R2')
# p2r2 <- LDchroms(data = ld.RILs.r2, parents = pg.2, type = 'R2')
# p3r2 <- LDchroms(data = ld.RIAILs.r2, parents = pg.2, type = 'R2')
# 
# # save them
# cowplot::ggsave2(p1, filename = "plots/LD_all_r_unimputed.png", width = 27, height = 12)
# cowplot::ggsave2(p2, filename = "plots/LD_RILs_r_unimputed.png", width = 27, height = 12)
# cowplot::ggsave2(p3, filename = "plots/LD_RIAILs_r_unimputed.png", width = 27, height = 12)
# 
# cowplot::ggsave2(p1r2, filename = "plots/LD_all_r2_unimputed.png", width = 27, height = 12)
# cowplot::ggsave2(p2r2, filename = "plots/LD_RILs_r2_unimputed.png", width = 27, height = 12)
# cowplot::ggsave2(p3r2, filename = "plots/LD_RIAILs_r2_unimputed.png", width = 27, height = 12)
# 
# cowplot::ggsave2(p1.D, filename = "plots/LD_all_D_unimputed.pdf", width = 27, height = 12)
# cowplot::ggsave2(p2.D, filename = "plots/LD_RILs_D_unimputed.png", width = 27, height = 12)
# cowplot::ggsave2(p3.D, filename = "plots/LD_RIAILs_D_unimputed.png", width = 27, height = 12)
