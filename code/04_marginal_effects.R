library(tidyverse)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# get the joined phenotype and imputed genotype data
d1 <- data.table::fread("data/processed/05_DFE_joined_imputed_prob_geno_pheno.csv")

# filter the observations to RILs with phenotypes and non-missing genotypes
d2 <- d1 %>%
  dplyr::mutate(type = stringr::str_extract(full_id, pattern = "^[^_]+")) %>% # make type from full_id
  dplyr::mutate(row.na.n = rowSums(is.na(.))) %>% # count missing genotypes per line
  dplyr::group_by(full_id, block) %>%
  dplyr::mutate(line_mean_fitness = mean(ln_ci)) %>% # calculate the mean fitness within all lines log_ci
  dplyr::select(-rep) %>% # drop rep because no longer relevant
  dplyr::ungroup() %>%
  dplyr::distinct(full_id, block, .keep_all = T) %>% # just keep line means
  dplyr::mutate(reg.line_mean_fitness = residuals(lm(line_mean_fitness ~ block))) %>%
  dplyr::mutate(overall_mean_ln_ci = mean(line_mean_fitness),
                reg.overall_mean_ln_ci = mean(reg.line_mean_fitness)) %>% # get the mean of all lines
  dplyr::filter(stringr::str_detect(full_id, pattern = "N2_G0|MA563|MA530") == F) %>% # remove parents and ancestor
  dplyr::filter(row.na.n != 169) %>% # remove lines with all missing genotypes
  dplyr::ungroup() %>%
  dplyr::select(-row.na.n) %>%
  dplyr::select(full_id, type, block:ln_ci, line_mean_fitness, reg.line_mean_fitness, overall_mean_ln_ci, reg.overall_mean_ln_ci, everything())

# plot assay effect
rawp <- ggplot(d2) +
  aes(x = block, y = line_mean_fitness, group = block, fill = type) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25, shape = 21)
rawp

regp <- ggplot(d2) +
  aes(x = block, y = reg.line_mean_fitness, group = block, fill = type) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25, shape = 21)
regp

# Check the block effect
m1 <- lm(data = d2, line_mean_fitness ~ block)
anova(m1)
summary(m1)

# reshape to long format
d3 <- d2 %>%
  tidyr::gather(variant, genotype, !(full_id:reg.overall_mean_ln_ci)) %>% # reshape the variant data to long format
  dplyr::mutate(genotype = dplyr::case_when(genotype > 0.9 ~ 1,
                                            genotype < 0.1 ~ 0,
                                            TRUE ~ genotype)) # set high confidence genotypes

# calculate variant effects
# Although I suppose we could calculate a weighted average, with each fitness value 
#weighted by its imputed probability.  Then, each variant is entered for its 1 value AND its 0 value,
# so if a locus is a “1” and its fitness is 0.75, its weighted “1” average is (1)(0.75) and
# its “0” average is (0)(0.75).  

# But then, for a variant with imputed probability 0.8 and fitness 0.75, 
# its “1” value is (0.8)(0.75) and its “0” value is (0.2)(0.75).  
# So the marginal effect at a locus is the difference between its “1” average and its “0” average.
mar_effects <- d3 %>%
  dplyr::mutate(reg.line_mean_fitness_mean_scaled = reg.line_mean_fitness - reg.overall_mean_ln_ci) %>% # calculate the mean scaled fitness for each line
  dplyr::mutate(r_genotype = 1-genotype,
                mut = genotype * reg.line_mean_fitness_mean_scaled,
                wt = r_genotype * reg.line_mean_fitness_mean_scaled) %>%
  dplyr::group_by(variant) %>%
  dplyr::mutate(reg.mar_eff = mean(mut) - mean(wt)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(variant, .keep_all = T) %>%
  dplyr::mutate(reg.overall_mean_mar_eff = mean(reg.mar_eff, na.rm =T),
                reg.overall_sd_mar_eff = sd(reg.mar_eff))

# export these data.
rio::export(mut_effects, file = 'data/processed/marginal_effects.csv')

#------------------------------------------------------------------------------#
# Plot both
#------------------------------------------------------------------------------#
method2 <- data.table::fread("data/processed/summary_bayesian_bootstrap.csv")

join <- method2 %>%
  dplyr::bind_cols(variant = mar_effects$variant) %>%
  dplyr::left_join(mar_effects)

# look at correlation between approaches
corr <- round(cor(x=join$Posterior_mean, y=join$reg.mar_eff), digits = 3)

# plot
p <- ggplot(join) +
  aes(x = Posterior_mean, y = reg.mar_eff) +
  geom_point(shape = 21) +
  xlim(c(-0.15, 0.15)) +
  ylim(c(-0.15, 0.15)) +
  #ggpmisc::stat_poly_line() +
  #ggpmisc::stat_poly_eq(ggpmisc::use_label(c("eq", "R2"))) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  annotate(geom = "text", x = -0.1, y = 0.125, label = glue::glue("r = {corr}")) +
  #geom_smooth(method = "lm", formula= y~x) +
  theme_classic()
p

p2 <- ggExtra::ggMarginal(p, type = "histogram")
p2

#summary(corr)
#lm.1 <- lm(data = join, formula = reg.mar_eff ~ Posterior_mean)
#summary(lm.1)
ggsave(p2, file = "figures/figure3.png", width = 5, height = 5)
