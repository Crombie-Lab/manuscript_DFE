library(tidyverse)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# get the joined phenotype and imputed genotype data
d1 <- data.table::fread("data/processed/05_DFE_joined_imputed_prob_geno_pheno.csv")

# get the list of genotyped lines
genotyped_lines <- data.table::fread("data/processed/genotyped_lines.csv")

#------------------------------------------------------------------------------#
# Bootstrap across lines
#------------------------------------------------------------------------------#
# # hold the output
# out <- NULL
# 
# # set the seed
# set.seed(99)
# for(i in 1:1000) {
# # resample lines from the genotyped lines
# line.samp <- sample(genotyped_lines$line, replace = TRUE)
# 
# # build full dataset with these lines
# dat.samp.list <- NULL
# for(j in 1:length(line.samp)){
#   # get the line name
#   line <- line.samp[j]
#   # get rep data
#   dat.samp <- d1 %>%
#     dplyr::filter(full_id == line)
#   # add to list
#   dat.samp.list[[j]] <-  dat.samp
# }
# 
# # bind sampled line data together and assign new id
# dat.samp.df <- data.table::rbindlist(dat.samp.list, idcol = "line.samp.id") %>%
#   dplyr::mutate(full_id2 = paste(full_id, line.samp.id, sep = "_"), .after = full_id)
# 
# # filter the observations to RILs with phenotypes and non-missing genotypes
# d2 <- dat.samp.df %>%
#   dplyr::mutate(type = stringr::str_extract(full_id, pattern = "^[^_]+")) %>% # make type from full_id
#   dplyr::group_by(full_id2, block) %>%
#   dplyr::mutate(line_mean_fitness = mean(ln_ci)) %>% # calculate the mean fitness within all lines ln_ci
#   dplyr::select(-rep) %>% # drop rep because no longer relevant
#   dplyr::ungroup() %>%
#   dplyr::distinct(full_id2, block, .keep_all = T) %>% # just keep line means
#   dplyr::mutate(reg.line_mean_fitness = residuals(lm(line_mean_fitness ~ as.factor(block)))) %>%
#   dplyr::mutate(overall_mean_ln_ci = mean(line_mean_fitness),
#                 reg.overall_mean_ln_ci = mean(reg.line_mean_fitness)) %>% # get the mean of all lines
#   dplyr::ungroup() %>%
#   dplyr::select(full_id2, type, full_id, block:ln_ci, line_mean_fitness, reg.line_mean_fitness, overall_mean_ln_ci, reg.overall_mean_ln_ci, everything(), -line.samp.id)
# 
# # reshape to long format
# d3 <- d2 %>%
#   tidyr::gather(variant, genotype, !(full_id2:reg.overall_mean_ln_ci)) %>% # reshape the variant data to long format
#   dplyr::mutate(genotype = dplyr::case_when(genotype > 0.9 ~ 1,
#                                             genotype < 0.1 ~ 0,
#                                             TRUE ~ genotype)) # set high confidence genotypes
# 
# # calculate variant effects
# mar_effects <- d3 %>%
#   dplyr::mutate(reg.line_mean_fitness_mean_scaled = reg.line_mean_fitness - reg.overall_mean_ln_ci) %>% # calculate the mean scaled fitness for each line
#   dplyr::mutate(r_genotype = 1-genotype,
#                 mut = genotype * reg.line_mean_fitness_mean_scaled,
#                 wt = r_genotype * reg.line_mean_fitness_mean_scaled) %>%
#   dplyr::group_by(variant) %>%
#   dplyr::mutate(reg.mar_eff = mean(mut) - mean(wt)) %>%
#   dplyr::ungroup() %>%
#   dplyr::distinct(variant, .keep_all = T) %>%
#   dplyr::select(variant, reg.mar_eff) %>%
#   dplyr::mutate(boot = glue::glue("{i}"))
# 
#   out[[i]] <-  mar_effects
#   message(glue::glue("{i}"))
# }
# 
# # bind these up
# mar_eff_boots <- data.table::rbindlist(out)
# 
# # export these data.
# rio::export(mar_eff_boots, file = 'data/processed/marginal_effects_boots.csv')

# load the marginal effects run 20240215
mar_eff_boots <- data.table::fread("data/processed/marginal_effects_boots.csv")

# find the mean, upper and lower 2.5% quantile
d4 <- mar_eff_boots %>%
  dplyr::group_by(variant) %>%
  dplyr::mutate(bs_mar_eff_mean = mean(reg.mar_eff),
                bs_mar_eff_lci = quantile(reg.mar_eff, probs = 0.025),
                bs_mar_eff_hci = quantile(reg.mar_eff, probs = 0.975)) %>%
  dplyr::distinct(variant, .keep_all = T) %>%
  dplyr::select(-boot, -reg.mar_eff)

#------------------------------------------------------------------------------#
# Plot both
#------------------------------------------------------------------------------#
method2 <- data.table::fread("data/processed/summary_bayesian_bootstrap.csv")

join <- method2 %>%
  dplyr::bind_cols(variant = d4$variant) %>%
  dplyr::left_join(d4) %>%
  tidyr::separate(variant, remove = F, into = c("chrom", "pos", "type"))

# look at correlation between approaches
corr <- round(cor(x=join$Posterior_mean, y=join$bs_mar_eff_mean), digits = 3)

# plot
pxy <- ggplot(join) +
  aes(x = Posterior_mean, y = bs_mar_eff_mean, fill = type) +
  #geom_point(shape = 21) +
  geom_linerange(aes(x = Posterior_mean, y = bs_mar_eff_mean,
                     xmin = `Posterior_hdi_3%`, xmax = `Posterior_hdi_97%`),
                 linewidth = 0.25) +
  geom_pointrange(aes(x = Posterior_mean, y = bs_mar_eff_mean,
                      ymin = bs_mar_eff_lci, ymax = bs_mar_eff_hci),
                  linewidth = 0.25,
                  shape = 21,
                  size = 0.75,
                  stroke = 0.25) +
  xlim(c(-0.42, 0.42)) +
  ylim(c(-0.42, 0.42)) +
  #ggpmisc::stat_poly_line() +
  #ggpmisc::stat_poly_eq(ggpmisc::use_label(c("eq", "R2"))) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  annotate(geom = "text", x = -0.28, y = 0.325, label = glue::glue("r = {corr}")) +
  #geom_smooth(method = "lm", formula= y~x) +
  labs(x = "Posterior effect", y = "Marginal effect", fill = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.825, 0.225),
        legend.background = element_rect(fill = "transparent")) 
pxy

# get the histogram for Posterior
options(scipen = 999)
pe <- ggplot(join) +
  aes(x = Posterior_mean) +
  geom_histogram(color = "black", binwidth = 0.01) +
  xlim(c(-0.165, 0.165)) +
  geom_vline(xintercept = mean(join$Posterior_mean), linetype = 2, color = "darkred") +
  annotate(geom = "label", x = mean(join$Posterior_mean), y = 3.6,
           label = glue::glue("{round(mean(join$Posterior_mean), digits=4)}"),
           color = "darkred") +
  labs(x = "Posterior effect", y = "") +
  theme_bw() +
  theme(panel.grid = element_blank())
pe

me <- ggplot(join) +
  aes(x = bs_mar_eff_mean, binwidth = 0.01) +
  geom_histogram(color = "black") +
  xlim(c(-0.165, 0.165)) +
  geom_vline(xintercept = mean(join$bs_mar_eff_mean), linetype = 2, color = "darkred") +
  annotate(geom = "label", x = mean(join$bs_mar_eff_mean), y = 2.5,
           label = glue::glue("{round(mean(join$bs_mar_eff_mean), digits=4)}"),
           color = "darkred") +
  labs(x = "Posterior effect", y = "") +
  labs(x = "Marginal effect", y = "") +
  theme_bw() +
  theme(panel.grid = element_blank())
me

# add them together
fig3ab <- cowplot::plot_grid(pe, me, nrow = 2, labels = c("a", "b"))
fig3ab
# add more
fig3 <- cowplot::plot_grid(fig3ab, pxy, ncol = 2, labels = c("", "c"), rel_widths = c(1,2))
fig3

# save it
ggsave(fig3, file = "figures/figure3.png", width = 7.5, height = 5)

# report the mean bootstrap marginal effect
mean_bs_mar_effect <- mean(join$bs_mar_eff_mean)
sd_bs_mar_effect <- sd(join$bs_mar_eff_mean)




