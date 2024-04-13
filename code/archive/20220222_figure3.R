library(tidyverse)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load the marginal effects from the marginal_effects_boot.R script
mar_eff_boots <- data.table::fread("data/processed/marginal_effects_boots.csv")

# find the mean, upper 97.5% and lower 2.5% quantiles among the marginal effect boots
mar <- mar_eff_boots %>%
  dplyr::group_by(variant) %>%
  dplyr::mutate(bs_mar_eff_mean = mean(reg.mar_eff),
                bs_mar_eff_lci = quantile(reg.mar_eff, probs = 0.025),
                bs_mar_eff_hci = quantile(reg.mar_eff, probs = 0.975)) %>%
  dplyr::distinct(variant, .keep_all = T) %>%
  dplyr::select(-boot, -reg.mar_eff)

# read in the bayesian method with posterior effects
post <- data.table::fread("data/processed/summary_bayesian_20240221.csv")

# join the the two methods together by variant - assuming same order.
join <- post %>%
  dplyr::bind_cols(variant = mar$variant) %>%
  dplyr::left_join(mar) %>%
  tidyr::separate(variant, remove = F, into = c("chrom", "pos", "type"))

# look at correlation between approaches
corr <- round(cor(x=join$Posterior_mean, y=join$bs_mar_eff_mean), digits = 3)

#------------------------------------------------------------------------------#
# Plots
#------------------------------------------------------------------------------#
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
  xlim(c(-0.525, 0.525)) +
  ylim(c(-0.525, 0.525)) +
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
fig3ab <- cowplot::plot_grid(me, pe, nrow = 2, labels = c("a", "b"))
fig3ab
# add more
fig3 <- cowplot::plot_grid(fig3ab, pxy, ncol = 2, labels = c("", "c"), rel_widths = c(1,2))
fig3

# save it
ggsave(fig3, file = "figures/figure3.png", width = 7.5, height = 5)

# report the mean bootstrap marginal effect
mean_bs_mar_effect <- mean(join$bs_mar_eff_mean)
sd_bs_mar_effect <- sd(join$bs_mar_eff_mean)
