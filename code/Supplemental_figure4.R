library(tidyverse)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

#------------------------------------------------------------------------------#
# Get allele freq
#------------------------------------------------------------------------------#
# calculate the allele freq from P-matrix
af <- data.table::fread("data/processed/05_DFE_joined_imputed_prob_geno_pheno.csv") %>%
  dplyr::distinct(full_id, .keep_all = T) %>%
  dplyr::filter(grepl(full_id, pattern = "RIAIL_|RIL_")) %>%
  dplyr::filter(!is.na(I_113752_snp)) %>%
  dplyr::summarise(across(7:175, ~ sum(.x)/n())) %>%
  tidyr::pivot_longer(cols = everything(), names_to = "variant", values_to = "mut_af")

#------------------------------------------------------------------------------#
# Get other stuff
#------------------------------------------------------------------------------#
# get the variant ids and the parents for plotting using hard coded parent_genos.csv file
vars <- data.table::fread("data/raw/parent_genos.csv") %>%
  tidyr::separate(loci, into = c("chrom", "pos", "type"), sep = "_", remove = T) %>%
  dplyr::mutate(parent = ifelse(MA530 == 1, "MA530", "MA563"),
                type = ifelse(type == "snp", "snp", "indel"),
                pos = as.numeric(pos))

# read in the bayesian method with posterior effects
post.all <- data.table::fread("data/processed/Table_posterior_mutation_effects_RIL+RIAIL.csv")
post.riails <- data.table::fread("data/processed/Table_posterior_mutation_effects_RIAIL.csv")
post.rils <- data.table::fread("data/processed/Table_posterior_mutation_effects_RIL.csv")

# read in the marginal effects
mar.all <- data.table::fread("data/processed/Table_marg_effects_RIL+RIAIL.csv") %>%
  dplyr::rename_with(~paste0(.,"_mar"))
mar.riails <- data.table::fread("data/processed/Table_marg_effects_RIAIL.csv") %>%
  dplyr::rename_with(~paste0(.,"_mar"))
mar.rils <- data.table::fread("data/processed/Table_marg_effects_RIL.csv") %>%
  dplyr::rename_with(~paste0(.,"_mar"))

# join the the two methods together by variant - assuming same order.
join.all <- post.all %>%
  dplyr::bind_cols(mar.all) %>%
  dplyr::bind_cols(vars) %>%
  dplyr::bind_cols(af)

join.riails <- post.riails %>%
  dplyr::bind_cols(mar.riails) %>%
  dplyr::bind_cols(vars) %>%
  dplyr::bind_cols(af)

join.rils <- post.rils %>%
  dplyr::bind_cols(mar.rils) %>%
  dplyr::bind_cols(vars) %>%
  dplyr::bind_cols(af)

# calcualte corr for marginal and mutational effects
mut_cor_test <- cor.test(x=join.all$mean, y=join.all$mut_af)
mar_cor_test <- cor.test(x=join.all$mean_mar, y=join.all$mut_af)
mut_corr <- round(cor(x=join.all$mean, y=join.all$mut_af), digits = 2)
mar_corr <- round(cor(x=join.all$mean_mar, y=join.all$mut_af), digits = 3)
#------------------------------------------------------------------------------#
# Plot mutant allele frequency by marginal and Bayesian posterior mutational effect
#------------------------------------------------------------------------------#
leg.df <- tibble::tibble(type = rep(c("indel", "snp"), 2),
                         parent = rep(c("MA530", "MA563"), 2),
                         x = c(0,2,5,5),
                         y = c(0,2,5,5))
legend.p <- ggplot(leg.df) +
  aes(x = x, y = y, color = parent, shape = type) +
  geom_point(size = 3) +
  scale_color_manual(values = c("MA530" = "#0072B2", "MA563" = "#D55E00")) +
  scale_shape_manual(values = c("indel" = 16, "snp" = 17)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.text = element_text(size = 9),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent"))
legend <- cowplot::get_legend(legend.p)

sfig4a <- ggplot(join.all) +
  geom_point(aes(x = mean_mar, y = mut_af, fill = parent, shape = type), stroke = 0.25, size = 2) +
  scale_fill_manual(values = c("MA530" = "#0072B2", "MA563" = "#D55E00")) +
  scale_shape_manual(values = c("indel" = 21, "snp" = 24)) +
  annotate(geom = "text", x = -0.475, y = 0.6, label = glue::glue("r = {mar_corr}"), size = 3) +
  #labs(x = expression(paste('Marginal effect (',italic('u'),')')), y = 'Mutant allele freq.') +
  labs(x="Raw difference (<i>u</i><sub>RAW</sub>)", y = 'Mutant allele freq.') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.title.x = ggtext::element_markdown(),
        legend.position = "none") # legend.background = element_rect(fill = "transparent")
sfig4a

sfig4b <- ggplot(join.all) +
  geom_point(aes(x = mean, y = mut_af, fill = parent, shape = type), stroke = 0.25, size = 2) +
  scale_fill_manual(values = c("MA530" = "#0072B2", "MA563" = "#D55E00")) +
  scale_shape_manual(values = c("indel" = 21, "snp" = 24)) +
  annotate(geom = "text", x = -0.125, y = 0.6, label = glue::glue("r = {mut_corr}"), size = 3) +
  labs(x = expression(paste('Mutational effect (',italic('u'),')')), y = 'Mutant allele freq.') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.position = "none") # legend.background = element_rect(fill = "transparent")
sfig4b

sfig4ab <- cowplot::plot_grid(sfig4a, sfig4b, labels = c("a", "b"), label_size = 12, ncol = 2, align = "vh", axis = "tblr")

sfig4 <- cowplot::plot_grid(sfig4ab, legend, labels = c("", ""), ncol = 2, rel_widths = c(1, 0.2))
sfig4

cowplot::ggsave2(sfig4, filename = "figures/supplemental_figure4.png", width = 6.25, height = 2.25, dpi = 350)
cowplot::ggsave2(sfig4, filename = "figures/supplemental_figure4.pdf", width = 6.25, height = 2.25)
