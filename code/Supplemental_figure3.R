library(tidyverse)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# source functions
source("code/LD_functions.R")

#------------------------------------------------------------------------------#
# Supplemental Figure 3 -Average pairwise LD (r2, all lines, RIAILs, RILs).
# (A) Intrachromosomal (B) Interchromosomal
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Calculate intrachromosomal LD from PLINK for all sets
#------------------------------------------------------------------------------#
r2_all_intra <- intra_ld(ld = "R2", ld_file = "data/processed/mut_all_allele_ld_r2.ld", plot = F) %>%
  dplyr::mutate(type = "RIAIL + RIL")
r2_ril_intra <- intra_ld(ld = "R2", ld_file = "data/processed/mut_ril_allele_ld_r2.ld", plot = F) %>%
  dplyr::mutate(type = "RIL")
r2_riail_intra <- intra_ld(ld = "R2", ld_file = "data/processed/mut_riail_allele_ld_r2.ld", plot = F) %>%
  dplyr::mutate(type = "RIAIL")

# merge them
r2_intra <- rbind(r2_all_intra, r2_ril_intra, r2_riail_intra)

# plot them
sfig3a <- ggplot(r2_intra) +
  aes(x = factor(chrom, levels = c("I", "II", "III", "IV", "V", "X")), y = mean_ld, fill = factor(type, levels = c("RIAIL + RIL", "RIL", "RIAIL"))) +
  geom_bar(stat = "identity", width = 0.5, position = "dodge") +
  scale_fill_manual(values = c("grey40", "#00A875", "#00B9F1")) +
  geom_errorbar(aes(ymin=mean_ld - sem_ld, ymax=mean_ld + sem_ld), width=.25, position=position_dodge(.5), size = 0.25) +
  #geom_col(fill = "dark blue", gr) +
  #geom_hline(yintercept = unique(inter_chr_ld$exp_ld), linetype = "dashed", color = "light blue") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9)) +
  labs(x = "Chromosome", y = expression(paste('Mean intra-chromosomal LD (', r^2, ')')), fill = "") +
  theme(legend.position = "none")
sfig3a

#------------------------------------------------------------------------------#
# Calculate interchromosomal LD from PLINK for all sets
#------------------------------------------------------------------------------#
r2_all_inter <- inter_ld(ld = "R2", ld_file = "data/processed/mut_all_allele_ld_r2.ld") %>%
  dplyr::mutate(type = "RIAIL + RIL")
r2_ril_inter <- inter_ld(ld = "R2", ld_file = "data/processed/mut_ril_allele_ld_r2.ld") %>%
  dplyr::mutate(type = "RIL")
r2_riail_inter <- inter_ld(ld = "R2", ld_file = "data/processed/mut_riail_allele_ld_r2.ld") %>%
  dplyr::mutate(type = "RIAIL")


# merge them
r2_inter <- rbind(r2_all_inter, r2_ril_inter, r2_riail_inter)

# plot them
sfig3b<- ggplot(r2_inter) +
  aes(x = factor(chrom_pair, levels = c(r2_inter$chrom_pair[1:15])), y = mean_ic_ld, fill = factor(type, levels = c("RIAIL + RIL", "RIL", "RIAIL"))) +
  geom_bar(stat = "identity", width = 0.5, position = "dodge") +
  scale_fill_manual(values = c("grey40", "#00A875", "#00B9F1")) +
  geom_errorbar(aes(ymin=mean_ic_ld - sem_ic_ld, ymax=mean_ic_ld + sem_ic_ld), width=.25, position=position_dodge(.5), size = 0.25) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9)) +
  labs(x = "Chromosome comparison", y = expression(paste('Mean inter-chromosomal LD (', r^2, ')')), fill = "") # expression('Mean inter-chromosomal\nlinkage disequillibrium r'^'2')
sfig3b

#------------------------------------------------------------------------------#
# Make Supplemental Figure 3 -  
# Average pairwise LD (r2, all lines). (A) Intrachromosomal (B) Interchromosomal
#------------------------------------------------------------------------------#
sfig3 <- cowplot::plot_grid(sfig3a, sfig3b, ncol = 2, rel_widths = c(.8, 2), labels = c("a", "b"), align = "vh", axis = "bl")
cowplot::ggsave2(sfig3, filename = "figures/Supplemental_figure3.png", width = 10, height = 4)
cowplot::ggsave2(sfig3, filename = "figures/Supplemental_figure3.pdf", width = 10, height = 4)

