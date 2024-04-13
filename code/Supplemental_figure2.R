library(tidyverse)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# source functions
source("code/LD_functions.R")

#------------------------------------------------------------------------------#
# Supplemental Figure 2 - Pairwise LD (D' signed). 
# (A) RIAILs + RILs; (B) RIAILs; (C) RILs
#------------------------------------------------------------------------------#
# make elements of fig 1
sfig2a <- LDchroms2(ld = "DP", ld_file = "data/processed/mut_all_allele_ld_r.ld")
sfig2b <- LDchroms2(ld = "DP", ld_file = "data/processed/mut_ril_allele_ld_r.ld")
sfig2c <- LDchroms2(ld = "DP", ld_file = "data/processed/mut_riail_allele_ld_r.ld")

# put them together
sfig2 <- cowplot::plot_grid(sfig2a, sfig2b, sfig2c, labels = c("a", "b", "c"), ncol = 1, align = "lb")

# save figure 1 
cowplot::ggsave2(sfig2, filename = "figures/supplemental_figure2.png", width = 27, height = 36)
cowplot::ggsave2(sfig2, filename = "figures/supplemental_figure2.pdf", width = 27, height = 36)
