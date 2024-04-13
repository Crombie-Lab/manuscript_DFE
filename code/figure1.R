library(tidyverse)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# source functions
source("code/LD_functions.R")

#------------------------------------------------------------------------------#
# Figure 1 - Pairwise LD (r2). Show mutations from different parents with
# different connecting lines (colors or dashes)
# (A) Correlation heat map of all lines; (B) RILs: (C) RIAILs
#------------------------------------------------------------------------------#
# make elements of fig 1
fig1a <- LDchroms2(ld = "R2", ld_file = "data/processed/mut_all_allele_ld_r2.ld")
fig1b <- LDchroms2(ld = "R2", ld_file = "data/processed/mut_ril_allele_ld_r2.ld")
fig1c <- LDchroms2(ld = "R2", ld_file = "data/processed/mut_riail_allele_ld_r2.ld")

# put them together
fig1 <- cowplot::plot_grid(fig1a, fig1b, fig1c, labels = c("a", "b", "c"), ncol = 1, align = "lb")

# save figure 1 
cowplot::ggsave2(fig1, filename = "figures/figure1.png", width = 27, height = 36)
cowplot::ggsave2(fig1, filename = "figures/figure1.pdf", width = 27, height = 36)

