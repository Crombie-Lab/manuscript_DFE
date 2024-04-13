library(tidyverse)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# assign genotypes based on parents for unimputed data
# load data for ploting heatmaps - calculated initially in LD.R
load("data/processed/ld.data.for.heatmaps.rda")

#==============================================================================#
# Plot heatmaps for intra chromosomal LD in r, R^2, and D with ggplot
#==============================================================================#
# get parents genotypes
parent_genos <- data.table::fread("data/raw/parent_genos.csv")
pg.2 <- parent_genos %>%
  dplyr::mutate(loci = stringr::str_replace(loci, pattern = regex("_[^_]*$"), replacement = ""),
                label = ifelse(MA530 == 1, "MA530", "MA563"))

# load the function to make the LD plots
source("code/functions.R")

# plot the intra chrom LD in D
p1.D <- LDchroms(data = ld.data.for.heatmaps$ld.all.D, parents = pg.2, type = 'D')
p2.D <- LDchroms(data = ld.data.for.heatmaps$ld.RILs.D, parents = pg.2, type = 'D')
p3.D <- LDchroms(data = ld.data.for.heatmaps$ld.RIAILs.D, parents = pg.2, type = 'D')

# plot the intra chrom LD in r
p1 <- LDchroms(data = ld.data.for.heatmaps$ld.all.r, parents = pg.2, type = 'r')
p2 <- LDchroms(data = ld.data.for.heatmaps$ld.RILs.r, parents = pg.2, type = 'r')
p3 <- LDchroms(data = ld.data.for.heatmaps$ld.RIAILs.r, parents = pg.2, type = 'r')

# plot the intra chrom LD in r2
p1r2 <- LDchroms(data = ld.data.for.heatmaps$ld.all.r2, parents = pg.2, type = 'R2')
p2r2 <- LDchroms(data = ld.data.for.heatmaps$ld.RILs.r2, parents = pg.2, type = 'R2')
p3r2 <- LDchroms(data = ld.data.for.heatmaps$ld.RIAILs.r2, parents = pg.2, type = 'R2')

# put them together
figure1 <- cowplot::plot_grid(p1.D, p1r2, ncol = 1, labels = c("a", "b"))

# save figure 1 
cowplot::ggsave2(figure1, filename = "figures/figure1.png", width = 27, height = 24)

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
# 
