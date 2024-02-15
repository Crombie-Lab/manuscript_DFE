library(tidyverse)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# Process output from python code
# get the marginal effects estimated from differences between genotypes (method 1)
# and the estimated effects from the MCML Hierarchical Bayesian Analysis (method 2)
eff <- data.table::fread("data/processed/summary_bayesian_bootstrap.csv")

# get the variant list to rename V1 above
vars <- data.table::fread("data/processed/04_DFE_prob_genotypes.csv")



