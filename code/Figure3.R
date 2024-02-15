library(tidyverse)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# get the p matrix data
p <- data.table::fread("data/processed/04_DFE_prob_genotypes.csv")

