library(tidyverse)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# source functions
source("code/LD_functions.R")

#------------------------------------------------------------------------------#
# Make files to run PLINK
#------------------------------------------------------------------------------#
make_plink_in(geno_path = "data/processed/02_DFE_unimputed_genotypes.csv",
              strain = "all",
              out_prefix = "dfe_all_mut") 

make_plink_in(geno_path = "data/processed/02_DFE_unimputed_genotypes.csv",
              strain = "ril",
              out_prefix = "dfe_ril_mut") 

make_plink_in(geno_path = "data/processed/02_DFE_unimputed_genotypes.csv",
              strain = "riail",
              out_prefix = "dfe_riail_mut")
#------------------------------------------------------------------------------#
# Run PLINK on QUEST with PLINK commands below
#------------------------------------------------------------------------------#
# # run it pairwise with a snp.txt list
# plink --file dfe_all_mut --allow-extra-chr --ld-window-r2 0 --r inter-chr dprime-signed --ld-snp-list snps.txt --out mut_all_allele_ld_r
# plink --file dfe_all_mut --allow-extra-chr --ld-window-r2 0 --r2 inter-chr --ld-snp-list snps.txt --out mut_all_allele_ld_r2
# 
# # run it for riails
# plink --file dfe_riail_mut --allow-extra-chr --ld-window-r2 0 --r inter-chr dprime-signed --ld-snp-list snps.txt --out mut_riail_allele_ld_r
# plink --file dfe_riail_mut --allow-extra-chr --ld-window-r2 0 --r2 inter-chr --ld-snp-list snps.txt --out mut_riail_allele_ld_r2
# 
# # run for rils
# plink --file dfe_ril_mut --allow-extra-chr --ld-window-r2 0 --r inter-chr dprime-signed --ld-snp-list snps.txt --out mut_ril_allele_ld_r
# plink --file dfe_ril_mut --allow-extra-chr --ld-window-r2 0 --r2 inter-chr --ld-snp-list snps.txt --out mut_ril_allele_ld_r2
