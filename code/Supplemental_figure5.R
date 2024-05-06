library(tidyverse)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# get the haplotype data
haplo <- data.table::fread("data/processed/mutatnt_haplotype_supfig5.csv")

# pull out the genotypes for the large effect locus on chr III
locus_pos <- haplo %>%
  dplyr::filter(hap.locus.id == 38) %>%
  dplyr::mutate(var_id = paste(chrom, pos, sep = "_")) %>%
  dplyr::pull(var_id)

# get imputed genotypes
geno <- data.table::fread("data/processed/03_DFE_prob_genotypes.csv")

geno_proc <- geno %>%
  dplyr::mutate(across(8:524, \(x) round(x, digits = 0))) %>% # weird new syntax, I hate it.
  dplyr::mutate(locus = paste(chrom, pos, sep = "_"), .before = full_var_id) %>%
  dplyr::filter(locus %in% locus_pos) %>%
  data.table::transpose(., keep.names = "name") %>%
  dplyr::slice(7:n()) %>%
  tidyr::unite(col = geno, V1:V13, sep = "")

neg.563 <- geno_proc %>%
  dplyr::filter(name == "MA563") %>%
  dplyr::pull(geno)

pos.530 <- geno_proc %>%
  dplyr::filter(name == "MA530") %>%
  dplyr::pull(geno)

# count number of diffs relative to each parental haplotype
geno_proc2 <- geno_proc %>%
  dplyr::rowwise() %>%
  dplyr::mutate(diff.from.563 = sum(unlist(strsplit(geno, "")) != unlist(strsplit(neg.563, ""))),
                diff.from.530 = sum(unlist(strsplit(geno, "")) != unlist(strsplit(pos.530, "")))) %>%
  dplyr::ungroup()

# get lines that match haplotypes
neg.lines <- geno_proc %>%
  dplyr::filter(!(name %in% c("MA530", "MA563")) & geno == neg.563) %>%
  dplyr::pull(name)

pos.lines <- geno_proc %>%
  dplyr::filter(!(name %in% c("MA530", "MA563")) & geno == pos.530) %>%
  dplyr::pull(name)

# get the full replicate phenotype data
pheno <- data.table::fread("data/processed/01_DFE_phenotypes.csv")

# subset the phneo data to the replicates we want to plot and plot
ancestor <- pheno %>%
  dplyr::filter(grepl(full_id, pattern = "N2_G0_ANC"))

supfig5a <- ggplot(ancestor) +
  aes(x = ln_ci) +
  geom_histogram(bins = 30, color = "black") + 
  theme_bw() +
  labs(y = "Count",
       x = "") + #expression(paste('Competitive index (',italic('W'),')'))
  scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 1)) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8))
supfig5a

# Do the neg.haplotype group
neg.pheno <- pheno %>%
  dplyr::filter(full_id %in% neg.lines)

supfig5b <- ggplot(neg.pheno) +
  aes(x = ln_ci) +
  geom_histogram(bins = 30, fill = "#D55E00", color = "black") +
  theme_bw() +
  labs(y = "",
       x = expression(paste('Competitive index (',italic('W'),')'))) +
  scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 1)) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8))
supfig5b

# Now do positive 530
pos.pheno <- pheno %>%
  dplyr::filter(full_id %in% pos.lines)

supfig5c <- ggplot(pos.pheno) +
  aes(x = ln_ci) +
  geom_histogram(bins = 30, fill = "#0072B2", color = "black") +
  theme_bw() +
  labs(y = "",
       x = "") + #expression(paste('Competitive index (',italic('W'),')'))
  scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 1)) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8))
supfig5c

# Put them together
supfig5 <- cowplot::plot_grid(supfig5a, supfig5b, supfig5c, labels = c("a", "b", "c"), 
                              label_size = 12, vjust = 1, ncol = 3, align = "vh", axis = "tblr")
supfig5

cowplot::ggsave2(supfig5, filename = "figures/supplemental_figure5.png", width = 6.25, height = 2.25)
cowplot::ggsave2(supfig5, filename = "figures/supplemental_figure5.pdf", width = 6.25, height = 2.25)
