library(tidyverse)
library(scales)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# get the variant ids and the parents for plotting using hard coded parent_genos.csv file
vars <- data.table::fread("data/raw/parent_genos.csv") %>%
  tidyr::separate(loci, into = c("chrom", "pos", "type"), sep = "_", remove = T) %>%
  dplyr::mutate(parent = ifelse(MA530 == 1, "MA530", "MA563"),
                type = ifelse(type == "snp", "snp", "indel"),
                pos = as.numeric(pos))

# ws276 gff3 genome size and bins
edge <- tibble::tibble(chrom = c("I", "II", "III", "IV", "V", "X", "I", "II", "III", "IV", "V", "X"),
                       center = c(0, 0, 0, 0, 0, 0, 15072434, 15279421 ,13783801 ,17493829, 20924180, 17718942),
                       alpha = "clear")  

# load the mutant haplotype table to lookup variants in haplotypes
# shape this table to estimate halpotype width and center
haplo.table1 <- data.table::fread("/projects/b1059/projects/Tim/manuscript_DFE/data/processed/Mutant_haplotypes_table.csv") %>%
  dplyr::mutate(n.var = stringr::str_count(Loci, pattern = " ")+1,
                hap.locus.id = stringr::str_extract(`Mutant Haplotype`, pattern = ".*(?=_)"))

haplo.table2 <- haplo.table1 %>%
  dplyr::group_by(hap.locus.id) %>%
  dplyr::mutate(n.mut.hap = n()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(hap.locus.id, .keep_all = T) %>%
  dplyr::select(3:5) 

haplo.table3 <- tibble::tibble(hap.locus.id = rep(haplo.table2$hap.locus.id, times = haplo.table2$n.var)) %>%
  dplyr::left_join(haplo.table2) %>%
  dplyr::bind_cols(vars) %>%
  dplyr::group_by(hap.locus.id) %>%
  dplyr::mutate(start = min(pos),
                stop = max(pos),
                center = (stop - start)/2 + start,
                type2 = case_when(n.var == 1 ~ type,
                                 n.var > 1 ~ "multi")) %>%
  dplyr::ungroup()

haplo.table4 <- haplo.table3 %>%
  dplyr::distinct(hap.locus.id, .keep_all = T) %>%
  dplyr::mutate(mut.hap.parent = case_when(n.mut.hap == 1 ~ parent,
                                           n.mut.hap == 2 ~ "multi",
                                           TRUE ~ "wtf?")) %>%
  dplyr::select(-5:-9) %>%
  dplyr::rename(type = type2) %>%
  dplyr::group_by(n.mut.hap) %>%
  slice(rep(1:n(), if_else(first(n.mut.hap) == 2, 2, 1))) %>% # add a row for each 2 mut hap
  dplyr::ungroup() %>%
  dplyr::arrange(chrom, start) %>%
  dplyr::group_by(hap.locus.id) %>%
  dplyr::mutate(mut.hap.id = 1:n(),
                parent = case_when(mut.hap.parent == "multi" & mut.hap.id == 1 ~ "MA530",
                                   mut.hap.parent == "multi" & mut.hap.id == 2 ~ "MA563",
                                   mut.hap.parent != "multi" ~ mut.hap.parent))

# export the haplotype data for supplmental figure 5
rio::export(haplo.table3, file = "data/processed/mutatnt_haplotype_supfig5.csv")
  
# load the data from JZ and bind var ids
me.all <- data.table::fread("data/processed/Table_posterior_haplotype_effects_RIL+RIAIL.csv") %>%
  dplyr::bind_cols(haplo.table4) %>%
  dplyr::full_join(edge) %>%
  dplyr::mutate(alpha = ifelse(is.na(alpha), "dark", alpha),
                mean = ifelse(is.na(mean), 0, mean),
                `CI2.5%` = ifelse(is.na(`CI2.5%`), 0, `CI2.5%`),
                `CI97.5%` = ifelse(is.na(`CI97.5%`), 0, `CI97.5%`),
                type = ifelse(is.na(type), "snp", type),
                parent = ifelse(is.na(parent), "MA530", parent))
  

me.ril <- data.table::fread("data/processed/Table_posterior_haplotype_effects_RIL.csv")  %>%
  dplyr::bind_cols(haplo.table4) %>%
  dplyr::full_join(edge) %>%
  dplyr::mutate(alpha = ifelse(is.na(alpha), "dark", alpha),
                mean = ifelse(is.na(mean), 0, mean),
                `CI2.5%` = ifelse(is.na(`CI2.5%`), 0, `CI2.5%`),
                `CI97.5%` = ifelse(is.na(`CI97.5%`), 0, `CI97.5%`),
                type = ifelse(is.na(type), "snp", type),
                parent = ifelse(is.na(parent), "MA530", parent))

me.riail <- data.table::fread("data/processed/Table_posterior_haplotype_effects_RIAIL.csv")  %>%
  dplyr::bind_cols(haplo.table4) %>%
  dplyr::full_join(edge) %>%
  dplyr::mutate(alpha = ifelse(is.na(alpha), "dark", alpha),
                mean = ifelse(is.na(mean), 0, mean),
                `CI2.5%` = ifelse(is.na(`CI2.5%`), 0, `CI2.5%`),
                `CI97.5%` = ifelse(is.na(`CI97.5%`), 0, `CI97.5%`),
                type = ifelse(is.na(type), "snp", type),
                parent = ifelse(is.na(parent), "MA530", parent))


#==============================================================================#
# Plot additive effect by genomic position
#==============================================================================#
leg.df <- tibble::tibble(type = rep(c("indel", "snp", "multi"), 2),
                         parent = rep(c("MA530", "MA563", "MA563"), 2),
                         x = c(0,2,5,5,3,3),
                         y = c(0,2,5,5,3,3))
legend.p <- ggplot(leg.df) +
  aes(x = x, y = y, color = parent, shape = factor(type, levels = c("indel", "snp", "multi"))) +
  geom_point(size = 3) +
  scale_color_manual(values = c("MA530" = "#0072B2", "MA563" = "#D55E00")) +
  scale_fill_discrete(na.translate=FALSE) +
  scale_fill_manual(values = c("MA530" = "#0072B2", "MA563" = "#D55E00")) +
  scale_shape_manual(values = c("multi" = 15, "indel" = 16, "snp" = 17)) +
  theme_bw() +
  labs(shape = "Type", color = "Parent") +
  theme(panel.grid = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "top",
        legend.background = element_rect(fill = "transparent"))
legend <- cowplot::get_legend(legend.p)

fig6a <- ggplot(me.all) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = 2, color = "darkgrey") +
  geom_linerange(aes(x = center/(10^6), y = mean,
                     ymin = `CI2.5%`, ymax = `CI97.5%`, color = parent),
                 linewidth = 0.5, alpha = 0.25) +
  geom_linerange(aes(x = center/(10^6), y = mean,
                     xmin = start/(10^6), xmax = stop/(10^6), color = parent),
                 linewidth = 3.25, alpha = 0.25) +
  geom_point(aes(x = center/(10^6), y = mean, fill = parent, alpha = alpha, shape = type), size = 1.5, stroke = 0.25) +
  scale_alpha_discrete(range = c(0,1), guide = "none") +
  scale_color_discrete(na.translate=FALSE) +
  scale_color_manual(values = c("MA530" = "#0072B2", "MA563" = "#D55E00")) +
  scale_fill_discrete(na.translate=FALSE) +
  scale_fill_manual(values = c("MA530" = "#0072B2", "MA563" = "#D55E00")) +
  scale_shape_manual(values = c("multi" = 22, "indel" = 21, "snp" = 24)) +
  theme_bw() +
  facet_grid(~chrom, scales = "free_x") +
  labs(y = expression(paste('RI(AI)Ls (',italic('\u0394W'),')')),
       x = "Genome postion (Mb)",
       color = "Parent",
       shape = "Type") +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10, face = "bold"))
fig6a

fig6b <- ggplot(me.ril) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = 2, color = "darkgrey") +
  geom_linerange(aes(x = center/(10^6), y = mean,
                     ymin = `CI2.5%`, ymax = `CI97.5%`, color = parent),
                 linewidth = 0.5, alpha = 0.25) +
  geom_linerange(aes(x = center/(10^6), y = mean,
                     xmin = start/(10^6), xmax = stop/(10^6), color = parent),
                 linewidth = 3.25, alpha = 0.25) +
  geom_point(aes(x = center/(10^6), y = mean, fill = parent, alpha = alpha, shape = type), size = 1.5, stroke = 0.25) +
  scale_alpha_discrete(range = c(0,1), guide = "none") +
  scale_color_discrete(na.translate=FALSE) +
  scale_color_manual(values = c("MA530" = "#0072B2", "MA563" = "#D55E00")) +
  scale_fill_discrete(na.translate=FALSE) +
  scale_fill_manual(values = c("MA530" = "#0072B2", "MA563" = "#D55E00")) +
  scale_shape_manual(values = c("multi" = 22, "indel" = 21, "snp" = 24)) +
  theme_bw() +
  facet_grid(~chrom, scales = "free_x") +
  labs(y = expression(paste('RILs (',italic('\u0394W'),')')),
       x = "Genome postion (Mb)",
       color = "Parent",
       shape = "Type") +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())
fig6b

fig6c <- ggplot(me.riail) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = 2, color = "darkgrey") +
  geom_linerange(aes(x = center/(10^6), y = mean,
                     ymin = `CI2.5%`, ymax = `CI97.5%`, color = parent),
                 linewidth = 0.5, alpha = 0.25) +
  geom_linerange(aes(x = center/(10^6), y = mean,
                     xmin = start/(10^6), xmax = stop/(10^6), color = parent),
                 linewidth = 3.25, alpha = 0.25) +
  geom_point(aes(x = center/(10^6), y = mean, fill = parent, alpha = alpha, shape = type), size = 1.5, stroke = 0.25) +
  scale_alpha_discrete(range = c(0,1), guide = "none") +
  scale_color_discrete(na.translate=FALSE) +
  scale_color_manual(values = c("MA530" = "#0072B2", "MA563" = "#D55E00")) +
  scale_fill_discrete(na.translate=FALSE) +
  scale_fill_manual(values = c("MA530" = "#0072B2", "MA563" = "#D55E00")) +
  scale_shape_manual(values = c("multi" = 22, "indel" = 21, "snp" = 24)) +
  theme_bw() +
  facet_grid(~chrom, scales = "free_x") +
  labs(y = expression(paste('RIAILs (',italic('\u0394W'),')')),
       x = "Genome postion (Mb)",
       color = "Parent",
       shape = "Type") +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        axis.text = element_text(size = 8),
        strip.background = element_blank(),
        strip.text.x = element_blank())
fig6c

fig6abc <- cowplot::plot_grid(fig6a, fig6b, fig6c, labels = c("a", "b", "c"), ncol = 1, align = "hv", axis = "tblr")
fig6 <- cowplot::plot_grid(fig6abc, legend, labels = c("", ""), ncol = 1, rel_heights = c(1, .05))

# save it
cowplot::ggsave2(fig6, filename = "figures/figure6.png", width = 7, height = 5)
cowplot::ggsave2(fig6, filename = "figures/figure6.pdf", width = 7, height = 5)


