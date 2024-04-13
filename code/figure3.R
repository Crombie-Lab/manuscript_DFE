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
                       pos = c(0, 0, 0, 0, 0, 0, 15072434, 15279421 ,13783801 ,17493829, 20924180, 17718942),
                       alpha = "clear")  

# load the data from JZ and bind var ids
me.all <- data.table::fread("data/processed/Table_posterior_mutation_effects_RIL+RIAIL.csv") %>%
  dplyr::bind_cols(vars) %>%
  dplyr::full_join(edge) %>%
  dplyr::mutate(alpha = ifelse(is.na(alpha), "dark", alpha),
                mean = ifelse(is.na(mean), 0, mean),
                `CI2.5%` = ifelse(is.na(`CI2.5%`), 0, `CI2.5%`),
                CI97.5 = ifelse(is.na(CI97.5), 0, CI97.5),
                type = ifelse(is.na(type), "snp", type),
                parent = ifelse(is.na(parent), "MA530", parent))

me.riail <- data.table::fread("data/processed/Table_posterior_mutation_effects_RIAIL.csv") %>%
  dplyr::bind_cols(vars) %>%
  dplyr::full_join(edge) %>%
  dplyr::mutate(alpha = ifelse(is.na(alpha), "dark", alpha),
                mean = ifelse(is.na(mean), 0, mean),
                `CI2.5%` = ifelse(is.na(`CI2.5%`), 0, `CI2.5%`),
                CI97.5 = ifelse(is.na(CI97.5), 0, CI97.5),
                type = ifelse(is.na(type), "snp", type),
                parent = ifelse(is.na(parent), "MA530", parent))

me.ril <- data.table::fread("data/processed/Table_posterior_mutation_effects_RIL.csv") %>%
  dplyr::bind_cols(vars) %>%
  dplyr::full_join(edge) %>%
  dplyr::mutate(alpha = ifelse(is.na(alpha), "dark", alpha),
                mean = ifelse(is.na(mean), 0, mean),
                `CI2.5%` = ifelse(is.na(`CI2.5%`), 0, `CI2.5%`),
                CI97.5 = ifelse(is.na(CI97.5), 0, CI97.5),
                type = ifelse(is.na(type), "snp", type),
                parent = ifelse(is.na(parent), "MA530", parent))

#==============================================================================#
# Plot additive effect by genomic position
#==============================================================================#
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
  labs(shape = "Type", color = "Parent") +
  theme(panel.grid = element_blank(),
        legend.text = element_text(size = 9),
        legend.position = "bottom",
        legend.background = element_rect(fill = "transparent"))
legend <- cowplot::get_legend(legend.p)

fig3a <- ggplot(me.all) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = 2, color = "darkgrey") +
  geom_linerange(aes(x = pos/(10^6), y = mean,
                     ymin = `CI2.5%`, ymax = `CI97.5`, color = parent),
                 linewidth = 0.5, alpha = 0.25) +
  geom_point(aes(x = pos/(10^6), y = mean, fill = parent, alpha = alpha, shape = type), size = 1.5, stroke = 0.25) +
  scale_alpha_discrete(range = c(0,1), guide = "none") +
  scale_color_discrete(na.translate=FALSE) +
  scale_color_manual(values = c("MA530" = "#0072B2", "MA563" = "#D55E00")) +
  scale_fill_discrete(na.translate=FALSE) +
  scale_fill_manual(values = c("MA530" = "#0072B2", "MA563" = "#D55E00")) +
  scale_shape_manual(values = c("indel" = 21, "snp" = 24)) +
  theme_bw() +
  facet_grid(~chrom, scales = "free_x") +
  labs(y = expression(paste('RI(AI)Ls (',italic('u'),')')),
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
fig3a
fig3b <- ggplot(me.ril) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = 2, color = "darkgrey") +
  geom_linerange(aes(x = pos/(10^6), y = mean,
                     ymin = `CI2.5%`, ymax = `CI97.5`, color = parent),
                 linewidth = 0.5, alpha = 0.25) +
  geom_point(aes(x = pos/(10^6), y = mean, fill = parent, alpha = alpha, shape = type), size = 1.5, stroke = 0.25) +
  scale_alpha_discrete(range = c(0,1), guide = "none") +
  scale_color_discrete(na.translate=FALSE) +
  scale_color_manual(values = c("MA530" = "#0072B2", "MA563" = "#D55E00")) +
  scale_fill_discrete(na.translate=FALSE) +
  scale_fill_manual(values = c("MA530" = "#0072B2", "MA563" = "#D55E00")) +
  scale_shape_manual(values = c("indel" = 21, "snp" = 24)) +
  theme_bw() +
  facet_grid(~chrom, scales = "free_x") +
  labs(y = expression(paste('RILs (',italic('u'),')')),
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

fig3c <- ggplot(me.riail) +
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = 2, color = "darkgrey") +
  geom_linerange(aes(x = pos/(10^6), y = mean,
                     ymin = `CI2.5%`, ymax = `CI97.5`, color = parent),
                 linewidth = 0.5, alpha = 0.25) +
  geom_point(aes(x = pos/(10^6), y = mean, fill = parent, alpha = alpha, shape = type), size = 1.5,stroke = 0.25) +
  scale_alpha_discrete(range = c(0,1), guide = "none") +
  scale_color_discrete(na.translate=FALSE) +
  scale_color_manual(values = c("MA530" = "#0072B2", "MA563" = "#D55E00")) +
  scale_fill_discrete(na.translate=FALSE) +
  scale_fill_manual(values = c("MA530" = "#0072B2", "MA563" = "#D55E00")) +
  scale_shape_manual(values = c("indel" = 21, "snp" = 24)) +
  theme_bw() +
  facet_grid(~chrom, scales = "free_x") +
  labs(y = expression(paste('RIAILs (',italic('u'),')')),
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

fig3abc <- cowplot::plot_grid(fig3a, fig3b, fig3c, labels = c("a", "b", "c"), ncol = 1, align = "hv", axis = "tblr")
fig3 <- cowplot::plot_grid(fig3abc, legend, labels = c("", ""), ncol = 1, rel_heights = c(1, .05))
fig3
# save it
cowplot::ggsave2(fig3, filename = "figures/figure3.png", width = 7, height = 5) #8 x 5 (wh)
cowplot::ggsave2(fig3, filename = "figures/figure3.pdf", width = 7, height = 5)

