library(tidyverse)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# get the variant ids and the parents for plotting using hard coded parent_genos.csv file
vars <- data.table::fread("data/raw/parent_genos.csv") %>%
  tidyr::separate(loci, into = c("chrom", "pos", "type"), sep = "_", remove = T) %>%
  dplyr::mutate(parent = ifelse(MA530 == 1, "MA530", "MA563"),
                type = ifelse(type == "snp", "snp", "indel"),
                pos = as.numeric(pos))

# read in the bayesian method with posterior effects
post.all <- data.table::fread("data/processed/Table_posterior_mutation_effects_RIL+RIAIL.csv") 
post.riails <- data.table::fread("data/processed/Table_posterior_mutation_effects_RIAIL.csv") 
post.rils <- data.table::fread("data/processed/Table_posterior_mutation_effects_RIL.csv")

# read in the marginal effects
mar.all <- data.table::fread("data/processed/Table_marg_effects_RIL+RIAIL.csv") %>%
  dplyr::rename_with(~paste0(.,"_mar"))
mar.riails <- data.table::fread("data/processed/Table_marg_effects_RIAIL.csv") %>%
  dplyr::rename_with(~paste0(.,"_mar"))
mar.rils <- data.table::fread("data/processed/Table_marg_effects_RIL.csv") %>%
  dplyr::rename_with(~paste0(.,"_mar"))

# join the the two methods together by variant - assuming same order.
join.all <- post.all %>%
  dplyr::bind_cols(mar.all) %>%
  dplyr::bind_cols(vars)

join.riails <- post.riails %>%
  dplyr::bind_cols(mar.riails) %>%
  dplyr::bind_cols(vars)

join.rils <- post.rils %>%
  dplyr::bind_cols(mar.rils) %>%
  dplyr::bind_cols(vars)

# look at correlation between approaches - pearson
corr.all2 <- round(cor(x=join.all$mean_mar, y=join.all$mean), digits = 3)
corr.riails <- round(cor(y=join.riails$mean, x=join.riails$mean_mar), digits = 3)
corr.rils <- round(cor(y=join.rils$mean, x=join.rils$mean_mar), digits = 3)

#------------------------------------------------------------------------------#
# Plots
#------------------------------------------------------------------------------#
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
        legend.position = "top")
legend <- cowplot::get_legend(legend.p)

fig4a <- ggplot(join.all) +
  scale_color_manual(values = c("MA530" = "#0072B2", "MA563" = "#D55E00")) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_linerange(aes(x = mean_mar, y = mean,
                     xmin = `CI2.5%_mar`, xmax = `CI97.5%_mar`, color = parent),
                 linewidth = 0.25, alpha = 0.25) +
  geom_linerange(aes(x = mean_mar, y = mean,
                     ymin = `CI2.5%`, ymax = `CI97.5`, color = parent),
                 linewidth = 0.25, alpha = 0.25) +
  ggnewscale::new_scale_color() +
  scale_color_manual(values = c("MA530" = "black", "MA563" = "black")) +
  geom_point(aes(x = mean_mar, y = mean, color = parent, fill = parent, shape = type), size = 1.5, stroke = 0.25) +
  #geom_point(size = 1.5) +
  scale_fill_manual(values = c("MA530" = "#0072B2", "MA563" = "#D55E00")) +
  scale_shape_manual(values = c("indel" = 21, "snp" = 24)) +
  annotate(geom = "text", x = -0.45, y = 0.25, label = glue::glue("r = {corr.all}"), size = 3) +
  labs(x = "", y = expression(paste('Mutational effect (',italic('u'),')'))) +
  scale_x_continuous(breaks = c(-0.8, -0.4, 0, 0.4)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        axis.text = element_text(size = 8)) # legend.background = element_rect(fill = "transparent")
fig4a
fig4b <- ggplot(join.rils) +
  scale_color_manual(values = c("MA530" = "#0072B2", "MA563" = "#D55E00")) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_linerange(aes(x = mean_mar, y = mean,
                     xmin = `CI2.5%_mar`, xmax = `CI97.5%_mar`, color = parent),
                 linewidth = 0.25, alpha = 0.25) +
  geom_linerange(aes(x = mean_mar, y = mean,
                      ymin = `CI2.5%`, ymax = `CI97.5`, color = parent),
                  linewidth = 0.25, alpha = 0.25) +
  ggnewscale::new_scale_color() +
  scale_color_manual(values = c("MA530" = "black", "MA563" = "black")) +
  geom_point(aes(x = mean_mar, y = mean, color = parent, fill = parent, shape = type), size = 1.5, stroke = 0.25) +
  #geom_point(size = 1.5) +
  scale_fill_manual(values = c("MA530" = "#0072B2", "MA563" = "#D55E00")) +
  scale_shape_manual(values = c("indel" = 21, "snp" = 24)) +
  annotate(geom = "text", x = -0.375, y = 0.275, label = glue::glue("r = {corr.rils}"), size = 3) +
  labs(x = expression(paste('Marginal effect (',italic('u'),')')), y = "") +
  scale_x_continuous(breaks = c(-0.6, -0.3, 0, 0.3, 0.6)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        axis.text = element_text(size = 8))
fig4b

fig4c <- ggplot(join.riails) +
  scale_color_manual(values = c("MA530" = "#0072B2", "MA563" = "#D55E00")) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_linerange(aes(x = mean_mar, y = mean,
                     xmin = `CI2.5%_mar`, xmax = `CI97.5%_mar`, color = parent),
                 linewidth = 0.25, alpha = 0.25) +
  geom_linerange(aes(x = mean_mar, y = mean,
                     ymin = `CI2.5%`, ymax = `CI97.5`, color = parent),
                 linewidth = 0.25, alpha = 0.25) +
  ggnewscale::new_scale_color() +
  scale_color_manual(values = c("MA530" = "black", "MA563" = "black")) +
  geom_point(aes(x = mean_mar, y = mean, color = parent, fill = parent, shape = type), size = 1.5, stroke = 0.25) +
  #geom_point(size = 1.5) +
  scale_fill_manual(values = c("MA530" = "#0072B2", "MA563" = "#D55E00")) +
  scale_shape_manual(values = c("indel" = 21, "snp" = 24)) +
  annotate(geom = "text", x = -0.5, y = 0.525, label = glue::glue("r = {corr.riails}"), size = 3) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        axis.text = element_text(size = 8)) # legend.background = element_rect(fill = "transparent")

fig4abc <- cowplot::plot_grid(fig4a, fig4b, fig4c, labels = c("a", "b", "c"),
                              ncol = 3, align = "vh", axis = "tbrl", label_size = 12)
fig4 <- cowplot::plot_grid(fig4abc, legend, labels = c("", ""), ncol = 1, rel_heights = c(1,0.1), label_size = 12)
fig4

cowplot::ggsave2(fig4, filename = "figures/figure4.png", width = 6.25, height = 2.25)
cowplot::ggsave2(fig4, filename = "figures/figure4.pdf", width = 6.25, height = 2.25)
