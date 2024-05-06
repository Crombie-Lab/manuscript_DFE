library(tidyverse)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load data
df_all <- data.table::fread("data/processed/Table_posterior_haplotype_effects_RIL+RIAIL.csv")
df_ril <- data.table::fread("data/processed/Table_posterior_haplotype_effects_RIL.csv")
df_riails <- data.table::fread("data/processed/Table_posterior_haplotype_effects_RIAIL.csv")

# get max x range for all dfs
ax <- c(max(df_all$mean), min(df_all$mean))
bx <- c(max(df_ril$mean), min(df_ril$mean))
cx <- c(max(df_riails$mean), min(df_riails$mean))

# plot the histogram with freq on y and Delta (w) on x
a <- ggplot(df_all) +
  aes(x = mean, y=stat(count)/sum(stat(count))) +
  geom_histogram(binwidth = 0.009, color = "black", linewidth = 0.25) +
  geom_vline(xintercept = mean(df_all$mean), linetype = 2, color = "darkred", linewidth = 0.25) +
  annotate(geom = "label", x = mean(df_all$mean), y = 0.35, size = 3,
           label = glue::glue("{round(mean(df_all$mean), digits=3)}"),
           color = "darkred", size = 3) +
  labs(y = "Freq.", x = expression(paste('RI(AI)Ls (',italic('u'),')'))) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1),
                     limits = c(-0.825, 0.225),
                     breaks = c(-0.8, -0.6, -0.4, -0.2, 0, 0.2)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1),
                     limits = c(0, 0.425)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8))
a

b <- ggplot(df_ril) +
  aes(x = mean, y=stat(count)/sum(stat(count))) +
  geom_histogram(binwidth = 0.009, color = "black", linewidth = 0.25) +
  geom_vline(xintercept = mean(df_ril$mean), linetype = 2, color = "darkred", linewidth = 0.25) +
  annotate(geom = "label", x = mean(df_ril$mean), y = 0.35, size = 3,
           label = glue::glue("{round(mean(df_ril$mean), digits=3)}"),
           color = "darkred", size = 3) +
  labs(y = "", x = expression(paste('RILs (',italic('u'),')'))) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1),
                     limits = c(-0.825, 0.225),
                     breaks = c(-0.8, -0.6, -0.4, -0.2, 0, 0.2)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1),
                     limits = c(0, 0.425)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8))
b

c <- ggplot(df_riails) +
  aes(x = mean, y=stat(count)/sum(stat(count))) +
  geom_histogram(binwidth = 0.009, color = "black", linewidth = 0.25) +
  geom_vline(xintercept = mean(df_riails$mean), linetype = 2, color = "darkred", linewidth = 0.25) +
  annotate(geom = "label", x = mean(df_riails$mean), y = 0.35, size = 3,
           label = glue::glue("{round(mean(df_riails$mean), digits=3)}"),
           color = "darkred", size = 3) +
  labs(y = "", x = expression(paste('RIAILs (',italic('u'),')'))) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1),
                     limits = c(-0.825, 0.225),
                     breaks = c(-0.8, -0.6, -0.4, -0.2, 0, 0.2)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1),
                     limits = c(0, 0.425)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8))
c

# put them together
#fig5 <- cowplot::plot_grid(a, b, c, labels = c("a", "b", "c"), ncol = 1, align = "vh", axis = "tblr", label_size = 12)
fig5 <- cowplot::plot_grid(a, b, c, labels = c("a", "b", "c"), ncol = 3, align = "vh", axis = "tblr", label_size = 12)
fig5

# save it
ggsave(fig5, filename = "figures/figure5.png", width = 6.25, height = 2.25)
ggsave(fig5, filename = "figures/figure5.pdf", width = 6.25, height = 2.25)
