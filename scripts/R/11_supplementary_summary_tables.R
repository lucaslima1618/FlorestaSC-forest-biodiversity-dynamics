source('scripts/R/00_config.R')

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

plot_table <- read_csv(file.path(processed_dir, 'plot_level', 'plot_level_taxonomic_phylogenetic_metrics.csv'), show_col_types = FALSE)

alpha_summary <- tibble(
  metric = c('Species richness (q = 0)', 'Shannon diversity (q = 1)', 'Simpson diversity (q = 2)'),
  cycle1_mean = c(mean(plot_table$hill_q0_cycle1, na.rm = TRUE), mean(plot_table$hill_q1_cycle1, na.rm = TRUE), mean(plot_table$hill_q2_cycle1, na.rm = TRUE)),
  cycle1_sd = c(sd(plot_table$hill_q0_cycle1, na.rm = TRUE), sd(plot_table$hill_q1_cycle1, na.rm = TRUE), sd(plot_table$hill_q2_cycle1, na.rm = TRUE)),
  cycle2_mean = c(mean(plot_table$hill_q0_cycle2, na.rm = TRUE), mean(plot_table$hill_q1_cycle2, na.rm = TRUE), mean(plot_table$hill_q2_cycle2, na.rm = TRUE)),
  cycle2_sd = c(sd(plot_table$hill_q0_cycle2, na.rm = TRUE), sd(plot_table$hill_q1_cycle2, na.rm = TRUE), sd(plot_table$hill_q2_cycle2, na.rm = TRUE)),
  mean_change = c(mean(plot_table$delta_hill_q0, na.rm = TRUE), mean(plot_table$delta_hill_q1, na.rm = TRUE), mean(plot_table$delta_hill_q2, na.rm = TRUE)),
  change_sd = c(sd(plot_table$delta_hill_q0, na.rm = TRUE), sd(plot_table$delta_hill_q1, na.rm = TRUE), sd(plot_table$delta_hill_q2, na.rm = TRUE))
)
safe_write_csv(alpha_summary, file.path(results_dir, 'tables', 'summary', 'alpha_diversity_summary_recalculated.csv'))

summarize_metric <- function(x) {
  tibble(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), min = min(x, na.rm = TRUE), q1 = quantile(x, 0.25, na.rm = TRUE), median = median(x, na.rm = TRUE), q3 = quantile(x, 0.75, na.rm = TRUE), max = max(x, na.rm = TRUE))
}

beta_metrics <- c('bray_curtis_dissimilarity', 'sorensen_dissimilarity', 'legendre_component', 'podani_component', 'species_recruitment', 'species_mortality', 'abundance_recruitment', 'abundance_mortality', 'relative_recruitment_percent', 'relative_mortality_percent', 'net_species_gain', 'net_abundance_gain', 'net_abundance_change_percent')
beta_summary <- bind_rows(lapply(beta_metrics, function(m) summarize_metric(plot_table[[m]]) %>% mutate(metric = m)), .id = NULL) %>% select(metric, everything())
safe_write_csv(beta_summary, file.path(results_dir, 'tables', 'summary', 'temporal_beta_demographic_summary_recalculated.csv'))

phylo_metrics <- c('PD_cycle1', 'PD_cycle2', 'PD_delta', 'MPD_cycle1', 'MPD_cycle2', 'MPD_delta', 'MNTD_cycle1', 'MNTD_cycle2', 'MNTD_delta')
phylo_summary <- bind_rows(lapply(phylo_metrics, function(m) summarize_metric(plot_table[[m]]) %>% mutate(metric = m)), .id = NULL) %>% select(metric, everything())
safe_write_csv(phylo_summary, file.path(results_dir, 'tables', 'summary', 'phylogenetic_diversity_summary_recalculated.csv'))

trajectory_summary <- plot_table %>% count(trajectory, name = 'plots') %>% mutate(percent = 100 * plots / sum(plots))
safe_write_csv(trajectory_summary, file.path(results_dir, 'tables', 'summary', 'community_trajectory_class_summary_recalculated.csv'))
message('Supplementary summary tables recalculated.')
