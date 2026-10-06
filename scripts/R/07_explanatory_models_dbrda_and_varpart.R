source('scripts/R/00_config.R')

suppressPackageStartupMessages({
  library(vegan)
  library(broom)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
})

env <- read_csv(file.path(processed_dir, 'environmental', 'environmental_predictors_scaled_vif5.csv'), show_col_types = FALSE) %>% mutate(UA = as.character(UA))
tax_change <- read_csv(file.path(processed_dir, 'plot_level', 'taxonomic_turnover_demography_by_plot.csv'), show_col_types = FALSE) %>% mutate(UA = as.character(UA))
phylo_alpha <- read_csv(file.path(processed_dir, 'plot_level', 'phylogenetic_alpha_diversity_by_plot.csv'), show_col_types = FALSE) %>% mutate(UA = as.character(UA))
phylo_beta <- read_csv(file.path(processed_dir, 'plot_level', 'temporal_phylogenetic_beta_by_plot.csv'), show_col_types = FALSE) %>% mutate(UA = as.character(UA))
tax_alpha <- read_csv(file.path(processed_dir, 'plot_level', 'taxonomic_alpha_diversity_hill_numbers.csv'), show_col_types = FALSE) %>% mutate(UA = as.character(UA))

analysis <- env %>%
  inner_join(tax_change, by = 'UA') %>%
  inner_join(phylo_alpha, by = 'UA') %>%
  inner_join(phylo_beta, by = 'UA') %>%
  inner_join(tax_alpha, by = 'UA') %>%
  drop_na()
safe_write_csv(analysis, file.path(processed_dir, 'plot_level', 'integrated_analysis_table.csv'))

response_vars <- c('bray_curtis_dissimilarity', 'PD_delta', 'temporal_phylo_beta_total')
env_vars <- setdiff(names(env), 'UA')
structural_vars <- c('species_recruitment', 'species_mortality', 'abundance_recruitment', 'abundance_mortality', 'hill_q0_cycle1', 'hill_q1_cycle1', 'hill_q2_cycle1')
phylo_vars <- c('PD_cycle1', 'MPD_cycle1', 'MNTD_cycle1')
structural_vars <- intersect(structural_vars, names(analysis))
phylo_vars <- intersect(phylo_vars, names(analysis))

fit_dbrda <- function(response, predictors) {
  d <- dist(response, method = 'euclidean')
  full <- vegan::capscale(d ~ ., data = predictors, add = TRUE)
  null <- vegan::capscale(d ~ 1, data = predictors, add = TRUE)
  step <- vegan::ordiR2step(null, scope = formula(full), direction = 'forward', R2scope = TRUE, permutations = 999, trace = FALSE)
  list(model = step, global = anova(step, permutations = 999), axes = anova(step, permutations = 999, by = 'axis'), terms = anova(step, permutations = 999, by = 'term'), adj_r2 = vegan::RsquareAdj(step)$adj.r.squared)
}

predictors <- analysis %>% select(all_of(c(env_vars, structural_vars, phylo_vars)))
model_summary <- list()

for (resp in response_vars) {
  fit <- fit_dbrda(analysis[[resp]], predictors)
  saveRDS(fit$model, file.path(model_dir, paste0('dbrda_', resp, '.rds')))
  safe_write_csv(anova_to_tibble(fit$global), file.path(results_dir, 'tables', 'models', 'dbRDA', paste0('dbrda_global_', resp, '.csv')))
  safe_write_csv(anova_to_tibble(fit$axes), file.path(results_dir, 'tables', 'models', 'dbRDA', paste0('dbrda_axes_', resp, '.csv')))
  safe_write_csv(anova_to_tibble(fit$terms), file.path(results_dir, 'tables', 'models', 'dbRDA', paste0('dbrda_terms_', resp, '.csv')))
  model_summary[[resp]] <- tibble(response = resp, adjusted_R2 = fit$adj_r2)
}
safe_write_csv(bind_rows(model_summary), file.path(results_dir, 'tables', 'models', 'dbRDA', 'dbrda_adjusted_R2_summary.csv'))

X_env <- analysis %>% select(all_of(env_vars))
X_structure <- analysis %>% select(all_of(structural_vars))
X_phylogeny <- analysis %>% select(all_of(phylo_vars))

for (resp in response_vars) {
  y <- dist(analysis[[resp]], method = 'euclidean')
  vp <- vegan::varpart(y, X_env, X_structure, X_phylogeny)
  saveRDS(vp, file.path(model_dir, paste0('variance_partitioning_environment_structure_phylogeny_', resp, '.rds')))
  png(file.path(fig_dir, 'exploratory', paste0('variance_partitioning_environment_structure_phylogeny_', resp, '.png')), width = 1000, height = 800, res = 140)
  plot(vp, Xnames = c('Environment', 'Structure', 'Phylogeny'))
  dev.off()
}

message('db-RDA models and environment-structure-phylogeny variance partitioning written.')
