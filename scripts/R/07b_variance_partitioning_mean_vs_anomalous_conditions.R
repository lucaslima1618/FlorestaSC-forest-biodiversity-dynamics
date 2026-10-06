source('scripts/R/00_config.R')

suppressPackageStartupMessages({
  library(vegan)
  library(readr)
  library(dplyr)
  library(tidyr)
})

env <- read_csv(file.path(processed_dir, 'environmental', 'environmental_predictors_scaled_vif5.csv'), show_col_types = FALSE) %>% mutate(UA = as.character(UA))
analysis <- read_csv(file.path(processed_dir, 'plot_level', 'integrated_analysis_table.csv'), show_col_types = FALSE) %>% mutate(UA = as.character(UA))

env_vars <- setdiff(names(env), 'UA')
anomaly_vars <- env_vars[grepl('^anom_|^delta_', env_vars)]
mean_vars <- setdiff(env_vars, anomaly_vars)

if (length(mean_vars) == 0 || length(anomaly_vars) == 0) stop('Both mean-condition and anomaly predictor sets must contain at least one variable.')

responses <- c('bray_curtis_dissimilarity', 'PD_delta')
summary_rows <- list()

for (resp in responses) {
  dat <- analysis %>% select(UA, all_of(resp), all_of(mean_vars), all_of(anomaly_vars)) %>% drop_na()
  y <- dist(dat[[resp]], method = 'euclidean')
  X_mean <- dat %>% select(all_of(mean_vars))
  X_anom <- dat %>% select(all_of(anomaly_vars))
  vp <- vegan::varpart(y, X_mean, X_anom)
  saveRDS(vp, file.path(model_dir, paste0('variance_partitioning_mean_vs_anomaly_', resp, '.rds')))

  total <- vegan::RsquareAdj(vegan::capscale(y ~ ., data = dat %>% select(all_of(c(mean_vars, anomaly_vars))), add = TRUE))$adj.r.squared
  mean_total <- vegan::RsquareAdj(vegan::capscale(y ~ ., data = X_mean, add = TRUE))$adj.r.squared
  anomaly_total <- vegan::RsquareAdj(vegan::capscale(y ~ ., data = X_anom, add = TRUE))$adj.r.squared
  shared <- mean_total + anomaly_total - total
  pure_mean <- total - anomaly_total
  pure_anom <- total - mean_total

  summary_rows[[resp]] <- tibble(
    response = resp,
    pure_mean_conditions_fraction = pure_mean,
    pure_anomaly_conditions_fraction = pure_anom,
    shared_fraction = shared,
    total_adjusted_R2 = total,
    mean_conditions_model_fraction = mean_total,
    anomaly_conditions_model_fraction = anomaly_total,
    pure_mean_conditions_percent = 100 * pure_mean,
    pure_anomaly_conditions_percent = 100 * pure_anom,
    shared_percent = 100 * shared,
    total_percent = 100 * total
  )

  png(file.path(fig_dir, 'exploratory', paste0('variance_partitioning_mean_vs_anomaly_', resp, '.png')), width = 900, height = 750, res = 140)
  plot(vp, Xnames = c('Mean conditions', 'Anomalies'))
  dev.off()
}

safe_write_csv(bind_rows(summary_rows), file.path(results_dir, 'tables', 'models', 'variance_partitioning', 'variance_partitioning_mean_vs_anomalous_conditions_recalculated.csv'))
message('Mean versus anomalous-condition variance partitioning written.')
