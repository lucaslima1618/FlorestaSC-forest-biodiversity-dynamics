source('scripts/R/00_config.R')

suppressPackageStartupMessages({
  library(spdep)
  library(readr)
  library(dplyr)
})

analysis <- read_csv(file.path(processed_dir, 'plot_level', 'integrated_analysis_table.csv'), show_col_types = FALSE) %>% mutate(UA = as.character(UA))
coords <- read_csv(paths$coordinates, show_col_types = FALSE) %>% mutate(UA = as.character(UA))
coord_cols <- get_coordinate_columns(coords)
coords <- coords %>% rename(longitude = all_of(coord_cols['longitude']), latitude = all_of(coord_cols['latitude']))
spatial_df <- analysis %>% inner_join(coords, by = 'UA') %>% arrange(UA)
knn <- spdep::knearneigh(as.matrix(spatial_df[, c('longitude', 'latitude')]), k = 8)
weights <- spdep::nb2listw(spdep::knn2nb(knn), style = 'W')
responses <- c('bray_curtis_dissimilarity', 'PD_delta', 'temporal_phylo_beta_total')
rows <- list()
for (resp in responses) {
  model_file <- file.path(model_dir, paste0('dbrda_', resp, '.rds'))
  if (!file.exists(model_file)) next
  model <- readRDS(model_file)
  residual_vector <- as.numeric(residuals(model))
  test <- spdep::moran.test(residual_vector, weights)
  rows[[resp]] <- tibble(
    response = resp,
    moran_I = unname(test$estimate[['Moran I statistic']]),
    expectation = unname(test$estimate[['Expectation']]),
    variance = unname(test$estimate[['Variance']]),
    p_value = test$p.value
  )
}
safe_write_csv(bind_rows(rows), file.path(results_dir, 'tables', 'models', 'morans_I_residuals_dbrda_models.csv'))
message('Moran I tests for db-RDA residuals written.')
