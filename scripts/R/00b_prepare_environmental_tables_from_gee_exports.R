source('scripts/R/00_config.R')

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

standardize_cycle_label <- function(x) {
  x <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    x %in% c('ciclo1', 'cycle1') ~ 'cycle1',
    x %in% c('ciclo2', 'cycle2') ~ 'cycle2',
    x %in% c('diferenca', 'difference') ~ 'difference',
    TRUE ~ x
  )
}

prepare_environmental_tables <- function() {
  raw_monthly <- file.path(raw_dir, 'environmental', 'gee_export_cycle_mean_conditions_original.csv')
  raw_complete <- file.path(raw_dir, 'environmental', 'gee_export_cycle_summaries_original.csv')
  if (!file.exists(raw_monthly) || !file.exists(raw_complete)) {
    stop('The original GEE exports are required: gee_export_cycle_mean_conditions_original.csv and gee_export_cycle_summaries_original.csv.')
  }
  monthly <- read_csv(raw_monthly, show_col_types = FALSE)
  complete <- read_csv(raw_complete, show_col_types = FALSE)
  if ('ciclo' %in% names(monthly) && !'cycle' %in% names(monthly)) monthly <- rename(monthly, cycle = ciclo)
  if ('ciclo' %in% names(complete) && !'cycle' %in% names(complete)) complete <- rename(complete, cycle = ciclo)
  monthly <- monthly %>% mutate(cycle = standardize_cycle_label(.data$cycle)) %>% filter(.data$cycle %in% c('cycle1', 'cycle2'))
  complete <- complete %>% mutate(cycle = standardize_cycle_label(.data$cycle)) %>% filter(.data$cycle %in% c('cycle1', 'cycle2'))

  monthly_value_cols <- monthly %>% select(where(is.numeric)) %>% names() %>% setdiff(c('UA', 'Lat', 'Long', 'altitude'))
  complete_value_cols <- complete %>% select(where(is.numeric)) %>% names() %>% setdiff(c('UA', 'Lat', 'Long'))

  elevation <- monthly %>% select(UA, elevation_m = altitude) %>% distinct(.data$UA, .keep_all = TRUE)

  long_complete <- complete %>% select(UA, any_of(c('Lat', 'Long')), cycle, all_of(complete_value_cols)) %>% left_join(elevation, by = 'UA')
  safe_write_csv(long_complete, paths$environmental_cycle_means)

  c1 <- complete %>% filter(.data$cycle == 'cycle1') %>% select(UA, all_of(complete_value_cols))
  c2 <- complete %>% filter(.data$cycle == 'cycle2') %>% select(UA, all_of(complete_value_cols))
  anomalies <- c1 %>% inner_join(c2, by = 'UA', suffix = c('_cycle1', '_cycle2'))
  for (v in complete_value_cols) {
    anomalies[[paste0('anom_', v)]] <- anomalies[[paste0(v, '_cycle2')]] - anomalies[[paste0(v, '_cycle1')]]
  }
  coords <- complete %>% select(UA, any_of(c('Lat', 'Long'))) %>% distinct(.data$UA, .keep_all = TRUE)
  anomalies <- coords %>% inner_join(anomalies %>% select(UA, starts_with('anom_')), by = 'UA')
  safe_write_csv(anomalies, paths$environmental_anomalies)

  monthly_wide <- monthly %>%
    select(UA, cycle, all_of(monthly_value_cols)) %>%
    pivot_wider(names_from = cycle, values_from = all_of(monthly_value_cols), names_glue = 'monthly_{.value}_{cycle}')

  complete_wide <- complete %>%
    select(UA, cycle, all_of(complete_value_cols)) %>%
    pivot_wider(names_from = cycle, values_from = all_of(complete_value_cols), names_glue = 'cycle_{.value}_{cycle}')

  all_variables <- monthly_wide %>%
    full_join(complete_wide, by = 'UA') %>%
    full_join(anomalies %>% select(UA, starts_with('anom_')), by = 'UA') %>%
    left_join(elevation, by = 'UA') %>%
    select(UA, sort(setdiff(names(.), 'UA')))

  safe_write_csv(all_variables, paths$environmental_all_variables)
  safe_write_csv(all_variables, file.path(processed_dir, 'environmental', 'all_variables_compiled_from_gee_exports.csv'))
  invisible(all_variables)
}

if (!file.exists(paths$environmental_all_variables)) {
  prepare_environmental_tables()
}
