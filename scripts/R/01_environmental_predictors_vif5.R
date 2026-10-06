source('scripts/R/00_config.R')

suppressPackageStartupMessages({
  library(usdm)
  library(readr)
  library(dplyr)
  library(tidyr)
})

if (!file.exists(paths$environmental_all_variables)) {
  source('scripts/R/00b_prepare_environmental_tables_from_gee_exports.R')
  prepare_environmental_tables()
}

env_raw <- read_csv(paths$environmental_all_variables, show_col_types = FALSE) %>%
  select(-any_of(c('Unnamed: 0', 'X', '...1', 'system:index', '.geo')))

if (!'UA' %in% names(env_raw)) stop('The environmental predictor table must contain a UA column.')

env_raw <- env_raw %>% mutate(UA = as.character(.data$UA))
predictor_names <- setdiff(names(env_raw), 'UA')
predictor_names <- predictor_names[vapply(env_raw[predictor_names], is.numeric, logical(1))]
env_numeric <- env_raw %>% select(UA, all_of(predictor_names))

complete_predictors <- env_numeric %>%
  select(-UA) %>%
  select(where(~ sum(!is.na(.x)) > 2 && stats::sd(.x, na.rm = TRUE) > 0))

cor_model <- usdm::vifcor(complete_predictors, th = 0.7)
cor_selected <- cor_model@results$Variables
vif_model <- usdm::vifstep(complete_predictors %>% select(all_of(cor_selected)), th = 5)
vif_selected <- vif_model@results$Variables

env_selected <- env_numeric %>% select(UA, all_of(vif_selected))
env_scaled <- zscore_numeric_columns(env_selected, exclude = 'UA')

selection_table <- tibble(variable = vif_selected) %>%
  mutate(variable_group = case_when(
    grepl('^anom_', variable) ~ 'inter-cycle anomaly',
    grepl('elevation|altitude', variable, ignore.case = TRUE) ~ 'topography',
    grepl('^monthly_', variable) ~ 'cycle-level monthly summary',
    grepl('^cycle_', variable) ~ 'cycle-level complete summary',
    TRUE ~ 'other environmental predictor'
  ))

safe_write_csv(env_raw, file.path(processed_dir, 'environmental', 'all_variables_compiled.csv'))
safe_write_csv(env_selected, file.path(processed_dir, 'environmental', 'environmental_predictors_unscaled_vif5.csv'))
safe_write_csv(env_scaled, file.path(processed_dir, 'environmental', 'environmental_predictors_scaled_vif5.csv'))
safe_write_csv(selection_table, file.path(results_dir, 'tables', 'summary', 'selected_environmental_variables_vif5.csv'))
safe_write_csv(as_tibble(cor_model@results), file.path(results_dir, 'tables', 'models', 'environmental_filtering', 'environmental_correlation_filter_vifcor.csv'))
safe_write_csv(as_tibble(vif_model@results), file.path(results_dir, 'tables', 'models', 'environmental_filtering', 'environmental_vif_filter_vif5.csv'))

message('Environmental predictors retained after correlation and VIF filtering: ', length(vif_selected))
