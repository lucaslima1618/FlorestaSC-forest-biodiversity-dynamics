source('scripts/R/00_config.R')

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(broom)
  library(randomForest)
})

analysis <- read_csv(file.path(processed_dir, 'plot_level', 'integrated_analysis_table.csv'), show_col_types = FALSE)
env <- read_csv(file.path(processed_dir, 'environmental', 'environmental_predictors_scaled_vif5.csv'), show_col_types = FALSE)
env_vars <- setdiff(names(env), 'UA')
base_predictors <- intersect(c('PD_cycle1'), names(analysis))
predictors <- c(base_predictors, env_vars)
model_df <- analysis %>% select(PD_delta, all_of(predictors)) %>% tidyr::drop_na()

full_formula <- as.formula(paste('PD_delta ~', paste(predictors, collapse = ' + ')))
glm_full <- lm(full_formula, data = model_df)
glm_step <- step(glm_full, direction = 'both', trace = FALSE)

safe_write_csv(broom::tidy(glm_step), file.path(results_dir, 'tables', 'models', 'GLM_PD_delta', 'glm_pd_delta_coefficients_recalculated.csv'))
safe_write_csv(broom::glance(glm_step), file.path(results_dir, 'tables', 'models', 'GLM_PD_delta', 'glm_pd_delta_model_summary_recalculated.csv'))
saveRDS(glm_step, file.path(model_dir, 'glm_pd_delta_stepwise.rds'))

rf <- randomForest::randomForest(full_formula, data = model_df, importance = TRUE, ntree = 1000, na.action = na.omit)
importance_df <- as.data.frame(randomForest::importance(rf)) %>% tibble::rownames_to_column('variable')
safe_write_csv(importance_df, file.path(results_dir, 'tables', 'models', 'random_forest', 'random_forest_pd_delta_importance_recalculated.csv'))
saveRDS(rf, file.path(model_dir, 'random_forest_pd_delta.rds'))
message('GLM and random forest models for PD_delta written.')
