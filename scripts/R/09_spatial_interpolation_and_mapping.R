source('scripts/R/00_config.R')

suppressPackageStartupMessages({
  library(sf)
  library(gstat)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

if (!file.exists(paths$state_boundary)) stop('Provide the Santa Catarina boundary shapefile at data/raw/spatial/santa_catarina_boundary.shp.')
analysis <- read_csv(file.path(processed_dir, 'plot_level', 'integrated_analysis_table.csv'), show_col_types = FALSE) %>% mutate(UA = as.character(UA))
coords <- read_csv(paths$coordinates, show_col_types = FALSE) %>% mutate(UA = as.character(UA))
coord_cols <- get_coordinate_columns(coords)
coords <- coords %>% rename(longitude = all_of(coord_cols['longitude']), latitude = all_of(coord_cols['latitude']))
state <- sf::st_read(paths$state_boundary, quiet = TRUE) %>% sf::st_make_valid()

points <- analysis %>% inner_join(coords, by = 'UA') %>% st_as_sf(coords = c('longitude', 'latitude'), crs = 4326) %>% st_transform(st_crs(state))
grid <- st_make_grid(state, cellsize = 10000, what = 'centers') %>% st_as_sf()
grid <- sf::st_filter(grid, state, .predicate = sf::st_within)

responses <- intersect(c('bray_curtis_dissimilarity', 'species_recruitment', 'species_mortality', 'PD_delta', 'temporal_phylo_beta_turnover', 'temporal_phylo_beta_nestedness'), names(points))

fit_best_variogram <- function(formula_obj, points_sf) {
  empirical <- gstat::variogram(formula_obj, data = as(points_sf, 'Spatial'))
  candidates <- list(
    Spherical = try(gstat::fit.variogram(empirical, gstat::vgm(model = 'Sph')), silent = TRUE),
    Exponential = try(gstat::fit.variogram(empirical, gstat::vgm(model = 'Exp')), silent = TRUE),
    Gaussian = try(gstat::fit.variogram(empirical, gstat::vgm(model = 'Gau')), silent = TRUE)
  )
  candidates <- candidates[!vapply(candidates, inherits, logical(1), 'try-error')]
  if (length(candidates) == 0) stop('No semivariogram model could be fitted.')
  sse <- vapply(candidates, function(x) attr(x, 'SSErr'), numeric(1))
  best <- names(which.min(sse))
  list(empirical = empirical, fitted = candidates[[best]], model = best, sse = sse[[best]])
}

variogram_rows <- list()

for (resp in responses) {
  f <- as.formula(paste(resp, '~ 1'))
  p <- points[!is.na(points[[resp]]), ]
  vg <- fit_best_variogram(f, p)
  idw <- gstat::idw(f, locations = as(p, 'Spatial'), newdata = as(grid, 'Spatial'), idp = 2)
  krig <- gstat::krige(f, locations = as(p, 'Spatial'), newdata = as(grid, 'Spatial'), model = vg$fitted)
  idw_sf <- st_as_sf(idw)
  krig_sf <- st_as_sf(krig)
  st_write(idw_sf, file.path(results_dir, 'figures', 'exploratory', paste0('idw_', resp, '.gpkg')), delete_dsn = TRUE, quiet = TRUE)
  st_write(krig_sf, file.path(results_dir, 'figures', 'exploratory', paste0('ordinary_kriging_', resp, '.gpkg')), delete_dsn = TRUE, quiet = TRUE)

  variogram_rows[[resp]] <- tibble(
    response = resp,
    interpolation_idw_power = 2,
    kriging_type = 'ordinary kriging',
    selected_semivariogram = vg$model,
    nugget = vg$fitted$psill[1],
    partial_sill = ifelse(nrow(vg$fitted) > 1, vg$fitted$psill[2], NA_real_),
    total_sill = sum(vg$fitted$psill, na.rm = TRUE),
    range = ifelse(nrow(vg$fitted) > 1, vg$fitted$range[2], vg$fitted$range[1]),
    sse = vg$sse
  )
}

safe_write_csv(bind_rows(variogram_rows), file.path(results_dir, 'tables', 'models', 'semivariogram_parameters_by_response.csv'))
message('IDW, ordinary kriging, and semivariogram parameters written.')
