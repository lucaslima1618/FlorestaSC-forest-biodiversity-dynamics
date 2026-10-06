required_packages <- c(
  'tidyverse', 'vegan', 'picante', 'ape', 'phytools', 'betapart',
  'V.PhyloMaker2', 'usdm', 'sf', 'terra', 'gstat', 'spdep',
  'broom', 'ggrepel', 'patchwork', 'hillR', 'randomForest'
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  message('Missing packages: ', paste(missing_packages, collapse = ', '))
}

suppressPackageStartupMessages({
  library(tidyverse)
})

project_root <- normalizePath('.', mustWork = TRUE)

data_dir <- file.path(project_root, 'data')
raw_dir <- file.path(data_dir, 'raw')
processed_dir <- file.path(data_dir, 'processed')
results_dir <- file.path(project_root, 'results')
fig_dir <- file.path(results_dir, 'figures')
tab_dir <- file.path(results_dir, 'tables')
model_dir <- file.path(results_dir, 'models')

for (x in c(
  data_dir, raw_dir, processed_dir, results_dir, fig_dir, tab_dir, model_dir,
  file.path(fig_dir, 'manuscript'),
  file.path(fig_dir, 'exploratory'),
  file.path(tab_dir, 'summary'),
  file.path(tab_dir, 'models'),
  file.path(tab_dir, 'models', 'dbRDA'),
  file.path(tab_dir, 'models', 'variance_partitioning'),
  file.path(tab_dir, 'models', 'environmental_filtering')
)) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE)
}

paths <- list(
  first_cycle_matrix = file.path(processed_dir, 'forest_inventory', 'first_cycle_species_by_plot.csv'),
  second_cycle_matrix = file.path(processed_dir, 'forest_inventory', 'second_cycle_species_by_plot.csv'),
  coordinates = file.path(raw_dir, 'spatial', 'sampling_unit_coordinates.csv'),
  environmental_all_variables = file.path(raw_dir, 'environmental', 'all_variables.csv'),
  environmental_cycle_means = file.path(raw_dir, 'environmental', 'florestasc_cycle_mean_environmental_conditions.csv'),
  environmental_anomalies = file.path(raw_dir, 'environmental', 'florestasc_intercycle_environmental_anomalies.csv'),
  taxon_list = file.path(raw_dir, 'phylogeny', 'sample_species_list.csv'),
  phylo_tree = file.path(raw_dir, 'phylogeny', 'arvore_filogenetica_S3.tre'),
  state_boundary = file.path(raw_dir, 'spatial', 'santa_catarina_boundary.shp')
)

normalize_species_names <- function(x) {
  x <- gsub('’|‘|`|\'|´', '', x)
  x <- iconv(x, to = 'ASCII//TRANSLIT')
  x <- gsub('\\b(subsp\\.|ssp\\.|var\\.|forma|f\\.|fo\\.|cv\\.|cultivar)\\b.*$', '', x, ignore.case = TRUE)
  x <- gsub('\\b(cf\\.|aff\\.|sp\\.|spp\\.)\\b', '', x, ignore.case = TRUE)
  x <- gsub('_', ' ', x)
  x <- gsub('\\s+', ' ', trimws(x))
  parts <- strsplit(x, ' ')
  out <- vapply(parts, function(p) {
    if (length(p) == 0 || is.na(p[1]) || p[1] == '') return(NA_character_)
    genus <- paste0(toupper(substr(p[1], 1, 1)), tolower(substr(p[1], 2, nchar(p[1]))))
    if (length(p) == 1) return(genus)
    epithet <- tolower(p[2])
    paste(genus, epithet, sep = '_')
  }, character(1))
  out
}

read_species_by_plot_matrix <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  if (!all(c('Nome', 'Familia') %in% names(df))) {
    names(df)[1] <- 'Nome'
    if (!'Familia' %in% names(df)) df$Familia <- NA_character_
  }
  plot_cols <- setdiff(names(df), c('Nome', 'Familia'))
  df <- df %>%
    mutate(species = normalize_species_names(.data$Nome)) %>%
    select(species, all_of(plot_cols)) %>%
    filter(!is.na(species), species != '') %>%
    mutate(across(all_of(plot_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
    group_by(species) %>%
    summarise(across(all_of(plot_cols), ~ sum(.x, na.rm = TRUE)), .groups = 'drop')
  mat <- df %>% tibble::column_to_rownames('species') %>% as.matrix()
  mat <- t(mat)
  rownames(mat) <- as.character(rownames(mat))
  mat[is.na(mat)] <- 0
  mat
}

expand_matrix <- function(mat, plots, species) {
  out <- matrix(0, nrow = length(plots), ncol = length(species), dimnames = list(plots, species))
  out[rownames(mat), colnames(mat)] <- mat
  out
}

align_community_matrices <- function(comm1, comm2) {
  all_plots <- sort(unique(c(rownames(comm1), rownames(comm2))))
  all_species <- sort(unique(c(colnames(comm1), colnames(comm2))))
  list(
    cycle1 = expand_matrix(comm1, all_plots, all_species),
    cycle2 = expand_matrix(comm2, all_plots, all_species)
  )
}

zscore_numeric_columns <- function(df, exclude = NULL) {
  out <- df
  numeric_cols <- setdiff(names(df)[vapply(df, is.numeric, logical(1))], exclude)
  out[numeric_cols] <- lapply(out[numeric_cols], function(x) as.numeric(scale(x)))
  out
}

get_coordinate_columns <- function(df) {
  lon <- names(df)[grepl('^(long|lon|x)$|longitude', names(df), ignore.case = TRUE)][1]
  lat <- names(df)[grepl('^(lat|y)$|latitude', names(df), ignore.case = TRUE)][1]
  if (is.na(lon) || is.na(lat)) stop('Coordinate table must contain longitude and latitude columns.')
  c(longitude = lon, latitude = lat)
}

safe_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(x, path)
}


anova_to_tibble <- function(x) {
  as.data.frame(x) %>% tibble::rownames_to_column('term')
}

set.seed(1234)
