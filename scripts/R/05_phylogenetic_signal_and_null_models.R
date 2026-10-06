source('scripts/R/00_config.R')

suppressPackageStartupMessages({
  library(ape)
  library(picante)
  library(phytools)
  library(readr)
  library(dplyr)
  library(purrr)
  library(tibble)
})

comm1_raw <- read_species_by_plot_matrix(paths$first_cycle_matrix)
comm2_raw <- read_species_by_plot_matrix(paths$second_cycle_matrix)
comm <- align_community_matrices(comm1_raw, comm2_raw)
comm1_abund <- comm$cycle1
comm2_abund <- comm$cycle2

tree <- ape::read.tree(paths$phylo_tree)
tree$tip.label <- normalize_species_names(tree$tip.label)
shared_species <- intersect(colnames(comm1_abund), tree$tip.label)
tree <- ape::keep.tip(tree, shared_species)
comm1 <- ifelse(comm1_abund[, shared_species, drop = FALSE] > 0, 1, 0)
comm2 <- ifelse(comm2_abund[, shared_species, drop = FALSE] > 0, 1, 0)
comm1 <- comm1[, tree$tip.label, drop = FALSE]
comm2 <- comm2[, tree$tip.label, drop = FALSE]

mortality_rate <- colMeans(comm1 == 1 & comm2 == 0)
recruitment_rate <- colMeans(comm1 == 0 & comm2 == 1)
mortality_rate <- mortality_rate[tree$tip.label]
recruitment_rate <- recruitment_rate[tree$tip.label]

rates <- tibble(species = tree$tip.label, mortality_rate = as.numeric(mortality_rate), recruitment_rate = as.numeric(recruitment_rate))
safe_write_csv(rates, file.path(results_dir, 'tables', 'summary', 'species_level_demographic_rates.csv'))

signal_test <- function(x, trait, method) {
  fit <- phytools::phylosig(tree, x, method = method, test = TRUE)
  if (method == 'K') tibble(trait = trait, statistic = 'Blomberg_K', estimate = fit$K, p_value = fit$P)
  else tibble(trait = trait, statistic = 'Pagel_lambda', estimate = fit$lambda, p_value = fit$P)
}

signal_summary <- bind_rows(
  signal_test(mortality_rate, 'mortality_rate', 'K'),
  signal_test(recruitment_rate, 'recruitment_rate', 'K'),
  signal_test(mortality_rate, 'mortality_rate', 'lambda'),
  signal_test(recruitment_rate, 'recruitment_rate', 'lambda')
)
safe_write_csv(signal_summary, file.path(results_dir, 'tables', 'models', 'phylogenetic_signal_mortality_recruitment.csv'))

pd1 <- picante::pd(comm1, tree, include.root = TRUE)$PD
pd2 <- picante::pd(comm2, tree, include.root = TRUE)$PD
pd_obs <- pd2 - pd1
names(pd_obs) <- rownames(comm1)

nreps <- 999
plot_ids <- rownames(comm1)

ses_value <- function(obs, null) {
  s <- sd(null, na.rm = TRUE)
  if (is.na(s) || s == 0) return(NA_real_)
  (obs - mean(null, na.rm = TRUE)) / s
}

null_results <- map_dfr(plot_ids, function(plot_id) {
  c1 <- comm1[plot_id, ]
  c2 <- comm2[plot_id, ]
  present1 <- names(c1)[c1 == 1]
  absent1 <- names(c1)[c1 == 0]
  lost <- names(c1)[c1 == 1 & c2 == 0]
  gained <- names(c1)[c1 == 0 & c2 == 1]
  obs <- pd_obs[plot_id]

  null_A <- replicate(nreps, {
    sampled_lost <- if (length(lost) > 0) sample(present1, length(lost), replace = FALSE) else character(0)
    retained <- setdiff(present1, sampled_lost)
    sampled_gains <- if (length(gained) > 0) sample(absent1, length(gained), replace = FALSE) else character(0)
    simulated_c2 <- integer(length(c1)); names(simulated_c2) <- names(c1)
    simulated_c2[unique(c(retained, sampled_gains))] <- 1
    picante::pd(matrix(simulated_c2, nrow = 1, dimnames = list(plot_id, names(simulated_c2))), tree, include.root = TRUE)$PD -
      picante::pd(matrix(c1, nrow = 1, dimnames = list(plot_id, names(c1))), tree, include.root = TRUE)$PD
  })

  null_B <- replicate(nreps, {
    shuffled_tree <- tree
    shuffled_tree$tip.label <- sample(shuffled_tree$tip.label)
    picante::pd(matrix(c2, nrow = 1, dimnames = list(plot_id, names(c2))), shuffled_tree, include.root = TRUE)$PD -
      picante::pd(matrix(c1, nrow = 1, dimnames = list(plot_id, names(c1))), shuffled_tree, include.root = TRUE)$PD
  })

  bind_rows(
    tibble(UA = plot_id, null_model = 'fixed_loss_gain', observed_PD_delta = obs, null_mean = mean(null_A), null_sd = sd(null_A), SES = ses_value(obs, null_A)),
    tibble(UA = plot_id, null_model = 'tip_shuffling', observed_PD_delta = obs, null_mean = mean(null_B), null_sd = sd(null_B), SES = ses_value(obs, null_B))
  )
})

safe_write_csv(null_results, file.path(results_dir, 'tables', 'models', 'pd_delta_null_models_999_randomizations.csv'))
message('Phylogenetic signal and null model analyses written.')
