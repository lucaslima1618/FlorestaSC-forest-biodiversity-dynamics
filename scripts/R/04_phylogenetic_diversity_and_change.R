source('scripts/R/00_config.R')

suppressPackageStartupMessages({
  library(ape)
  library(picante)
  library(betapart)
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
if (length(shared_species) == 0) stop('No species are shared between the community matrices and the phylogeny.')

tree <- ape::keep.tip(tree, shared_species)
comm1_abund <- comm1_abund[, shared_species, drop = FALSE]
comm2_abund <- comm2_abund[, shared_species, drop = FALSE]
comm1_pa <- ifelse(comm1_abund > 0, 1, 0)
comm2_pa <- ifelse(comm2_abund > 0, 1, 0)

pd1 <- picante::pd(comm1_pa, tree, include.root = TRUE)
pd2 <- picante::pd(comm2_pa, tree, include.root = TRUE)
mpd1 <- picante::mpd(comm1_pa, cophenetic(tree))
mpd2 <- picante::mpd(comm2_pa, cophenetic(tree))
mntd1 <- picante::mntd(comm1_pa, cophenetic(tree))
mntd2 <- picante::mntd(comm2_pa, cophenetic(tree))

alpha <- tibble(
  UA = rownames(comm1_pa),
  PD_cycle1 = pd1$PD,
  PD_cycle2 = pd2$PD,
  PD_delta = pd2$PD - pd1$PD,
  MPD_cycle1 = mpd1,
  MPD_cycle2 = mpd2,
  MPD_delta = mpd2 - mpd1,
  MNTD_cycle1 = mntd1,
  MNTD_cycle2 = mntd2,
  MNTD_delta = mntd2 - mntd1
)

plot_level_beta <- map_dfr(rownames(comm1_pa), function(plot_id) {
  pa_pair <- rbind(cycle1 = comm1_pa[plot_id, ], cycle2 = comm2_pa[plot_id, ])
  abundance_pair <- rbind(cycle1 = comm1_abund[plot_id, ], cycle2 = comm2_abund[plot_id, ])
  beta_pair <- betapart::phylo.beta.pair(pa_pair, tree)
  phylo_sor <- as.matrix(beta_pair$phylo.beta.sor)[1, 2]
  phylo_sim <- as.matrix(beta_pair$phylo.beta.sim)[1, 2]
  phylo_sne <- as.matrix(beta_pair$phylo.beta.sne)[1, 2]
  unifrac <- tryCatch(as.matrix(picante::unifrac(abundance_pair, tree))[1, 2], error = function(e) NA_real_)
  tibble(
    UA = plot_id,
    temporal_phylo_beta_total = phylo_sor,
    temporal_phylo_beta_turnover = phylo_sim,
    temporal_phylo_beta_nestedness = phylo_sne,
    temporal_phylosor_similarity = 1 - phylo_sor,
    temporal_unifrac_distance = unifrac
  )
})

compute_cycle_phylo_beta <- function(comm_pa, label) {
  beta <- betapart::phylo.beta.pair(comm_pa, tree)
  list(
    cycle = label,
    beta_sor = beta$phylo.beta.sor,
    beta_sim = beta$phylo.beta.sim,
    beta_sne = beta$phylo.beta.sne,
    summary = tibble(
      cycle = label,
      beta_sor_mean = mean(as.matrix(beta$phylo.beta.sor)[upper.tri(as.matrix(beta$phylo.beta.sor))], na.rm = TRUE),
      beta_sim_mean = mean(as.matrix(beta$phylo.beta.sim)[upper.tri(as.matrix(beta$phylo.beta.sim))], na.rm = TRUE),
      beta_sne_mean = mean(as.matrix(beta$phylo.beta.sne)[upper.tri(as.matrix(beta$phylo.beta.sne))], na.rm = TRUE)
    )
  )
}

cycle1_beta <- compute_cycle_phylo_beta(comm1_pa, 'cycle1')
cycle2_beta <- compute_cycle_phylo_beta(comm2_pa, 'cycle2')

saveRDS(cycle1_beta, file.path(model_dir, 'phylogenetic_beta_diversity_cycle1.rds'))
saveRDS(cycle2_beta, file.path(model_dir, 'phylogenetic_beta_diversity_cycle2.rds'))
safe_write_csv(alpha, file.path(processed_dir, 'plot_level', 'phylogenetic_alpha_diversity_by_plot.csv'))
safe_write_csv(plot_level_beta, file.path(processed_dir, 'plot_level', 'temporal_phylogenetic_beta_by_plot.csv'))
safe_write_csv(bind_rows(cycle1_beta$summary, cycle2_beta$summary), file.path(results_dir, 'tables', 'summary', 'phylogenetic_beta_diversity_by_cycle_summary.csv'))
message('Phylogenetic diversity and beta-diversity outputs written.')
