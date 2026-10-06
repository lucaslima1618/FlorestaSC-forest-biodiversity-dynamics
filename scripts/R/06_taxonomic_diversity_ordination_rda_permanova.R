source('scripts/R/00_config.R')

suppressPackageStartupMessages({
  library(vegan)
  library(readr)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(broom)
})

comm1_raw <- read_species_by_plot_matrix(paths$first_cycle_matrix)
comm2_raw <- read_species_by_plot_matrix(paths$second_cycle_matrix)
comm <- align_community_matrices(comm1_raw, comm2_raw)
comm1_abund <- comm$cycle1
comm2_abund <- comm$cycle2
comm1_pa <- ifelse(comm1_abund > 0, 1, 0)
comm2_pa <- ifelse(comm2_abund > 0, 1, 0)

hill_numbers <- function(mat) {
  tibble(
    UA = rownames(mat),
    hill_q0 = rowSums(mat > 0),
    hill_q1 = exp(vegan::diversity(mat, index = 'shannon')),
    hill_q2 = vegan::diversity(mat, index = 'invsimpson')
  )
}

alpha <- hill_numbers(comm1_abund) %>%
  rename(hill_q0_cycle1 = hill_q0, hill_q1_cycle1 = hill_q1, hill_q2_cycle1 = hill_q2) %>%
  inner_join(hill_numbers(comm2_abund) %>% rename(hill_q0_cycle2 = hill_q0, hill_q1_cycle2 = hill_q1, hill_q2_cycle2 = hill_q2), by = 'UA') %>%
  mutate(delta_hill_q0 = hill_q0_cycle2 - hill_q0_cycle1, delta_hill_q1 = hill_q1_cycle2 - hill_q1_cycle1, delta_hill_q2 = hill_q2_cycle2 - hill_q2_cycle1)
safe_write_csv(alpha, file.path(processed_dir, 'plot_level', 'taxonomic_alpha_diversity_hill_numbers.csv'))

make_beta <- function(abund, pa) {
  jaccard <- vegan::vegdist(pa, method = 'jaccard', binary = TRUE)
  bray <- vegan::vegdist(abund, method = 'bray')
  sorensen <- as.dist(as.matrix(jaccard) / (2 - as.matrix(jaccard)))
  list(jaccard = jaccard, sorensen = sorensen, bray_curtis = bray)
}

beta1 <- make_beta(comm1_abund, comm1_pa)
beta2 <- make_beta(comm2_abund, comm2_pa)
saveRDS(beta1, file.path(model_dir, 'taxonomic_beta_diversity_matrices_cycle1.rds'))
saveRDS(beta2, file.path(model_dir, 'taxonomic_beta_diversity_matrices_cycle2.rds'))

run_pcoa <- function(dist_obj, label) {
  fit <- cmdscale(dist_obj, eig = TRUE, k = nrow(as.matrix(dist_obj)) - 1)
  positive <- which(fit$eig > 0)
  scores <- as.data.frame(fit$points[, positive, drop = FALSE]) %>% rownames_to_column('UA')
  names(scores)[-1] <- paste0(label, '_PCoA', seq_len(length(positive)))
  eigen <- tibble(matrix = label, axis = paste0('PCoA', seq_len(length(positive))), eigenvalue = fit$eig[positive], variance_explained = fit$eig[positive] / sum(fit$eig[positive]))
  list(scores = scores, eigenvalues = eigen)
}

pcoa <- list(
  cycle1_jaccard = run_pcoa(beta1$jaccard, 'cycle1_jaccard'),
  cycle1_sorensen = run_pcoa(beta1$sorensen, 'cycle1_sorensen'),
  cycle1_bray = run_pcoa(beta1$bray_curtis, 'cycle1_bray'),
  cycle2_jaccard = run_pcoa(beta2$jaccard, 'cycle2_jaccard'),
  cycle2_sorensen = run_pcoa(beta2$sorensen, 'cycle2_sorensen'),
  cycle2_bray = run_pcoa(beta2$bray_curtis, 'cycle2_bray')
)

safe_write_csv(bind_rows(lapply(pcoa, `[[`, 'eigenvalues')), file.path(results_dir, 'tables', 'models', 'pcoa_positive_eigenvalues_taxonomic_beta.csv'))
safe_write_csv(bind_rows(lapply(pcoa, `[[`, 'scores'), .id = 'matrix'), file.path(processed_dir, 'plot_level', 'pcoa_scores_taxonomic_beta.csv'))

tax_change <- read_csv(file.path(processed_dir, 'plot_level', 'taxonomic_turnover_demography_by_plot.csv'), show_col_types = FALSE)
turnover_dist <- dist(tax_change$bray_curtis_dissimilarity)

mantel_one <- function(x, label) {
  test <- vegan::mantel(x, turnover_dist, permutations = 999)
  tibble(matrix = label, mantel_r = unname(test$statistic), p_value = test$signif)
}

mantel_results <- bind_rows(
  mantel_one(beta1$jaccard, 'cycle1_jaccard'),
  mantel_one(beta1$sorensen, 'cycle1_sorensen'),
  mantel_one(beta1$bray_curtis, 'cycle1_bray_curtis'),
  mantel_one(beta2$jaccard, 'cycle2_jaccard'),
  mantel_one(beta2$sorensen, 'cycle2_sorensen'),
  mantel_one(beta2$bray_curtis, 'cycle2_bray_curtis')
)
safe_write_csv(mantel_results, file.path(results_dir, 'tables', 'models', 'mantel_beta_diversity_vs_temporal_turnover.csv'))

rda_predictors <- alpha %>%
  left_join(pcoa$cycle1_bray$scores %>% select(UA, cycle1_bray_PCoA1, cycle1_bray_PCoA2), by = 'UA') %>%
  left_join(pcoa$cycle2_bray$scores %>% select(UA, cycle2_bray_PCoA1, cycle2_bray_PCoA2), by = 'UA') %>%
  inner_join(tax_change %>% select(UA, bray_curtis_dissimilarity), by = 'UA') %>%
  drop_na()

rda_model <- vegan::rda(bray_curtis_dissimilarity ~ hill_q0_cycle1 + hill_q1_cycle1 + hill_q2_cycle1 + cycle1_bray_PCoA1 + cycle1_bray_PCoA2 + cycle2_bray_PCoA1 + cycle2_bray_PCoA2, data = rda_predictors)
saveRDS(rda_model, file.path(model_dir, 'rda_alpha_beta_descriptors_vs_bray_curtis.rds'))
safe_write_csv(anova_to_tibble(anova(rda_model, permutations = 999)), file.path(results_dir, 'tables', 'models', 'rda_global_alpha_beta_descriptors_vs_bray_curtis.csv'))
safe_write_csv(anova_to_tibble(anova(rda_model, permutations = 999, by = 'axis')), file.path(results_dir, 'tables', 'models', 'rda_axes_alpha_beta_descriptors_vs_bray_curtis.csv'))

phylo_cycle1 <- readRDS(file.path(model_dir, 'phylogenetic_beta_diversity_cycle1.rds'))
phylo_pcoa <- run_pcoa(phylo_cycle1$beta_sor, 'phylo_cycle1_beta_sor')
phylo_scores <- phylo_pcoa$scores %>% select(UA, phylo_cycle1_beta_sor_PCoA1, phylo_cycle1_beta_sor_PCoA2)
perm_data <- phylo_scores
rownames(perm_data) <- perm_data$UA
perm_data <- perm_data %>% select(-UA)
permanova <- vegan::adonis2(beta1$bray_curtis ~ phylo_cycle1_beta_sor_PCoA1 + phylo_cycle1_beta_sor_PCoA2, data = perm_data, permutations = 999)
safe_write_csv(anova_to_tibble(permanova), file.path(results_dir, 'tables', 'models', 'permanova_taxonomic_bray_by_phylogenetic_beta_axes.csv'))

message('Hill numbers, PCoA, RDA, Mantel tests, and PERMANOVA outputs written.')
