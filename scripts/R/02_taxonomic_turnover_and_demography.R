source('scripts/R/00_config.R')

suppressPackageStartupMessages({
  library(vegan)
  library(readr)
  library(dplyr)
  library(tibble)
  library(purrr)
})

comm1_raw <- read_species_by_plot_matrix(paths$first_cycle_matrix)
comm2_raw <- read_species_by_plot_matrix(paths$second_cycle_matrix)
comm <- align_community_matrices(comm1_raw, comm2_raw)
comm1 <- comm$cycle1
comm2 <- comm$cycle2
plot_ids <- rownames(comm1)

metrics <- map_dfr(plot_ids, function(plot_id) {
  a1 <- comm1[plot_id, ]
  a2 <- comm2[plot_id, ]
  pa1 <- as.integer(a1 > 0)
  pa2 <- as.integer(a2 > 0)
  shared <- sum(pa1 == 1 & pa2 == 1)
  lost <- sum(pa1 == 1 & pa2 == 0)
  gained <- sum(pa1 == 0 & pa2 == 1)
  sorensen <- ifelse((2 * shared + lost + gained) == 0, NA_real_, (lost + gained) / (2 * shared + lost + gained))
  bray <- as.numeric(vegan::vegdist(rbind(a1, a2), method = 'bray'))
  jaccard <- as.numeric(vegan::vegdist(rbind(pa1, pa2), method = 'jaccard', binary = TRUE))
  abundance_recruitment <- sum(pmax(a2 - a1, 0))
  abundance_mortality <- sum(pmax(a1 - a2, 0))
  total1 <- sum(a1)
  total2 <- sum(a2)
  tibble(
    UA = plot_id,
    sorensen_dissimilarity = sorensen,
    jaccard_dissimilarity = jaccard,
    bray_curtis_dissimilarity = bray,
    species_recruitment = gained,
    species_mortality = lost,
    abundance_recruitment = abundance_recruitment,
    abundance_mortality = abundance_mortality,
    total_abundance_cycle1 = total1,
    total_abundance_cycle2 = total2,
    relative_recruitment = ifelse(total2 == 0, NA_real_, abundance_recruitment / total2),
    relative_mortality = ifelse(total1 == 0, NA_real_, abundance_mortality / total1),
    relative_recruitment_percent = 100 * ifelse(total2 == 0, NA_real_, abundance_recruitment / total2),
    relative_mortality_percent = 100 * ifelse(total1 == 0, NA_real_, abundance_mortality / total1),
    net_species_gain = gained - lost,
    net_abundance_gain = abundance_recruitment - abundance_mortality,
    net_abundance_change_percent = 100 * ((ifelse(total1 == 0, NA_real_, total2 / total1)) - 1),
    trajectory = case_when(
      gained > lost ~ 'Regeneration',
      gained < lost ~ 'Decline',
      TRUE ~ 'Stagnation'
    )
  )
})

safe_write_csv(as.data.frame(comm1) %>% rownames_to_column('UA'), file.path(processed_dir, 'forest_inventory', 'community_matrix_cycle1_aligned.csv'))
safe_write_csv(as.data.frame(comm2) %>% rownames_to_column('UA'), file.path(processed_dir, 'forest_inventory', 'community_matrix_cycle2_aligned.csv'))
safe_write_csv(metrics, file.path(processed_dir, 'plot_level', 'taxonomic_turnover_demography_by_plot.csv'))

summary_table <- metrics %>% summarise(
  plots = n(),
  mean_bray_curtis = mean(bray_curtis_dissimilarity, na.rm = TRUE),
  median_bray_curtis = median(bray_curtis_dissimilarity, na.rm = TRUE),
  mean_species_recruitment = mean(species_recruitment, na.rm = TRUE),
  mean_species_mortality = mean(species_mortality, na.rm = TRUE),
  mean_net_species_gain = mean(net_species_gain, na.rm = TRUE)
)
safe_write_csv(summary_table, file.path(results_dir, 'tables', 'summary', 'taxonomic_turnover_demography_summary.csv'))
message('Taxonomic turnover and demographic metrics written for ', nrow(metrics), ' plots.')
