# Methods-to-script crosswalk

This document maps each Methods section to the repository files that reproduce or support it.

| Methods component | Repository implementation | Key files |
|---|---|---|
| Study area and sampling design | Inventory matrices and sampling-unit coordinates (restricted, request from SinFlor) and Santa Catarina boundary | `data/raw/forest_inventory/` (restricted), `data/processed/forest_inventory/` (restricted), `data/raw/spatial/sampling_unit_coordinates.csv` (restricted), `data/raw/spatial/santa_catarina_boundary.shp` |
| Environmental variables from GEE | Original GEE scripts and original exported CSVs | `scripts/gee/`, `data/raw/environmental/gee_export_cycle_mean_conditions_original.csv`, `data/raw/environmental/gee_export_cycle_summaries_original.csv` |
| Environmental preprocessing and anomalies | Recomputes anomalies as cycle 2 minus cycle 1 and creates the full candidate predictor table | `scripts/R/00b_prepare_environmental_tables_from_gee_exports.R`, `data/raw/environmental/all_variables.csv` |
| VIF 5 environmental selection | Applies pairwise correlation filtering followed by VIF <= 5 | `scripts/R/01_environmental_predictors_vif5.R` |
| Species recruitment, mortality, and Bray-Curtis change | Calculates temporal species and abundance change per plot | `scripts/R/02_taxonomic_turnover_and_demography.R` |
| V.PhyloMaker2 phylogeny | Builds the scenario S3 phylogeny from the species list | `scripts/R/03_build_phylogeny_vphylomaker2.R`, `data/raw/phylogeny/sample_species_list.csv`, `data/raw/phylogeny/arvore_filogenetica_S3.tre` |
| PD, MPD, MNTD, ΔPD, and phylogenetic beta diversity | Calculates phylogenetic metrics and temporal changes | `scripts/R/04_phylogenetic_diversity_and_change.R` |
| Phylogenetic signal and null models | Runs Blomberg's K, Pagel's lambda, fixed-loss/gain null model, and tip-shuffling null model | `scripts/R/05_phylogenetic_signal_and_null_models.R` |
| Hill numbers, beta diversity, PCoA, RDA, Mantel, and PERMANOVA | Reproduces diversity and ordination analyses | `scripts/R/06_taxonomic_diversity_ordination_rda_permanova.R` |
| db-RDA and variance partitioning | Separates the environmental/structural/phylogenetic model from the mean/anomaly model | `scripts/R/07_explanatory_models_dbrda_and_varpart.R`, `scripts/R/07b_variance_partitioning_mean_vs_anomalous_conditions.R` |
| Spatial interpolation and Moran's I | Performs IDW with p = 2, ordinary kriging, semivariogram selection, and residual autocorrelation tests | `scripts/R/09_spatial_interpolation_and_mapping.R`, `scripts/R/10_spatial_autocorrelation_of_residuals.R` |
| Supplementary descriptive tables | Recomputes the descriptive tables used in supplementary material | `scripts/R/11_supplementary_summary_tables.R` |
