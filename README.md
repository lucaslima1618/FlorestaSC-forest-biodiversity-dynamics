# FlorestaSC forest biodiversity dynamics


This repository contains the reproducibility package for the manuscript **Spatiotemporal patterns and environmental drivers of tropical forest biodiversity dynamics in the subtropical Atlantic Forest**.

The project evaluates temporal forest reorganization in Santa Catarina, southern Brazil, using permanent sampling units measured in two inventory cycles: 2007–2010 and 2014–2020. The workflow integrates taxonomic turnover, demographic change, environmental predictors, phylogenetic diversity, phylogenetic signal, null models, ordination, db-RDA, variance partitioning, and spatial interpolation.

## Repository status

This version includes the Santa Catarina boundary shapefile, Google Earth Engine environmental exports (without sampling-unit coordinates), the V.PhyloMaker2 species list, the final S3 phylogeny, plot-level metric tables, analysis scripts, supplementary files, and documentation.

The environmental data are stored in two forms:

1. original Google Earth Engine exports in `data/raw/environmental/`;
2. standardized analysis-ready tables generated from those exports, including `all_variables.csv`.

The original anomaly export contained rows labelled as `diferenca`, but these rows did not preserve the intended `delta_` fields in the CSV schema. Therefore, the repository includes a reproducible preparation script that recalculates anomalies directly from cycle 1 and cycle 2 values as cycle 2 minus cycle 1.

## Data availability and restricted data

The forest inventory data used in this study belong to the Santa Catarina Forest and Floristic Inventory (IFFSC) and are managed through **SinFlor**. The following data are **not distributed in this repository** and must be requested from SinFlor:

- species-by-plot abundance matrices for both inventory cycles (expected in `data/raw/forest_inventory/` and `data/processed/forest_inventory/`);
- geographic coordinates of the sampling units (expected in `data/raw/spatial/sampling_unit_coordinates.csv`). Coordinate columns (`Lat`, `Long`, `.geo`) were also removed from the environmental tables and Google Earth Engine exports.

Requests: SinFlor (<https://www.sinflor.com.br/>)

After obtaining the data, place the files at the paths listed in `scripts/R/00_config.R` to run the complete workflow. Scripts that depend on the restricted data are `02`–`06` (species-by-plot matrices) and `09`–`10` (coordinates). The environmental tables, plot-level metrics, phylogeny, and result tables included here are keyed by the sampling-unit identifier (`UA`) and can be used without the restricted data.

## Repository structure

```text
FlorestaSC-biodiversity-dynamics/
├── data/
│   ├── raw/
│   │   ├── environmental/
│   │   ├── forest_inventory/      (restricted — request from SinFlor)
│   │   ├── phylogeny/
│   │   └── spatial/
│   └── processed/
│       ├── environmental/
│       ├── forest_inventory/      (restricted — request from SinFlor)
│       └── plot_level/
├── scripts/
│   ├── gee/
│   └── R/
├── results/
│   ├── figures/        (created when the scripts are run)
│   └── tables/
└── docs/
    ├── metadata/
    ├── review_support/
    └── supplementary/
```

## Analytical workflow

Run the scripts from the repository root.

```r
source('scripts/R/run_all.R')
```

For the environmental section, the first R script after configuration is:

```r
source('scripts/R/00b_prepare_environmental_tables_from_gee_exports.R')
```

This script converts the original GEE exports into analysis-ready tables and creates `data/raw/environmental/all_variables.csv`.

## Main scripts

| Script | Purpose |
|---|---|
| `scripts/gee/00_extract_cycle_mean_conditions.js` | Extract cycle-level environmental summaries from MODIS, CHIRPS, ERA5-Land, and SRTM. |
| `scripts/gee/01_extract_cycle_anomalies.js` | Extract complete cycle summaries and support the calculation of inter-cycle anomalies. |
| `scripts/R/00b_prepare_environmental_tables_from_gee_exports.R` | Convert original GEE exports into standardized environmental tables and recompute anomalies from cycle 2 minus cycle 1. |
| `scripts/R/01_environmental_predictors_vif5.R` | Compile environmental predictors, standardize variables, and retain predictors using correlation filtering and VIF <= 5. |
| `scripts/R/02_taxonomic_turnover_and_demography.R` | Calculate recruitment, mortality, Bray-Curtis, Jaccard, and Sørensen temporal dissimilarities. |
| `scripts/R/03_build_phylogeny_vphylomaker2.R` | Generate the V.PhyloMaker2 S3 phylogeny from the species list. |
| `scripts/R/04_phylogenetic_diversity_and_change.R` | Calculate PD, MPD, MNTD, Delta PD, temporal phylogenetic beta diversity, PhyloSor, and UniFrac. |
| `scripts/R/05_phylogenetic_signal_and_null_models.R` | Estimate Blomberg's K and Pagel's lambda for mortality and recruitment rates and run two null models for Delta PD. |
| `scripts/R/06_taxonomic_diversity_ordination_rda_permanova.R` | Calculate Hill numbers, PCoA, RDA, Mantel tests, and PERMANOVA. |
| `scripts/R/07_explanatory_models_dbrda_and_varpart.R` | Fit db-RDA models and partition variation among environmental, structural, and phylogenetic predictor groups. |
| `scripts/R/07b_variance_partitioning_mean_vs_anomalous_conditions.R` | Partition environmental variation between long-term mean conditions and inter-cycle anomalies. |
| `scripts/R/08_glm_and_random_forest_pd_delta.R` | Fit complementary GLM and random forest models for Delta PD. |
| `scripts/R/09_spatial_interpolation_and_mapping.R` | Perform IDW interpolation with p = 2 and ordinary kriging with fitted semivariograms. |
| `scripts/R/10_spatial_autocorrelation_of_residuals.R` | Test Moran's I on residuals from db-RDA models. |
| `scripts/R/11_supplementary_summary_tables.R` | Recalculate manuscript and supplementary descriptive tables. |

## Environmental predictor accounting

The GEE exports support 72 remote-sensing/reanalysis predictor columns when both exported environmental products are retained. Elevation is included as an additional static topographic covariate (`elevation_m`). Thus, the analysis-ready `all_variables.csv` contains 72 satellite/reanalysis predictors plus elevation.

## Spatial data

The Santa Catarina boundary is available as:

```text
data/raw/spatial/santa_catarina_boundary.shp
```

with its `.shx`, `.dbf`, `.prj`, and auxiliary files. The projection is WGS 84 / UTM zone 22S.

## Citation

Please cite the manuscript and this repository when using these data or scripts.

## License

- **Code** (`scripts/`): MIT License — see [`LICENSE`](LICENSE).
- **Data, results, and documentation** (`data/`, `results/`, `docs/`): [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/). You may share and adapt these materials for any purpose, provided appropriate credit is given by citing the manuscript and this repository.
