source('scripts/R/00_config.R')

suppressPackageStartupMessages({
  library(V.PhyloMaker2)
  library(ape)
  library(readr)
  library(dplyr)
})

taxa <- read_csv(paths$taxon_list, show_col_types = FALSE)
required_cols <- c('species', 'genus', 'family')
missing_cols <- setdiff(required_cols, names(taxa))
if (length(missing_cols) > 0) stop('The taxon list must contain species, genus, and family columns.')

taxa <- taxa %>% mutate(
  species = gsub('_', ' ', species),
  genus = ifelse(is.na(genus) | genus == '', sub(' .*$', '', species), genus)
) %>% select(species, genus, family, any_of(c('species.relative', 'genus.relative')))

phylo_object <- V.PhyloMaker2::phylo.maker(
  sp.list = taxa,
  tree = GBOTB.extended.TPL,
  nodes = nodes.info.1.TPL,
  scenarios = 'S3'
)

phylo_tree <- phylo_object$scenario.3
ape::write.tree(phylo_tree, file = paths$phylo_tree)
saveRDS(phylo_tree, file.path(model_dir, 'phylogeny_s3_vphylomaker2.rds'))

if (!is.null(phylo_object$species.list)) {
  safe_write_csv(as_tibble(phylo_object$species.list), file.path(results_dir, 'tables', 'summary', 'vphylomaker2_species_insertion_log.csv'))
}
message('V.PhyloMaker2 S3 phylogeny written to ', paths$phylo_tree)
