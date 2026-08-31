# Regenerate only when intentionally reviewing a backend compatibility change.
source('tools/load-local.R')
invisible(load_local_package())
if (!requireNamespace('PhylogeneticEM', quietly=TRUE)) stop('Run make setup first.')
set.seed(20260831)
phy <- ape::stree(16L, type='balanced')
phy$edge.length <- rep(1, nrow(phy$edge))
phy$tip.label <- paste0('S', seq_len(16L), '_sp')
traits <- matrix(rnorm(32), 2, dimnames=list(c('trait one', 'trait_two'), phy$tip.label))
traits[1,3] <- NA_real_
traits[2,7] <- NA_real_
fit_model <- function(data) {
    PhylogeneticEM::PhyloEM(phy, data, process='scOU', K_max=1L, alpha=1,
        method.selection='LINselect', progress.bar=FALSE, parallel_alpha=FALSE)
}
shifted_traits <- traits
# A non-root clade avoids the equivalent solutions at the two root edges.
shifted_traits[,1:4] <- shifted_traits[,1:4] + c(10, -10)
models <- list(zero_shift=fit_model(traits), one_shift=fit_model(shifted_traits),
    univariate=fit_model(traits[1,,drop=FALSE]))
fixture <- list(models=models, seed=20260831L,
    backend_version=as.character(utils::packageVersion('PhylogeneticEM')))
dir.create('tests/testthat/fixtures', showWarnings=FALSE)
saveRDS(fixture, 'tests/testthat/fixtures/phyloem-fit.rds', compress='xz')
message('Wrote a real fitted PhyloEM object; review adapter tests before committing it.')
