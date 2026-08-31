# The optional backend stays behind a small, testable adapter boundary.
.phylogeneticem_trait_vector = function(values, traits, param, diagonal=FALSE) {
    if (diagonal) {
        values = as.matrix(values)
        if (nrow(values) != ncol(values)) stop(param, ' must be a square matrix.')
        if (nrow(values) > 1L && !is.null(rownames(values))) {
            if (anyDuplicated(rownames(values)) || !setequal(rownames(values), traits) ||
                    is.null(colnames(values)) || anyDuplicated(colnames(values)) ||
                    !setequal(colnames(values), traits)) {
                stop(param, ' matrix names must match Y_data traits.')
            }
            values = values[traits,traits,drop=FALSE]
        }
        values = diag(values)
    }
    if (!is.numeric(values) || any(!is.finite(values))) {
        stop(param, ' must contain finite numeric values.')
    }
    if (length(values) == 1L && length(traits) > 1L) {
        values = rep(unname(values), length(traits))
    }
    if (length(values) != length(traits)) stop(param, ' must match the number of traits.')
    if (!is.null(names(values))) {
        if (anyDuplicated(names(values)) || !setequal(names(values), traits)) {
            stop(param, ' names must match Y_data traits.')
        }
        values = values[traits]
    }
    unname(values)
}

.phylogeneticem_regime_table = function(pcm_out) {
    pp = .phylogeneticem_params_process(pcm_out)
    tree = fill_node_labels(pcm_out[['phylo']])
    .validate_phylo_input(tree, 'pcm_out$phylo', unique_tips=TRUE)
    traits = rownames(pcm_out[['Y_data']])
    if (!length(traits) || anyNA(traits) || any(traits == '') || anyDuplicated(traits) ||
            any(traits %in% c('regime', 'node_name', 'param'))) {
        stop('Y_data must have unique, non-empty trait row names that are not metadata columns.')
    }
    edges = pp[['shifts']][['edges']]
    if (is.null(edges)) edges = integer(0)
    edges = .normalize_integerish(edges, 'shift edges',
        min_value=1L, max_value=nrow(tree[['edge']]), allow_empty=TRUE)
    n_shift = length(edges)
    values = matrix(NA_real_, nrow=n_shift + 3L, ncol=length(traits))
    if (n_shift) {
        shifts = t(as.matrix(pp[['shifts']][['values']]))
        if (!is.numeric(shifts) || any(!is.finite(shifts)) ||
                nrow(shifts) != n_shift || ncol(shifts) != length(traits)) {
            stop('Shift values must be a numeric trait-by-shift matrix matching shift edges.')
        }
        if (!is.null(colnames(shifts))) {
            if (anyDuplicated(colnames(shifts)) || !setequal(colnames(shifts), traits)) {
                stop('Shift value names must match Y_data traits.')
            }
            shifts = shifts[,traits,drop=FALSE]
        }
        values[seq_len(n_shift),] = shifts
    }
    values[n_shift + 1L,] = .phylogeneticem_trait_vector(
        pp[['selection.strength']], traits, 'selection.strength', diagonal=TRUE)
    values[n_shift + 2L,] = .phylogeneticem_trait_vector(
        pp[['variance']], traits, 'variance', diagonal=TRUE)
    values[n_shift + 3L,] = .phylogeneticem_trait_vector(
        pp[['optimal.value']], traits, 'optimal.value')
    colnames(values) = traits
    data.frame(
        regime=c(seq_len(n_shift), rep(NA_integer_, 3L)),
        node_name=c(get_node_name_by_num(tree, tree[['edge']][edges,2]), rep(NA_character_, 3L)),
        param=c(rep('shift_value', n_shift), 'alpha', 'sigma2', 'intercept'),
        values, check.names=FALSE, stringsAsFactors=FALSE)
}

.phylogeneticem_leaf_table = function(pcm_out) {
    tree = pcm_out[['phylo']]
    .validate_phylo_input(tree, 'pcm_out$phylo', unique_tips=TRUE)
    traits = rownames(pcm_out[['Y_data']])
    if (!length(traits) || anyNA(traits) || any(traits == '') || anyDuplicated(traits)) {
        stop('Y_data must have unique, non-empty trait row names.')
    }
    regimes = get_leaf_regimes(pcm_out, mode='PhylogeneticEM')[['regime']]
    rows = lapply(c('imputed', 'expectations'), function(param) {
        values = t(.phylogeneticem_imputed_traits(pcm_out,
            trait=seq_along(traits), where='tips', what=param))
        if (!is.numeric(values) || nrow(values) != length(tree[['tip.label']]) ||
                ncol(values) != length(traits)) {
            stop('PhylogeneticEM returned incompatible trait dimensions.')
        }
        if (!is.null(rownames(values))) {
            if (anyDuplicated(rownames(values)) ||
                    !setequal(rownames(values), tree[['tip.label']])) {
                stop('PhylogeneticEM trait rows must match tree tip labels.')
            }
            values = values[tree[['tip.label']],,drop=FALSE]
        }
        if (!is.null(colnames(values))) {
            if (anyDuplicated(colnames(values)) || !setequal(colnames(values), traits)) {
                stop('PhylogeneticEM trait names must match Y_data.')
            }
            values = values[,traits,drop=FALSE]
        }
        colnames(values) = traits
        data.frame(regime=regimes, node_name=tree[['tip.label']], param=param,
            values, check.names=FALSE, stringsAsFactors=FALSE)
    })
    out = do.call(rbind, rows)
    rownames(out) = NULL
    out
}
# Comparative-model and trait-table utilities.

.phylogeneticem_params_process = function(...) {
    do.call(.get_optional_pkg_fun('PhylogeneticEM', 'params_process'), args=list(...))
}


.phylogeneticem_imputed_traits = function(...) {
    do.call(.get_optional_pkg_fun('PhylogeneticEM', 'imputed_traits'), args=list(...))
}
