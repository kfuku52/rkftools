
.rphylopars_phylopars = function(...) {
    do.call(.get_optional_pkg_fun('Rphylopars', 'phylopars'), args=list(...))
}


#' Restore observed values in imputed leaf rows
#'
#' @param leaf_table A leaf table with `param` and node-name columns.
#' @param original_trait_table Original traits named by leaf rows.
#' @return A leaf table with observed values restored in imputed rows.
#' @export
restore_imputed_leaves = function(leaf_table, original_trait_table) {
    if (!is.data.frame(leaf_table)) {
        stop('leaf_table must be a data.frame.')
    }
    if (!('param' %in% colnames(leaf_table))) {
        stop('leaf_table must contain a "param" column.')
    }
    if (is.null(rownames(original_trait_table))) {
        stop('original_trait_table must have row names.')
    }
    .validate_numeric_trait_table(original_trait_table, 'original_trait_table')
    if (anyDuplicated(rownames(original_trait_table))) {
        stop('original_trait_table must have unique row names.')
    }
    traits = colnames(original_trait_table)
    if (is.null(traits) || anyNA(traits) || any(traits == '') || anyDuplicated(traits)) {
        stop('original_trait_table must have unique, non-empty trait names.')
    }
    out_leaf_table = leaf_table
    if ('node_name' %in% colnames(out_leaf_table)) {
        leaf_col = 'node_name'
    } else if ('label' %in% colnames(out_leaf_table)) {
        leaf_col = 'label'
    } else {
        stop('leaf_table must contain either "node_name" or "label" column.')
    }
    missing_trait_cols = setdiff(traits, colnames(out_leaf_table))
    if (length(missing_trait_cols)) {
        stop(
            'leaf_table is missing trait column(s): ',
            paste(missing_trait_cols, collapse=', ')
        )
    }
    .validate_numeric_trait_table(out_leaf_table[,traits,drop=FALSE], 'leaf_table trait columns')
    leaf_names = as.character(out_leaf_table[[leaf_col]])
    conditions = !is.na(leaf_names) & leaf_names %in% rownames(original_trait_table) &
        !is.na(out_leaf_table$param) & out_leaf_table$param == 'imputed'
    target_rows = which(conditions)
    for (trait in traits) {
        observed = original_trait_table[leaf_names[target_rows], trait]
        has_value = !is.na(observed)
        out_leaf_table[target_rows[has_value], trait] = observed[has_value]
    }
    return(out_leaf_table)
}


#' Impute missing traits with Rphylopars
#'
#' @param tree A `phylo` tree.
#' @param trait_table A numeric trait table named by tree tips.
#' @param verbose Whether to report the number of imputed values.
#' @return A data frame of reconstructed tip traits.
#' @export
phylogenetic_imputation = function(tree, trait_table, verbose=FALSE) {
    verbose = .normalize_single_logical_arg(verbose, 'verbose')
    .validate_numeric_trait_table(trait_table)
    if (is.null(rownames(trait_table)) || any(is.na(rownames(trait_table)) | rownames(trait_table) == '')) {
        stop('trait_table must have non-empty row names matching tree tip labels.')
    }
    .validate_phylo_input(tree, context='tree', unique_tips=TRUE, finite_lengths=TRUE)
    if (anyDuplicated(rownames(trait_table))) stop('trait_table must have unique row names.')
    trait_table = as.data.frame(trait_table, check.names=FALSE)
    if ('species' %in% colnames(trait_table)) {
        stop('Trait name \"species\" is reserved for tip labels by Rphylopars.')
    }
    trait_table2 = data.frame(species=rownames(trait_table), trait_table, check.names=FALSE)
    trait_table2 = sort_exp(trait_table2, tree, col='species')
    num_missing = sum(is.na(trait_table2[setdiff(colnames(trait_table2), 'species')]))
    num_all = nrow(trait_table2) * (ncol(trait_table2) - 1L)
    if (verbose) {
        message('Phylogenetic imputation with phylopars: ', num_missing, '/',
            num_all, ' traits will be imputed.')
    }
    rp_out = .rphylopars_phylopars(tree=tree, trait_data=trait_table2, phylo_correlated=TRUE, pheno_correlated=TRUE)
    imputed_matrix = data.frame(rp_out[['anc_recon']], check.names=FALSE)
    imputed_matrix = imputed_matrix[tree[['tip.label']],,drop=FALSE]
    return(imputed_matrix)
}
