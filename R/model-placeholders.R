
#' Build a placeholder leaf table
#'
#' @param tree A `phylo` tree.
#' @param original_trait_table Original traits named by tree tips.
#' @return A long-format placeholder leaf table.
#' @export
get_placeholder_leaf = function(tree, original_trait_table) {
    if (missing(tree) || !inherits(tree, 'phylo')) {
        stop('tree must be an object of class "phylo" in get_placeholder_leaf().')
    }
    if (is.null(rownames(original_trait_table))) {
        stop('original_trait_table must have row names in get_placeholder_leaf().')
    }
    missing_tips = setdiff(tree[['tip.label']], rownames(original_trait_table))
    if (length(missing_tips)) {
        stop('original_trait_table is missing tree tip row(s): ', paste(missing_tips, collapse=', '))
    }
    params = c('Y', 'optima', 'mu', 'residuals')
    rows = lapply(params, function(param) {
        tmp = data.frame(
            regime=rep(0,nrow(original_trait_table)),
            node_name=rownames(original_trait_table),
            param=param
        )
        tmp = cbind(tmp, original_trait_table)
        rownames(tmp) = NULL
        tmp
    })
    out = do.call(rbind, rows)
    rownames(out) = NULL
    return(out)
}


#' Build a placeholder regime table
#'
#' @param tree A `phylo` tree.
#' @param original_trait_table Original numeric trait table.
#' @return A placeholder regime table.
#' @export
get_placeholder_regime = function(tree, original_trait_table) {
    if (missing(tree) || !inherits(tree, 'phylo')) {
        stop('tree must be an object of class "phylo" in get_placeholder_regime().')
    }
    .validate_numeric_trait_table(original_trait_table, 'original_trait_table')
    params = c('alpha', 'sigma2', 'intercept', 'log_likelihood')
    trait_cols = colnames(original_trait_table)
    if (is.null(trait_cols)) {
        trait_cols = paste0("trait", seq_len(ncol(original_trait_table)))
    }
    out = do.call(rbind, lapply(params, function(param) {
        c(NA, NA, param, rep(NA, ncol(original_trait_table)))
    }))
    colnames(out) = c('regime', 'node_name', 'param', trait_cols)
    return(out)
}


#' Build a placeholder tree-level model summary
#'
#' @param tree A `phylo` tree.
#' @param original_trait_table Original numeric trait table.
#' @return A one-row placeholder model summary.
#' @export
get_placeholder_tree = function(tree, original_trait_table) {
    if (missing(tree) || !inherits(tree, 'phylo')) {
        stop('tree must be an object of class "phylo" in get_placeholder_tree().')
    }
    if (missing(original_trait_table)) {
        stop('original_trait_table is required in get_placeholder_tree().')
    }
    .validate_numeric_trait_table(original_trait_table, 'original_trait_table')
    data.frame(
        num_shift=0L,
        num_regime=1L,
        num_conv_regime=0L,
        num_uniq_regime=1L,
        num_species=length(unique(tree[['tip.label']])),
        num_leaf=length(tree[['tip.label']]),
        model_score=NA_real_,
        stringsAsFactors=FALSE
    )
}
