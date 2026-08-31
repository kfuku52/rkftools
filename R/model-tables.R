
#' Summarize a comparative-model fit
#'
#' @param pcm_out An l1ou or PhylogeneticEM result object.
#' @param mode Either `"l1ou"` or `"PhylogeneticEM"`.
#' @param species_parser Species-label convention. Species-only labels are
#'   accepted; labels without a genus/species pair are counted verbatim.
#' @param sep Literal separator in leaf labels.
#' @return A one-row model summary data frame.
#' @export
get_tree_table = function(pcm_out, mode, species_parser='legacy', sep='_') {
    mode = .normalize_choice_arg(
        value=mode,
        arg_name='mode',
        choices=c('l1ou', 'PhylogeneticEM')
    )
    species_parser = .normalize_species_parser_arg(
        value=species_parser,
        arg_name='species_parser'
    )
    sep = .normalize_single_string_arg(
        value=sep,
        arg_name='sep',
        allow_empty=FALSE
    )
    model_tree = if (mode == 'l1ou') pcm_out[['tree']] else pcm_out[['phylo']]
    .validate_phylo_input(model_tree, 'model tree', unique_tips=TRUE)
    spp = unique(.parse_species_labels(
        model_tree[['tip.label']], species_parser=species_parser,
        sep=sep, output_sep=sep, require_gene=FALSE, fallback_label=TRUE
    )[['species_labels']])
    tree_table = data.frame()
    if (mode=='l1ou') {
        leaves = pcm_out$tree$tip.label
        if (is.null(names(pcm_out$shift.configuration))) {
            shift_conf = c("0", as.character(seq_along(pcm_out$shift.configuration)))
        } else {
            shift_conf = c("0", names(pcm_out$shift.configuration))
        }
        num_regime = length(unique(shift_conf))
        num_conv_regime = length(unique(shift_conf[duplicated(shift_conf)]))
        tree_table = data.frame(
            num_shift = pcm_out$nShifts,
            num_regime = num_regime,
            num_conv_regime = num_conv_regime,
            num_uniq_regime = num_regime - num_conv_regime,
            num_species = length(spp),
            num_leaf = length(leaves),
            model_score = pcm_out$score,
            stringsAsFactors = FALSE
        )
    } else if (mode=='PhylogeneticEM') {
        pp = .phylogeneticem_params_process(pcm_out)
        df = data.frame(matrix(NA,1,4))
        colnames(df) = c('num_shift','log_likelihood','num_species','num_leaf')
        df[['num_shift']] = ncol(pp[['shifts']][['values']])
        df[['log_likelihood']] = attr(pp, 'log_likelihood')
        df[['num_species']] = length(spp)
        df[['num_leaf']] = length(pcm_out[['phylo']][['tip.label']])
        tree_table = df
    } else {
        cat('mode "', mode, '" is not supported.\n')
    }
    return(tree_table)
}


.order_shift_indices_ancestor_first = function(tree, shift_idx) {
    if (length(shift_idx) == 0L) {
        return(integer(0))
    }
    node_depths = .get_node_depths_from_root(tree)
    child_nodes = tree[['edge']][shift_idx,2]
    shift_depths = node_depths[child_nodes]
    if (any(is.na(shift_depths))) {
        stop('Unable to order shift.configuration by tree depth.')
    }
    order(shift_depths, seq_along(shift_idx), na.last=TRUE)
}


#' Build a regime-level comparative-model table
#'
#' @param pcm_out An l1ou or PhylogeneticEM result object.
#' @param mode Either `"l1ou"` or `"PhylogeneticEM"`.
#' @return A long-format regime and parameter data frame.
#' @export
get_regime_table = function(pcm_out, mode) {
    mode = .normalize_choice_arg(mode, 'mode', choices=c('l1ou', 'PhylogeneticEM'))
    switch(mode, l1ou=.l1ou_regime_table(pcm_out),
        PhylogeneticEM=.phylogeneticem_regime_table(pcm_out))
}


.apply_shift_conf_to_leaf_regimes = function(leaf_regimes, tree, shift_conf) {
    if (length(shift_conf) == 0) {
        return(leaf_regimes)
    }
    shift_idx = tryCatch(
        .normalize_integerish(
            shift_conf,
            'shift.configuration',
            allow_empty=TRUE
        ),
        error=function(e) integer(0)
    )
    if (length(shift_idx) != length(shift_conf)) {
        stop('shift.configuration must contain integer edge indices.')
    }
    invalid_idx = shift_idx[(shift_idx < 1L) | (shift_idx > nrow(tree[['edge']]))]
    if (length(invalid_idx)) {
        stop(
            'shift.configuration contains invalid edge index/indices: ',
            paste(unique(invalid_idx), collapse=', ')
        )
    }
    shift_names = names(shift_conf)
    if (is.null(shift_names)) {
        shift_names = rep(NA_character_, length(shift_idx))
    }
    shift_order = .order_shift_indices_ancestor_first(tree, shift_idx)
    shift_idx = shift_idx[shift_order]
    shift_names = shift_names[shift_order]
    out_leaf_regimes = leaf_regimes
    table_leaves = as.character(out_leaf_regimes[['label']])
    for (i in seq_along(shift_idx)) {
        node_index = shift_idx[[i]]
        node_num = tree[['edge']][node_index,2]
        regime_name = shift_names[[i]]
        regime = ifelse(is.na(regime_name) || regime_name == '', 1, regime_name)
        subtree_leaves = as.character(get_tip_labels(tree, node_num))
        out_leaf_regimes[table_leaves %in% subtree_leaves, "regime"] = regime
    }
    return(out_leaf_regimes)
}


#' Assign comparative-model regimes to leaves
#'
#' @param pcm_out An l1ou or PhylogeneticEM result object.
#' @param mode Either `"l1ou"` or `"PhylogeneticEM"`.
#' @return A data frame mapping leaf labels to regimes.
#' @export
get_leaf_regimes = function(pcm_out, mode) {
    mode = .normalize_choice_arg(
        value=mode,
        arg_name='mode',
        choices=c('l1ou', 'PhylogeneticEM')
    )
    if (mode=='l1ou') {
        tree = pcm_out[['tree']]
        y_labels = rownames(pcm_out[['Y']])
        if (is.null(y_labels)) {
            stop('pcm_out$Y must have row names matching tree tip labels in get_leaf_regimes(mode=\"l1ou\").')
        }
        missing_labels = setdiff(tree[['tip.label']], y_labels)
        extra_labels = setdiff(y_labels, tree[['tip.label']])
        if (length(missing_labels) || length(extra_labels)) {
            stop(
                'Row names of pcm_out$Y must match tree tip labels. Missing: ',
                ifelse(length(missing_labels), paste(missing_labels, collapse=', '), 'none'),
                '. Extra: ',
                ifelse(length(extra_labels), paste(extra_labels, collapse=', '), 'none'),
                '.'
            )
        }
        leaf_regimes = data.frame(regime=0, label=tree[['tip.label']])
        shift_conf = pcm_out[['shift.configuration']]
        if ((length(shift_conf)>0)&(is.null(names(shift_conf)))) {
            names(shift_conf) = seq_along(shift_conf)
        }
        leaf_regimes = .apply_shift_conf_to_leaf_regimes(leaf_regimes=leaf_regimes, tree=tree, shift_conf=shift_conf)
    } else if (mode=='PhylogeneticEM') {
        tree = pcm_out[['phylo']]
        pp = .phylogeneticem_params_process(pcm_out)
        leaf_regimes = data.frame(regime=0, label=tree[['tip.label']])
        shift_conf = pp[['shifts']][['edges']]
        if (!is.null(shift_conf)) {
            if (is.null(names(shift_conf))) {
                names(shift_conf) = seq_along(shift_conf)
            }
        }
        shift_conf = shift_conf[order(tree[['edge']][shift_conf,1], decreasing=FALSE)]
        leaf_regimes = .apply_shift_conf_to_leaf_regimes(leaf_regimes=leaf_regimes, tree=tree, shift_conf=shift_conf)
    } else {
        cat('mode "', mode, '" is not supported.\n')
        leaf_regimes = data.frame()
    }
    return(leaf_regimes)
}


#' Build a leaf-level comparative-model table
#'
#' @param pcm_out An l1ou or PhylogeneticEM result object.
#' @param mode Either `"l1ou"` or `"PhylogeneticEM"`.
#' @return A long-format leaf, regime, and parameter data frame.
#' @export
get_leaf_table = function(pcm_out, mode) {
    mode = .normalize_choice_arg(mode, 'mode', choices=c('l1ou', 'PhylogeneticEM'))
    switch(mode, l1ou=.l1ou_leaf_table(pcm_out),
        PhylogeneticEM=.phylogeneticem_leaf_table(pcm_out))
}


#' Build a bootstrap-support table
#'
#' @param pcm_out An l1ou result object.
#' @param bootstrap_result An object containing `detection.rate`.
#' @param mode Processing mode; currently only `"l1ou"`.
#' @return A node-name and bootstrap-support data frame.
#' @export
get_bootstrap_table = function(pcm_out, bootstrap_result, mode='l1ou') {
    mode = .normalize_choice_arg(
        value=mode,
        arg_name='mode',
        choices='l1ou'
    )
    if (mode=='l1ou') {
        tree = fill_node_labels(pcm_out$tree)
        node_names = c(tree$tip.label, tree$node.label)
        nleaf = length(tree$tip.label)
        nnode = length(tree$node.label)
        column_names = c("node_name", "bootstrap_support")
        bp_table = data.frame(matrix(0,length(node_names),length(column_names)))
        colnames(bp_table) = column_names
        bp_table["node_name"] = node_names
        detection_rate_raw = bootstrap_result$detection.rate
        if (is.factor(detection_rate_raw)) {
            detection_rate_raw = as.character(detection_rate_raw)
        }
        detection_rate = suppressWarnings(as.numeric(detection_rate_raw))
        invalid_numeric_idx = which(is.na(detection_rate) & !is.na(detection_rate_raw))
        if (length(invalid_numeric_idx)) {
            stop(
                'bootstrap_result$detection.rate must be numeric or coercible to numeric. ',
                'Invalid value(s): ',
                paste(unique(as.character(detection_rate_raw[invalid_numeric_idx])), collapse=', ')
            )
        }
        expected_rate_length = nleaf + nnode - 1L
        if (length(detection_rate) != expected_rate_length) {
            stop(
                'Length mismatch in bootstrap_result$detection.rate: expected ',
                expected_rate_length, ', got ', length(detection_rate), '.'
            )
        }
        tail_rate = numeric(0)
        if (length(detection_rate) > nleaf) {
            tail_rate = detection_rate[(nleaf + 1L):length(detection_rate)]
        }
        bp_table["bootstrap_support"] = c(
            detection_rate[seq_len(nleaf)],
            NA,
            tail_rate
        )
    } else {
        cat('mode "', mode, '" is not supported.\n')
        bp_table = data.frame(matrix(NA, 0, 2))
        colnames(bp_table) = c("node_name", "bootstrap_support")
    }
    return(bp_table)
}
