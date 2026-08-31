
.collapse_clade_root_state = function(subtree_trait_values, subtree) {
    trait_values = suppressWarnings(as.numeric(subtree_trait_values))
    non_na_values = trait_values[!is.na(trait_values)]
    if (!length(non_na_values)) {
        return(NA_real_)
    }
    if (length(unique(non_na_values)) == 1L) {
        return(non_na_values[[1]])
    }
    child_counts = table(subtree[['edge']][,1])
    if (any(child_counts != 2L)) {
        gls_state = .nonbinary_root_gls(trait_values, subtree)
        if (length(gls_state) == 1L && !is.na(gls_state) && is.finite(gls_state)) {
            return(gls_state)
        }
        warning(
            'Could not estimate a non-binary collapsed-clade root state ',
            'with Brownian GLS; using the arithmetic mean.'
        )
        return(mean(non_na_values))
    }
    ace_error = NULL
    ace_result = tryCatch(
        suppressWarnings(
            ape::ace(
                trait_values,
                subtree,
                type='continuous',
                method='REML',
                CI=FALSE,
                model='BM'
            )
        ),
        error=function(e) {
            ace_error <<- conditionMessage(e)
            NULL
        }
    )
    if (is.null(ace_result) || is.null(ace_result$ace) || !length(ace_result$ace)) {
        warning(
            'ape::ace could not estimate a collapsed-clade root state; using the arithmetic mean. ',
            if (is.null(ace_error)) '' else paste0('Cause: ', ace_error)
        )
        return(mean(non_na_values))
    }
    root_num = get_root_num(subtree)
    root_state = ace_result$ace[names(ace_result$ace) == as.character(root_num)]
    if (!length(root_state)) {
        root_state = ace_result$ace[[1]]
    } else {
        root_state = root_state[[1]]
    }
    root_state = suppressWarnings(as.numeric(root_state))
    if (length(root_state) != 1L || is.na(root_state) || !is.finite(root_state)) {
        return(mean(non_na_values))
    }
    root_state
}


.nonbinary_root_gls = function(trait_values, subtree) {
    complete_indices = which(!is.na(trait_values))
    if (!length(complete_indices)) {
        return(NA_real_)
    }
    if (length(complete_indices) == 1L) {
        return(trait_values[[complete_indices]])
    }
    if (is.null(subtree[['edge.length']]) ||
            length(subtree[['edge.length']]) != nrow(subtree[['edge']]) ||
            anyNA(subtree[['edge.length']]) ||
            any(!is.finite(subtree[['edge.length']])) ||
            any(subtree[['edge.length']] < 0)) {
        return(NA_real_)
    }
    covariance = tryCatch(
        ape::vcv.phylo(subtree, corr=FALSE),
        error=function(e) NULL
    )
    if (is.null(covariance)) {
        return(NA_real_)
    }
    complete_tips = subtree[['tip.label']][complete_indices]
    covariance = covariance[complete_tips, complete_tips, drop=FALSE]
    covariance = (covariance + t(covariance)) / 2
    covariance_scale = max(abs(covariance))
    if (!is.finite(covariance_scale)) {
        return(NA_real_)
    }
    if (covariance_scale == 0) {
        return(mean(trait_values[complete_indices]))
    }
    decomposition = tryCatch(
        eigen(covariance / covariance_scale, symmetric=TRUE),
        error=function(e) NULL
    )
    if (is.null(decomposition)) {
        return(NA_real_)
    }
    max_value = max(abs(decomposition[['values']]))
    tolerance = max(dim(covariance)) * max_value * .Machine$double.eps
    if (any(decomposition[['values']] < -tolerance)) {
        return(NA_real_)
    }
    vectors = decomposition[['vectors']]
    values = decomposition[['values']]
    values[values <= 0] = tolerance
    one = rep(1, length(complete_indices))
    inverse_one = drop(vectors %*% (drop(crossprod(vectors, one)) / values))
    denominator = sum(inverse_one)
    if (!is.finite(denominator) || denominator == 0) {
        return(NA_real_)
    }
    estimate = sum(inverse_one * trait_values[complete_indices]) / denominator
    suppressWarnings(as.numeric(estimate))
}


#' Collapse selected clades and estimate their root traits
#'
#' @param tree A rooted `phylo` tree.
#' @param trait_table A numeric table named by tree tips.
#' @param collapse_node_nums Integer internal-node numbers to collapse.
#' @param verbose Whether to emit progress messages.
#' @return A list containing the collapsed tree, aligned traits, and collapse map.
#' @export
collapse_clades = function(tree, trait_table, collapse_node_nums, verbose=FALSE) {
    verbose = .normalize_single_logical_arg(verbose, 'verbose')
    .validate_phylo_input(tree, context='tree', rooted=TRUE, unique_tips=TRUE)
    if (verbose) {
        message('Number of leaves before collapse: ', length(tree$tip.label))
    }
    .validate_numeric_trait_table(trait_table)
    if (is.null(rownames(trait_table))) {
        stop('trait_table must have row names that match tree tip labels.')
    }
    if (anyDuplicated(rownames(trait_table))) {
        duplicated_rows = unique(rownames(trait_table)[duplicated(rownames(trait_table))])
        stop(
            'trait_table has duplicated row name(s): ',
            paste(duplicated_rows, collapse=', ')
        )
    }
    missing_tip_rows = setdiff(tree$tip.label, rownames(trait_table))
    if (length(missing_tip_rows)) {
        stop(
            'trait_table is missing rows for tree tip label(s): ',
            paste(missing_tip_rows, collapse=', ')
        )
    }
    if (length(collapse_node_nums)) {
        collapse_node_nums = tryCatch(
            .normalize_integerish(
                collapse_node_nums,
                'collapse_node_nums',
                min_value=1L,
                max_value=max(tree[['edge']])
            ),
            error=function(e) {
                stop(
                    'collapse_node_nums must be integer node numbers. ',
                    conditionMessage(e),
                    call.=FALSE
                )
            }
        )
    }
    trait_table_aligned = trait_table[tree$tip.label,,drop=FALSE]
    if (anyDuplicated(collapse_node_nums)) {
        warning('Duplicate node id(s) found in collapse_node_nums; duplicates will be ignored.')
        collapse_node_nums = unique(collapse_node_nums)
    }
    if (length(collapse_node_nums)) {
        invalid_node_nums = collapse_node_nums[
            collapse_node_nums <= length(tree$tip.label) |
            collapse_node_nums > max(tree$edge)
        ]
        if (length(invalid_node_nums)) {
            stop(
                'collapse_node_nums must contain only internal node numbers in [',
                length(tree$tip.label) + 1L, ', ', max(tree$edge), ']. Invalid node(s): ',
                paste(invalid_node_nums, collapse=', ')
            )
        }
    }
    if (length(collapse_node_nums) > 1) {
        is_redundant_descendant = vapply(collapse_node_nums, function(nn) {
            ancestors = get_ancestor_num(tree, nn)
            any(ancestors %in% collapse_node_nums)
        }, logical(1))
        if (any(is_redundant_descendant)) {
            redundant_nodes = collapse_node_nums[is_redundant_descendant]
            warning(
                'Descendant node id(s) found in collapse_node_nums and ignored because an ancestor is also selected: ',
                paste(redundant_nodes, collapse=', ')
            )
            collapse_node_nums = collapse_node_nums[!is_redundant_descendant]
        }
    }
    collapse_leaf_names = list()
    retain_tips = character(0)
    all_drop_tips = character(0)
    collapsed_trait_rows = list()
    for (cnn in collapse_node_nums) {
        collapse_key = as.character(cnn)
        collapse_leaf_names[[collapse_key]] = get_tip_labels(tree, cnn)
        retain_tip = collapse_leaf_names[[collapse_key]][[1]]
        drop_tips = collapse_leaf_names[[collapse_key]][-1]
        if (collapse_key %in% setdiff(tree$tip.label, collapse_leaf_names[[collapse_key]])) {
            stop(
                'Collapsed node label "', collapse_key,
                '" conflicts with an existing tip label outside the clade.'
            )
        }
        retain_tips[[collapse_key]] = retain_tip
        all_drop_tips = c(all_drop_tips, drop_tips)
        subtree = ape::extract.clade(tree, cnn)
        subtree_traits = trait_table_aligned[subtree$tip.label,,drop=FALSE]
        subtree_root_states = vapply(seq_len(ncol(trait_table_aligned)), function(i) {
            .collapse_clade_root_state(subtree_traits[,i], subtree)
        }, numeric(1))
        subtree_root_row = as.data.frame(t(subtree_root_states), stringsAsFactors=FALSE)
        colnames(subtree_root_row) = colnames(trait_table_aligned)
        rownames(subtree_root_row) = collapse_key
        collapsed_trait_rows[[collapse_key]] = subtree_root_row
    }
    out_tree = tree
    out_trait = trait_table_aligned
    if (length(collapse_node_nums)) {
        out_tree = ape::drop.tip(out_tree, tip=all_drop_tips, trim.internal=TRUE)
        for (collapse_key in names(retain_tips)) {
            out_tree$tip.label[out_tree$tip.label == retain_tips[[collapse_key]]] = collapse_key
        }
        collapsed_leaves = unique(unlist(collapse_leaf_names, use.names=FALSE))
        out_trait = out_trait[!(rownames(out_trait) %in% collapsed_leaves),,drop=FALSE]
        out_trait = rbind(out_trait, do.call(rbind, collapsed_trait_rows))
    }
    out_trait = out_trait[out_tree$tip.label,,drop=FALSE]
    if (verbose) {
        message('Number of leaves after collapse: ', length(out_tree$tip.label))
    }
    out = list(tree=out_tree, trait=out_trait, collapse_leaf_names=collapse_leaf_names)
    return(out)
}


#' Map original nodes to collapsed-tree nodes
#'
#' @param tree_original Original `phylo` tree.
#' @param tree_collapsed Collapsed `phylo` tree.
#' @param collapse_leaf_names Named list mapping placeholders to original tips.
#' @param verbose Whether to emit a summary message.
#' @return A two-column node-number mapping data frame.
#' @export
map_node_num = function(tree_original, tree_collapsed, collapse_leaf_names=list(), verbose=FALSE) {
    original_index = .validate_phylo_input(tree_original, 'tree_original', unique_tips=TRUE)
    collapsed_index = .validate_phylo_input(tree_collapsed, 'tree_collapsed', unique_tips=TRUE)
    # tree_collapsed should be a collapsed tree
    if (length(tree_original$tip.label) < length(tree_collapsed$tip.label)) {
        stop('tree_original must have at least as many tips as tree_collapsed.')
    }
    if (is.null(collapse_leaf_names)) {
        collapse_leaf_names = list()
    }
    if (!is.list(collapse_leaf_names)) {
        stop('collapse_leaf_names must be a named list in map_node_num().')
    }
    if (length(collapse_leaf_names)) {
        cln = names(collapse_leaf_names)
        if (is.null(cln) || any(is.na(cln) | trimws(cln) == '') || anyDuplicated(cln)) {
            stop('collapse_leaf_names must be a named list in map_node_num().')
        }
    }
    verbose = .normalize_single_logical_arg(
        value=verbose,
        arg_name='verbose'
    )
    expand_collapsed_tips = function(tip_labels, collapse_map) {
        out = character(0)
        for (tip_label in as.character(tip_labels)) {
            if (!is.null(collapse_map[[tip_label]])) {
                out = c(out, as.character(collapse_map[[tip_label]]))
            } else {
                out = c(out, tip_label)
            }
        }
        out
    }
    expanded_collapsed_tips = expand_collapsed_tips(tree_collapsed$tip.label, collapse_leaf_names)
    if (anyDuplicated(expanded_collapsed_tips)) {
        duplicated_tips = unique(expanded_collapsed_tips[duplicated(expanded_collapsed_tips)])
        stop(
            'Expanded collapsed tip labels are duplicated in map_node_num(): ',
            paste(duplicated_tips, collapse=', ')
        )
    }
    missing_tips = setdiff(tree_original$tip.label, expanded_collapsed_tips)
    extra_tips = setdiff(expanded_collapsed_tips, tree_original$tip.label)
    if (length(missing_tips) || length(extra_tips)) {
        stop(
            'Tip labels are inconsistent between tree_original and tree_collapsed in map_node_num(). Missing: ',
            ifelse(length(missing_tips), paste(missing_tips, collapse=', '), 'none'),
            '. Extra: ',
            ifelse(length(extra_tips), paste(extra_tips, collapse=', '), 'none'),
            '.'
        )
    }
    original_node_nums = seq_len(original_index[['max_node']])
    original_tip_sets = .get_node_tip_sets(tree_original, original_index)
    original_signatures = vapply(original_tip_sets, .encode_tip_signature, character(1))
    collapsed_tip_sets = .get_node_tip_sets(tree_collapsed, collapsed_index)
    expanded_sets = lapply(collapsed_tip_sets, expand_collapsed_tips, collapse_map=collapse_leaf_names)
    collapsed_signatures = vapply(expanded_sets, .encode_tip_signature, character(1))
    mapped_collapsed_num = rep(NA_integer_, length(original_node_nums))
    original_depth = integer(original_index[['max_node']])
    for (node in original_index[['preorder']][-1L]) {
        original_depth[[node]] = original_depth[[original_index[['parent']][[node]]]] + 1L
    }
    # Tips and internal nodes are different roles even on unary chains.
    for (tip in seq_len(collapsed_index[['num_tip']])) {
        label = tree_collapsed[['tip.label']][[tip]]
        if (label %in% names(collapse_leaf_names)) {
            candidates = which(original_signatures == collapsed_signatures[[tip]])
            if (!length(candidates)) stop('Collapsed tips must represent monophyletic clades.')
            # collapse_clades() records the original node number as the key.
            keyed_node = suppressWarnings(as.integer(label))
            if (!is.na(keyed_node) && keyed_node %in% candidates) {
                clade_root = keyed_node
            } else {
                clade_root = candidates[which.min(original_depth[candidates])]
            }
            inside = logical(original_index[['max_node']])
            inside[[clade_root]] = TRUE
            for (node in original_index[['preorder']]) {
                if (inside[[node]]) inside[original_index[['children']][[node]]] = TRUE
            }
            if (any(!is.na(mapped_collapsed_num[inside]))) {
                stop('Collapsed clades must not overlap in map_node_num().')
            }
            mapped_collapsed_num[inside] = tip
        } else {
            mapped_collapsed_num[match(label, tree_original[['tip.label']])] = tip
        }
    }
    internal_nodes = collapsed_index[['preorder']][
        collapsed_index[['preorder']] > collapsed_index[['num_tip']]]
    for (signature in unique(collapsed_signatures[internal_nodes])) {
        targets = internal_nodes[collapsed_signatures[internal_nodes] == signature]
        candidates = original_node_nums[is.na(mapped_collapsed_num) &
            original_node_nums > original_index[['num_tip']] & original_signatures == signature]
        candidates = candidates[order(original_depth[candidates])]
        if (!length(candidates) && collapsed_index[['num_tip']] == 1L &&
                identical(as.integer(targets), collapsed_index[['root']])) {
            next # The synthetic root above a completely collapsed tree.
        }
        if (length(targets) == 1L && length(candidates)) {
            mapped_collapsed_num[candidates] = targets
        } else if (length(candidates) == length(targets)) {
            mapped_collapsed_num[candidates] = targets
        } else {
            stop('Ambiguous node mapping: incompatible unary chains in map_node_num().')
        }
    }

    if (anyNA(mapped_collapsed_num)) {
        stop(
            'Cannot map original node number(s) in map_node_num(): ',
            paste(original_node_nums[is.na(mapped_collapsed_num)], collapse=', ')
        )
    }
    if (verbose) {
        message('Mapped ', length(original_node_nums), ' original nodes to ',
            length(unique(mapped_collapsed_num)), ' collapsed nodes.')
    }
    data.frame(
        tree_original=original_node_nums,
        tree_collapsed=mapped_collapsed_num,
        stringsAsFactors=FALSE
    )
}
