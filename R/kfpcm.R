# Title     : TODO
# Objective : TODO
# Created by: kf
# Created on: 5/19/18

.phylogeneticem_params_process = function(...) {
    do.call(.get_optional_pkg_fun('PhylogeneticEM', 'params_process'), args=list(...))
}

.phylogeneticem_imputed_traits = function(...) {
    do.call(.get_optional_pkg_fun('PhylogeneticEM', 'imputed_traits'), args=list(...))
}

.rphylopars_phylopars = function(...) {
    do.call(.get_optional_pkg_fun('Rphylopars', 'phylopars'), args=list(...))
}

.validate_numeric_trait_table = function(trait_table, context='trait_table') {
    if (is.null(dim(trait_table)) || length(dim(trait_table)) != 2L) {
        stop(context, ' must be a matrix or data.frame.')
    }
    if (ncol(trait_table) == 0L) {
        stop(context, ' must contain at least one trait column.')
    }
    valid_columns = vapply(seq_len(ncol(trait_table)), function(i) {
        values = trait_table[,i]
        is.numeric(values) || is.integer(values) || is.logical(values)
    }, logical(1))
    if (any(!valid_columns)) {
        stop(context, ' must contain only numeric or logical trait columns.')
    }
    finite_columns = vapply(seq_len(ncol(trait_table)), function(i) {
        values = as.numeric(trait_table[,i])
        all(is.finite(values[!is.na(values)]))
    }, logical(1))
    if (any(!finite_columns)) {
        trait_names = colnames(trait_table)
        if (is.null(trait_names)) {
            trait_names = as.character(seq_len(ncol(trait_table)))
        }
        stop(
            context, ' contains non-finite trait values in column(s): ',
            paste(trait_names[!finite_columns], collapse=', ')
        )
    }
    invisible(TRUE)
}

remove_invariant_traits = function(trait_table, small_dif=0.001) {
    .validate_numeric_trait_table(trait_table)
    small_dif = suppressWarnings(as.numeric(small_dif))
    if (length(small_dif) != 1L || is.na(small_dif) || !is.finite(small_dif) || small_dif < 0) {
        stop('small_dif must be a single finite non-negative numeric value.')
    }
    out_trait_table = trait_table
    num_traits = ncol(trait_table)
    trait_names = colnames(trait_table)
    if (is.null(trait_names)) {
        trait_names = as.character(seq_len(num_traits))
    }

    is_small_dif = logical(num_traits)
    for (i in seq_len(num_traits)) {
        trait_values = trait_table[,i]
        trait_values = as.numeric(trait_values)
        trait_values = trait_values[!is.na(trait_values)]
        if (length(trait_values) == 0) {
            is_small_dif[i] = TRUE
            next
        }
        trait_min = min(trait_values)
        trait_max = max(trait_values)
        is_small_dif[i] = (trait_max - trait_min < small_dif)
    }

    removed_traits = trait_names[is_small_dif]
    if (length(removed_traits)) {
        cat("Trait removed due to small difference (<", small_dif, '):',  removed_traits, '\n')
    } else {
        cat('All traits passed small difference check.\n')
    }
    out_trait_table = out_trait_table[,!is_small_dif, drop=FALSE]
    out = list(trait_table=out_trait_table, removed_traits=removed_traits)
    return(out)
}

.replicate_base_name = function(column_name, replicate_sep) {
    split_name = strsplit(column_name, replicate_sep, fixed=TRUE)[[1]]
    if (length(split_name) <= 1L) {
        return(column_name)
    }
    paste(split_name[-length(split_name)], collapse=replicate_sep)
}

merge_replicates = function(trait_table, replicate_sep) {
    .validate_numeric_trait_table(trait_table)
    replicate_sep = .normalize_single_string_arg(
        value=replicate_sep,
        arg_name='replicate_sep',
        allow_empty=TRUE
    )
    if (replicate_sep=='') {
        return(trait_table)
    }
    if (is.null(colnames(trait_table)) || any(is.na(colnames(trait_table)) | colnames(trait_table) == '')) {
        stop('trait_table must have non-empty column names in merge_replicates().')
    }
    without_reps = vapply(
        colnames(trait_table),
        .replicate_base_name,
        character(1),
        replicate_sep=replicate_sep
    )
    if (length(unique(without_reps))==ncol(trait_table)) {
        cat(paste0('No replicate was found with --replicate_sep="', replicate_sep, '"\n'))
        return(trait_table)
    } else {
        cat(paste0('Replicates were found with --replicate_sep="', replicate_sep, '". Mean values will be used.\n'))
    }
    new_cols = unique(without_reps)
    out = data.frame(matrix(ncol=length(new_cols), nrow=nrow(trait_table)))
    colnames(out) = new_cols
    rownames(out) = rownames(trait_table)
    for (new_col in colnames(out)) {
        is_col = without_reps == new_col
        if (sum(is_col)==1) {
            values = trait_table[,is_col]
        } else {
            values = apply(trait_table[,is_col], 1, function(x){
                out = mean(x, na.rm=TRUE)
                if (is.nan(out)) {
                    return(NA_real_)
                }
                out
            })
        }
        out[,new_col] = values
    }
    return(out)
}

.trait_similarity = function(trait1, trait2, method) {
    trait1 = suppressWarnings(as.numeric(trait1))
    trait2 = suppressWarnings(as.numeric(trait2))
    if (length(trait1) != length(trait2)) {
        stop('Trait vectors must have equal length in get_high_similarity_clades().')
    }
    complete = stats::complete.cases(trait1, trait2) & is.finite(trait1) & is.finite(trait2)
    if (method == 'complementarity') {
        if (!any(complete)) {
            return(NA_real_)
        }
        return(1 - calc_complementarity(
            trait1[complete],
            trait2[complete],
            method='weighted'
        ))
    }
    if (sum(complete) < 2L) {
        return(NA_real_)
    }
    suppressWarnings(stats::cor(
        trait1[complete],
        trait2[complete],
        method=method
    ))
}

get_high_similarity_clades = function(tree, trait_table, method, threshold, verbose=FALSE, num_test=0) {
    if (length(method) != 1 || is.na(method)) {
        stop('method must be a single non-missing value in get_high_similarity_clades().')
    }
    method_name = as.character(method)
    supported_methods = c('complementarity', 'pearson', 'spearman', 'kendall')
    if (!(method_name %in% supported_methods)) {
        stop(
            'Unsupported method in get_high_similarity_clades(): ',
            method_name,
            '. Supported methods: ',
            paste(supported_methods, collapse=', ')
        )
    }
    threshold_num = suppressWarnings(as.numeric(threshold))
    if (length(threshold_num) != 1 || is.na(threshold_num) || !is.finite(threshold_num)) {
        stop('threshold must be a single finite numeric value in get_high_similarity_clades().')
    }
    num_test_numeric = suppressWarnings(as.numeric(num_test))
    if (length(num_test_numeric) != 1 || is.na(num_test_numeric) ||
            !is.finite(num_test_numeric) || num_test_numeric < 0 ||
            num_test_numeric != as.integer(num_test_numeric)) {
        stop('num_test must be a single non-negative integer in get_high_similarity_clades().')
    }
    num_test_num = as.integer(num_test_numeric)
    verbose = .normalize_single_logical_arg(
        value=verbose,
        arg_name='verbose'
    )
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
    .validate_numeric_trait_table(trait_table)
    num_all_leaves = length(tree$tip.label)
    collapse_node_nums = c()
    root_num = get_root_num(tree)
    subroot_nums = get_children_num(tree, root_num)
    next_node_nums = subroot_nums
    for (nnn in next_node_nums) {
        if (length(get_tip_labels(tree, nnn))==1) {
            next_node_nums = next_node_nums[next_node_nums!=nnn]
        }
    }
    num_processed_node = 0
    while (length(next_node_nums) > 0) {
        current_node_num = next_node_nums[1]
        tip_labels = get_tip_labels(tree, current_node_num)
        num_leaves = length(tip_labels)
        num_combinations = choose(num_leaves, 2)
        if (verbose) {
            cat('number of next_node_nums =', length(next_node_nums), 'number of leaf combinations =', num_combinations, '\n')
        }
        min_similarity = 1
        do_collapse = TRUE
        combination_index = 0L
        if (num_leaves >= 2L) {
            for (c1 in seq_len(num_leaves - 1L)) {
                for (c2 in seq.int(c1 + 1L, num_leaves)) {
                    combination_index = combination_index + 1L
                    leaf_c1 = tip_labels[[c1]]
                    leaf_c2 = tip_labels[[c2]]
                    current_similarity = .trait_similarity(
                        trait1=unlist(trait_table[leaf_c1,,drop=FALSE], use.names=FALSE),
                        trait2=unlist(trait_table[leaf_c2,,drop=FALSE], use.names=FALSE),
                        method=method_name
                    )
                    if (is.na(current_similarity)) {
                        do_collapse = FALSE
                        min_similarity = NA_real_
                    } else {
                        min_similarity = min(min_similarity, current_similarity)
                        do_collapse = min_similarity >= threshold_num
                    }
                    if (!do_collapse) {
                        if (verbose) {
                            cat(
                                'A low or undefined similarity (', min_similarity,
                                ') found at combination ', combination_index,
                                ' of ', num_combinations, '\n'
                            )
                        }
                        break
                    }
                }
                if (!do_collapse) {
                    break
                }
            }
        }
        if (length(next_node_nums)==1) {
            next_node_nums = c()
        } else {
            next_node_nums = next_node_nums[-1]
        }
        if (do_collapse) {
            cat('node_num =', current_node_num, 'size =', num_leaves, 'Min similarity =', min_similarity, '\n')
            collapse_node_nums = c(collapse_node_nums, current_node_num)
       } else {
            children_nums = get_children_num(tree, current_node_num)
            children_nums = stats::na.omit(children_nums)
            children_nums = children_nums[children_nums > num_all_leaves]
            next_node_nums = c(next_node_nums, children_nums)
            if (verbose) {
                cat('node_num =', current_node_num, 'size =', num_leaves, 'Min similarity =', min_similarity, '\n')
            }
        }
        num_processed_node = num_processed_node + 1
        if (num_processed_node%%100==0) {
            cat('processed', num_processed_node, 'nodes\n')
        }
        if (num_test_num != 0 && num_processed_node == num_test_num) {
            cat('Reaching num_test, Exiting at the ', num_test_num, ' th processed node.\n')
            break
        }
    }
    return(collapse_node_nums)

}

.collapse_clade_root_state = function(subtree_trait_values, subtree) {
    trait_values = suppressWarnings(as.numeric(subtree_trait_values))
    non_na_values = trait_values[!is.na(trait_values)]
    if (!length(non_na_values)) {
        return(NA_real_)
    }
    if (length(unique(non_na_values)) == 1L) {
        return(non_na_values[[1]])
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

collapse_clades = function(tree, trait_table, collapse_node_nums) {
    cat('number of leaves before collapse =', length(tree$tip.label), '\n')
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
    collapse_node_nums_input = collapse_node_nums
    if (length(collapse_node_nums)) {
        collapse_node_nums = suppressWarnings(as.integer(collapse_node_nums))
        if (any(is.na(collapse_node_nums))) {
            invalid_values = unique(as.character(collapse_node_nums_input[is.na(collapse_node_nums)]))
            stop(
                'collapse_node_nums must be integer node numbers. Invalid value(s): ',
                paste(invalid_values, collapse=', ')
            )
        }
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
    cat('number of leaves after collapse =', length(out_tree$tip.label), '\n')
    out = list(tree=out_tree, trait=out_trait, collapse_leaf_names=collapse_leaf_names)
    return(out)
}

map_node_num = function(tree_original, tree_collapsed, collapse_leaf_names=list(), verbose=FALSE) {
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
        if (is.null(cln) || any(is.na(cln) | trimws(cln) == '')) {
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
    original_node_nums = seq_len(max(tree_original[['edge']]))
    collapsed_node_nums = seq_len(max(tree_collapsed[['edge']]))
    original_tip_sets = .get_node_tip_sets(tree_original)
    original_signatures = vapply(
        original_tip_sets, .encode_tip_signature, character(1)
    )
    collapsed_tip_sets = .get_node_tip_sets(tree_collapsed)
    mapped_collapsed_num = rep(NA_integer_, length(original_node_nums))

    for (nn2 in collapsed_node_nums) {
        collapsed_tips = as.character(collapsed_tip_sets[[nn2]])
        collapsed_placeholders = intersect(collapsed_tips, names(collapse_leaf_names))
        expanded_tips = setdiff(collapsed_tips, collapsed_placeholders)
        for (placeholder in collapsed_placeholders) {
            expanded_tips = c(expanded_tips, as.character(collapse_leaf_names[[placeholder]]))
        }
        expanded_tips = sort(unique(expanded_tips))
        is_collapsed_tip = length(collapsed_tips) == 1L &&
            collapsed_tips[[1]] %in% names(collapse_leaf_names)

        if (is_collapsed_tip) {
            matched_original = vapply(original_tip_sets, function(tips) {
                length(tips) > 0L && all(tips %in% expanded_tips)
            }, logical(1))
        } else {
            expanded_signature = .encode_tip_signature(expanded_tips)
            matched_original = original_signatures == expanded_signature
        }
        conflicts = matched_original & !is.na(mapped_collapsed_num) & mapped_collapsed_num != nn2
        if (any(conflicts)) {
            stop(
                'Ambiguous node mapping in map_node_num() for original node(s): ',
                paste(original_node_nums[conflicts], collapse=', ')
            )
        }
        mapped_collapsed_num[matched_original] = nn2
    }

    if (anyNA(mapped_collapsed_num)) {
        stop(
            'Cannot map original node number(s) in map_node_num(): ',
            paste(original_node_nums[is.na(mapped_collapsed_num)], collapse=', ')
        )
    }
    if (verbose) {
        cat('Mapped', length(original_node_nums), 'original nodes to',
            length(unique(mapped_collapsed_num)), 'collapsed nodes.\n')
    }
    data.frame(
        tree_original=original_node_nums,
        tree_collapsed=mapped_collapsed_num,
        stringsAsFactors=FALSE
    )
}

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
    tree_table = data.frame()
    if (mode=='l1ou') {
        leaves = pcm_out$tree$tip.label
        spp = unique(.parse_species_labels(
            labels=leaves,
            species_parser=species_parser,
            sep=sep,
            output_sep=sep,
            require_gene=FALSE,
            fallback_label=TRUE
        )[['species_labels']])
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
        df[['num_species']] = length(unique(suppressWarnings(leaf2species(
            leaf_names=pcm_out[['phylo']][['tip.label']],
            species_parser=species_parser,
            sep=sep
        ))))
        df[['num_leaf']] = length(pcm_out[['phylo']][['tip.label']])
        tree_table = df
    } else {
        cat('mode "', mode, '" is not supported.\n')
    }
    return(tree_table)
}

.as_l1ou_shift_matrix = function(values, param_name, expected_rows, expected_cols) {
    mat = as.matrix(values)
    if (expected_rows == 0L) {
        return(matrix(numeric(0), nrow=0, ncol=expected_cols))
    }
    if (is.null(dim(mat)) || length(dim(mat)) != 2L) {
        stop(
            'pcm_out$',
            param_name,
            ' must be a matrix/data.frame with ',
            expected_rows,
            ' row(s) and ',
            expected_cols,
            ' column(s) in get_regime_table(mode="l1ou").'
        )
    }
    if (nrow(mat) != expected_rows || ncol(mat) != expected_cols) {
        stop(
            'pcm_out$',
            param_name,
            ' has invalid dimensions in get_regime_table(mode="l1ou"): expected ',
            expected_rows,
            'x',
            expected_cols,
            ', got ',
            nrow(mat),
            'x',
            ncol(mat),
            '.'
        )
    }
    mat
}

.as_l1ou_param_vector = function(values, param_name, expected_cols) {
    original_names = names(values)
    vec = as.vector(values)
    if (!is.null(original_names) && length(original_names) == length(vec)) {
        names(vec) = original_names
    }
    if (expected_cols == 0L) {
        return(numeric(0))
    }
    if (length(vec) == expected_cols) {
        return(vec)
    }
    if (length(vec) == 1L) {
        return(unname(rep(vec, expected_cols)))
    }
    stop(
        'pcm_out$',
        param_name,
        ' must have length 1 or ',
        expected_cols,
        ' in get_regime_table(mode="l1ou"), got ',
        length(vec),
        '.'
    )
}

.align_l1ou_trait_columns = function(param_table, param_name, trait_cols) {
    if (ncol(param_table) != length(trait_cols)) {
        stop(
            'pcm_out$',
            param_name,
            ' must have ',
            length(trait_cols),
            ' column(s) to match pcm_out$Y.'
        )
    }

    param_cols = colnames(param_table)
    if (is.null(param_cols)) {
        colnames(param_table) = trait_cols
        return(param_table)
    }
    if (anyDuplicated(param_cols)) {
        duplicated_cols = unique(param_cols[duplicated(param_cols)])
        stop(
            'Duplicate column name(s) in pcm_out$',
            param_name,
            ': ',
            paste(duplicated_cols, collapse=', ')
        )
    }

    missing_cols = setdiff(trait_cols, param_cols)
    extra_cols = setdiff(param_cols, trait_cols)
    if (length(missing_cols) || length(extra_cols)) {
        stop(
            'Column names of pcm_out$',
            param_name,
            ' must match pcm_out$Y. Missing: ',
            ifelse(length(missing_cols), paste(missing_cols, collapse=', '), 'none'),
            '. Extra: ',
            ifelse(length(extra_cols), paste(extra_cols, collapse=', '), 'none'),
            '.'
        )
    }

    param_table[, trait_cols, drop=FALSE]
}

.get_node_depths_from_root = function(tree) {
    edge = tree[['edge']]
    root_candidates = setdiff(unique(edge[,1]), unique(edge[,2]))
    if (length(root_candidates) != 1L) {
        stop('Unable to identify a unique root node while ordering shift.configuration.')
    }

    depths = rep(NA_integer_, max(edge))
    root_node = root_candidates[[1]]
    depths[root_node] = 0L
    children_by_parent = split(edge[,2], edge[,1])
    current_nodes = root_node
    while (length(current_nodes)) {
        next_nodes = integer(0)
        for (parent_node in current_nodes) {
            children = children_by_parent[[as.character(parent_node)]]
            if (!is.null(children)) {
                depths[children] = depths[parent_node] + 1L
                next_nodes = c(next_nodes, children)
            }
        }
        current_nodes = next_nodes
    }
    depths
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

get_regime_table = function(pcm_out, mode) {
    mode = .normalize_choice_arg(
        value=mode,
        arg_name='mode',
        choices=c('l1ou', 'PhylogeneticEM')
    )
    regime_table = data.frame()
    if (mode=='l1ou') {
        if (is.null(colnames(pcm_out$Y))) {
            cols = 'trait1'
        } else {
            cols = colnames(pcm_out$Y)
        }
        column_names = c("regime", "node_name", "param", cols)
        regime_table = data.frame(matrix(0,0,ncol(pcm_out$Y)+3))
        colnames(regime_table) = column_names
        shift_conf_raw = pcm_out$shift.configuration
        shift_conf = suppressWarnings(as.integer(shift_conf_raw))
        if (length(shift_conf) != length(shift_conf_raw) || any(is.na(shift_conf))) {
            stop('shift.configuration must contain integer edge indices in get_regime_table(mode="l1ou").')
        }
        if (is.null(names(shift_conf_raw))) {
            regimes = as.character(seq_along(shift_conf_raw))
        } else {
            regimes = as.character(names(shift_conf_raw))
        }
        invalid_shift_idx = shift_conf[(shift_conf < 1L) | (shift_conf > nrow(pcm_out[['tree']][['edge']]))]
        if (length(invalid_shift_idx)) {
            stop(
                'shift.configuration contains invalid edge index/indices in get_regime_table(mode="l1ou"): ',
                paste(unique(invalid_shift_idx), collapse=', ')
            )
        }
        branch_index = shift_conf
        expected_shift_count = length(branch_index)
        n_shifts_field = suppressWarnings(as.integer(pcm_out$nShifts))
        if (length(n_shifts_field) != 1 || is.na(n_shifts_field) || n_shifts_field < 0) {
            stop('pcm_out$nShifts must be a single non-negative integer in get_regime_table(mode="l1ou").')
        }
        if (n_shifts_field != expected_shift_count) {
            stop(
                'Mismatch between pcm_out$nShifts (',
                n_shifts_field,
                ') and length(pcm_out$shift.configuration) (',
                expected_shift_count,
                ') in get_regime_table(mode="l1ou").'
            )
        }
        tree = fill_node_labels(pcm_out[['tree']])
        node_names = vapply(branch_index, function(ind) {
            if (is.na(ind) || ind < 1 || ind > nrow(tree$edge)) {
                return(NA_character_)
            }
            node_values = get_node_name_by_num(phy=tree, node_num=tree$edge[ind,2])
            if (length(node_values) != 1) {
                return(NA_character_)
            }
            as.character(node_values)
        }, character(1))

        if (expected_shift_count > 0) {
            shift_values = .as_l1ou_shift_matrix(
                values=pcm_out$shift.values,
                param_name='shift.values',
                expected_rows=expected_shift_count,
                expected_cols=ncol(pcm_out$Y)
            )
            shift_means = .as_l1ou_shift_matrix(
                values=pcm_out$shift.means,
                param_name='shift.means',
                expected_rows=expected_shift_count,
                expected_cols=ncol(pcm_out$Y)
            )
            shift_values = .align_l1ou_trait_columns(
                shift_values, 'shift.values', cols
            )
            shift_means = .align_l1ou_trait_columns(
                shift_means, 'shift.means', cols
            )
            tmp_table = data.frame(
                regime=regimes,
                node_name=node_names,
                param='shift_value',
                shift_values,
                stringsAsFactors=FALSE,
                check.names=FALSE
            )
            colnames(tmp_table) = column_names
            regime_table = rbind(regime_table, tmp_table)

            tmp_table = data.frame(
                regime=regimes,
                node_name=node_names,
                param='shift_mean',
                shift_means,
                stringsAsFactors=FALSE,
                check.names=FALSE
            )
            colnames(tmp_table) = column_names
            regime_table = rbind(regime_table, tmp_table)
        }

        param_values = list(
            alpha=pcm_out$alpha,
            sigma2=pcm_out$sigma2,
            intercept=pcm_out$intercept,
            log_likelihood=pcm_out$logLik
        )
        for (param_name in names(param_values)) {
            values = .as_l1ou_param_vector(
                values=param_values[[param_name]],
                param_name=param_name,
                expected_cols=ncol(pcm_out$Y)
            )
            value_names = names(values)
            if (!is.null(value_names)) {
                if (anyDuplicated(value_names) || !setequal(value_names, cols)) {
                    stop(
                        'Names of pcm_out$', param_name,
                        ' must match pcm_out$Y columns.'
                    )
                }
                values = values[cols]
            }
            value_table = as.data.frame(
                matrix(as.numeric(values), 1, ncol(pcm_out$Y)),
                stringsAsFactors=FALSE
            )
            colnames(value_table) = cols
            tmp_table = data.frame(
                regime=NA_character_,
                node_name=NA_character_,
                param=param_name,
                value_table,
                stringsAsFactors=FALSE,
                check.names=FALSE
            )
            colnames(tmp_table) = column_names
            regime_table = rbind(regime_table, tmp_table)
        }
    } else if (mode=='PhylogeneticEM') {
        pp = .phylogeneticem_params_process(pcm_out)
        traits = rownames(pcm_out[['Y_data']])
        num_trait = length(traits)
        shift_node_indices = pp[['shifts']][['edges']]
        shift_node_nums = pcm_out[['phylo']]$edge[shift_node_indices,2]
        num_shift = length(shift_node_nums)
        df = data.frame(matrix(NA, num_shift+3, num_trait+3))
        colnames(df) = c('regime','node_name','param',traits)
        row=1
        if (num_shift) {
            df[row:(row+num_shift-1),'regime'] = 1:num_shift
            df[row:(row+num_shift-1),'node_name'] = get_node_name_by_num(pcm_out[['phylo']], shift_node_nums)
            df[row:(row+num_shift-1),'param'] = 'shift_value'
            df[row:(row+num_shift-1),traits] = t(pp[['shifts']][['values']])
        }
        row = row+num_shift - 1
        row = row+1; df[row,] = c(NA, NA, 'alpha', diag(as.matrix(pp$selection.strength)))
        row = row+1; df[row,] = c(NA, NA, 'sigma2', diag(as.matrix(pp$variance)))
        row = row+1; df[row,] = c(NA, NA, 'intercept', pp$optimal.value)
        regime_table = df
    } else {
        cat('mode "', mode, '" is not supported.\n')
    }
    rownames(regime_table) = NULL
    return(regime_table)
}

.apply_shift_conf_to_leaf_regimes = function(leaf_regimes, tree, shift_conf) {
    if (length(shift_conf) == 0) {
        return(leaf_regimes)
    }
    shift_idx = suppressWarnings(as.integer(shift_conf))
    if (any(is.na(shift_idx))) {
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

get_leaf_table = function(pcm_out, mode) {
    mode = .normalize_choice_arg(
        value=mode,
        arg_name='mode',
        choices=c('l1ou', 'PhylogeneticEM')
    )
    leaf_table = data.frame()
    if (mode=='l1ou') {
        trait_cols = colnames(pcm_out$Y)
        if (is.null(trait_cols)) {
            trait_cols = paste0("trait", seq_len(ncol(pcm_out$Y)))
        }
        column_names = c("regime", 'node_name', "param", trait_cols)
        leaf_table = data.frame(matrix(0,0,ncol(pcm_out$Y)+3))
        colnames(leaf_table) = column_names
        leaf_regimes = get_leaf_regimes(pcm_out, mode)
        param2table = list(
            Y=pcm_out$Y,
            optima=pcm_out$optima,
            mu=pcm_out$mu,
            residuals=pcm_out$residuals
        )
        for (param_name in names(param2table)) {
            param_table_raw = param2table[[param_name]]
            param_col_names = colnames(param_table_raw)
            param_table = as.data.frame(param_table_raw, stringsAsFactors=FALSE)
            if (!is.null(param_col_names)) {
                colnames(param_table) = param_col_names
            } else {
                colnames(param_table) = NULL
            }
            param_row_names = rownames(param_table)
            if (!is.null(param_row_names)) {
                if (anyDuplicated(param_row_names)) {
                    duplicated_rows = unique(param_row_names[duplicated(param_row_names)])
                    stop(
                        'Duplicate row name(s) in pcm_out$',
                        param_name,
                        ': ',
                        paste(duplicated_rows, collapse=', ')
                    )
                }
                leaf_labels = as.character(leaf_regimes[['label']])
                missing_rows = setdiff(leaf_labels, param_row_names)
                if (length(missing_rows)) {
                    stop(
                        'pcm_out$',
                        param_name,
                        ' is missing row(s) for tree tip label(s): ',
                        paste(missing_rows, collapse=', ')
                    )
                }
                param_table = param_table[leaf_labels,,drop=FALSE]
            } else if (nrow(param_table) != nrow(leaf_regimes)) {
                stop(
                    'pcm_out$',
                    param_name,
                    ' must have ',
                    nrow(leaf_regimes),
                    ' row(s) to match tree tips.'
                )
            }
            param_table = .align_l1ou_trait_columns(param_table, param_name, trait_cols)
            tmp_table = cbind(leaf_regimes, param_name, param_table)
            colnames(tmp_table) = column_names
            leaf_table = rbind(leaf_table, tmp_table)
        }
    } else if (mode=='PhylogeneticEM') {
        num_leaf = length(pcm_out[['phylo']][['tip.label']])
        traits = rownames(pcm_out[['Y_data']])
        num_trait = length(traits)
        params = c('imputed','expectations')
        num_param = length(params)
        df_leaf_regimes = get_leaf_regimes(pcm_out, mode='PhylogeneticEM')
        regimes = df_leaf_regimes[['regime']]
        df = data.frame(matrix(NA, num_leaf*num_param, num_trait+3))
        colnames(df) = c('regime','node_name','param', rownames(pcm_out[['Y_data']]))
        ind_start = 1
        for (param in params) {
            ind_end = ind_start + num_leaf - 1
            df[ind_start:ind_end,'regime'] = regimes
            df[ind_start:ind_end,'node_name'] = pcm_out[['phylo']][['tip.label']]
            df[ind_start:ind_end,'param'] = param
            df[ind_start:ind_end,traits] = t(.phylogeneticem_imputed_traits(pcm_out, trait=1:num_trait, where='tips', what=param))
            ind_start = ind_end + 1
        }
        leaf_table = df
    } else {
        cat('mode "', mode, '" is not supported.\n')
    }
    rownames(leaf_table) = NULL
    return(leaf_table)
}

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

tree_table_collapse2original = function(tree_table, tree_original, species_parser='legacy', sep='_') {
    num_leaf = length(tree_original$tip.label)
    species_names = suppressWarnings(leaf2species(
        tree_original$tip.label,
        species_parser=species_parser,
        sep=sep
    ))
    if (length(species_names) != num_leaf || !length(species_names)) {
        species_names = tree_original$tip.label
    } else if (anyNA(species_names)) {
        species_names[is.na(species_names)] = tree_original$tip.label[is.na(species_names)]
    }
    num_species = length(unique(stats::na.omit(species_names)))
    if (num_species == 0L) {
        num_species = num_leaf
    }
    out_tree_table = tree_table
    out_tree_table[,'num_leaf'] = num_leaf
    out_tree_table[,'num_species'] = num_species
    return(out_tree_table)
}

get_deepest_node_num = function(tree, node_nums) {
    if (!length(node_nums)) {
        return(NA_integer_)
    }
    node_nums = suppressWarnings(as.integer(node_nums))
    node_nums = stats::na.omit(node_nums)
    if (!length(node_nums)) {
        return(NA_integer_)
    }
    depth_values = ape::node.depth(tree)
    node_nums = node_nums[(node_nums >= 1) & (node_nums <= length(depth_values))]
    if (!length(node_nums)) {
        return(NA_integer_)
    }
    target_depth_values = depth_values[node_nums]
    deepest_node_num = node_nums[which.max(target_depth_values)]
    return(deepest_node_num)
}

regime_table_collapse2original = function(regime_table, tree_original, tree_collapsed, node_num_mapping) {
    out_regime_table = regime_table
    node_names_collapsed = unique(stats::na.omit(as.character(out_regime_table[['node_name']])))
    for (node_name_collapsed in node_names_collapsed) {
        if (node_name_collapsed %in% c(tree_collapsed$tip.label, tree_collapsed$node.label)) {
            node_num_collapsed = get_node_num_by_name(tree_collapsed, node_name_collapsed)
            if (length(node_num_collapsed) != 1) {
                next
            }
            node_num_originals = node_num_mapping[(node_num_mapping['tree_collapsed']==node_num_collapsed),'tree_original']
            node_num_original = get_deepest_node_num(tree_original, node_num_originals)
            if (is.na(node_num_original)) {
                next
            }
            node_name_original = get_node_name_by_num(tree_original, node_num_original)
            if (length(node_name_original) != 1) {
                next
            }
            is_target = (as.character(out_regime_table[['node_name']]) == node_name_collapsed)
            is_target[is.na(is_target)] = FALSE
            out_regime_table[is_target,'node_name'] = as.character(node_name_original)
        }
    }
    rownames(out_regime_table) = NULL
    return(out_regime_table)
}

leaf_table_collapse2original = function(leaf_table, tree_original, tree_collapsed, node_num_mapping) {
    num_leaf_original = length(tree_original$tip.label)
    node_names_collapsed = unique(stats::na.omit(as.character(leaf_table[['node_name']])))
    params = unique(stats::na.omit(as.character(leaf_table[['param']])))
    df = data.frame()
    for (param in params) {
        for (node_name_collapsed in node_names_collapsed) {
            if (node_name_collapsed %in% c(tree_collapsed$tip.label, tree_collapsed$node.label)) {
                node_num_collapsed = get_node_num_by_name(tree_collapsed, node_name_collapsed)
                if (length(node_num_collapsed) != 1) {
                    next
                }
                node_num_originals = node_num_mapping[(node_num_mapping['tree_collapsed']==node_num_collapsed),'tree_original']
                for (node_num_original in node_num_originals) {
                    if (node_num_original <= num_leaf_original) {
                        node_name_original = get_node_name_by_num(tree_original, node_num_original)
                        conditions = (as.character(leaf_table[['param']]) == param) &
                            (as.character(leaf_table[['node_name']]) == node_name_collapsed)
                        conditions[is.na(conditions)] = FALSE
                        row = leaf_table[conditions,]
                        if (!nrow(row)) {
                            next
                        }
                        row[,'node_name'] = node_name_original
                        df = rbind(df, row)
                    }
                }
            }
        }
    }
    rownames(df) = NULL
    return(df)
}

restore_imputed_leaves = function(leaf_table, original_trait_table) {
    traits = colnames(original_trait_table)
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
    leaf_names = as.character(out_leaf_table[[leaf_col]])
    conditions = !is.na(leaf_names) & leaf_names %in% rownames(original_trait_table) &
        !is.na(out_leaf_table$param) & out_leaf_table$param == 'imputed'
    target_rows = which(conditions)
    if (length(target_rows)) {
        out_leaf_table[target_rows, traits] = original_trait_table[
            leaf_names[target_rows], traits, drop=FALSE
        ]
    }
    return(out_leaf_table)
}

get_placeholder_leaf = function(tree, original_trait_table) {
    if (missing(tree) || !inherits(tree, 'phylo')) {
        stop('tree must be an object of class "phylo" in get_placeholder_leaf().')
    }
    if (is.null(rownames(original_trait_table))) {
        stop('original_trait_table must have row names in get_placeholder_leaf().')
    }
    out = data.frame()
    params = c('Y', 'optima', 'mu', 'residuals')
    for (param in params) {
        tmp = data.frame(
            regime=rep(0,nrow(original_trait_table)),
            node_name=rownames(original_trait_table),
            param=param
        )
        tmp = cbind(tmp, original_trait_table)
        rownames(tmp) = NULL
        out = rbind(out, tmp)
    }
    return(out)
}

get_placeholder_regime = function(tree, original_trait_table) {
    if (missing(tree) || !inherits(tree, 'phylo')) {
        stop('tree must be an object of class "phylo" in get_placeholder_regime().')
    }
    out = data.frame()
    params = c('alpha', 'sigma2', 'intercept', 'log_likelihood')
    trait_cols = colnames(original_trait_table)
    if (is.null(trait_cols)) {
        trait_cols = paste0("trait", seq_len(ncol(original_trait_table)))
    }
    for (param in params) {
        tmp = c(NA, NA, param, rep(NA, ncol(original_trait_table)))
        out = rbind(out, tmp)
    }
    colnames(out) = c('regime', 'node_name', 'param', trait_cols)
    return(out)
}

get_placeholder_tree = function(tree, original_trait_table) {
    if (missing(tree) || !inherits(tree, 'phylo')) {
        stop('tree must be an object of class "phylo" in get_placeholder_tree().')
    }
    if (missing(original_trait_table)) {
        stop('original_trait_table is required in get_placeholder_tree().')
    }
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

sort_exp = function(exp, tree, col='gene_id') {
    out_exp = exp
    if (!(col %in% colnames(out_exp))) {
        stop('Column "', col, '" is not present in exp.')
    }
    key_values = as.character(out_exp[,col])
    if (any(is.na(key_values) | trimws(key_values) == '')) {
        stop('exp contains missing/blank values in key column "', col, '".')
    }
    if (anyDuplicated(key_values)) {
        duplicated_keys = unique(key_values[duplicated(key_values)])
        stop(
            'exp contains duplicated values in key column "', col, '": ',
            paste(duplicated_keys, collapse=', ')
        )
    }
    missing_tips = setdiff(tree[['tip.label']], key_values)
    if (length(missing_tips)) {
        stop(
            'exp is missing rows for tree tip label(s): ',
            paste(missing_tips, collapse=', ')
        )
    }
    rownames(out_exp) = key_values
    out_exp = out_exp[tree[['tip.label']],,drop=FALSE]
    out_exp[,col] = tree[['tip.label']]
    rownames(out_exp) = NULL
    return(out_exp)
}

phylogenetic_imputation = function(tree, trait_table) {
    .validate_numeric_trait_table(trait_table)
    if (is.null(rownames(trait_table)) || any(is.na(rownames(trait_table)) | rownames(trait_table) == '')) {
        stop('trait_table must have non-empty row names matching tree tip labels.')
    }
    trait_table2 = cbind(species=rownames(trait_table), trait_table)
    trait_table2 = sort_exp(trait_table2, tree, col='species')
    num_missing = sum(is.na(trait_table2[setdiff(colnames(trait_table2), 'species')]))
    num_all = nrow(trait_table2) * (ncol(trait_table2) - 1L)
    cat('Phylogenetic imputation with phylopars:', num_missing, '/', num_all, 'traits will be imputed.\n')
    rp_out = .rphylopars_phylopars(tree=tree, trait_data=trait_table2, phylo_correlated=TRUE, pheno_correlated=TRUE)
    imputed_matrix = data.frame(rp_out[['anc_recon']])
    imputed_matrix = imputed_matrix[tree[['tip.label']],,drop=FALSE]
    return(imputed_matrix)
}

get_expression_bases = function(trait_table, replicate_sep) {
    replicate_sep = .normalize_single_string_arg(
        value=replicate_sep,
        arg_name='replicate_sep',
        allow_empty=TRUE
    )
    if (replicate_sep=='') {
        out = unique(colnames(trait_table))
        out = out[out!='gene']
        return(out)
    }
    if (is.null(colnames(trait_table)) || any(is.na(colnames(trait_table)) | colnames(trait_table) == '')) {
        stop('trait_table must have non-empty column names in get_expression_bases().')
    }
    without_reps = vapply(
        colnames(trait_table),
        .replicate_base_name,
        character(1),
        replicate_sep=replicate_sep
    )
    out = unique(without_reps)
    out = out[out!='gene']
    return(out)
}

count_foreground_lineage = function(tree, trait) {
    num_fg_lineage = list()
    if (!inherits(tree, 'phylo') || is.null(tree[['edge']]) || !nrow(tree[['edge']])) {
        stop('tree must be a non-empty object of class "phylo".')
    }
    if (anyDuplicated(tree[['tip.label']])) {
        stop('tree must have unique tip labels in count_foreground_lineage().')
    }
    if (!("species" %in% colnames(trait))) {
        stop('trait must contain a "species" column.')
    }
    trait_cols = setdiff(colnames(trait), "species")
    if (!length(trait_cols)) {
        stop('trait must contain at least one foreground trait column.')
    }
    species_names = as.character(trait[['species']])
    if (any(is.na(species_names) | trimws(species_names) == '')) {
        stop('trait$species must contain only non-missing, non-empty names.')
    }
    if (anyDuplicated(species_names)) {
        stop(
            'trait$species contains duplicated names: ',
            paste(unique(species_names[duplicated(species_names)]), collapse=', ')
        )
    }
    unknown_species = setdiff(stats::na.omit(species_names), tree[['tip.label']])
    if (length(unknown_species)) {
        stop(
            'trait$species contains names not present in tree tip labels: ',
            paste(unique(unknown_species), collapse=', ')
        )
    }
    missing_species = setdiff(tree[['tip.label']], species_names)
    if (length(missing_species)) {
        stop(
            'trait is missing rows for tree tip label(s): ',
            paste(missing_species, collapse=', ')
        )
    }
    node_nums = seq_len(max(tree[['edge']]))
    num_tip = length(tree[['tip.label']])
    children_by_parent = split(tree[['edge']][,2], tree[['edge']][,1])
    parent_by_child = stats::setNames(tree[['edge']][,1], tree[['edge']][,2])
    root_num = get_root_num(tree)
    if (length(root_num) != 1L) {
        stop('tree must have exactly one root in count_foreground_lineage().')
    }
    traversal_order = integer(length(node_nums))
    traversal_order[[1]] = root_num
    queue_size = 1L
    next_index = 1L
    while (next_index <= queue_size) {
        children = children_by_parent[[as.character(traversal_order[[next_index]])]]
        if (length(children)) {
            target_indices = queue_size + seq_along(children)
            traversal_order[target_indices] = children
            queue_size = queue_size + length(children)
        }
        next_index = next_index + 1L
    }
    traversal_order = traversal_order[seq_len(queue_size)]
    if (queue_size != length(node_nums) || anyDuplicated(traversal_order) ||
            !setequal(traversal_order, node_nums)) {
        stop('tree has an invalid or disconnected topology in count_foreground_lineage().')
    }
    for (trait_col in trait_cols) {
        trait_values = trait[[trait_col]]
        if (is.factor(trait_values)) {
            trait_values = as.character(trait_values)
        }
        if (is.logical(trait_values)) {
            trait_values = as.numeric(trait_values)
        } else {
            suppressWarnings(trait_values <- as.numeric(trait_values))
        }
        is_binary = all(is.na(trait_values) | (trait_values %in% c(0, 1)))
        if (!is_binary) {
            cat(trait_col, ': count_foreground_lineage() supports 0/1 traits. Returning NA.\n')
            num_fg_lineage[[trait_col]] = NA
            next
        }
        if (anyNA(trait_values)) {
            warning(
                trait_col,
                ': missing foreground states make the lineage count undefined. Returning NA.'
            )
            num_fg_lineage[[trait_col]] = NA
            next
        }
        fg_spp = species_names[(!is.na(trait_values)) & (trait_values == 1)]
        is_fg_only_clade = logical(length(node_nums))
        is_fg_only_clade[seq_len(num_tip)] = tree[['tip.label']] %in% fg_spp
        for (node_num in rev(traversal_order)) {
            if (node_num <= num_tip) {
                next
            }
            children = children_by_parent[[as.character(node_num)]]
            is_fg_only_clade[[node_num]] = length(children) > 0L &&
                all(is_fg_only_clade[children])
        }
        is_fg_stem = logical(length(node_nums))
        is_fg_stem[[root_num]] = is_fg_only_clade[[root_num]]
        nonroot_nodes = setdiff(node_nums, root_num)
        nonroot_parents = unname(parent_by_child[as.character(nonroot_nodes)])
        if (anyNA(nonroot_parents)) {
            stop('Invalid tree topology in count_foreground_lineage().')
        }
        is_fg_stem[nonroot_nodes] = is_fg_only_clade[nonroot_nodes] &
            !is_fg_only_clade[nonroot_parents]
        num_fg_lineage[[trait_col]] = sum(is_fg_stem)
    }
    return(num_fg_lineage)
}
