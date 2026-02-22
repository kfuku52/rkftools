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

remove_invariant_traits = function(trait_table, small_dif=0.001) {
    stopifnot(!is.null(dim(trait_table)))
    out_trait_table = trait_table
    num_traits = ncol(trait_table)
    trait_names = colnames(trait_table)
    if (is.null(trait_names)) {
        trait_names = as.character(seq_len(num_traits))
    }

    is_small_dif = logical(num_traits)
    for (i in seq_len(num_traits)) {
        trait_values = trait_table[,i]
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

merge_replicates = function(trait_table, replicate_sep) {
    replicate_sep = .normalize_single_string_arg(
        value=replicate_sep,
        arg_name='replicate_sep',
        allow_empty=TRUE
    )
    if (replicate_sep=='') {
        return(trait_table)
    }
    without_reps = c()
    for (col in colnames(trait_table)) {
        split = strsplit(col, replicate_sep, fixed=TRUE)[[1]]
        num_split = length(split)
        if (num_split == 1) {
            without_rep = split
        } else if (num_split > 1) {
            without_rep = split[1:(num_split-1)]
        }
        if (length(without_rep)>1) {
            cat(col, without_rep, '\n')
        }
        without_reps = c(without_reps, without_rep)
    }
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
        is_col = (colnames(trait_table) == new_col) |
            startsWith(colnames(trait_table), paste0(new_col, replicate_sep))
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
    if (length(threshold_num) != 1 || is.na(threshold_num)) {
        stop('threshold must be a single numeric value in get_high_similarity_clades().')
    }
    num_test_num = suppressWarnings(as.integer(num_test))
    if (length(num_test_num) != 1 || is.na(num_test_num) || num_test_num < 0) {
        stop('num_test must be a single non-negative integer in get_high_similarity_clades().')
    }
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
        leaf_combinations = expand.grid(
            c1=seq_len(num_leaves), c2=seq_len(num_leaves),
            KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE
        )
        leaf_combinations = leaf_combinations[leaf_combinations[['c1']] < leaf_combinations[['c2']],,drop=FALSE]
        rownames(leaf_combinations) = NULL
        num_combinations = nrow(leaf_combinations)
        if (verbose) {
            cat('number of next_node_nums =', length(next_node_nums), 'number of leaf combinations =', num_combinations, '\n')
        }
        min_similarity = 1
        do_collapse = TRUE
        for (i in seq_len(num_combinations)) {
            current_similarity = 1
            leaf_c1 = tip_labels[leaf_combinations[i,'c1']]
            leaf_c2 = tip_labels[leaf_combinations[i,'c2']]
            trait_c1 = unlist(trait_table[leaf_c1,])
            trait_c2 = unlist(trait_table[leaf_c2,])
            if (method_name=='complementarity') {
                current_similarity = 1 - calc_complementarity(trait_c1, trait_c2, method='weighted')
            } else {
                current_similarity = suppressWarnings(stats::cor(trait_c1, trait_c2, method=method_name))
                current_similarity = ifelse(is.na(current_similarity), 1, current_similarity)
            }
            min_similarity = min(min_similarity, current_similarity)
            if (min_similarity < threshold_num) {
                do_collapse = FALSE
                if (verbose) {
                    cat('A low similarity (', min_similarity, ') found at', i, 'th leaf combinations (', num_combinations, ')\n')
                }
                break
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
        error=function(e) NULL
    )
    if (is.null(ace_result) || is.null(ace_result$ace) || !length(ace_result$ace)) {
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
    retain_tips = list()
    drop_tips = list()
    collapse_leaf_names = list()
    for (cnn in collapse_node_nums) {
        collapse_leaf_names[[as.character(cnn)]] = get_tip_labels(tree, cnn)
        retain_tips[[cnn]] = collapse_leaf_names[[as.character(cnn)]][1]
        drop_tips[[cnn]] = collapse_leaf_names[[as.character(cnn)]][2:length(collapse_leaf_names[[as.character(cnn)]])]
    }
    out_tree = tree
    out_trait = trait_table_aligned
    for (cnn in collapse_node_nums) {
        out_tree = ape::drop.tip(out_tree, tip=drop_tips[[cnn]], trim.internal=TRUE)
        out_tree$tip.label[out_tree$tip.label==retain_tips[[cnn]]] = cnn
        subtree = ape::extract.clade(tree, cnn)
        subtree_traits = trait_table_aligned[subtree$tip.label,,drop=FALSE]
        subtree_root_states = c()
        for (i in seq_len(ncol(trait_table_aligned))) {
            subtree_root_state = .collapse_clade_root_state(subtree_traits[,i], subtree)
            subtree_root_states = c(subtree_root_states, subtree_root_state)
        }
        out_trait = out_trait[!(rownames(out_trait) %in% subtree$tip.label),,drop=FALSE]
        subtree_root_row = as.data.frame(t(subtree_root_states), stringsAsFactors=FALSE)
        colnames(subtree_root_row) = colnames(out_trait)
        out_trait = rbind(out_trait, subtree_root_row)
        rownames(out_trait)[nrow(out_trait)] = cnn
    }
    out_trait = out_trait[out_tree$tip.label,,drop=FALSE]
    cat('number of leaves after collapse =', length(out_tree$tip.label), '\n')
    out = list(tree=out_tree, trait=out_trait, collapse_leaf_names=collapse_leaf_names)
    return(out)
}

map_node_num = function(tree_original, tree_collapsed, collapse_leaf_names=list(), verbose=FALSE) {
    # tree_collapsed should be a collapsed tree
    stopifnot(length(tree_original$tip.label)>=length(tree_collapsed$tip.label))
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
    df = data.frame(matrix(NA, max(tree_original$edge[,1]), 2))
    colnames(df) = c('tree_original', 'tree_collapsed')
    tree_original_node_nums = 1:max(tree_original$edge[,1])
    tree_collapsed_node_nums = 1:max(tree_collapsed$edge[,1])
    while (length(tree_original_node_nums)>0) {
        nn1 = tree_original_node_nums[1]
        nn1_leaf_names = sort(get_tip_labels(tree_original, nn1))
        matched = FALSE
        for (i in seq_along(tree_collapsed_node_nums)) {
            nn2 = tree_collapsed_node_nums[i]
            nn2_leaf_names = sort(get_tip_labels(tree_collapsed, nn2))
            collapsed_node_nums = nn2_leaf_names[as.character(nn2_leaf_names) %in% names(collapse_leaf_names)]
            any_collapse = (length(collapsed_node_nums)>0)
            is_collapse_nn2 = (length(nn2_leaf_names)==1)&(any_collapse)
            if (any_collapse) {
                nn2_leaf_names = nn2_leaf_names[!nn2_leaf_names %in% collapsed_node_nums]
                for (cnn in collapsed_node_nums) {
                    nn2_leaf_names = c(nn2_leaf_names, collapse_leaf_names[[as.character(cnn)]])
                }
                nn2_leaf_names = sort(nn2_leaf_names)
            }
            if (is_collapse_nn2) {
                flag = all(nn1_leaf_names %in% nn2_leaf_names)
            } else if (length(nn1_leaf_names)==length(nn2_leaf_names)) {
                flag = identical(nn1_leaf_names, nn2_leaf_names)
            } else {
                flag = FALSE
            }
            if (verbose) {
                noncollapse_match = identical(nn1_leaf_names, nn2_leaf_names)
                cat('nn1:', nn1, 'nn2:', nn2,
                    'collapse match:', all(nn1_leaf_names %in% nn2_leaf_names),
                    'noncollapse match:', noncollapse_match,
                    'is_collapse_nn2:', is_collapse_nn2, '\n')
            }
            if (flag) {
                matched = TRUE
                df[nn1,] = c(nn1,nn2)
                tree_original_node_nums = tree_original_node_nums[-1]
                if (!is_collapse_nn2) {
                    tree_collapsed_node_nums = tree_collapsed_node_nums[-i]
                    tree_collapsed_node_nums = tree_collapsed_node_nums[tree_collapsed_node_nums!=nn2]
                }
                break
            }
        }
        if (!matched) {
            stop(
                'Cannot map node number in map_node_num(). First unmatched original node: ',
                nn1,
                '. Remaining original nodes: ',
                paste(tree_original_node_nums, collapse=', '),
                '. Remaining collapsed nodes: ',
                paste(tree_collapsed_node_nums, collapse=', '),
                '.'
            )
        }
    }
    num_match_before = nrow(df)
    df = df[(!is.na(df[['tree_original']])),]
    df = df[(!is.na(df[['tree_collapsed']])),]
    num_match_after = nrow(df)
    if (num_match_before!=num_match_after) {
        num_na = num_match_before-num_match_after
        cat('Warning:', num_na, 'NA are produced in map_node_num().\n')
    }
    return(df)
}

get_tree_table = function(pcm_out, mode) {
    mode = .normalize_single_string_arg(
        value=mode,
        arg_name='mode',
        allow_empty=FALSE
    )
    tree_table = data.frame()
    if (mode=='l1ou') {
        leaves = pcm_out$tree$tip.label
        leaf2sp = function(leaf_name) {
            leaf_split = strsplit(leaf_name, "_")
            sp = paste(leaf_split[[1]][1], leaf_split[[1]][2], sep="_")
            return(sp)
        }
        spp = unique(sapply(leaves, leaf2sp))
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
        df[['num_species']] = length(unique(suppressWarnings(leaf2species(leaf_names=pcm_out[['phylo']][['tip.label']]))))
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
    vec = as.vector(values)
    if (expected_cols == 0L) {
        return(numeric(0))
    }
    if (length(vec) == expected_cols) {
        return(vec)
    }
    if (length(vec) == 1L) {
        return(rep(vec, expected_cols))
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

get_regime_table = function(pcm_out, mode) {
    mode = .normalize_single_string_arg(
        value=mode,
        arg_name='mode',
        allow_empty=FALSE
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
        shift_conf_raw = sort(pcm_out$shift.configuration, decreasing=TRUE, na.last=TRUE)
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
            tmp_table = cbind(regimes, node_names, "shift_value", shift_values)
            colnames(tmp_table) = column_names
            regime_table = rbind(regime_table, tmp_table)

            tmp_table = cbind(regimes, node_names, "shift_mean", shift_means)
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
            tmp_table = cbind(NA, NA, param_name, matrix(values, 1, ncol(pcm_out$Y)))
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
    mode = .normalize_single_string_arg(
        value=mode,
        arg_name='mode',
        allow_empty=FALSE
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
        shift_conf = sort(pcm_out[['shift.configuration']], decreasing=TRUE, na.last=TRUE)
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
    mode = .normalize_single_string_arg(
        value=mode,
        arg_name='mode',
        allow_empty=FALSE
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
            param_table = param2table[[param_name]]
            param_table = as.data.frame(param_table, stringsAsFactors=FALSE)
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
            if (is.null(colnames(param_table))) {
                colnames(param_table) = trait_cols
            }
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
    mode = .normalize_single_string_arg(
        value=mode,
        arg_name='mode',
        allow_empty=FALSE
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

tree_table_collapse2original = function(tree_table, tree_original) {
    num_leaf = length(tree_original$tip.label)
    species_names = suppressWarnings(leaf2species(tree_original$tip.label))
    if (length(species_names) != num_leaf || !length(species_names) || all(is.na(species_names))) {
        species_names = tree_original$tip.label
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
    leaf_names = as.character(out_leaf_table[[leaf_col]])
    for (leaf_name in leaf_names) {
        if (leaf_name %in% rownames(original_trait_table)) {
            conditions = ((as.character(out_leaf_table[[leaf_col]]) == leaf_name) &
                (out_leaf_table$param == 'imputed'))
            conditions[is.na(conditions)] = FALSE
            out_leaf_table[conditions, traits] = original_trait_table[leaf_name, traits]
        }
    }
    return(out_leaf_table)
}

get_placeholder_leaf = function(tree, original_trait_table) {
    stopifnot(!missing(tree))
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
    stopifnot(!missing(tree))
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
    stopifnot(!missing(tree))
    stopifnot(!missing(original_trait_table))
    out = data.frame()
    params = c('num_shift','num_regime','num_conv_regime','num_uniq_regime','num_species','num_leaf','model_score')
    out = rbind(out, rep(0, length(params)))
    colnames(out) = params
    return(out)
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
    trait_table2 = cbind(species=rownames(trait_table), trait_table)
    trait_table2 = sort_exp(trait_table2, tree, col='species')
    num_missing = sum(is.na(trait_table2))
    num_all = nrow(trait_table2) * ncol(trait_table2)
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
    without_reps = c()
    for (col in colnames(trait_table)) {
        split = strsplit(col, replicate_sep, fixed=TRUE)[[1]]
        num_split = length(split)
        if (num_split == 1) {
            without_rep = split
        } else if (num_split > 1) {
            without_rep = split[1:(num_split-1)]
        }
        without_reps = c(without_reps, without_rep)
    }
    out = unique(without_reps)
    out = out[out!='gene']
    return(out)
}

count_foreground_lineage = function(tree, trait) {
    num_fg_lineage = list()
    if (!("species" %in% colnames(trait))) {
        stop('trait must contain a "species" column.')
    }
    trait_cols = setdiff(colnames(trait), "species")
    if (!length(trait_cols)) {
        stop('trait must contain at least one foreground trait column.')
    }
    species_names = as.character(trait[['species']])
    unknown_species = setdiff(stats::na.omit(species_names), tree[['tip.label']])
    if (length(unknown_species)) {
        stop(
            'trait$species contains names not present in tree tip labels: ',
            paste(unique(unknown_species), collapse=', ')
        )
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
        fg_spp = species_names[(!is.na(trait_values)) & (trait_values == 1)]
        fg_tip_nums = get_node_num_by_name(phy=tree, node_name=fg_spp)
        node_nums = 1:max(tree[['edge']])
        is_fg_only_clade = rep(NA, length(node_nums))
        for (node_num in node_nums) {
            if (node_num <= length(tree[['tip.label']])) {
                clade_leaf_nums = node_num
            } else {
                clade_leaf_nums = get_descendent_num(phy=tree, node_num=node_num, leaf_only=TRUE)
            }
            is_all_leaf_fg = all(clade_leaf_nums %in% fg_tip_nums)
            is_fg_only_clade[node_num] = is_all_leaf_fg
        }
        stopifnot(!any(is.na(is_fg_only_clade)))
        is_fg_stem = rep(NA, length(node_nums))
        for (node_num in node_nums) {
            is_node_fg = is_fg_only_clade[node_num]
            if (is_root(phy=tree, node_num=node_num)) {
                if (is_node_fg) {
                    is_fg_stem[node_num] = TRUE
                } else {
                    is_fg_stem[node_num] = FALSE
                }
                next
            }
            parent_node_num = get_parent_num(phy=tree, node_num=node_num)
            is_parent_bg = !is_fg_only_clade[parent_node_num]
            if (is_node_fg & is_parent_bg) {
                is_fg_stem[node_num] = TRUE
            } else {
                is_fg_stem[node_num] = FALSE
            }
        }
        stopifnot(!any(is.na(is_fg_stem)))
        num_fg_lineage[[trait_col]] = sum(is_fg_stem)
    }
    return(num_fg_lineage)
}
