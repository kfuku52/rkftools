
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


#' Find clades with highly similar trait profiles
#'
#' Pairwise similarities are cached across nested clades.
#'
#' @param tree A rooted `phylo` tree with unique tips.
#' @param trait_table A numeric trait table named by tree tips.
#' @param method One of `"complementarity"`, `"pearson"`, `"spearman"`, or
#'   `"kendall"`.
#' @param threshold Similarity threshold in the method's valid range.
#' @param verbose Whether to emit progress messages.
#' @param num_test Optional maximum number of processed nodes; zero is unlimited.
#' @return Integer node numbers whose clades meet the threshold.
#' @export
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
    threshold_num = .normalize_finite_numeric_scalar(
        threshold,
        'threshold',
        min_value=if (method_name == 'complementarity') 0 else -1,
        max_value=1
    )
    num_test_num = tryCatch(
        .normalize_integerish(num_test, 'num_test', min_value=0L),
        error=function(e) integer(0)
    )
    if (length(num_test_num) != 1L) {
        stop('num_test must be a single non-negative integer in get_high_similarity_clades().')
    }
    verbose = .normalize_single_logical_arg(
        value=verbose,
        arg_name='verbose'
    )
    .validate_phylo_input(tree, context='tree', rooted=TRUE, unique_tips=TRUE)
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
    trait_table = trait_table[tree[['tip.label']],,drop=FALSE]
    num_all_leaves = length(tree$tip.label)
    collapse_node_nums = c()
    root_num = get_root_num(tree)
    subroot_nums = get_children_num(tree, root_num)
    next_node_nums = subroot_nums
    tip_sets = .get_node_tip_sets(tree)
    similarity_cache = new.env(parent=emptyenv(), hash=TRUE)
    get_pair_similarity = function(leaf1, leaf2) {
        key = .encode_tip_signature(c(leaf1, leaf2))
        if (!exists(key, envir=similarity_cache, inherits=FALSE)) {
            value = .trait_similarity(
                trait1=unlist(trait_table[leaf1,,drop=FALSE], use.names=FALSE),
                trait2=unlist(trait_table[leaf2,,drop=FALSE], use.names=FALSE),
                method=method_name
            )
            assign(key, value, envir=similarity_cache)
        }
        get(key, envir=similarity_cache, inherits=FALSE)
    }
    for (nnn in next_node_nums) {
        if (length(get_tip_labels(tree, nnn))==1) {
            next_node_nums = next_node_nums[next_node_nums!=nnn]
        }
    }
    num_processed_node = 0
    while (length(next_node_nums) > 0) {
        current_node_num = next_node_nums[1]
        tip_labels = as.character(tip_sets[[current_node_num]])
        num_leaves = length(tip_labels)
        num_combinations = choose(num_leaves, 2)
        if (verbose) {
            message('Queued nodes: ', length(next_node_nums),
                '; leaf combinations: ', num_combinations)
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
                    current_similarity = get_pair_similarity(leaf_c1, leaf_c2)
                    if (is.na(current_similarity)) {
                        do_collapse = FALSE
                        min_similarity = NA_real_
                    } else {
                        min_similarity = min(min_similarity, current_similarity)
                        do_collapse = min_similarity >= threshold_num
                    }
                    if (!do_collapse) {
                        if (verbose) {
                            message(
                                'A low or undefined similarity (', min_similarity,
                                ') was found at combination ', combination_index,
                                ' of ', num_combinations, '.'
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
            if (verbose) {
                message('node_num = ', current_node_num, '; size = ', num_leaves,
                    '; min similarity = ', min_similarity)
            }
            collapse_node_nums = c(collapse_node_nums, current_node_num)
       } else {
            children_nums = get_children_num(tree, current_node_num)
            children_nums = stats::na.omit(children_nums)
            children_nums = children_nums[children_nums > num_all_leaves]
            next_node_nums = c(next_node_nums, children_nums)
            if (verbose) {
                message('node_num = ', current_node_num, '; size = ', num_leaves,
                    '; min similarity = ', min_similarity)
            }
        }
        num_processed_node = num_processed_node + 1
        if (verbose && num_processed_node%%100==0) {
            message('Processed ', num_processed_node, ' nodes.')
        }
        if (num_test_num != 0 && num_processed_node == num_test_num) {
            if (verbose) {
                message('Reached num_test after ', num_test_num, ' processed nodes.')
            }
            break
        }
    }
    return(collapse_node_nums)

}
