
#' Restore original-tree counts in a model summary
#'
#' @param tree_table A comparative-model summary table.
#' @param tree_original The original `phylo` tree.
#' @param species_parser Species-label convention. Species-only labels are
#'   accepted; labels without a genus/species pair are counted verbatim.
#' @param sep Literal separator in leaf labels.
#' @return `tree_table` with original leaf and species counts.
#' @export
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


#' Select the deepest candidate node
#'
#' @param tree A `phylo` tree.
#' @param node_nums Candidate integer node numbers.
#' @return A single integer node number or `NA`.
#' @export
get_deepest_node_num = function(tree, node_nums) {
    if (!length(node_nums)) {
        return(NA_integer_)
    }
    node_nums = .normalize_integerish(
        node_nums, 'node_nums', allow_na=TRUE, allow_empty=TRUE
    )
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


#' Restore original node names in a regime table
#'
#' @param regime_table A collapsed-tree regime table.
#' @param tree_original Original `phylo` tree.
#' @param tree_collapsed Collapsed `phylo` tree.
#' @param node_num_mapping Mapping returned by [map_node_num()].
#' @return A regime table using original-tree node names.
#' @export
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


#' Expand a collapsed-tree leaf table
#'
#' @param leaf_table A collapsed-tree leaf table.
#' @param tree_original Original `phylo` tree.
#' @param tree_collapsed Collapsed `phylo` tree.
#' @param node_num_mapping Mapping returned by [map_node_num()].
#' @return A leaf table expanded to original tips.
#' @export
leaf_table_collapse2original = function(leaf_table, tree_original, tree_collapsed, node_num_mapping) {
    num_leaf_original = length(tree_original$tip.label)
    node_names_collapsed = unique(stats::na.omit(as.character(leaf_table[['node_name']])))
    params = unique(stats::na.omit(as.character(leaf_table[['param']])))
    rows = list()
    row_index = 0L
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
                        row_index = row_index + 1L
                        rows[[row_index]] = row
                    }
                }
            }
        }
    }
    df = if (length(rows)) do.call(rbind, rows) else data.frame()
    rownames(df) = NULL
    return(df)
}
