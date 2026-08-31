
#' Transfer internal-node labels between congruent trees
#'
#' @param phy_from A labeled `phylo` tree.
#' @param phy_to A `phylo` tree with the same tips.
#' @return `phy_to` with labels transferred by clade signature.
#' @export
transfer_node_labels = function(phy_from, phy_to) {
    out_phy_to = phy_to
    if (!setequal(phy_from$tip.label, out_phy_to$tip.label)) {
        stop('phy_from and phy_to must contain the same tip labels.')
    }
    if (anyDuplicated(phy_from$tip.label) || anyDuplicated(out_phy_to$tip.label)) {
        stop('phy_from and phy_to must have unique tip labels.')
    }
    if (is.null(phy_from$node.label)) {
        stop('phy_from has no node labels to transfer.')
    }
    num_tip_from = length(phy_from$tip.label)
    num_tip_to = length(out_phy_to$tip.label)
    num_int_from = as.integer(phy_from$Nnode)
    num_int_to = as.integer(out_phy_to$Nnode)
    if (is.null(out_phy_to$node.label)) {
        out_phy_to$node.label = rep(NA_character_, num_int_to)
    } else {
        out_phy_to$node.label = as.character(out_phy_to$node.label)
        if (length(out_phy_to$node.label) < num_int_to) {
            out_phy_to$node.label = c(
                out_phy_to$node.label,
                rep(NA_character_, num_int_to - length(out_phy_to$node.label))
            )
        } else if (length(out_phy_to$node.label) > num_int_to) {
            out_phy_to$node.label = out_phy_to$node.label[seq_len(num_int_to)]
        }
    }
    from_node_nums = num_tip_from + seq_len(num_int_from)
    to_node_nums = num_tip_to + seq_len(num_int_to)
    from_signatures = .get_node_tip_signatures(phy_from)[as.character(from_node_nums)]
    to_signatures = .get_node_tip_signatures(out_phy_to)[as.character(to_node_nums)]
    if (anyDuplicated(from_signatures)) {
        stop('phy_from contains duplicated internal clade signatures.')
    }
    from_labels_by_signature = stats::setNames(
        as.character(phy_from$node.label[seq_len(num_int_from)]),
        from_signatures
    )
    matched_labels = unname(from_labels_by_signature[to_signatures])
    has_match = !is.na(matched_labels)
    if (any(has_match)) {
        out_phy_to$node.label[has_match] = matched_labels[has_match]
    }
    return(out_phy_to)
}


#' Fill missing internal-node labels
#'
#' @param phy A `phylo` tree.
#' @param verbose Whether to emit a progress message.
#' @return A `phylo` tree with unique generated labels for missing nodes.
#' @export
fill_node_labels = function(phy, verbose=FALSE) {
    out_phy = phy
    verbose = .normalize_single_logical_arg(verbose, 'verbose')
    num_internal = as.integer(out_phy[['Nnode']])
    if (length(num_internal) != 1 || is.na(num_internal) || num_internal < 0) {
        num_internal = max(out_phy[['edge']]) - length(out_phy[['tip.label']])
    }
    if (is.null(out_phy[['node.label']])) {
        nl = rep(NA_character_, num_internal)
    } else {
        nl = as.character(out_phy[['node.label']])
        if (length(nl) < num_internal) {
            nl = c(nl, rep(NA_character_, num_internal - length(nl)))
        } else if (length(nl) > num_internal) {
            nl = nl[seq_len(num_internal)]
        }
    }
    out_phy[['node.label']] = nl
    is_missing = is.na(nl) | (trimws(nl) == '')
    if (sum(is_missing)==0) {
        return(out_phy)
    }
    missing_index = (seq_along(nl))[is_missing]
    if (verbose) {
        message('Filling ', length(missing_index), ' node names.')
    }
    used_names = unique(c(out_phy[['tip.label']], nl[!is_missing]))
    candidates = paste0('n', seq.int(0L, sum(is_missing) + length(used_names)))
    candidates = candidates[!candidates %in% used_names]
    out_phy[['node.label']][missing_index] = candidates[seq_along(missing_index)]
    return(out_phy)
}
