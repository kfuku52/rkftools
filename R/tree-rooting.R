
#' Get the smallest root-child clade
#'
#' @param phy A rooted `phylo` tree.
#' @return Tip labels in the smallest root-child clade.
#' @export
get_outgroup = function(phy) {
    if (!ape::is.rooted(phy)) {
        stop('phy is unrooted. get_outgroup() requires a rooted tree.')
    }
    root_num = get_root_num(phy)
    if (length(root_num) != 1) {
        stop('Unable to identify a unique root node in phy.')
    }
    children_nums = get_children_num(phy, root_num)
    if (length(children_nums) < 2) {
        stop('Root node must have at least two child clades.')
    }
    outgroup_labels = phy[['tip.label']]
    for (cn in children_nums) {
        og_labels = get_tip_labels(phy, cn)
        if (length(outgroup_labels)>length(og_labels)) {
            outgroup_labels = og_labels
        }
    }
    return(outgroup_labels)
}


#' Compare root splits between trees
#'
#' @param phy1,phy2 Rooted `phylo` trees containing the same tips.
#' @return A single logical value.
#' @export
is_same_root = function(phy1, phy2) {
    if (! ape::is.rooted(phy1)) {
        stop('phy1 is unrooted.')
    }
    if (! ape::is.rooted(phy2)) {
        stop('phy2 is unrooted.')
    }
    if (! identical(sort(phy1$tip.label), sort(phy2$tip.label))) {
        stop('phy1 and phy2 have different sets of leaves.')
    }
    root1 = get_root_num(phy1)
    root2 = get_root_num(phy2)
    if (length(root1) != 1) {
        stop('Unable to identify a unique root node in phy1.')
    }
    if (length(root2) != 1) {
        stop('Unable to identify a unique root node in phy2.')
    }
    get_root_child_signatures = function(phy, root_num) {
        children = get_children_num(phy, root_num)
        if (!length(children)) {
            return(character(0))
        }
        signatures = vapply(children, function(cn) {
            paste(sort(get_tip_labels(phy, cn)), collapse='\r')
        }, character(1))
        sort(signatures)
    }
    sig1 = get_root_child_signatures(phy1, root1)
    sig2 = get_root_child_signatures(phy2, root2)
    return(identical(sig1, sig2))
}


#' Locate one tree's root split in another tree
#'
#' For a binary target root, counts descendant tips and target-side membership
#' in one traversal of `phy1`. For a multifurcating target root, compares the
#' component tip sets at candidate nodes. Neither path repeatedly reroots `phy1`
#' or creates a worker pool.
#'
#' @param phy1 A `phylo` tree in which to locate the split.
#' @param phy2 A rooted `phylo` tree providing the target root partition.
#' @param nslots Legacy worker-count argument; accepted but no workers are used.
#' @param mode Return `"node_num"` or edge `"index"`.
#' @return An integer node number or edge index, or `NA` when unmatched. When a
#'   multifurcating target root matches the current root of `phy1`, `mode =
#'   "index"` returns `NA` because the root has no incoming edge.
#' @export
get_phy2_root_in_phy1 = function(phy1, phy2, nslots=NULL, mode=c("node_num", "index")) {
    mode_name = match.arg(mode)
    .validate_phylo_input(phy1, context='phy1', unique_tips=TRUE)
    .validate_phylo_input(phy2, context='phy2', unique_tips=TRUE)
    # nslots is retained for API compatibility. Both root-matching paths run
    # serially without rerooting phy1 or creating a worker pool.
    if (!is.null(nslots)) {
        .resolve_parallel_cores(
            requested=nslots,
            max_tasks=nrow(phy1[['edge']]),
            auto_when_missing=FALSE
        )
    }
    if (! identical(sort(phy1$tip.label), sort(phy2$tip.label))) {
        stop('phy1 and phy2 have different sets of leaves.')
    }

    phy2_root = get_root_num(phy2)
    if (length(phy2_root) != 1L) {
        stop('phy2 must have exactly one root node.')
    }
    phy2_root_children = get_children_num(phy2, phy2_root)
    if (length(phy2_root_children) < 2L) {
        stop('phy2 root must have at least two children.')
    }

    if (length(phy2_root_children) > 2L) {
        phy2_tip_sets = .get_node_tip_sets(phy2)
        target_partition = .get_node_component_signatures(
            phy2,
            phy2_root,
            include_complement=FALSE,
            tip_sets=phy2_tip_sets
        )
        phy1_index = .build_phy_index(phy1, context='phy1')
        phy1_tip_sets = .get_node_tip_sets(phy1)
        candidate_nodes = unique(c(
            phy1_index[['root']],
            phy1[['edge']][,1]
        ))
        matched_node = NA_integer_
        for (candidate_node in candidate_nodes) {
            candidate_partition = .get_node_component_signatures(
                phy1,
                candidate_node,
                include_complement=(candidate_node != phy1_index[['root']]),
                tip_sets=phy1_tip_sets
            )
            if (identical(candidate_partition, target_partition)) {
                matched_node = as.integer(candidate_node)
                break
            }
        }
        if (is.na(matched_node)) {
            return(NA)
        }
        if (mode_name == 'node_num') {
            return(matched_node)
        }
        if (matched_node == phy1_index[['root']]) {
            return(NA_integer_)
        }
        return(which(phy1[['edge']][,2] == matched_node)[[1]])
    }

    reference_index = .build_phy_index(phy2, context='phy2')
    in_reference = logical(reference_index[['max_node']])
    in_reference[phy2_root_children[[1L]]] = TRUE
    for (node in reference_index[['preorder']]) {
        if (in_reference[[node]]) {
            in_reference[reference_index[['children']][[node]]] = TRUE
        }
    }
    reference_tips = phy2[['tip.label']][in_reference[seq_along(phy2[['tip.label']])]]
    index = .build_phy_index(phy1, context='phy1')
    total = integer(index[['max_node']])
    matches = integer(index[['max_node']])
    total[seq_len(index[['num_tip']])] = 1L
    matches[seq_len(index[['num_tip']])] = phy1[['tip.label']] %in% reference_tips
    for (node in index[['postorder']]) {
        children = index[['children']][[node]]
        if (length(children)) {
            total[[node]] = sum(total[children])
            matches[[node]] = sum(matches[children])
        }
    }
    k = length(reference_tips)
    child_nodes = phy1[['edge']][,2]
    is_match = (total[child_nodes] == k & matches[child_nodes] == k) |
        (index[['num_tip']] - total[child_nodes] == k & matches[child_nodes] == 0L)
    matched_indices = which(is_match)
    matched_index = if (length(matched_indices)) matched_indices[[1]] else NA_integer_

    if (is.na(matched_index)) {
        return(NA)
    }
    if (mode_name=="node_num") {
        root_pos = phy1$edge[matched_index,2]
    } else {
        root_pos = matched_index
    }
    return(root_pos)
}
