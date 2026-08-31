# Validate topology once, then use integer-indexed adjacency and traversal data.
.build_phy_index = function(phy, context='phy') {
    edge = phy[['edge']]
    if (!inherits(phy, 'phylo') || !is.matrix(edge) ||
            ncol(edge) != 2L || !nrow(edge)) {
        stop(context, ' must be a non-empty object of class "phylo".')
    }
    tips = phy[['tip.label']]
    if (!is.character(tips) || !length(tips) || anyNA(tips) ||
            any(trimws(tips) == '')) {
        stop(context, ' must contain non-missing, non-empty tip labels.')
    }
    num_tip = length(tips)
    num_internal = phy[['Nnode']]
    if (!is.numeric(num_internal) || length(num_internal) != 1L ||
            is.na(num_internal) || !is.finite(num_internal) ||
            num_internal < 1 || num_internal != floor(num_internal)) {
        stop(context, ' must have a positive integer Nnode.')
    }
    max_node = num_tip + num_internal
    if (!is.numeric(edge) || anyNA(edge) || any(!is.finite(edge)) ||
            any(edge != floor(edge)) || any(edge < 1 | edge > max_node)) {
        stop(context, ' must use integer node numbers from 1 to Ntip + Nnode.')
    }
    if (nrow(edge) != max_node - 1L ||
            any(tabulate(edge, nbins=max_node) == 0L)) {
        stop(context, ' has edge/node counts inconsistent with Ntip and Nnode.')
    }
    if (any(edge[,1] <= num_tip)) {
        stop(context, ' has a tip node used as a parent.')
    }
    if (anyDuplicated(edge[,2])) {
        stop(context, ' contains a multiply-parented node or duplicate edge.')
    }
    root = setdiff(seq_len(max_node), edge[,2])
    if (length(root) != 1L || root != num_tip + 1L) {
        stop(context, ' must have exactly one root numbered Ntip + 1.')
    }
    internal = seq.int(num_tip + 1L, max_node)
    if (any(tabulate(edge[,1], nbins=max_node)[internal] == 0L)) {
        stop(context, ' contains an internal node without children.')
    }
    children = rep(list(integer(0)), max_node)
    grouped = split(as.integer(edge[,2]), edge[,1])
    children[as.integer(names(grouped))] = unname(grouped)
    parent = rep(NA_integer_, max_node)
    parent[edge[,2]] = as.integer(edge[,1])
    edge_index = rep(NA_integer_, max_node)
    edge_index[edge[,2]] = seq_len(nrow(edge))
    # A breadth-first order also places every parent before its children.
    preorder = integer(max_node)
    preorder[[1L]] = root
    size = 1L
    next_index = 1L
    while (next_index <= size) {
        child_nodes = children[[preorder[[next_index]]]]
        if (length(child_nodes)) {
            preorder[size + seq_along(child_nodes)] = child_nodes
            size = size + length(child_nodes)
        }
        next_index = next_index + 1L
    }
    if (size != max_node) {
        stop(context, ' contains a disconnected or cyclic topology.')
    }
    list(root=as.integer(root), children=children, parent=parent,
        edge_index=edge_index, preorder=preorder, postorder=rev(preorder),
        num_tip=num_tip, max_node=as.integer(max_node))
}

.validate_phylo_input = function(
    phy, context='phy', rooted=NULL, binary=NULL, require_lengths=FALSE,
    finite_lengths=FALSE, unique_tips=FALSE
) {
    index = .build_phy_index(phy, context)
    if (unique_tips && anyDuplicated(phy[['tip.label']])) {
        stop(context, ' must contain unique tip labels.')
    }
    if (!is.null(rooted) && !identical(isTRUE(ape::is.rooted(phy)), rooted)) {
        stop(context, if (rooted) ' must be rooted.' else ' must be unrooted.')
    }
    if (!is.null(binary) && !identical(isTRUE(ape::is.binary(phy)), binary)) {
        stop(context, if (binary) ' must be binary.' else ' must be non-binary.')
    }
    if (require_lengths || finite_lengths) {
        lengths = phy[['edge.length']]
        if (!is.numeric(lengths) || length(lengths) != nrow(phy[['edge']])) {
            stop(context, ' must contain one numeric branch length per edge.')
        }
        if (finite_lengths && any(!is.finite(lengths))) {
            stop(context, ' must contain only finite, non-missing branch lengths.')
        }
    }
    invisible(index)
}

.node_tip_cache = function(index, tip_values=seq_len(index[['num_tip']])) {
    cache = vector('list', index[['max_node']])
    cache[seq_len(index[['num_tip']])] = as.list(tip_values)
    for (node in index[['postorder']]) {
        if (node > index[['num_tip']]) {
            cache[[node]] = unlist(cache[index[['children']][[node]]], use.names=FALSE)
        }
    }
    cache
}

.get_node_tip_sets = function(phy, index=.build_phy_index(phy)) {
    cache = .node_tip_cache(index, as.character(phy[['tip.label']]))
    cache = lapply(cache, function(tips) sort(unique(tips)))
    names(cache) = as.character(seq_len(index[['max_node']]))
    cache
}

.get_node_tip_signatures = function(phy) {
    vapply(.get_node_tip_sets(phy), .encode_tip_signature, character(1))
}

.descendant_nodes = function(index, nodes, leaf_only=FALSE) {
    active = logical(index[['max_node']])
    descendants = active
    active[nodes] = TRUE
    for (node in index[['preorder']]) {
        if (active[[node]]) {
            children = index[['children']][[node]]
            active[children] = TRUE
            descendants[children] = TRUE
        }
    }
    if (leaf_only) descendants = descendants[seq_len(index[['num_tip']])]
    which(descendants)
}

.encode_tip_signature = function(tip_labels) {
    tip_labels = sort(as.character(tip_labels))
    paste0(nchar(tip_labels), ':', tip_labels, collapse='|')
}


.get_node_component_signatures = function(
    phy,
    node_num,
    include_complement=FALSE,
    tip_sets=NULL
) {
    if (is.null(tip_sets)) {
        tip_sets = .get_node_tip_sets(phy)
    }
    children = get_children_num(phy, node_num)
    component_sets = tip_sets[as.character(children)]
    if (include_complement) {
        component_sets = c(
            component_sets,
            list(setdiff(phy[['tip.label']], tip_sets[[as.character(node_num)]]))
        )
    }
    sort(unname(vapply(component_sets, .encode_tip_signature, character(1))))
}


.get_node_depths_from_root = function(tree) {
    index = .build_phy_index(tree, context='tree')
    depths = rep(NA_integer_, index[['max_node']])
    depths[index[['root']]] = 0L
    for (node_num in index[['preorder']]) {
        children = index[['children']][[node_num]]
        if (length(children)) {
            depths[children] = depths[[node_num]] + 1L
        }
    }
    depths
}
