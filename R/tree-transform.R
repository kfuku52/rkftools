
# alias
#' Pad short external edges
#'
#' Convenience wrapper around [pad_short_edges()] restricted to terminal edges.
#'
#' @param tree A `phylo` tree with branch lengths.
#' @param threshold Finite non-negative minimum edge length.
#' @param verbose Whether to emit progress messages.
#' @return A modified `phylo` tree.
#' @export
collapse_short_external_edges = function(tree, threshold=1e-6, verbose=FALSE) {
    return(pad_short_edges(
        tree,
        threshold=threshold,
        external_only=TRUE,
        verbose=verbose
    ))
}


#' Pad short phylogenetic edges
#'
#' Transfers length through ancestors while preserving every root-to-tip
#' distance when possible. If the root must move, all finite root-to-tip
#' distances increase by the same minimum amount. Missing lengths remain
#' missing and impose no constraint across that edge. Unary, binary, and
#' multifurcating trees use the same algorithm.
#'
#' @param tree A `phylo` tree with branch lengths.
#' @param threshold Finite non-negative minimum edge length.
#' @param external_only Whether the minimum applies only to terminal edges.
#'   Internal edges can still transfer length, but never become negative.
#' @param verbose Whether to emit progress messages.
#' @return A modified `phylo` tree.
#' @export
pad_short_edges = function(tree, threshold=1e-6, external_only=FALSE, verbose=FALSE) {
    external_only = .normalize_single_logical_arg(external_only, 'external_only')
    verbose = .normalize_single_logical_arg(verbose, 'verbose')
    threshold = .normalize_finite_numeric_scalar(threshold, 'threshold', min_value=0)
    index = .validate_phylo_input(tree, context='tree', require_lengths=TRUE)
    lengths = tree[['edge.length']]
    if (any(!is.finite(lengths[!is.na(lengths)]))) {
        stop('tree must not contain infinite branch lengths in pad_short_edges().')
    }
    minimum = rep(threshold, length(lengths))
    if (external_only) {
        minimum[tree[['edge']][,2] > index[['num_tip']]] = 0
    }
    # Lower bounds on node shifts come from fixed tip positions. Negative
    # bounds allow a short internal edge to borrow from longer child edges,
    # avoiding an unnecessary change to the root height.
    required = rep(-Inf, index[['max_node']])
    required[seq_len(index[['num_tip']])] = 0
    for (node in index[['postorder']]) {
        children = index[['children']][[node]]
        if (!length(children)) next
        edges = index[['edge_index']][children]
        known = !is.na(lengths[edges])
        if (any(known)) {
            required[[node]] = max(minimum[edges[known]] +
                required[children[known]] - lengths[edges[known]])
        }
    }
    # Choose the smallest non-negative root extension, then keep each other
    # internal node as close to its original position as the bounds permit.
    shift = numeric(index[['max_node']])
    shift[[index[['root']]]] = max(0, required[[index[['root']]]])
    for (node in index[['preorder']][-1L]) {
        if (node <= index[['num_tip']]) next
        incoming = index[['edge_index']][[node]]
        upper = if (is.na(lengths[[incoming]])) Inf else
            shift[[index[['parent']][[node]]]] + lengths[[incoming]] - minimum[[incoming]]
        shift[[node]] = min(max(0, required[[node]]), upper)
    }
    delta = shift[tree[['edge']][,1]] - shift[tree[['edge']][,2]]
    tree[['edge.length']] = lengths + delta
    # Round-off at a transferred edge must not leave a negative/short length.
    known = !is.na(lengths)
    tree[['edge.length']][known] = pmax(tree[['edge.length']][known], minimum[known])
    if (verbose) {
        message('Padded ', sum(known & lengths < minimum), ' short edges.')
        if (anyNA(lengths)) message('NA edge lengths were left unchanged.')
        if (shift[[index[['root']]]] > 0) {
            message('Root-to-tip distances increased by ', shift[[index[['root']]]], '.')
        }
    }
    tree
}


#' Test whether a tree contains a polytomy
#'
#' @param phy A `phylo` tree.
#' @return A single logical value.
#' @export
contains_polytomy = function(phy) {
    .validate_phylo_input(phy, context='phy')
    child_counts = table(phy[['edge']][,1])
    if (ape::is.rooted(phy)) {
        return(any(child_counts > 2L))
    }
    root_num = get_root_num(phy)
    if (length(root_num) != 1L) {
        stop('phy must have exactly one root node in contains_polytomy().')
    }
    node_nums = as.integer(names(child_counts))
    underlying_degrees = as.integer(child_counts)
    underlying_degrees[node_nums != root_num] =
        underlying_degrees[node_nums != root_num] + 1L
    any(underlying_degrees > 3L)
}


#' Compare descendant leaves at two nodes
#'
#' @param phy1,phy2 `phylo` trees with identical tip-label multisets.
#' @param phy1_node,phy2_node Integer node numbers to compare.
#' @return A single logical value.
#' @export
has_same_leaves = function(phy1, phy1_node, phy2, phy2_node) {
    if (!identical(sort(phy1$tip.label), sort(phy2$tip.label))) {
        stop('phy1 and phy2 must contain identical tip-label multisets.')
    }
    phy1_leaves = sort(get_tip_labels(phy1, phy1_node))
    phy2_leaves = sort(get_tip_labels(phy2, phy2_node))
    is_same_leaves = identical(phy1_leaves, phy2_leaves)
    return(is_same_leaves)
}


#' Map polytomy nodes to a binary resolution
#'
#' @param multifurcated_tree A multifurcating `phylo` tree.
#' @param bifurcated_tree A binary resolution with identical ordered tips.
#' @param verbose Whether to emit progress messages.
#' @return A data frame mapping multifurcated to binary node numbers.
#' @export
multi2bi_node_number_transfer = function(multifurcated_tree, bifurcated_tree, verbose=FALSE) {
    verbose = .normalize_single_logical_arg(verbose, 'verbose')
    mtree = multifurcated_tree
    btree = bifurcated_tree
    if (!identical(mtree[['tip.label']], btree[['tip.label']])) {
        stop('multifurcated_tree and bifurcated_tree must have identical tip-label order.')
    }
    .validate_phylo_input(mtree, context='multifurcated_tree', unique_tips=TRUE)
    .validate_phylo_input(
        btree, context='bifurcated_tree', binary=TRUE, unique_tips=TRUE
    )
    mtree_signatures = .get_node_tip_signatures(mtree)
    btree_signatures = .get_node_tip_signatures(btree)
    mtree_internal_nodes = unique(mtree[['edge']][,1])
    missing_signatures = setdiff(
        unname(mtree_signatures[as.character(mtree_internal_nodes)]),
        unname(btree_signatures)
    )
    if (length(missing_signatures)) {
        stop(
            'bifurcated_tree must be a rooted binary refinement of ',
            'multifurcated_tree.'
        )
    }
    internal_node_counts = table(mtree[['edge']][,1])
    polytomy_parents = as.integer(names(internal_node_counts)[internal_node_counts > 2])
    if (!length(polytomy_parents)) {
        return(data.frame(mtree_node=integer(0), btree_node=integer(0)))
    }
    if (verbose) {
        message('Polytomy parent nodes: ', paste(polytomy_parents, collapse=', '))
    }
    df = data.frame(mtree_node=integer(0), btree_node=integer(0))
    for (mtree_pp in polytomy_parents) {
        matched_nodes = as.integer(names(btree_signatures)[
            btree_signatures == mtree_signatures[[as.character(mtree_pp)]]
        ])
        if (length(matched_nodes)) {
            df = rbind(df, data.frame(
                mtree_node=mtree_pp,
                btree_node=matched_nodes[[1]]
            ))
        } else {
            warning('Failed to map polytomy parent node to bifurcated tree: ', mtree_pp)
        }
    }
    return(df)
}


#' Collapse extremely short internal branches
#'
#' @param tree A `phylo` tree with branch lengths and unique tips.
#' @param tol Finite non-negative collapse tolerance.
#' @param verbose Whether to emit progress messages.
#' @return A modified `phylo` tree with unaffected edge lengths restored.
#' @export
collapse_short_branches = function(tree, tol=1e-8, verbose=FALSE) {
    verbose = .normalize_single_logical_arg(verbose, 'verbose')
    tol = .normalize_finite_numeric_scalar(tol, 'tol', min_value=0)
    .validate_phylo_input(
        tree, context='tree', require_lengths=TRUE, unique_tips=TRUE
    )
    if (is.null(tree[['edge.length']]) || length(tree[['edge.length']]) != nrow(tree[['edge']])) {
        stop('tree must contain one branch length per edge in collapse_short_branches().')
    }
    if (any(!is.finite(tree[['edge.length']][!is.na(tree[['edge.length']])]))) {
        stop('tree must not contain infinite branch lengths in collapse_short_branches().')
    }
    if (anyDuplicated(tree[['tip.label']])) {
        stop('tree must have unique tip labels in collapse_short_branches().')
    }
    out_tree = tree
    edge_lengths = as.numeric(out_tree$edge.length)
    num_tip = length(out_tree[['tip.label']])
    is_internal_edge = out_tree[['edge']][,2] > num_tip
    is_short_edge = is_internal_edge & !is.na(edge_lengths) & (abs(edge_lengths) < tol)
    if (verbose && any(is.na(edge_lengths))) {
        message('NA edge lengths were ignored when searching extremely short internal branches.')
    }
    if (!any(is_short_edge)) {
        if (verbose) {
            message('No extremely short internal branch was detected. tol = ', tol)
        }
        return(out_tree)
    }

    original_signatures = .get_node_tip_signatures(out_tree)
    original_length_by_signature = stats::setNames(
        edge_lengths[!is_short_edge],
        original_signatures[as.character(out_tree[['edge']][!is_short_edge,2])]
    )
    working_tree = out_tree
    protected_internal = is_internal_edge & !is_short_edge &
        !is.na(edge_lengths) & edge_lengths <= tol
    working_tree[['edge.length']][protected_internal] = max(1, tol * 2)
    working_tree[['edge.length']][is_short_edge] = 0
    if (verbose) {
        message(
            'Extremely short internal branches (n = ', sum(is_short_edge),
            ') were collapsed. tol = ', tol
        )
    }
    collapsed = ape::di2multi(working_tree, tol=tol)
    collapsed_signatures = .get_node_tip_signatures(collapsed)
    child_signatures = collapsed_signatures[as.character(collapsed[['edge']][,2])]
    restored_index = match(child_signatures, names(original_length_by_signature))
    if (anyNA(restored_index)) {
        stop('Failed to restore branch lengths after collapsing short branches.')
    }
    restored_lengths = unname(original_length_by_signature[restored_index])
    collapsed[['edge.length']] = as.numeric(restored_lengths)
    collapsed
}


#' Adjust a tree to be ultrametric
#'
#' Binary trees are adjusted with [ape::chronoMPL()]. For multifurcating trees,
#' terminal branches are extended to the maximum root-to-tip distance so the
#' original topology is preserved.
#'
#' @param tree A `phylo` tree with finite branch lengths.
#' @param stop_if_larger_change Maximum fractional total adjustment.
#' @param verbose Whether to emit progress messages.
#' @return An ultrametric `phylo` tree.
#' @export
force_ultrametric = function(tree, stop_if_larger_change=0.01, verbose=FALSE) {
    out_tree = tree
    verbose = .normalize_single_logical_arg(verbose, 'verbose')
    stop_if_larger_change = .normalize_finite_numeric_scalar(
        stop_if_larger_change, 'stop_if_larger_change', min_value=0
    )
    .validate_phylo_input(out_tree, context='tree', require_lengths=TRUE)
    if (is.null(out_tree[['edge.length']]) || length(out_tree[['edge.length']]) != nrow(out_tree[['edge']])) {
        stop('tree must contain one branch length per edge in force_ultrametric().')
    }
    if (any(is.na(out_tree[['edge.length']]))) {
        stop('tree contains NA edge lengths. Please resolve NA values before force_ultrametric().')
    }
    if (any(!is.finite(out_tree[['edge.length']]))) {
        stop('tree contains non-finite edge lengths. Please resolve them before force_ultrametric().')
    }
    if (ape::is.ultrametric(out_tree)) {
        if (verbose) {
            message('The tree is ultrametric.')
        }
    } else {
        if (verbose) {
            message('The tree is not ultrametric. Adjusting the branch length.')
        }
        edge_length_before = out_tree[['edge.length']]
        if (ape::is.binary(out_tree)) {
            out_tree = ape::chronoMPL(out_tree)
        } else {
            num_tip = length(out_tree[['tip.label']])
            tip_depths = ape::node.depth.edgelength(out_tree)[seq_len(num_tip)]
            target_depth = max(tip_depths)
            tip_edge_indices = match(
                seq_len(num_tip),
                out_tree[['edge']][,2]
            )
            out_tree[['edge.length']][tip_edge_indices] =
                out_tree[['edge.length']][tip_edge_indices] +
                (target_depth - tip_depths)
        }
        edge_length_after = out_tree[['edge.length']]
        sum_adjustment = sum(abs(edge_length_after-edge_length_before))
        if (verbose) {
            message('Total branch length difference: ', sum_adjustment)
        }
        allowed_adjustment = sum(abs(out_tree[['edge.length']])) * stop_if_larger_change
        if (!(sum_adjustment < allowed_adjustment)) {
            stop(
                'Ultrametric adjustment exceeded the allowed fraction: ',
                sum_adjustment, ' >= ', allowed_adjustment, '.'
            )
        }
    }
    return(out_tree)
}


#' Construct a one-tip tree
#'
#' @param name A single non-empty tip label.
#' @param dist A finite branch length.
#' @return A one-tip `phylo` tree.
#' @export
get_single_branch_tree = function(name, dist) {
    if (length(name) != 1 || is.na(name) || trimws(as.character(name)) == '') {
        stop('name must be a single non-empty tip label in get_single_branch_tree().')
    }
    dist_num = .normalize_finite_numeric_scalar(dist, 'dist')
    phy = list(
      edge = matrix(c(2,1),1,2),
      tip.label = as.character(name),
      edge.length = dist_num,
      Nnode = 1
    )
    class(phy) = "phylo"
    return(phy)
}


#' Collapse redundant singleton root edges
#'
#' @param phy A `phylo` tree.
#' @return A `phylo` tree with singleton nodes collapsed where applicable.
#' @export
remove_redundant_root_edge = function(phy) {
    if (is.null(phy[['tip.label']]) || length(phy[['tip.label']]) <= 1) {
        return(phy)
    }
    if (is.null(phy[['edge']]) || nrow(phy[['edge']]) <= 1) {
        return(phy)
    }
    out_phy = ape::collapse.singles(phy, root.edge=TRUE)
    return(out_phy)
}
