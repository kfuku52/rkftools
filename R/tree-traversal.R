# Phylogenetic tree traversal, transformation, rooting, and conversion tools.

# Naming convention for phylogenetic identifiers in this package:
# * *_num is an ape::phylo node number. These are the integer node indices used
#   in phy[['edge']], with tips numbered before internal nodes.
# * *_name is a biological or display label, usually from phy[['tip.label']] or
#   phy[['node.label']]. In branch tables this role is often held by a column
#   such as node_name, label, or another user-specified name_col.
# * *_id is a table or algorithm identifier, not necessarily an ape node number.
#   In branch tables, branch_id/parent/sister are *_id values that identify rows
#   and refer to each other. Other local *_id values, such as species_ids in
#   reconciliation helpers or branch_ids in MAD helpers, are compact
#   algorithm-specific indices.
# * Current *_id assignment sites are phylo2table(), which assigns branch-table
#   branch_id/parent/sister values from a bare ape::phylo; .tip_species_id_map(),
#   which assigns internal species_ids for species-overlap calculations; and
#   .compute_mad_scores(), which creates local branch_ids for iterating over edge
#   rows. Most other functions consume *_id values supplied by input tables or
#   use them only as local algorithm indices.
# * table2phylo() consumes branch-table *_id values and resolves
#   *_id -> *_name -> *_num while constructing an ape::phylo.
# * Because a phylo does not preserve original input branch_id values,
#   phylo2table() generates numerical *_id values from clade signatures, matching
#   genegalleon's numerical_label convention.

#' Look up phylogenetic node numbers by label
#'
#' @param phy A `phylo` tree.
#' @param node_name One or more tip or internal-node labels.
#' @return An integer vector of matching node numbers.
#' @export
get_node_num_by_name = function(phy, node_name) {
    node_names = c(phy[['tip.label']], phy$node.label)
    out = integer(0)
    for (name in as.character(node_name)) {
        if (is.na(name)) {
            next
        }
        out = c(out, which(node_names == name))
    }
    return(out)
}


#' Look up phylogenetic node labels by number
#'
#' @param phy A `phylo` tree.
#' @param node_num One or more integer node numbers. Missing and out-of-range
#'   values are ignored.
#' @return A character vector of node labels.
#' @export
get_node_name_by_num = function(phy, node_num) {
    node_names = c(phy[['tip.label']], phy$node.label)
    out = character(0)
    max_node_num = length(node_names)
    node_nums = .normalize_integerish(
        node_num, 'node_num', allow_na=TRUE, allow_empty=TRUE
    )
    for (num in node_nums) {
        if (is.na(num) || num < 1L || num > max_node_num) {
            next
        }
        out = c(out, node_names[num])
    }
    return(out)
}


#' Get the root node number
#'
#' @param phy A `phylo` tree.
#' @return The integer root-node number, or multiple values for malformed trees.
#' @examples
#' tree <- ape::read.tree(text="((A:1,B:1):1,C:1);")
#' get_root_num(tree)
#' @export
get_root_num = function(phy) {
    root_num = setdiff(phy[['edge']][,1], phy[['edge']][,2])
    return(root_num)
}


#' Get child node numbers
#'
#' @param phy A `phylo` tree.
#' @param node_num One or more integer node numbers.
#' @return An integer vector of child node numbers.
#' @export
get_children_num = function(phy, node_num) {
    out = integer(0)
    max_node_num = max(phy[['edge']])
    node_nums = .normalize_integerish(
        node_num, 'node_num', allow_na=TRUE, allow_empty=TRUE
    )
    for (nn in node_nums) {
        if (is.na(nn) || nn < 1L || nn > max_node_num) {
            next
        }
        out = c(out, phy[['edge']][(phy[['edge']][,1]==nn),2])
    }
    return(out)
}


#' Test whether nodes are the root
#'
#' @param phy A `phylo` tree.
#' @param node_num One or more integer node numbers.
#' @return A logical vector.
#' @export
is_root = function(phy, node_num) {
    node_num = .normalize_integerish(
        node_num, 'node_num', allow_na=TRUE, allow_empty=TRUE
    )
    root_num = get_root_num(phy)
    return(node_num==root_num)
}


#' Test whether a node is a leaf
#'
#' @param phy A `phylo` tree.
#' @param node_num A single integer node number.
#' @return A single logical value.
#' @export
is_leaf = function(phy, node_num) {
    node_num = .normalize_integerish(
        node_num, 'node_num', allow_na=TRUE, allow_empty=TRUE
    )
    max_node_num = max(phy[['edge']])
    if (length(node_num) != 1L || is.na(node_num) || node_num < 1L || node_num > max_node_num) {
        return(FALSE)
    }
    return(node_num <= length(phy[['tip.label']]))
}


#' Get descendant node numbers
#'
#' @param phy A `phylo` tree.
#' @param node_num One or more integer node numbers.
#' @param leaf_only Whether to return only tip nodes.
#' @return A sorted integer vector of descendants.
#' @export
get_descendent_num = function(phy, node_num, leaf_only=FALSE) {
    leaf_only = .normalize_single_logical_arg(
        value=leaf_only,
        arg_name='leaf_only'
    )
    index = .build_phy_index(phy)
    nodes = .normalize_integerish(node_num, 'node_num', allow_na=TRUE, allow_empty=TRUE)
    nodes = nodes[!is.na(nodes) & nodes >= 1L & nodes <= index[['max_node']]]
    .descendant_nodes(index, nodes, leaf_only)
}


#' Get parent node numbers
#'
#' @param phy A `phylo` tree.
#' @param node_num One or more integer node numbers.
#' @return An integer vector of parent node numbers.
#' @export
get_parent_num = function(phy, node_num) {
    out = integer(0)
    max_node_num = max(phy[['edge']])
    node_nums = .normalize_integerish(
        node_num, 'node_num', allow_na=TRUE, allow_empty=TRUE
    )
    for (nn in node_nums) {
        if (is.na(nn) || nn < 1L || nn > max_node_num) {
            next
        }
        out = c(out, phy[['edge']][(phy[['edge']][,2]==nn),1])
    }
    return(out)
}


#' Get sister node numbers
#'
#' @param phy A `phylo` tree.
#' @param node_num One or more integer node numbers.
#' @return An integer vector of sibling node numbers.
#' @export
get_sister_num = function(phy, node_num) {
    out = integer(0)
    max_node_num = max(phy[['edge']])
    node_nums = .normalize_integerish(
        node_num, 'node_num', allow_na=TRUE, allow_empty=TRUE
    )
    for (nn in node_nums) {
        if (is.na(nn) || nn < 1L || nn > max_node_num) {
            next
        }
        parent_num = phy[['edge']][(phy[['edge']][,2]==nn),1]
        sibling_num = phy[['edge']][(phy[['edge']][,1]==parent_num),2]
        sister_num = sibling_num[sibling_num!=nn]
        out = c(out, sister_num)
    }
    return(out)
}


#' Get ancestor node numbers
#'
#' @param phy A `phylo` tree.
#' @param node_num A single integer node number.
#' @return Ancestor node numbers ordered from parent to root.
#' @export
get_ancestor_num = function(phy, node_num) {
    node_num = .normalize_integerish(
        node_num, 'node_num', allow_na=TRUE, allow_empty=TRUE
    )
    index = .build_phy_index(phy)
    if (length(node_num) != 1L || is.na(node_num) ||
            node_num < 1L || node_num > index[['max_node']]) return(integer(0))
    ancestors = integer(index[['max_node']])
    count = 0L
    parent = index[['parent']][[node_num]]
    while (!is.na(parent)) {
        count = count + 1L
        ancestors[[count]] = parent
        parent = index[['parent']][[parent]]
    }
    ancestors[seq_len(count)]
}


#' Get tip labels below nodes
#'
#' @param phy A `phylo` tree.
#' @param node_num One or more integer node numbers.
#' @param out Optional character values to prepend for compatibility.
#' @return A character vector of descendant tip labels.
#' @export
get_tip_labels = function(phy, node_num, out=NULL) {
    index = .build_phy_index(phy)
    num_leaf = index[['num_tip']]
    node_nums = .normalize_integerish(
        node_num, 'node_num', allow_na=TRUE, allow_empty=TRUE
    )
    max_node_num = max(phy[['edge']])
    invalid_node_nums = node_nums[!is.na(node_nums) & (node_nums > max_node_num)]
    if (length(invalid_node_nums)) {
        stop(
            'node_num contains value(s) outside valid node range [1, ',
            max_node_num,
            ']: ',
            paste(unique(invalid_node_nums), collapse=', ')
        )
    }
    tip_labels = character(0)
    for (nn in node_nums) {
        if (is.na(nn) || nn < 1) {
            next
        }
        if (nn > num_leaf) {
            tips = .descendant_nodes(index, nn, leaf_only=TRUE)
            tip_labels = c(tip_labels, phy[['tip.label']][tips])
        } else {
            tip_labels = c(tip_labels, phy[['tip.label']][nn])
        }
    }
    if (!is.null(out)) {
        tip_labels = c(out, tip_labels)
    }
    return(tip_labels)
}


#' Find nearest subject tips
#'
#' @param phy A `phylo` tree.
#' @param query A single tip label.
#' @param subjects Candidate tip labels.
#' @param mrca_matrix A named matrix returned by `ape::mrca()`.
#' @return A list containing nearest tips and their MRCA.
#' @export
get_nearest_tips = function(phy, query, subjects, mrca_matrix) {
    query = .normalize_single_string_arg(
        value=query,
        arg_name='query',
        allow_empty=FALSE
    )
    subjects = as.character(subjects)
    if (length(subjects) == 0 || any(is.na(subjects) | trimws(subjects) == '')) {
        stop('subjects must contain at least one non-missing tip label.')
    }
    if (!(query %in% rownames(mrca_matrix))) {
        stop('query must be present in mrca_matrix row names: ', query)
    }
    missing_subjects = setdiff(subjects, colnames(mrca_matrix))
    if (length(missing_subjects)) {
        stop(
            'subjects are missing in mrca_matrix column names: ',
            paste(missing_subjects, collapse=', ')
        )
    }
    query_num = get_node_num_by_name(phy, query)
    if (length(query_num) != 1) {
        stop('query must map to exactly one node in phy: ', query)
    }
    mrcas = mrca_matrix[query,subjects]
    if (length(subjects)==1) {
        names(mrcas) = subjects
    }
    uniq_mrcas = unique(mrcas)
    path_lens = rep(NA, length(uniq_mrcas))
    names(path_lens) = uniq_mrcas
    for (i in seq_along(uniq_mrcas)) {
        path_lens[i] = length(ape::nodepath(phy, from=query_num, to=uniq_mrcas[i]))
    }
    nearest_mrca = names(path_lens)[path_lens==min(path_lens)]
    nearest_tips = names(mrcas)[mrcas==nearest_mrca]
    return(list(nearests=nearest_tips, mrca=nearest_mrca))
}


#' Get the age of a node in an ultrametric tree
#'
#' @param phy An ultrametric `phylo` tree.
#' @param node_num A single integer node number.
#' @return Numeric time from the node to its descendant tips.
#' @export
get_node_age = function(phy, node_num) {
    if (!ape::is.ultrametric(phy)) {
        stop('phy must be ultrametric in get_node_age().')
    }
    node_num = .normalize_integerish(node_num, 'node_num')
    max_node_num = max(phy[['edge']])
    if (length(node_num) != 1 || is.na(node_num) || node_num < 1 || node_num > max_node_num) {
        stop('node_num must be a single integer in [1, ', max_node_num, '] in get_node_age().')
    }
    age = 0
    current_node_num = node_num
    while (!is.na(current_node_num)) {
        descendent_node_num = get_children_num(phy, current_node_num)[1]
        if (!is.na(descendent_node_num)) {
            edge_length = phy[['edge.length']][(phy[['edge']][,1]==current_node_num)&(phy[['edge']][,2]==descendent_node_num)]
            age = age + edge_length
        }
        current_node_num = descendent_node_num
    }
    return(age)
}
