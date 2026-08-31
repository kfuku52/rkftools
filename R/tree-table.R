
# See the file-level naming convention above for *_num, *_id, and *_name.
.table2phylo_is_parent_sentinel = function(values) {
    if (is.factor(values)) {
        values = as.character(values)
    }
    is_sentinel = is.na(values)
    values_chr = trimws(as.character(values))
    is_sentinel = is_sentinel | (values_chr == '')
    suppressWarnings(values_num <- as.numeric(values_chr))
    is_sentinel = is_sentinel | ((!is.na(values_num)) & (values_num == -999))
    return(is_sentinel)
}


.table2phylo_normalize_dist_column = function(values) {
    if (is.factor(values)) {
        values = as.character(values)
    }
    if (is.character(values)) {
        values = trimws(values)
        values[values == ''] = NA_character_
    }
    numeric_values = suppressWarnings(as.numeric(values))
    invalid = is.na(numeric_values) & !is.na(values)
    if (any(invalid)) {
        stop(
            'dist_col contains non-numeric value(s) in table2phylo(): ',
            paste(unique(as.character(values[invalid])), collapse=', ')
        )
    }
    if (any(!is.finite(numeric_values[!is.na(numeric_values)]))) {
        stop('dist_col must contain only finite numeric branch lengths in table2phylo().')
    }
    numeric_values
}


.table2phylo_detect_root_id = function(df, branch_col='branch_id', parent_col='parent') {
    branch_ids = df[[branch_col]]
    is_sentinel_parent = .table2phylo_is_parent_sentinel(df[[parent_col]])
    sentinel_rows = which(is_sentinel_parent)
    if (length(sentinel_rows) == 1) {
        return(branch_ids[[sentinel_rows]])
    }
    if (length(sentinel_rows) > 1) {
        sentinel_ids = unique(as.character(branch_ids[sentinel_rows]))
        stop(
            'Ambiguous root candidate(s) in table2phylo(): multiple sentinel-parent rows found for branch_id: ',
            paste(sentinel_ids, collapse=', ')
        )
    }

    stop('Unable to infer root in table2phylo(): exactly one sentinel-parent row is required.')
}


.table2phylo_validate_graph = function(df, root_id, name_col, dist_col) {
    branch_ids = as.character(df[['branch_id']])
    parent_ids = as.character(df[['parent']])
    sister_ids = as.character(df[['sister']])
    node_names = as.character(df[[name_col]])

    if (anyDuplicated(branch_ids)) {
        stop('Duplicate branch_id values detected in table2phylo().')
    }
    missing_names = is.na(node_names) | trimws(node_names) == ''
    if (any(missing_names)) {
        stop('name_col contains missing/blank value(s) in table2phylo().')
    }
    root_index = match(as.character(root_id), branch_ids)
    nonroot_index = setdiff(seq_along(branch_ids), root_index)
    if (!.table2phylo_is_parent_sentinel(parent_ids[[root_index]])) {
        stop('The root row must have a sentinel parent in table2phylo().')
    }
    if (!.table2phylo_is_parent_sentinel(sister_ids[[root_index]])) {
        stop('The root row must have a sentinel sister in table2phylo().')
    }

    if (length(nonroot_index)) {
        nonroot_parents = parent_ids[nonroot_index]
        sentinel_parent = .table2phylo_is_parent_sentinel(nonroot_parents)
        if (any(sentinel_parent)) {
            stop(
                'Only the root row may have a sentinel parent in table2phylo(). branch_id: ',
                paste(branch_ids[nonroot_index][sentinel_parent], collapse=', ')
            )
        }
        unknown_parents = setdiff(nonroot_parents, branch_ids)
        if (length(unknown_parents)) {
            stop(
                'Unknown parent branch_id reference(s) in table2phylo(): ',
                paste(unique(unknown_parents), collapse=', ')
            )
        }
        self_parent = branch_ids[nonroot_index] == nonroot_parents
        if (any(self_parent)) {
            stop(
                'Self-referential parent branch_id(s) in table2phylo(): ',
                paste(branch_ids[nonroot_index][self_parent], collapse=', ')
            )
        }
    }

    children_by_parent = split(branch_ids[nonroot_index], parent_ids[nonroot_index])
    reachable = as.character(root_id)
    frontier = reachable
    while (length(frontier)) {
        next_frontier = unique(unlist(children_by_parent[frontier], use.names=FALSE))
        next_frontier = setdiff(next_frontier, reachable)
        reachable = c(reachable, next_frontier)
        frontier = next_frontier
    }
    unreachable = setdiff(branch_ids, reachable)
    if (length(unreachable)) {
        stop(
            'Disconnected or cyclic branch_id row(s) in table2phylo(): ',
            paste(unreachable, collapse=', ')
        )
    }

    tip_ids = setdiff(branch_ids, names(children_by_parent))
    tip_names = node_names[match(tip_ids, branch_ids)]
    if (anyDuplicated(tip_names)) {
        stop(
            'table2phylo() requires unique tip labels. Duplicated label(s): ',
            paste(unique(tip_names[duplicated(tip_names)]), collapse=', ')
        )
    }

    for (parent_id in names(children_by_parent)) {
        children = children_by_parent[[parent_id]]
        if (length(children) > 2L) {
            child_indices = match(children, branch_ids)
            if (any(!.table2phylo_is_parent_sentinel(sister_ids[child_indices]))) {
                stop(
                    'Children of a multifurcating parent must use sentinel sister values ',
                    'in table2phylo(): parent branch_id ', parent_id, '.'
                )
            }
        } else if (length(children) == 1L) {
            child_index = match(children, branch_ids)
            if (!.table2phylo_is_parent_sentinel(sister_ids[[child_index]])) {
                stop(
                    'A child without a sister must use a sentinel sister in table2phylo(): ',
                    children
                )
            }
        } else if (length(children) == 2L) {
            first_index = match(children[[1]], branch_ids)
            second_index = match(children[[2]], branch_ids)
            if (!identical(sister_ids[[first_index]], children[[2]]) ||
                    !identical(sister_ids[[second_index]], children[[1]])) {
                stop(
                    'Non-reciprocal sister references under parent branch_id ',
                    parent_id, ' in table2phylo().'
                )
            }
        }
    }

    if (any(is.na(df[[dist_col]][nonroot_index]))) {
        stop('Non-root branches must have finite branch lengths in table2phylo().')
    }
    invisible(children_by_parent)
}


#' Convert a branch table to a phylogenetic tree
#'
#' The table must describe one connected rooted tree. `branch_id`, `parent`,
#' and `sister` are identifiers; `-999`, blank, or missing values denote root
#' and missing-sister sentinels. `parent` is the canonical topology field.
#' Binary children use reciprocal `sister` identifiers, while unary and
#' multifurcating children use sister sentinels.
#'
#' @param df A non-empty branch data frame.
#' @param name_col Name of the node-label column.
#' @param dist_col Name of the finite branch-length column.
#' @return A rooted `phylo` tree. Root-row distance is stored as `root.edge`.
#' @examples
#' tab <- data.frame(
#'     branch_id=c(2, 0, 1), parent=c(-999, 2, 2), sister=c(-999, 1, 0),
#'     label=c("Root", "A", "B"), dist=c(0, 0.5, 0.75)
#' )
#' table2phylo(tab, "label", "dist")
#' @export
table2phylo = function(df, name_col, dist_col) {
    name_col = .normalize_single_string_arg(name_col, 'name_col', allow_empty=FALSE)
    dist_col = .normalize_single_string_arg(dist_col, 'dist_col', allow_empty=FALSE)
    if (!is.data.frame(df) || nrow(df) == 0L) {
        stop('df must be a non-empty data.frame in table2phylo().')
    }
    df_local = df
    .validate_branch_columns(name_col, dist_col)
    required_cols = unique(c('branch_id', 'parent', 'sister', name_col, dist_col))
    missing_cols = required_cols[!(required_cols %in% colnames(df_local))]
    if (length(missing_cols) > 0) {
        stop('Missing required columns in table2phylo(): ', paste(missing_cols, collapse=', '))
    }
    id_cols = c('branch_id', 'parent', 'sister')
    for (id_col in id_cols) {
        df_local[[id_col]] = as.character(df_local[[id_col]])
    }
    branch_ids_chr = as.character(df_local[['branch_id']])
    missing_branch_ids = is.na(branch_ids_chr) | (trimws(branch_ids_chr) == '')
    if (any(missing_branch_ids)) {
        stop('branch_id contains missing/blank value(s) in table2phylo().')
    }

    df_local[[dist_col]] = .table2phylo_normalize_dist_column(df_local[[dist_col]])
    root_id = .table2phylo_detect_root_id(df_local, branch_col='branch_id', parent_col='parent')
    if (is.na(root_id) || trimws(as.character(root_id)) == '') {
        stop('Failed to infer a valid root branch_id in table2phylo().')
    }
    is_root_row = (as.character(df_local[,'branch_id']) == as.character(root_id))
    is_root_row[is.na(is_root_row)] = FALSE
    if (sum(is_root_row) != 1) {
        stop('Failed to identify a unique root row in table2phylo().')
    }
    .table2phylo_validate_graph(
        df=df_local,
        root_id=root_id,
        name_col=name_col,
        dist_col=dist_col
    )

    branch_ids = as.character(df_local[['branch_id']])
    parent_ids = as.character(df_local[['parent']])
    node_names = as.character(df_local[[name_col]])
    root_index = match(as.character(root_id), branch_ids)
    root_name = node_names[[root_index]]
    root_dist = df_local[[dist_col]][[root_index]]
    if (is.na(root_name) || trimws(as.character(root_name)) == '') {
        stop('Failed to resolve root metadata in table2phylo().')
    }
    if (is.na(root_dist)) {
        root_dist = 0
    }
    if (nrow(df_local) == 1L) {
        return(get_single_branch_tree(root_name, root_dist))
    }

    nonroot_index = setdiff(seq_len(nrow(df_local)), root_index)
    parent_node_ids = unique(parent_ids[nonroot_index])
    tip_ids = branch_ids[!(branch_ids %in% parent_node_ids)]
    internal_ids = c(as.character(root_id), setdiff(parent_node_ids, as.character(root_id)))
    num_tip = length(tip_ids)
    num_internal = length(internal_ids)
    node_num_by_id = c(
        stats::setNames(seq_len(num_tip), tip_ids),
        stats::setNames(num_tip + seq_len(num_internal), internal_ids)
    )

    edge = cbind(
        as.integer(node_num_by_id[parent_ids[nonroot_index]]),
        as.integer(node_num_by_id[branch_ids[nonroot_index]])
    )
    edge_lengths = as.numeric(df_local[[dist_col]][nonroot_index])
    names_by_id = stats::setNames(node_names, branch_ids)
    phy = list(
        edge=edge,
        tip.label=unname(names_by_id[tip_ids]),
        edge.length=edge_lengths,
        Nnode=as.integer(num_internal),
        node.label=unname(names_by_id[internal_ids])
    )
    phy[['root.edge']] = as.numeric(root_dist)
    class(phy) = 'phylo'
    phy = ape::reorder.phylo(phy, order='cladewise')
    if (length(phy[['tip.label']]) > 1L) {
        phy = ape::ladderize(phy, right=TRUE)
    }
    return(phy)
}


.phylo2table_make_branch_id_map = function(phy, node_nums, index=.build_phy_index(phy)) {
    num_tip = index[['num_tip']]
    tip_rank = integer(num_tip)
    tip_rank[order(as.character(phy[['tip.label']]))] = seq_len(num_tip)
    tip_cache = .node_tip_cache(index, tip_rank)
    # Preserve the historical numerical_label ordering, including unary ties.
    signatures = vapply(tip_cache[node_nums], function(ranks) {
        bits = rep('0', num_tip)
        bits[num_tip - ranks + 1L] = '1'
        paste0(bits, collapse='')
    }, character(1))
    branch_ids = integer(length(node_nums))
    branch_ids[order(signatures, node_nums)] = seq_along(node_nums) - 1L
    branch_ids
}


#' Convert a phylogenetic tree to a branch table
#'
#' Supports rooted trees with unary, binary, or multifurcating internal nodes
#' and assigns stable branch identifiers from clade signatures. Because the
#' legacy `sister` column is singular, children with multiple siblings receive
#' the sentinel `-999`; their shared `parent` fully represents the topology.
#'
#' @param phy A rooted `phylo` tree with finite branch lengths and unique tips.
#' @param name_col Output node-label column name.
#' @param dist_col Output branch-length column name.
#' @return A branch data frame with `branch_id`, `parent`, `sister`, labels,
#'   and distances.
#' @examples
#' tree <- ape::read.tree(text="((A:1,B:1):1,C:1);")
#' phylo2table(tree)
#' @export
phylo2table = function(phy, name_col='label', dist_col='dist') {
    index = .validate_phylo_input(phy, context='phy', rooted=TRUE,
        finite_lengths=TRUE, unique_tips=TRUE)
    name_col = .normalize_single_string_arg(name_col, 'name_col', allow_empty=FALSE)
    dist_col = .normalize_single_string_arg(dist_col, 'dist_col', allow_empty=FALSE)
    .validate_branch_columns(name_col, dist_col)
    root_edge = phy[['root.edge']]
    if (is.null(root_edge)) {
        root_edge = 0
    } else {
        root_edge = .normalize_finite_numeric_scalar(root_edge, 'phy$root.edge')
    }
    phy = fill_node_labels(phy)
    node_names = c(phy[['tip.label']], phy[['node.label']])
    branch_ids = .phylo2table_make_branch_id_map(
        phy, seq_len(index[['max_node']]), index=index)
    sisters = rep(-999L, index[['max_node']])
    for (children in index[['children']]) {
        if (length(children) == 2L) sisters[children] = branch_ids[rev(children)]
    }
    ordered_nodes = c(index[['root']], as.integer(phy[['edge']][,2]))
    out = data.frame(
        branch_id=branch_ids[ordered_nodes],
        parent=c(-999L, branch_ids[phy[['edge']][,1]]),
        sister=sisters[ordered_nodes],
        label=node_names[ordered_nodes],
        dist=c(root_edge, phy[['edge.length']]),
        stringsAsFactors=FALSE
    )
    names(out) = c('branch_id', 'parent', 'sister', name_col, dist_col)
    rownames(out) = NULL
    out
}


.validate_branch_columns = function(name_col, dist_col) {
    reserved = c('branch_id', 'parent', 'sister')
    if (name_col %in% reserved) {
        stop('name_col must not be one of: ', paste(reserved, collapse=', '))
    }
    if (dist_col %in% c(reserved, name_col)) {
        stop('dist_col must not duplicate branch_id, parent, sister, or name_col.')
    }
}
