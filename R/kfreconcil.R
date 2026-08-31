# Gene-tree reconciliation and NOTUNG parsing utilities.

#' Calculate a duplication confidence score
#'
#' @param phy A `phylo` tree with unique tip labels.
#' @param node_num A single integer internal-node number.
#' @param species_parser Species-label convention: `"legacy"` or
#'   `"taxonomic"`.
#' @param sep Literal separator in gene labels.
#' @return The maximum pairwise Jaccard overlap among child clades, or `NA`
#'   when the node has fewer than two children or is the trivalent display root
#'   of an unrooted tree. For binary nodes this is the original two-child score.
#' @export
get_duplication_confidence_score = function(phy, node_num, species_parser='legacy', sep='_') {
    .validate_phylo_input(phy, context='phy', unique_tips=TRUE)
    node_num = .normalize_integerish(
        node_num,
        'node_num',
        min_value=1L,
        max_value=max(phy[['edge']])
    )
    if (length(node_num) != 1L) {
        stop('node_num must be a single integer node number.')
    }
    species_parser = .normalize_species_parser_arg(
        value=species_parser,
        arg_name='species_parser'
    )
    sep = .normalize_single_string_arg(
        value=sep,
        arg_name='sep',
        allow_empty=FALSE
    )
    children_num = get_children_num(phy, node_num)
    root_num = get_root_num(phy)
    if (!ape::is.rooted(phy) && length(root_num) == 1L &&
            node_num == root_num && length(children_num) == 3L) {
        return(NA)
    }
    if (length(children_num)>=2) {
        child_leaves = vector(mode='list', length(children_num))
        for (j in seq_along(children_num)) {
            child_leaves[[j]] = get_tip_labels(phy, children_num[j])
            parsed_species = .parse_species_labels(
                labels=child_leaves[[j]],
                species_parser=species_parser,
                sep=sep,
                output_sep=' ',
                require_gene=FALSE,
                fallback_label=FALSE
            )
            if (any(!parsed_species[['parsed_ok']])) {
                stop(
                    'Unable to parse species from tip label(s): ',
                    paste(child_leaves[[j]][!parsed_species[['parsed_ok']]], collapse=', ')
                )
            }
            child_leaves[[j]] = parsed_species[['species_labels']]
        }
        dc_score = .max_pairwise_species_overlap(child_leaves)
    } else {
        dc_score = NA
    }
    return(dc_score)
}

.species_jaccard = function(species1, species2) {
    species1 = unique(species1)
    species2 = unique(species2)
    if (!length(species1) && !length(species2)) {
        return(0)
    }
    intersect_count = if (length(species1) <= length(species2)) {
        sum(species1 %in% species2)
    } else {
        sum(species2 %in% species1)
    }
    union_count = length(species1) + length(species2) - intersect_count
    if (union_count == 0L) {
        return(0)
    }
    intersect_count / union_count
}

.max_pairwise_species_overlap = function(child_species) {
    num_children = length(child_species)
    if (num_children < 2L) {
        return(NA_real_)
    }
    max_overlap = 0
    for (first_index in seq_len(num_children - 1L)) {
        for (second_index in seq.int(first_index + 1L, num_children)) {
            max_overlap = max(
                max_overlap,
                .species_jaccard(
                    child_species[[first_index]],
                    child_species[[second_index]]
                )
            )
        }
    }
    max_overlap
}

# Assigns compact internal species_id values for species-overlap algorithms.
# These IDs are algorithm identifiers, not ape::phylo node numbers.
.tip_species_id_map = function(tip_labels, species_parser='legacy', sep='_') {
    species_parser = .normalize_species_parser_arg(
        value=species_parser,
        arg_name='species_parser'
    )
    sep = .normalize_single_string_arg(
        value=sep,
        arg_name='sep',
        allow_empty=FALSE
    )
    parsed_species = .parse_species_labels(
        labels=tip_labels,
        species_parser=species_parser,
        sep=sep,
        output_sep=sep,
        require_gene=FALSE,
        fallback_label=FALSE
    )
    if (any(!parsed_species[['parsed_ok']])) {
        stop(
            'Unable to parse species from tip label(s): ',
            paste(tip_labels[!parsed_species[['parsed_ok']]], collapse=', '),
            '. Use labels compatible with species_parser or choose the correct sep.'
        )
    }
    species_labels = parsed_species[['species_labels']]
    species_levels = unique(species_labels)
    species_ids = match(species_labels, species_levels)
    names(species_ids) = tip_labels
    species_ids
}

.species_overlap_score_fast_impl = function(
    phy,
    dc_cutoff=0,
    species_id_by_label=NULL,
    species_parser='legacy',
    sep='_'
) {
    index = .build_phy_index(phy)
    tip_count = index[['num_tip']]
    if (is.null(species_id_by_label)) {
        tip_species_ids = unname(.tip_species_id_map(
            phy[['tip.label']], species_parser=species_parser, sep=sep
        ))
    } else {
        tip_species_ids = as.integer(species_id_by_label[phy[['tip.label']]])
        if (anyNA(tip_species_ids)) {
            stop('species_id_by_label must contain every tip label.')
        }
    }
    species_cache = vector('list', index[['max_node']])
    species_cache[seq_len(tip_count)] = as.list(tip_species_ids)
    skip_root = !ape::is.rooted(phy) &&
        length(index[['children']][[index[['root']]]]) == 3L
    overlap_count = 0L
    for (node in index[['postorder']]) {
        if (node <= tip_count) next
        children = index[['children']][[node]]
        child_species = species_cache[children]
        species_cache[[node]] = unique(unlist(child_species, use.names=FALSE))
        if (length(children) < 2L || (skip_root && node == index[['root']])) next
        dc_score = .max_pairwise_species_overlap(child_species)
        if (!is.na(dc_score) && dc_score > dc_cutoff) {
            overlap_count = overlap_count + 1L
        }
    }
    as.numeric(overlap_count)
}

#' Score species-overlap duplications in a gene tree
#'
#' @param phy A `phylo` tree with unique tip labels. Multifurcating nodes use
#'   their maximum pairwise child-clade overlap.
#' @param dc_cutoff Finite duplication-confidence cutoff in `[0, 1]`.
#' @param species_parser Species-label convention.
#' @param sep Literal separator in gene labels.
#' @return The number of internal nodes whose score exceeds `dc_cutoff`.
#' @export
get_species_overlap_score = function(phy, dc_cutoff=0, species_parser='legacy', sep='_') {
    # this function assumes that leaf names are: GENUS_SPECIES_GENEID (e.g. Bos_taurus_AF492351.1)
    species_parser = .normalize_species_parser_arg(
        value=species_parser,
        arg_name='species_parser'
    )
    sep = .normalize_single_string_arg(
        value=sep,
        arg_name='sep',
        allow_empty=FALSE
    )
    .validate_phylo_input(
        phy,
        context='phy',
        unique_tips=TRUE
    )
    dc_cutoff = .normalize_finite_numeric_scalar(
        dc_cutoff,
        'dc_cutoff',
        min_value=0,
        max_value=1
    )
    .species_overlap_score_fast_impl(
        phy=phy,
        dc_cutoff=dc_cutoff,
        species_id_by_label=NULL,
        species_parser=species_parser,
        sep=sep
    )
}

.directed_species_overlap_scores = function(phy, species_ids, dc_cutoff=0) {
    root_num = get_root_num(phy)
    root_children = get_children_num(phy, root_num)
    if (length(root_num) != 1L || length(root_children) < 2L) {
        stop('Root-position scoring requires a tree with one root and at least two root children.')
    }

    original_edges = phy[['edge']]
    if (length(root_children) == 2L) {
        nonroot_edges = original_edges[original_edges[,1] != root_num,,drop=FALSE]
        unrooted_edges = rbind(nonroot_edges, root_children)
    } else {
        unrooted_edges = original_edges
    }
    adjacency = split(
        c(unrooted_edges[,2], unrooted_edges[,1]),
        c(unrooted_edges[,1], unrooted_edges[,2])
    )
    nodes = sort(unique(as.integer(c(unrooted_edges))))
    if (nrow(unrooted_edges) != length(nodes) - 1L) {
        stop('Invalid cyclic topology in root-position scoring.')
    }
    traversal_root = nodes[[1]]
    parent_by_node = rep(NA_integer_, max(nodes))
    traversal_order = integer(length(nodes))
    traversal_order[[1]] = traversal_root
    visited = logical(max(nodes))
    visited[[traversal_root]] = TRUE
    queue_size = 1L
    next_index = 1L
    while (next_index <= queue_size) {
        node = traversal_order[[next_index]]
        next_nodes = setdiff(
            adjacency[[as.character(node)]],
            parent_by_node[[node]]
        )
        next_nodes = next_nodes[!visited[next_nodes]]
        if (length(next_nodes)) {
            parent_by_node[next_nodes] = node
            visited[next_nodes] = TRUE
            target_indices = queue_size + seq_along(next_nodes)
            traversal_order[target_indices] = next_nodes
            queue_size = queue_size + length(next_nodes)
        }
        next_index = next_index + 1L
    }
    traversal_order = traversal_order[seq_len(queue_size)]
    if (queue_size != length(nodes) || !setequal(traversal_order, nodes)) {
        stop('Invalid disconnected topology in root-position scoring.')
    }

    component_species = new.env(parent=emptyenv(), hash=TRUE)
    component_scores = new.env(parent=emptyenv(), hash=TRUE)
    local_overlap_by_excluded_neighbor = new.env(parent=emptyenv(), hash=TRUE)

    overlap_indicator = function(child_species) {
        score = .max_pairwise_species_overlap(child_species)
        if (is.na(score)) 0L else as.integer(score > dc_cutoff)
    }

    build_species_message = function(from, to) {
        key = paste0(from, '>', to)
        next_nodes = setdiff(adjacency[[as.character(to)]], from)
        next_keys = paste0(to, '>', next_nodes)
        if (to <= length(species_ids)) {
            species = species_ids[[to]]
        } else {
            if (any(!vapply(next_keys, exists, logical(1),
                    envir=component_species, inherits=FALSE))) {
                stop('Invalid species traversal state in root-position scoring.')
            }
            next_species = mget(
                next_keys, envir=component_species, inherits=FALSE
            )
            species = unique(unlist(next_species, use.names=FALSE))
        }
        assign(key, species, envir=component_species)
        invisible(NULL)
    }

    # Species sets do not depend on overlap scores, so calculate every directed
    # species message before scoring nodes. This permits each node's qualifying
    # child pairs to be counted once rather than recomputed for every excluded
    # neighbor.
    for (to in rev(traversal_order[-1L])) {
        build_species_message(parent_by_node[[to]], to)
    }
    for (from in traversal_order) {
        children = adjacency[[as.character(from)]]
        is_child = !is.na(parent_by_node[children]) &
            parent_by_node[children] == from
        children = children[is_child]
        for (to in children) {
            build_species_message(to, from)
        }
    }

    internal_nodes = nodes[nodes > length(species_ids)]
    for (node in internal_nodes) {
        neighbors = adjacency[[as.character(node)]]
        neighbor_species = mget(
            paste0(node, '>', neighbors),
            envir=component_species,
            inherits=FALSE
        )
        num_qualifying_pairs = 0L
        incident_qualifying_pairs = integer(length(neighbors))
        if (length(neighbors) >= 2L) {
            for (first_index in seq_len(length(neighbors) - 1L)) {
                for (second_index in seq.int(first_index + 1L, length(neighbors))) {
                    if (.species_jaccard(
                            neighbor_species[[first_index]],
                            neighbor_species[[second_index]]
                        ) > dc_cutoff) {
                        num_qualifying_pairs = num_qualifying_pairs + 1L
                        incident_qualifying_pairs[[first_index]] =
                            incident_qualifying_pairs[[first_index]] + 1L
                        incident_qualifying_pairs[[second_index]] =
                            incident_qualifying_pairs[[second_index]] + 1L
                    }
                }
            }
        }
        excluded_indicators = as.integer(
            num_qualifying_pairs - incident_qualifying_pairs > 0L
        )
        for (neighbor_index in seq_along(neighbors)) {
            assign(
                paste0(neighbors[[neighbor_index]], '>', node),
                excluded_indicators[[neighbor_index]],
                envir=local_overlap_by_excluded_neighbor
            )
        }
    }

    build_score_message = function(from, to) {
        key = paste0(from, '>', to)
        if (to <= length(species_ids)) {
            score = 0
        } else {
            next_nodes = setdiff(adjacency[[as.character(to)]], from)
            next_keys = paste0(to, '>', next_nodes)
            if (any(!vapply(next_keys, exists, logical(1),
                    envir=component_scores, inherits=FALSE))) {
                stop('Invalid score traversal state in root-position scoring.')
            }
            local_score = get(
                key,
                envir=local_overlap_by_excluded_neighbor,
                inherits=FALSE
            )
            score = as.numeric(local_score + sum(unlist(mget(
                next_keys, envir=component_scores, inherits=FALSE
            ), use.names=FALSE)))
        }
        assign(key, score, envir=component_scores)
        invisible(NULL)
    }

    for (to in rev(traversal_order[-1L])) {
        build_score_message(parent_by_node[[to]], to)
    }
    for (from in traversal_order) {
        children = adjacency[[as.character(from)]]
        is_child = !is.na(parent_by_node[children]) &
            parent_by_node[children] == from
        children = children[is_child]
        for (to in children) {
            build_score_message(to, from)
        }
    }

    score_unrooted_edge = function(node1, node2) {
        key12 = paste0(node1, '>', node2)
        key21 = paste0(node2, '>', node1)
        get(key12, envir=component_scores, inherits=FALSE) +
            get(key21, envir=component_scores, inherits=FALSE) +
            overlap_indicator(list(
                get(key12, envir=component_species, inherits=FALSE),
                get(key21, envir=component_species, inherits=FALSE)
            ))
    }

    vapply(seq_len(nrow(original_edges)), function(edge_index) {
        parent = original_edges[edge_index,1]
        child = original_edges[edge_index,2]
        if (length(root_children) == 2L && parent == root_num) {
            score_unrooted_edge(root_children[[1]], root_children[[2]])
        } else {
            score_unrooted_edge(parent, child)
        }
    }, numeric(1))
}

#' Score every candidate root position by species overlap
#'
#' Uses a bidirectional traversal; the legacy `nslots` argument is accepted but
#' no worker pool is created.
#'
#' @param phy A rooted or unrooted `phylo` tree with unique tip labels.
#'   Multifurcating nodes use their maximum pairwise child-clade overlap.
#' @param nslots Legacy requested worker count.
#' @param species_parser Species-label convention.
#' @param sep Literal separator in gene labels.
#' @return A numeric vector with one score for each row of `phy$edge`, in the
#'   same order. Unrooting or reordering a tree can change these edge indices;
#'   use the same tree object when labeling candidate branches and their scores.
#' @export
get_root_position_dependent_species_overlap_scores = function(
    phy,
    nslots=NULL,
    species_parser='legacy',
    sep='_'
) {
    species_parser = .normalize_species_parser_arg(
        value=species_parser,
        arg_name='species_parser'
    )
    sep = .normalize_single_string_arg(
        value=sep,
        arg_name='sep',
        allow_empty=FALSE
    )

    .validate_phylo_input(
        phy,
        context='phy',
        unique_tips=TRUE
    )
    root_num = get_root_num(phy)
    if (length(root_num) != 1L || length(get_children_num(phy, root_num)) < 2L) {
        stop('phy must contain one root with at least two children.')
    }
    num_edges = nrow(phy[['edge']])
    if (num_edges == 0) {
        return(numeric(0))
    }
    if (!is.null(nslots)) {
        .resolve_parallel_cores(
            requested=nslots,
            max_tasks=num_edges,
            auto_when_missing=FALSE
        )
    }
    species_ids = unname(.tip_species_id_map(
        tip_labels=phy[['tip.label']],
        species_parser=species_parser,
        sep=sep
    ))
    .directed_species_overlap_scores(
        phy=phy,
        species_ids=species_ids,
        dc_cutoff=0
    )
}

#' Read NOTUNG parsable duplication records
#'
#' @param file Path to a NOTUNG parsable-output file.
#' @param mode Event mode. Currently only `"D"` is supported.
#' @return A data frame of duplication-event fields.
#' @export
read_notung_parsable = function(file, mode='D') {
    file = .normalize_single_string_arg(file, 'file', allow_empty=FALSE)
    if (!file.exists(file)) {
        stop('file does not exist: ', file)
    }
    mode = .normalize_choice_arg(
        value=mode,
        arg_name='mode',
        choices='D'
    )
    cols = c('event', 'gn_node', 'lower_sp_node', 'upper_sp_node')
    empty_df = data.frame(matrix(NA_character_, 0, length(cols)), stringsAsFactors=FALSE)
    colnames(empty_df) = cols

    con = base::file(file, "r")
    on.exit(close(con), add=TRUE)
    dup_lines = character(0)
    repeat {
        chunk = readLines(con=con, n=10000L, warn=FALSE)
        if (!length(chunk)) {
            break
        }
        matched_lines = grep("^\\s*#D\\b", chunk, value=TRUE)
        if (length(matched_lines)) {
            dup_lines = c(dup_lines, matched_lines)
        }
    }
    if (length(dup_lines)==0) {
        return(empty_df)
    }
    dup_items = strsplit(dup_lines, "\\s+")
    parsed = lapply(dup_items, function(item_vec) {
        item_values = item_vec[nchar(item_vec)>0]
        if (length(item_values) && item_values[1] == '#D') {
            item_values = item_values[-1]
        }
        if (length(item_values) != length(cols) - 1L) {
            stop(
                'Malformed NOTUNG duplication line: expected 3 fields after #D, got ',
                length(item_values), '. Line: ', paste(item_vec, collapse=' ')
            )
        }
        out = rep(NA_character_, length(cols))
        out[1] = 'D'
        out[1 + seq_along(item_values)] = item_values
        out
    })
    df = data.frame(do.call(rbind, parsed), stringsAsFactors=FALSE)
    rownames(df) = NULL
    colnames(df) = cols
    return(df)
}
