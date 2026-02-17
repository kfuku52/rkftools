# Title     : TODO
# Objective : TODO
# Created by: kef74yk
# Created on: 2019-01-02

get_duplication_confidence_score = function(phy, node_num) {
    children_num = get_children_num(phy, node_num)
    if (length(children_num)==2) {
        child_leaves = vector(mode='list', length(children_num))
        for (j in seq_along(children_num)) {
            child_leaves[[j]] = get_tip_labels(phy, children_num[j])
            child_leaves[[j]] = sub('_', ' ', child_leaves[[j]])
            child_leaves[[j]] = sub('_.*', '', child_leaves[[j]])
        }
        sp_intersect = intersect(child_leaves[[1]], child_leaves[[2]])
        sp_union = union(child_leaves[[1]], child_leaves[[2]])
        dc_score = length(sp_intersect) / length(sp_union)
    } else {
        dc_score = NA
    }
    return(dc_score)
}

.species_overlap_score_core = function(phy, dc_cutoff=0, dc_score_fun) {
    tip_count = length(phy[['tip.label']])
    max_node_id = max(phy[['edge']])
    if (max_node_id <= tip_count) {
        return(0)
    }
    internal_nodes = seq.int(tip_count + 1, max_node_id)
    dc_scores = vapply(internal_nodes, function(int_node) {
        dc_score_fun(phy, int_node)
    }, numeric(1))
    sum(dc_scores > dc_cutoff)
}

.tip_species_id_map = function(tip_labels) {
    species_labels = sub('_', ' ', tip_labels)
    species_labels = sub('_.*', '', species_labels)
    species_levels = unique(species_labels)
    species_ids = match(species_labels, species_levels)
    names(species_ids) = tip_labels
    species_ids
}

.species_overlap_score_fast_impl = function(phy, dc_cutoff=0, species_id_by_label=NULL) {
    tip_count = length(phy[['tip.label']])
    if (tip_count == 0) {
        return(0)
    }

    edges = phy[['edge']]
    if (!nrow(edges)) {
        return(0)
    }
    children_by_parent = split(edges[,2], edges[,1])
    if (!length(children_by_parent)) {
        return(0)
    }

    if (is.null(species_id_by_label)) {
        tip_species_ids = unname(.tip_species_id_map(phy[['tip.label']]))
    } else {
        tip_species_ids = as.integer(species_id_by_label[phy[['tip.label']]])
        if (anyNA(tip_species_ids)) {
            tip_species_ids = unname(.tip_species_id_map(phy[['tip.label']]))
        }
    }

    max_node_id = max(edges)
    species_cache = vector(mode='list', length=max_node_id)
    for (tip_index in seq_len(tip_count)) {
        species_cache[[tip_index]] = tip_species_ids[tip_index]
    }

    resolve_species = function(node_num) {
        cached = species_cache[[node_num]]
        if (!is.null(cached)) {
            return(cached)
        }

        children_num = children_by_parent[[as.character(node_num)]]
        if (is.null(children_num) || !length(children_num)) {
            out = integer(0)
        } else {
            child_species = lapply(children_num, function(cn) {
                resolve_species(as.integer(cn))
            })
            out = unique(unlist(child_species, use.names=FALSE))
        }
        species_cache[[node_num]] <<- out
        out
    }

    internal_nodes = as.integer(names(children_by_parent))
    internal_nodes = internal_nodes[internal_nodes > tip_count]
    if (!length(internal_nodes)) {
        return(0)
    }

    overlap_count = 0L
    for (node_num in internal_nodes) {
        children_num = children_by_parent[[as.character(node_num)]]
        if (length(children_num) != 2) {
            next
        }

        child_sp1 = resolve_species(as.integer(children_num[1]))
        child_sp2 = resolve_species(as.integer(children_num[2]))
        if (length(child_sp1) <= length(child_sp2)) {
            sp_intersect = sum(child_sp1 %in% child_sp2)
        } else {
            sp_intersect = sum(child_sp2 %in% child_sp1)
        }
        sp_union = length(child_sp1) + length(child_sp2) - sp_intersect
        if (sp_union == 0) {
            next
        }
        dc_score = sp_intersect / sp_union
        if (!is.na(dc_score) && dc_score > dc_cutoff) {
            overlap_count = overlap_count + 1L
        }
    }
    as.numeric(overlap_count)
}

get_species_overlap_score = function(phy, dc_cutoff=0) {
    # this function assumes that leaf names are: GENUS_SPECIES_GENEID (e.g. Bos_taurus_AF492351.1)
    .species_overlap_score_fast_impl(
        phy=phy,
        dc_cutoff=dc_cutoff,
        species_id_by_label=NULL
    )
}

get_root_position_dependent_species_overlap_scores = function(phy, nslots) {
    if (!requireNamespace('phytools', quietly=TRUE)) {
        stop("'phytools' package not found, please install it.")
    }

    num_edges = nrow(phy[['edge']])
    if (num_edges == 0) {
        return(numeric(0))
    }
    edge_indices = seq_len(num_edges)
    num_parallel = .resolve_parallel_cores(
        requested=nslots,
        max_tasks=num_edges,
        auto_when_missing=FALSE
    )
    species_id_by_label = .tip_species_id_map(phy[['tip.label']])

    score_edge = function(i, phy_obj, dc_cutoff=0) {
        rt = phytools::reroot(tree=phy_obj, node.number=phy_obj[['edge']][i,2])
        .species_overlap_score_fast_impl(
            phy=rt,
            dc_cutoff=dc_cutoff,
            species_id_by_label=species_id_by_label
        )
    }

    if (num_parallel == 1) {
        species_overlap_scores = vapply(edge_indices, function(i) {
            score_edge(i=i, phy_obj=phy, dc_cutoff=0)
        }, numeric(1))
        return(species_overlap_scores)
    }

    if (.Platform$OS.type != "windows") {
        species_overlap_scores = unlist(parallel::mclapply(
            X=edge_indices, FUN=score_edge,
            phy_obj=phy, dc_cutoff=0,
            mc.cores=num_parallel
        ), use.names=FALSE)
        return(species_overlap_scores)
    }

    cluster = parallel::makeCluster(num_parallel, 'PSOCK', outfile='')
    on.exit(parallel::stopCluster(cluster), add=TRUE)
    species_overlap_scores = unlist(parallel::parLapply(
        cl=cluster, X=edge_indices, fun=score_edge,
        phy_obj=phy, dc_cutoff=0
    ), use.names=FALSE)
    return(species_overlap_scores)
}

read_notung_parsable = function(file, mode='D') {
    cols = c('event', 'gn_node', 'lower_sp_node', 'upper_sp_node')
    empty_df = data.frame(matrix(NA_character_, 0, length(cols)), stringsAsFactors=FALSE)
    colnames(empty_df) = cols

    if (mode!='D') {
        cat('mode', mode, 'is not supported.')
        return(empty_df)
    }

    con = base::file(file, "r")
    on.exit(close(con), add=TRUE)
    event_lines = readLines(con=con, warn=FALSE)
    dup_positions = grep("^#D", event_lines)
    if (length(dup_positions)<=1) {
        return(empty_df)
    }

    dup_lines = event_lines[dup_positions]
    dup_items = strsplit(dup_lines[2:length(dup_lines)], "\\s+")
    parsed = lapply(dup_items, function(item_vec) {
        item_values = item_vec[nchar(item_vec)>0]
        out = rep(NA_character_, length(cols))
        num_copy = min(length(item_values), length(cols))
        if (num_copy>0) {
            out[seq_len(num_copy)] = item_values[seq_len(num_copy)]
        }
        out
    })
    df = data.frame(do.call(rbind, parsed), stringsAsFactors=FALSE)
    rownames(df) = NULL
    colnames(df) = cols
    df$event = 'D'
    return(df)
}
