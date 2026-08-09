# Title     : TODO
# Objective : TODO
# Created by: kf
# Created on: 3/22/18

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

get_node_name_by_num = function(phy, node_num) {
    node_names = c(phy[['tip.label']], phy$node.label)
    out = character(0)
    max_node_num = length(node_names)
    for (num in suppressWarnings(as.integer(node_num))) {
        if (is.na(num) || num < 1L || num > max_node_num) {
            next
        }
        out = c(out, node_names[num])
    }
    return(out)
}

get_root_num = function(phy) {
    root_num = setdiff(phy[['edge']][,1], phy[['edge']][,2])
    return(root_num)
}

get_children_num = function(phy, node_num) {
    out = integer(0)
    max_node_num = max(phy[['edge']])
    for (nn in suppressWarnings(as.integer(node_num))) {
        if (is.na(nn) || nn < 1L || nn > max_node_num) {
            next
        }
        out = c(out, phy[['edge']][(phy[['edge']][,1]==nn),2])
    }
    return(out)
}

is_root = function(phy, node_num) {
    root_num = get_root_num(phy)
    return(node_num==root_num)
}

is_leaf = function(phy, node_num) {
    node_num = suppressWarnings(as.integer(node_num))
    max_node_num = max(phy[['edge']])
    if (length(node_num) != 1L || is.na(node_num) || node_num < 1L || node_num > max_node_num) {
        return(FALSE)
    }
    return(node_num <= length(phy[['tip.label']]))
} 

get_descendent_num = function(phy, node_num, leaf_only=FALSE) {
    leaf_only = .normalize_single_logical_arg(
        value=leaf_only,
        arg_name='leaf_only'
    )
    descendent_nums = c()
    children_nums = get_children_num(phy, node_num)
    while(!all(is.na(children_nums))) {
        children_nums = children_nums[!is.na(children_nums)]
        descendent_nums = c(descendent_nums, children_nums)
            childrens = c()
            for (nn in children_nums) {
                childrens = c(childrens, get_children_num(phy, nn))
            }
            children_nums = childrens
        }
    descendent_nums = sort(unique(descendent_nums))
    if (leaf_only) {
        descendent_nums = descendent_nums[descendent_nums <= length(phy[['tip.label']])]
    }
    return(descendent_nums)
}

get_parent_num = function(phy, node_num) {
    out = integer(0)
    max_node_num = max(phy[['edge']])
    for (nn in suppressWarnings(as.integer(node_num))) {
        if (is.na(nn) || nn < 1L || nn > max_node_num) {
            next
        }
        out = c(out, phy[['edge']][(phy[['edge']][,2]==nn),1])
    }
    return(out)
}

get_sister_num = function(phy, node_num) {
    out = integer(0)
    max_node_num = max(phy[['edge']])
    for (nn in suppressWarnings(as.integer(node_num))) {
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

get_ancestor_num = function(phy, node_num) {
    node_num = suppressWarnings(as.integer(node_num))
    max_node_num = max(phy[['edge']])
    if (length(node_num) != 1 || is.na(node_num) || node_num < 1 || node_num > max_node_num) {
        return(integer(0))
    }
    ancestor_num = c()
    root_num = get_root_num(phy)
    current_node_num = node_num
    for (i in seq_len(phy[['Nnode']])) {
        if (!length(current_node_num) || is.na(current_node_num)) {
            break
        }
        if (current_node_num==root_num) {
            break
        }
        parent_num = get_parent_num(phy, current_node_num)
        if (!length(parent_num)) {
            break
        }
        ancestor_num = c(ancestor_num, parent_num)
        current_node_num = parent_num
    }
    return(ancestor_num)
}

# alias
collapse_short_external_edges = function(tree, threshold=1e-6) {
    return(pad_short_edges(tree, threshold=threshold, external_only=TRUE))
}

pad_short_edges = function(tree, threshold=1e-6, external_only=FALSE) {
    external_only = .normalize_single_logical_arg(
        value=external_only,
        arg_name='external_only'
    )
    threshold = suppressWarnings(as.numeric(threshold))
    if (length(threshold) != 1L || is.na(threshold) || !is.finite(threshold) || threshold < 0) {
        stop('threshold must be a single finite non-negative numeric value.')
    }
    if (!ape::is.binary(tree)) {
        stop('tree must be binary in pad_short_edges().')
    }
    if (is.null(tree[['edge.length']]) || length(tree[['edge.length']]) != nrow(tree[['edge']])) {
        stop('tree must contain one branch length per edge in pad_short_edges().')
    }
    out_tree = tree
    edge_idx = seq_len(nrow(out_tree$edge))
    is_target_edge = rep(TRUE, nrow(out_tree$edge))
    if (external_only) {
        is_target_edge = is_target_edge & (out_tree$edge[,2]<=length(out_tree$tip.label))
    }
    edge_lengths = out_tree[['edge.length']][is_target_edge]
    non_na_edge_lengths = edge_lengths[!is.na(edge_lengths)]
    min_eel = if (length(non_na_edge_lengths)) min(non_na_edge_lengths) else NA_real_
    cat('Minimum edge length:', min_eel, '\n')
    if (any(is.na(edge_lengths))) {
        cat('NA edge lengths were ignored when searching short edges.\n')
    }
    is_short_eel = is_target_edge & (!is.na(out_tree$edge.length)) & (out_tree$edge.length < threshold)
    num_short_eel = sum(is_short_eel)
    cat('Number of short edges ( length <', threshold, '):', num_short_eel, '\n')
    if (num_short_eel>0) {
        short_eel_idx = edge_idx[is_short_eel]
        for (i in short_eel_idx) {
            if (!is.na(out_tree$edge.length[i]) && (out_tree$edge.length[i] < threshold)) {
                shift_value = threshold - out_tree$edge.length[i]
                sister_node_num = get_sister_num(out_tree, out_tree$edge[i,2])
                sister_edge_idx = edge_idx[out_tree$edge[,2]==sister_node_num]
                root_num = get_root_num(out_tree)
                flag = TRUE
                flag_root = FALSE
                current_idx = i
                while (flag==TRUE) {
                    parent_node_num = out_tree$edge[current_idx,1]
                    parent_edge_idx = edge_idx[out_tree$edge[,2]==parent_node_num]
                    parent_edge_length = out_tree$edge.length[parent_edge_idx]
                    if (parent_node_num==root_num) {
                        flag = FALSE
                        flag_root = TRUE
                    } else if (is.na(parent_edge_length)) {
                        flag = FALSE
                        flag_root = TRUE
                    } else if (parent_edge_length>=threshold+shift_value) {
                        flag = FALSE
                    } else {
                        current_idx = edge_idx[out_tree$edge[,2]==parent_node_num]
                    }
                }

                out_tree$edge.length[i] = out_tree$edge.length[i] +shift_value
                out_tree$edge.length[sister_edge_idx] = out_tree$edge.length[sister_edge_idx] + shift_value
                if (flag_root) {
                    cat('Adding branch length to subroot edges,', i, 'and', sister_edge_idx, '\n')
                } else {
                    cat('Transfering branch length from edge', parent_edge_idx, 'to', i, 'and', sister_edge_idx, '\n')
                    out_tree$edge.length[parent_edge_idx] = out_tree$edge.length[parent_edge_idx] - shift_value
                }
            }
        }
    }
    return(out_tree)
}

get_tip_labels = function(phy, node_num, out=NULL) {
    num_leaf = length(phy[['tip.label']])
    node_nums = suppressWarnings(as.integer(node_num))
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
            subtree = ape::extract.clade(phy, nn)
            tip_labels = c(tip_labels, subtree$tip.label)
        } else {
            tip_labels = c(tip_labels, phy[['tip.label']][nn])
        }
    }
    if (!is.null(out)) {
        tip_labels = c(out, tip_labels)
    }
    return(tip_labels)
}

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

get_node_age = function(phy, node_num) {
    if (!ape::is.ultrametric(phy)) {
        stop('phy must be ultrametric in get_node_age().')
    }
    node_num = suppressWarnings(as.integer(node_num))
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

get_phy2_root_in_phy1 = function(phy1, phy2, nslots=NULL, mode=c("node_num", "index")) {
    mode_name = match.arg(mode)
    num_slots = .resolve_parallel_cores(
        requested=nslots,
        max_tasks=nrow(phy1[['edge']]),
        auto_when_missing=FALSE
    )
    if (! ape::is.rooted(phy2)) {
        stop('phy2 is unrooted.')
    }
    if (! identical(sort(phy1$tip.label), sort(phy2$tip.label))) {
        stop('phy1 and phy2 have different sets of leaves.')
    }

    get_root_num_local = function(phy_obj) {
        setdiff(phy_obj[['edge']][,1], phy_obj[['edge']][,2])
    }
    get_children_num_local = function(phy_obj, node_num) {
        phy_obj[['edge']][(phy_obj[['edge']][,1]==node_num),2]
    }
    get_tip_labels_local = function(phy_obj, node_num) {
        num_leaf = length(phy_obj[['tip.label']])
        if (node_num > num_leaf) {
            return(ape::extract.clade(phy_obj, node_num)$tip.label)
        }
        phy_obj[['tip.label']][node_num]
    }
    get_root_split = function(phy_obj) {
        root_num = get_root_num_local(phy_obj)
        children = get_children_num_local(phy_obj, root_num)
        sort(get_tip_labels_local(phy_obj, children[1]))
    }

    all_tips_sorted = sort(phy1$tip.label)
    ref_root_split = get_root_split(phy2)
    is_matching_root = function(i, phy_obj, ref_split, all_tips) {
        rerooted = phytools::reroot(tree=phy_obj, node.number=phy_obj$edge[i,2])
        split_a = get_root_split(rerooted)
        if (identical(split_a, ref_split)) {
            return(TRUE)
        }
        split_b = sort(setdiff(all_tips, split_a))
        identical(split_b, ref_split)
    }

    edge_indices = seq_len(nrow(phy1$edge))
    matched_index = NA_integer_
    if (num_slots == 1L || length(edge_indices) <= 1L) {
        for (i in edge_indices) {
            if (is_matching_root(i, phy_obj=phy1, ref_split=ref_root_split, all_tips=all_tips_sorted)) {
                matched_index = i
                break
            }
        }
    } else {
        num_slots = min(num_slots, length(edge_indices))
        if (.Platform$OS.type != 'windows') {
            is_match = unlist(parallel::mclapply(
                X=edge_indices,
                FUN=is_matching_root,
                phy_obj=phy1,
                ref_split=ref_root_split,
                all_tips=all_tips_sorted,
                mc.cores=num_slots
            ), use.names=FALSE)
        } else {
            cl = parallel::makeCluster(num_slots)
            on.exit(parallel::stopCluster(cl), add=TRUE)
            is_match = unlist(parallel::parLapply(
                cl=cl, X=edge_indices, fun=is_matching_root,
                phy_obj=phy1, ref_split=ref_root_split, all_tips=all_tips_sorted
            ), use.names=FALSE)
        }
        hit = which(is_match)[1]
        if (length(hit) && !is.na(hit)) {
            matched_index = edge_indices[hit]
        }
    }

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

get_rooted_newick = function(t, madr, rho) {
    madr_numeric = suppressWarnings(as.numeric(madr))
    if (length(madr_numeric) != 1 || is.na(madr_numeric) || !is.finite(madr_numeric) ||
            madr_numeric != as.integer(madr_numeric) || madr_numeric < 1 ||
            madr_numeric > nrow(t$edge)) {
        stop('madr must be a single edge index in [1, ', nrow(t$edge), '] in get_rooted_newick().')
    }
    madr = as.integer(madr_numeric)
    if (length(rho) != nrow(t$edge)) {
        stop('rho must have length ', nrow(t$edge), ' in get_rooted_newick().')
    }
    if (!is.numeric(rho) || is.na(rho[madr]) || !is.finite(rho[madr])) {
        stop('rho[madr] must be finite in get_rooted_newick().')
    }
    notu <- length(t$tip.label)
    dis <- ape::dist.nodes(t)
    pp <- rho[madr]*t$edge.length[madr]
    nn <- t$edge[madr,]
    rt <- phytools::reroot(t,nn[2], pos = pp)
    rooted_newick <- ape::write.tree(rt)
    dd <- dis[1:notu,nn]
    sp <- dd[,1]<dd[,2]
    otu2root <- vector(mode="numeric",notu)
    otu2root[sp] <- dd[sp,1] + pp
    otu2root[!sp] <- dd[!sp,1] - pp
    ccv <- 100*stats::sd(otu2root)/mean(otu2root)
    return(list(rooted_newick, rt, ccv))
}

.format_mad_result = function(t, rho, bad, output_mode=NULL) {
    if (!any(is.finite(bad))) {
        stop('MAD could not score any branch; check that the tree has positive pairwise distances.')
    }
    jj = sort(bad, index.return=TRUE)
    tf = bad == jj$x[1]
    tf[is.na(tf)] = FALSE
    nroots = sum(tf)
    if (nroots > 1) {
        warning("More than one possible root position. Multiple newick strings printed")
    }
    madr = which(tf)
    rai = if (length(jj$x) >= 2L && is.finite(jj$x[2]) && jj$x[2] != 0) {
        jj$x[1] / jj$x[2]
    } else {
        NA_real_
    }
    badr = bad[tf]

    rt = vector("list", nroots)
    ccv = numeric(nroots)
    rooted_newick = character(nroots)
    for (i in seq_along(madr)) {
        out = get_rooted_newick(t, madr[i], rho)
        rooted_newick[i] = out[[1]]
        rt[[i]] = out[[2]]
        ccv[i] = out[[3]]
    }
    rooted_newick = sub(')Root;', ');', rooted_newick)

    if (is.null(output_mode) || output_mode == 'newick') {
        return(rooted_newick)
    } else if (output_mode == 'stats') {
        root_stats = data.frame(ambiguity_index=rai, clock_cv=ccv, ancestor_deviation=badr, n_roots=nroots)
        return(list(rooted_newick, root_stats))
    } else if (output_mode == 'full') {
        root_stats = data.frame(ambiguity_index=rai, clock_cv=ccv, ancestor_deviation=badr, n_roots=nroots)
        return(list(rooted_newick, root_stats, t, madr, bad, rt))
    } else if (output_mode == 'custom') {
        root_stats = data.frame(ambiguity_index=rai, clock_cv=ccv, ancestor_deviation=badr, n_roots=nroots)
        return(list(rooted_newick, root_stats, t, madr, bad, rt, rho))
    }
    return(rooted_newick)
}

.handle_mad_duplicate_tip = function(t, output_mode=NULL, rerun_fun) {
    notu = length(t$tip.label)
    dis = ape::dist.nodes(t)
    sdis = dis[seq_len(notu), seq_len(notu)]
    ii = which(sdis == 0, arr.ind=TRUE)
    k = which(ii[,1] != ii[,2])
    if (!length(k)) {
        return(NULL)
    }

    dup_row = ii[k[1],1]
    dup_col = ii[k[1],2]
    vv = c(
        paste('@#', t$tip.label[dup_row], '@#', sep=''),
        paste('(', t$tip.label[dup_row], ':0,', t$tip.label[dup_col], ':0)', sep='')
    )
    st = ape::drop.tip(t, dup_col)
    st$tip.label[st$tip.label == t$tip.label[dup_row]] = vv[1]
    res = rerun_fun(st, output_mode)
    if (is.list(res)) {
        res[[1]] = sub(vv[1], vv[2], res[[1]], fixed=TRUE)
    } else {
        res = sub(vv[1], vv[2], res, fixed=TRUE)
    }
    return(res)
}

.calc_mad_branch_stats = function(br, t, dis, sdis, disbr, nodeids, otuids, npairs, notu, nbranch) {
    dij = t$edge.length[br]
    if (dij == 0) {
        return(c(rho=NA_real_, bad=NA_real_))
    }

    rbca = numeric(npairs)
    i = t$edge[br,1]
    j = t$edge[br,2]
    sp = dis[seq_len(notu),i] < dis[seq_len(notu),j]
    dbc = matrix(sdis[sp,!sp], nrow=sum(sp), ncol=sum(!sp))
    dbi = replicate(ncol(dbc), dis[(seq_len(notu))[sp],i])
    rho_br = sum((dbc - 2 * dbi) * dbc^-2) / (2 * dij * sum(dbc^-2))
    rho_br = min(max(0, rho_br), 1)

    dab = dbi + (dij * rho_br)
    ndab = length(dab)
    rbca[seq_len(ndab)] = as.vector(2 * dab / dbc - 1)

    bcsp = rbind(sp, !sp)
    ij = c(i, j)
    counter = ndab
    i2p = matrix(FALSE, nrow=nbranch + 1, ncol=notu)
    for (w in c(1, 2)) {
        if (sum(bcsp[w,]) >= 2) {
            disbrw = disbr[,ij[w]]
            pairids = otuids[bcsp[w,]]
            for (z in pairids) {
                i2p[,z] = (disbr[z,] + disbrw == disbrw[z])
            }
            for (z_idx in seq_len(length(pairids) - 1)) {
                p1 = pairids[z_idx]
                disp1 = dis[p1,]
                pan = nodeids[i2p[,p1]]
                for (y_idx in (z_idx + 1):length(pairids)) {
                    p2 = pairids[y_idx]
                    pan1 = pan[i2p[pan,p2]]
                    an = pan1[which.max(disbrw[pan1])]
                    counter = counter + 1
                    rbca[counter] = 2 * disp1[an] / disp1[p2] - 1
                }
            }
        }
    }
    if (length(rbca) != npairs) {
        stop("Unexpected number of pairs.")
    }

    bad_br = sqrt(mean(rbca^2))
    c(rho=rho_br, bad=bad_br)
}

.prepare_mad_tree = function(unrooted_newick) {
    if (!requireNamespace('ape', quietly = TRUE)) {
        stop("'ape' package not found, please install it to run MAD")
    }
    if (!requireNamespace('phytools', quietly = TRUE)) {
        stop("'phytools' package not found, please install it to run MAD")
    }

    t = if (inherits(unrooted_newick, "phylo")) unrooted_newick else ape::read.tree(text=unrooted_newick)
    if (is.null(t) || !inherits(t, 'phylo')) {
        stop('unrooted_newick must be a valid Newick string or an object of class "phylo".')
    }
    if (anyDuplicated(t$tip.label)) {
        stop('Input tree tip labels must be unique for MAD.')
    }
    if (ape::is.rooted(t)) {
        t = ape::unroot(t)
    }
    if (!ape::is.binary(t)) {
        warning("Input tree is not binary! Internal multifurcations will be converted to branches of length zero and identical OTUs will be collapsed!")
        t = ape::multi2di(t)
    }
    if (is.null(t$edge.length)) {
        stop("Input tree has no branch lengths. MAD requires branch lengths.")
    }
    if (length(t$edge.length) != nrow(t$edge) || any(is.na(t$edge.length))) {
        stop("Input tree contains missing branch lengths. MAD requires complete branch lengths.")
    }
    if (any(!is.finite(t$edge.length))) {
        stop('Input tree contains non-finite branch lengths. MAD requires finite branch lengths.')
    }
    has_negative = (t$edge.length < 0)
    if (any(has_negative)) {
        warning("Input tree contains negative branch lengths. They will be converted to zeros!")
        t$edge.length[has_negative] = 0
    }
    if (all(t$edge.length == 0)) {
        stop('Input tree has no positive branch lengths. MAD cannot root an all-zero tree.')
    }
    return(t)
}

.compute_mad_scores = function(t, ncpu=1, use_parallel=FALSE) {
    notu = length(t$tip.label)
    nbranch = nrow(t$edge)
    dis = ape::dist.nodes(t)
    sdis = dis[seq_len(notu), seq_len(notu)]

    t2 = t
    t2$edge.length = rep(1, nbranch)
    disbr = ape::dist.nodes(t2)
    nodeids = seq_len(nbranch + 1)
    otuids = seq_len(notu)
    npairs = notu * (notu - 1) / 2

    mad_branch_stat_fun = .calc_mad_branch_stats
    calc_branch_stats = function(br) {
        do.call(
            what=mad_branch_stat_fun,
            args=list(
                br=br, t=t, dis=dis, sdis=sdis, disbr=disbr,
                nodeids=nodeids, otuids=otuids, npairs=npairs,
                notu=notu, nbranch=nbranch
            )
        )
    }

    num_parallel = .resolve_parallel_cores(
        requested=ncpu,
        max_tasks=nbranch,
        auto_when_missing=FALSE
    )
    branch_ids = seq_len(nbranch)
    if (!use_parallel || num_parallel == 1 || nbranch == 1) {
        result_list = lapply(branch_ids, calc_branch_stats)
    } else {
        num_parallel = min(num_parallel, nbranch)
        if (.Platform$OS.type != "windows") {
            result_list = parallel::mclapply(
                X=branch_ids, FUN=calc_branch_stats,
                mc.cores=num_parallel
            )
        } else {
            result_list = local({
                cluster = parallel::makeCluster(num_parallel)
                on.exit(parallel::stopCluster(cluster), add=TRUE)
                parallel::parLapply(cluster, branch_ids, calc_branch_stats)
            })
        }
    }
    results = do.call(rbind, result_list)
    if (is.null(dim(results))) {
        results = matrix(results, nrow=1)
    }
    list(rho=results[,1], bad=results[,2])
}

.run_mad_with_tree = function(t, output_mode, ncpu, use_parallel, rerun_fun) {
    dup_res = .handle_mad_duplicate_tip(
        t=t,
        output_mode=output_mode,
        rerun_fun=rerun_fun
    )
    if (!is.null(dup_res)) {
        return(dup_res)
    }

    if (!use_parallel) {
        gc()
    }
    scores = .compute_mad_scores(t=t, ncpu=ncpu, use_parallel=use_parallel)
    return(.format_mad_result(t=t, rho=scores[['rho']], bad=scores[['bad']], output_mode=output_mode))
}

MAD <- function(unrooted_newick,output_mode){
    # this function was modified from the original MAD function from:
    # https://www.mikrobio.uni-kiel.de/de/ag-dagan/ressourcen
    if(nargs()==0){ #print help message
        return(cat("Minimal Ancestor Deviation (MAD) rooting","","Usage: res <- MAD(unrooted_newick,output_mode)","",
        "unrooted_newick: Unrooted tree string in newick format or a tree object of class 'phylo'","",
        "output_mode: Amount of information to return.", "  If 'newick' (default) only the rooted newick string",
        "  If 'stats' also a structure with the ambiguity index, clock cv, the minimum ancestor deviation and the number of roots",
        "  If 'full' also an unrooted tree object, the index of the root branch, the branch ancestor deviations and a rooted tree object",
        "","res: a list with the results containing one ('newick'), two ('stats') or six elements ('full')","",
        "Dependencies: 'ape' and 'phytools'","","Version: 1.1, 03-May-2017",sep="\n"))
    }
    mode = if (missing(output_mode)) NULL else output_mode
    if (!is.null(mode)) {
        mode = .normalize_single_string_arg(
            value=mode,
            arg_name='output_mode',
            allow_empty=FALSE
        )
        mode = match.arg(mode, c('newick', 'stats', 'full', 'custom'))
    }
    t <- .prepare_mad_tree(unrooted_newick)
    return(.run_mad_with_tree(
        t=t, output_mode=mode, ncpu=1, use_parallel=FALSE,
        rerun_fun=function(tree_obj, mode) MAD(tree_obj, output_mode=mode)
    ))
}

MAD_parallel = function(unrooted_newick, output_mode, ncpu=NULL) {
    mode = if (missing(output_mode)) NULL else output_mode
    if (!is.null(mode)) {
        mode = .normalize_single_string_arg(
            value=mode,
            arg_name='output_mode',
            allow_empty=FALSE
        )
        mode = match.arg(mode, c('newick', 'stats', 'full', 'custom'))
    }

    t = .prepare_mad_tree(unrooted_newick)
    num_parallel = .resolve_parallel_cores(
        requested=ncpu,
        max_tasks=nrow(t$edge),
        auto_when_missing=TRUE
    )
    if (is.null(ncpu) && nrow(t$edge) < 256L) {
        num_parallel = 1L
    }
    return(.run_mad_with_tree(
        t=t, output_mode=mode, ncpu=num_parallel, use_parallel=TRUE,
        rerun_fun=function(tree_obj, mode) MAD_parallel(tree_obj, output_mode=mode, ncpu=num_parallel)
    ))
}

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

get_species_name = function(a, species_parser='legacy', sep='_') {
    species_parser = .normalize_species_parser_arg(
        value=species_parser,
        arg_name='species_parser'
    )
    sep = .normalize_single_string_arg(
        value=sep,
        arg_name='sep',
        allow_empty=FALSE
    )
    parsed = .parse_species_labels(
        labels=a,
        species_parser=species_parser,
        sep=sep,
        output_sep=' ',
        require_gene=FALSE,
        fallback_label=TRUE
    )
    return(parsed[['species_labels']])
}

get_species_names = function(phy, sep='_', species_parser='legacy') {
    sep = .normalize_single_string_arg(
        value=sep,
        arg_name='sep',
        allow_empty=FALSE
    )
    species_parser = .normalize_species_parser_arg(
        value=species_parser,
        arg_name='species_parser'
    )
    parsed = .parse_species_labels(
        labels=phy[['tip.label']],
        species_parser=species_parser,
        sep=sep,
        output_sep=sep,
        require_gene=FALSE,
        fallback_label=FALSE
    )
    species_names = parsed[['species_labels']]
    if (any(!parsed[['parsed_ok']])) {
        bad_labels = phy[['tip.label']][!parsed[['parsed_ok']]]
        for (i in seq_along(bad_labels)) {
            warning(
                'leaf name could not be interpreted with species_parser="',
                species_parser,
                '": ',
                bad_labels[[i]]
            )
        }
    }
    return(species_names)
}

leaf2species = function(leaf_names, use_underbar=FALSE, species_parser='legacy', sep='_') {
    use_underbar = .normalize_single_logical_arg(
        value=use_underbar,
        arg_name='use_underbar'
    )
    species_parser = .normalize_species_parser_arg(
        value=species_parser,
        arg_name='species_parser'
    )
    sep = .normalize_single_string_arg(
        value=sep,
        arg_name='sep',
        allow_empty=FALSE
    )
    parsed = .parse_species_labels(
        labels=leaf_names,
        species_parser=species_parser,
        sep=sep,
        output_sep=if (use_underbar) '_' else ' ',
        require_gene=TRUE,
        fallback_label=FALSE
    )
    species_names = parsed[['species_labels']]
    if (any(!parsed[['parsed_ok']])) {
        bad_labels = leaf_names[!parsed[['parsed_ok']]]
        for (i in seq_along(bad_labels)) {
            warning(
                'leaf name could not be interpreted as a species-bearing label with species_parser="',
                species_parser,
                '": ',
                bad_labels[[i]],
                '\n'
            )
        }
    }
    return(species_names)
}

contains_polytomy = function(phy) {
    if (max(table(phy[['edge']][,1]))>2) {
        is_polytomy = TRUE
    } else {
        is_polytomy = FALSE
    }
    return(is_polytomy)
}

has_same_leaves = function(phy1, phy1_node, phy2, phy2_node) {
    if (!identical(sort(phy1$tip.label), sort(phy2$tip.label))) {
        stop('phy1 and phy2 must contain identical tip-label multisets.')
    }
    phy1_leaves = sort(get_tip_labels(phy1, phy1_node))
    phy2_leaves = sort(get_tip_labels(phy2, phy2_node))
    is_same_leaves = identical(phy1_leaves, phy2_leaves)
    return(is_same_leaves)
}

multi2bi_node_number_transfer = function(multifurcated_tree, bifurcated_tree) {
    mtree = multifurcated_tree
    btree = bifurcated_tree
    if (!identical(mtree[['tip.label']], btree[['tip.label']])) {
        stop('multifurcated_tree and bifurcated_tree must have identical tip-label order.')
    }
    if (!is_same_root(mtree, btree)) {
        stop('multifurcated_tree and bifurcated_tree must have the same root split.')
    }
    topo_dist = ape::dist.topo(ape::unroot(mtree), ape::unroot(btree), method='PH85')
    if (is.na(topo_dist)) {
        stop('Unable to compute topological distance in multi2bi_node_number_transfer().')
    }
    internal_node_counts = table(mtree[['edge']][,1])
    polytomy_parents = as.integer(names(internal_node_counts)[internal_node_counts > 2])
    if (!length(polytomy_parents)) {
        return(data.frame(mtree_node=integer(0), btree_node=integer(0)))
    }
    cat('Polytomy parent nodes:', polytomy_parents, '\n')
    df = data.frame()
    for (mtree_pp in polytomy_parents) {
        mtree_pp_leaves = sort(get_tip_labels(mtree, mtree_pp))
        matched = FALSE
        for (btree_internal_node in btree$edge[,1]) {
            btree_in_leaves = sort(get_tip_labels(btree, btree_internal_node))
            if (length(mtree_pp_leaves)==length(btree_in_leaves)) {
                if (all(mtree_pp_leaves==btree_in_leaves)) {
                    df = rbind(df, data.frame(mtree_node=mtree_pp, btree_node=btree_internal_node))
                    matched = TRUE
                    break
                }
            }
        }
        if (!matched) {
            warning('Failed to map polytomy parent node to bifurcated tree: ', mtree_pp)
        }
    }
    return(df)
}

.encode_tip_signature = function(tip_labels) {
    tip_labels = sort(as.character(tip_labels))
    paste0(nchar(tip_labels), ':', tip_labels, collapse='|')
}

.get_node_tip_sets = function(phy) {
    if (!inherits(phy, 'phylo') || is.null(phy[['edge']]) || !nrow(phy[['edge']])) {
        stop('phy must be a non-empty object of class "phylo".')
    }
    num_tip = length(phy[['tip.label']])
    max_node = max(phy[['edge']])
    children_by_parent = split(phy[['edge']][,2], phy[['edge']][,1])
    tip_cache = vector('list', max_node)
    for (tip_num in seq_len(num_tip)) {
        tip_cache[[tip_num]] = phy[['tip.label']][[tip_num]]
    }
    root_num = get_root_num(phy)
    if (length(root_num) != 1L) {
        stop('Invalid phylo topology while computing clade signatures.')
    }
    traversal_order = integer(max_node)
    traversal_order[[1]] = root_num
    queue_size = 1L
    next_index = 1L
    while (next_index <= queue_size) {
        children = children_by_parent[[as.character(traversal_order[[next_index]])]]
        if (length(children)) {
            target_indices = queue_size + seq_along(children)
            traversal_order[target_indices] = children
            queue_size = queue_size + length(children)
        }
        next_index = next_index + 1L
    }
    traversal_order = traversal_order[seq_len(queue_size)]
    if (queue_size != max_node || anyDuplicated(traversal_order) ||
            !setequal(traversal_order, seq_len(max_node))) {
        stop('Invalid or disconnected phylo topology while computing clade signatures.')
    }
    for (node_num in rev(traversal_order)) {
        if (node_num <= num_tip) {
            next
        }
        children = children_by_parent[[as.character(node_num)]]
        if (is.null(children) || !length(children) || any(vapply(
                tip_cache[children], is.null, logical(1)))) {
            stop('Invalid phylo topology while computing clade signatures.')
        }
        tip_cache[[node_num]] = sort(unique(unlist(tip_cache[children], use.names=FALSE)))
    }
    names(tip_cache) = as.character(seq_len(max_node))
    tip_cache
}

.get_node_tip_signatures = function(phy) {
    tip_sets = .get_node_tip_sets(phy)
    signatures = vapply(tip_sets, .encode_tip_signature, character(1))
    names(signatures) = names(tip_sets)
    signatures
}

collapse_short_branches = function(tree, tol=1e-8) {
    tol = suppressWarnings(as.numeric(tol))
    if (length(tol) != 1L || is.na(tol) || !is.finite(tol) || tol < 0) {
        stop('tol must be a single finite non-negative numeric value.')
    }
    if (is.null(tree[['edge.length']]) || length(tree[['edge.length']]) != nrow(tree[['edge']])) {
        stop('tree must contain one branch length per edge in collapse_short_branches().')
    }
    if (anyDuplicated(tree[['tip.label']])) {
        stop('tree must have unique tip labels in collapse_short_branches().')
    }
    out_tree = tree
    edge_lengths = as.numeric(out_tree$edge.length)
    num_tip = length(out_tree[['tip.label']])
    is_internal_edge = out_tree[['edge']][,2] > num_tip
    is_short_edge = is_internal_edge & !is.na(edge_lengths) & (abs(edge_lengths) < tol)
    if (any(is.na(edge_lengths))) {
        cat('NA edge lengths were ignored when searching extremely short internal branches.\n')
    }
    if (!any(is_short_edge)) {
        cat('No extremely short internal branch was detected. tol =', tol, '\n')
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
    cat(
        'Extremely short internal branches ( n =', sum(is_short_edge),
        ') were collapsed. tol =', tol, '\n'
    )
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

force_ultrametric = function(tree, stop_if_larger_change=0.01) {
    out_tree = tree
    stop_if_larger_change = suppressWarnings(as.numeric(stop_if_larger_change))
    if (length(stop_if_larger_change) != 1L || is.na(stop_if_larger_change) ||
            !is.finite(stop_if_larger_change) || stop_if_larger_change < 0) {
        stop('stop_if_larger_change must be a single finite non-negative numeric value.')
    }
    if (is.null(out_tree[['edge.length']]) || length(out_tree[['edge.length']]) != nrow(out_tree[['edge']])) {
        stop('tree must contain one branch length per edge in force_ultrametric().')
    }
    if (any(is.na(out_tree[['edge.length']]))) {
        stop('tree contains NA edge lengths. Please resolve NA values before force_ultrametric().')
    }
    if (ape::is.ultrametric(out_tree)) {
        cat('The tree is ultrametric.\n')
    } else {
        cat('The tree is not ultrametric. Adjusting the branch length.\n')
        edge_length_before = out_tree[['edge.length']]
        out_tree = ape::chronoMPL(out_tree)
        edge_length_after = out_tree[['edge.length']]
        sum_adjustment = sum(abs(edge_length_after-edge_length_before))
        cat('Total branch length difference between before- and after-adjustment:', sum_adjustment, '\n')
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

get_single_branch_tree = function(name, dist) {
    if (length(name) != 1 || is.na(name) || trimws(as.character(name)) == '') {
        stop('name must be a single non-empty tip label in get_single_branch_tree().')
    }
    dist_num = suppressWarnings(as.numeric(dist))
    if (length(dist_num) != 1 || is.na(dist_num)) {
        stop('dist must be a single numeric branch length in get_single_branch_tree().')
    }
    phy = list(
      edge = matrix(c(2,1),1,2),
      tip.label = as.character(name),
      edge.length = dist_num,
      Nnode = 1
    )
    class(phy) = "phylo"
    return(phy)
}

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
            stop(
                'table2phylo() requires a binary tree; parent branch_id ',
                parent_id, ' has ', length(children), ' children.'
            )
        }
        if (length(children) == 1L) {
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

table2phylo = function(df, name_col, dist_col) {
    name_col = .normalize_single_string_arg(name_col, 'name_col', allow_empty=FALSE)
    dist_col = .normalize_single_string_arg(dist_col, 'dist_col', allow_empty=FALSE)
    if (!is.data.frame(df) || nrow(df) == 0L) {
        stop('df must be a non-empty data.frame in table2phylo().')
    }
    df_local = df
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
    class(phy) = 'phylo'
    phy = ape::reorder.phylo(phy, order='cladewise')
    if (length(phy[['tip.label']]) > 1L) {
        phy = ape::ladderize(phy, right=TRUE)
    }
    return(phy)
}

.phylo2table_get_descendant_tip_nums = function(phy, node_num, num_tip) {
    if (node_num <= num_tip) {
        return(node_num)
    }
    get_descendent_num(phy, node_num, leaf_only=TRUE)
}

.phylo2table_make_clade_signature = function(tip_nums, tip_rank_by_num, num_tip) {
    tip_ranks = as.integer(tip_rank_by_num[as.character(tip_nums)])
    if (any(is.na(tip_ranks))) {
        stop('phylo2table() could not rank descendant tips for a clade.')
    }
    clade_bits = rep('0', num_tip)
    clade_bits[num_tip - tip_ranks + 1L] = '1'
    paste0(clade_bits, collapse='')
}

.phylo2table_make_branch_id_map = function(phy, node_nums) {
    num_tip = length(phy[['tip.label']])
    tip_order = order(as.character(phy[['tip.label']]))
    tip_rank_by_num = integer(num_tip)
    tip_rank_by_num[tip_order] = seq_len(num_tip)
    names(tip_rank_by_num) = as.character(seq_len(num_tip))

    clade_signatures = vapply(
        X=node_nums,
        FUN=function(node_num) {
            tip_nums = .phylo2table_get_descendant_tip_nums(
                phy=phy,
                node_num=node_num,
                num_tip=num_tip
            )
            .phylo2table_make_clade_signature(
                tip_nums=tip_nums,
                tip_rank_by_num=tip_rank_by_num,
                num_tip=num_tip
            )
        },
        FUN.VALUE=character(1)
    )
    clade_order = order(clade_signatures, node_nums)
    branch_ids = integer(length(node_nums))
    branch_ids[clade_order] = seq_along(node_nums) - 1L
    names(branch_ids) = as.character(node_nums)
    branch_ids
}

phylo2table = function(phy, name_col='label', dist_col='dist') {
    if (!inherits(phy, 'phylo')) {
        stop('phy must be an object of class "phylo".')
    }
    name_col = .normalize_single_string_arg(
        value=name_col,
        arg_name='name_col',
        allow_empty=FALSE
    )
    dist_col = .normalize_single_string_arg(
        value=dist_col,
        arg_name='dist_col',
        allow_empty=FALSE
    )
    reserved_cols = c('branch_id', 'parent', 'sister')
    if (name_col %in% reserved_cols) {
        stop('name_col must not be one of: ', paste(reserved_cols, collapse=', '))
    }
    if (dist_col %in% c(reserved_cols, name_col)) {
        stop('dist_col must not duplicate branch_id, parent, sister, or name_col.')
    }
    if (is.null(phy[['edge']]) || nrow(phy[['edge']]) == 0) {
        stop('phylo2table() requires a tree with at least one edge.')
    }
    if (!ape::is.rooted(phy)) {
        stop('phylo2table() requires a rooted tree.')
    }
    if (anyDuplicated(phy[['tip.label']])) {
        stop('phylo2table() requires unique tip labels.')
    }
    if (is.null(phy[['edge.length']]) || length(phy[['edge.length']]) != nrow(phy[['edge']])) {
        stop('phylo2table() requires branch lengths for all edges.')
    }

    num_tip = length(phy[['tip.label']])
    num_internal = as.integer(phy[['Nnode']])
    internal_nodes = seq.int(num_tip + 1L, num_tip + num_internal)
    max_node = max(phy[['edge']])
    node_name_by_num = rep(NA_character_, max_node)
    node_name_by_num[seq_len(num_tip)] = as.character(phy[['tip.label']])

    internal_node_names = phy[['node.label']]
    if (is.null(internal_node_names)) {
        internal_node_names = rep(NA_character_, num_internal)
    } else {
        internal_node_names = as.character(internal_node_names)
        if (length(internal_node_names) < num_internal) {
            internal_node_names = c(internal_node_names, rep(NA_character_, num_internal - length(internal_node_names)))
        } else if (length(internal_node_names) > num_internal) {
            internal_node_names = internal_node_names[seq_len(num_internal)]
        }
    }
    missing_node_names = is.na(internal_node_names) | (trimws(internal_node_names) == '')
    if (any(missing_node_names)) {
        used_names = c(
            as.character(phy[['tip.label']]),
            internal_node_names[!missing_node_names]
        )
        counter = 0L
        for (missing_index in which(missing_node_names)) {
            candidate = paste0('n', counter)
            while (candidate %in% used_names) {
                counter = counter + 1L
                candidate = paste0('n', counter)
            }
            internal_node_names[[missing_index]] = candidate
            used_names = c(used_names, candidate)
            counter = counter + 1L
        }
    }
    node_name_by_num[internal_nodes] = internal_node_names

    node_nums = sort(unique(as.integer(c(phy[['edge']]))))
    node_names_for_table = node_name_by_num[node_nums]
    missing_node_names_for_table = is.na(node_names_for_table) | (trimws(node_names_for_table) == '')
    if (any(missing_node_names_for_table)) {
        stop('phylo2table() could not resolve labels for all nodes.')
    }
    root_num = get_root_num(phy)
    if (length(root_num) != 1) {
        stop('phylo2table() requires a tree with exactly one root.')
    }
    branch_id_by_node_num = .phylo2table_make_branch_id_map(
        phy=phy,
        node_nums=node_nums
    )
    child_nums_by_parent_num = split(phy[['edge']][,2], phy[['edge']][,1])
    parent_num_by_child_num = stats::setNames(phy[['edge']][,1], phy[['edge']][,2])
    edge_length_by_child_num = stats::setNames(phy[['edge.length']], phy[['edge']][,2])
    ordered_nodes = c(root_num, as.integer(phy[['edge']][,2]))

    table_rows = lapply(ordered_nodes, function(node_num) {
        if (node_num == root_num) {
            parent_id = -999L
            sister_id = -999L
            branch_dist = 0
        } else {
            parent_num = as.integer(parent_num_by_child_num[as.character(node_num)])
            sister_nums = setdiff(as.integer(child_nums_by_parent_num[[as.character(parent_num)]]), node_num)
            if (length(sister_nums) == 0 && parent_num == root_num) {
                sister_id = -999L
            } else if (length(sister_nums) == 1) {
                sister_id = unname(branch_id_by_node_num[as.character(sister_nums)])
            } else {
                stop('phylo2table() currently requires a binary rooted tree.')
            }
            parent_id = unname(branch_id_by_node_num[as.character(parent_num)])
            branch_dist = as.numeric(edge_length_by_child_num[as.character(node_num)])
        }
        data.frame(
            branch_id=unname(branch_id_by_node_num[as.character(node_num)]),
            parent=parent_id,
            sister=sister_id,
            label=unname(node_name_by_num[node_num]),
            dist=branch_dist,
            stringsAsFactors=FALSE
        )
    })

    out_table = do.call(rbind, table_rows)
    rownames(out_table) = NULL
    colnames(out_table)[colnames(out_table) == 'label'] = name_col
    colnames(out_table)[colnames(out_table) == 'dist'] = dist_col
    return(out_table)
}

fill_node_labels = function(phy) {
    out_phy = phy
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
    cat('Filling', length(missing_index), 'node names.\n')
    counter = 0
    for (i in missing_index) {
        lab = paste0('n', counter)
        while (any(lab == c(out_phy[['tip.label']], out_phy[['node.label']]), na.rm=TRUE)) {
            counter = counter + 1L
            lab = paste0('n', counter)
        }
        out_phy[['node.label']][i] = lab
        counter = counter + 1
    }
    return(out_phy)
}
