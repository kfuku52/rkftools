# Title     : TODO
# Objective : TODO
# Created by: kf
# Created on: 3/22/18


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
    stopifnot(ape::is.binary(tree))
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
    stopifnot(ape::is.ultrametric(phy))
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

get_phy2_root_in_phy1 = function(phy1, phy2, nslots, mode=c("node_num", "index")) {
    mode_name = match.arg(mode)
    num_slots = suppressWarnings(as.integer(nslots))
    if (!length(num_slots)) {
        num_slots = 1L
    } else {
        num_slots = num_slots[[1]]
    }
    if (is.na(num_slots) || num_slots < 1) {
        num_slots = 1L
    }
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
        cl = parallel::makeCluster(num_slots)
        on.exit(parallel::stopCluster(cl), add=TRUE)
        is_match = unlist(parallel::parLapply(
            cl=cl, X=edge_indices, fun=is_matching_root,
            phy_obj=phy1, ref_split=ref_root_split, all_tips=all_tips_sorted
        ), use.names=FALSE)
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
    if (length(madr) != 1 || is.na(madr) || madr < 1 || madr > nrow(t$edge)) {
        stop('madr must be a single edge index in [1, ', nrow(t$edge), '] in get_rooted_newick().')
    }
    if (length(rho) != nrow(t$edge)) {
        stop('rho must have length ', nrow(t$edge), ' in get_rooted_newick().')
    }
    if (is.na(rho[madr])) {
        stop('rho[madr] is NA in get_rooted_newick().')
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
    jj = sort(bad, index.return=TRUE)
    tf = bad == jj$x[1]
    tf[is.na(tf)] = FALSE
    nroots = sum(tf)
    if (nroots > 1) {
        warning("More than one possible root position. Multiple newick strings printed")
    }
    madr = which(tf)
    rai = jj$x[1] / jj$x[2]
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
        res[[1]] = sub(vv[1], vv[2], res[[1]])
    } else {
        res = sub(vv[1], vv[2], res)
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
    if (ape::is.rooted(t)) {
        t = ape::unroot(t)
    }
    if (!ape::is.binary(t)) {
        warning("Input tree is not binary! Internal multifurcations will be converted to branches of length zero and identical OTUs will be collapsed!")
        t = ape::multi2di(t)
    }
    has_negative = (t$edge.length < 0)
    if (any(has_negative)) {
        warning("Input tree contains negative branch lengths. They will be converted to zeros!")
        t$edge.length[has_negative] = 0
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
    }

    t = .prepare_mad_tree(unrooted_newick)
    num_parallel = .resolve_parallel_cores(
        requested=ncpu,
        max_tasks=nrow(t$edge),
        auto_when_missing=TRUE
    )
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
    for (t in seq_len(num_int_to)) {
        to_node_num = num_tip_to + t
        to_clade = ape::extract.clade(phy=out_phy_to, node=to_node_num, root.edge=0, interactive=FALSE)
        to_leaves = to_clade$tip.label
        for (f in seq_len(num_int_from)) {
            from_node_num = num_tip_from + f
            from_clade = ape::extract.clade(phy=phy_from, node=from_node_num, root.edge=0, interactive=FALSE)
            from_leaves = from_clade$tip.label
            if (setequal(to_leaves, from_leaves)) {
                out_phy_to$node.label[t] = as.character(phy_from$node.label[f])
                break
            }
        }
    }
    return(out_phy_to)
}

get_species_name = function(a) {
    out = sub('_',' ', a)
    out = sub('_.*','', out)
    return(out)
}

get_species_names = function(phy, sep='_') {
    sep = .normalize_single_string_arg(
        value=sep,
        arg_name='sep',
        allow_empty=FALSE
    )
    split_names = strsplit(phy[['tip.label']], sep, fixed=TRUE)
    species_names = character(length(split_names))
    for (i in seq_along(split_names)) {
        sn = split_names[[i]]
        if (length(sn) >= 2) {
            species_names[i] = paste0(sn[1], sep, sn[2])
        } else {
            warning(
                'leaf name could not be interpreted as genus', sep, 'species: ',
                phy[['tip.label']][i]
            )
            species_names[i] = NA_character_
        }
    }
    return(species_names)
}

leaf2species = function(leaf_names, use_underbar=FALSE) {
    use_underbar = .normalize_single_logical_arg(
        value=use_underbar,
        arg_name='use_underbar'
    )
    split = strsplit(leaf_names, '_')
    species_names = rep(NA_character_, length(split))
    for (i in seq_along(split)) {
        if (length(split[[i]])>=3) {
            species_names[i] = paste(split[[i]][[1]], split[[i]][[2]])
        } else {
            warning('leaf name could not be interpreted as genus_species_gene: ', leaf_names[i], '\n')
        }
    }
    if (use_underbar) {
        species_names = gsub(' ', '_', species_names)
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
    stopifnot(all(sort(phy1$tip.label)==sort(phy2$tip.label)))
    phy1_leaves = sort(get_tip_labels(phy1, phy1_node))
    phy2_leaves = sort(get_tip_labels(phy2, phy2_node))
    is_same_leaves = all(phy1_leaves==phy2_leaves)
    return(is_same_leaves)
}

multi2bi_node_number_transfer = function(multifurcated_tree, bifurcated_tree) {
    mtree = multifurcated_tree
    btree = bifurcated_tree
    stopifnot(all(mtree[['tip.label']]==btree[['tip.label']]))
    stopifnot(is_same_root(mtree, btree))
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

collapse_short_branches = function(tree, tol=1e-8) {
    out_tree = tree
    edge_lengths = out_tree$edge.length
    is_short_edge = (!is.na(edge_lengths)) & (abs(edge_lengths) < tol)
    if (any(is.na(edge_lengths))) {
        cat('NA edge lengths were ignored when searching extremely short branches.\n')
    }
    if (any(is_short_edge)) {
        cat('Extremely short branches ( n =', sum(is_short_edge), ') were collapsed. tol =', tol, '\n')
        out_tree = ape::di2multi(out_tree, tol=tol)
    } else {
        cat('No extremely short branch was detected. tol =', tol, '\n')
    }
    return(out_tree)
}

force_ultrametric = function(tree, stop_if_larger_change=0.01) {
    out_tree = tree
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
        stopifnot(sum_adjustment<(sum(out_tree[['edge.length']]) * stop_if_larger_change))
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
    suppressWarnings(as.numeric(values))
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

    branch_ids_chr = as.character(branch_ids)
    referenced_parent_ids = as.character(df[[parent_col]][!is_sentinel_parent])
    root_rows = which(!(branch_ids_chr %in% referenced_parent_ids))
    if (length(root_rows) == 1) {
        return(branch_ids[[root_rows]])
    }
    if (length(root_rows) == 0) {
        stop(
            'Unable to infer root in table2phylo(): no sentinel-parent row and no unreferenced parent candidate.'
        )
    }
    root_ids = unique(as.character(branch_ids[root_rows]))
    stop(
        'Ambiguous root candidate(s) in table2phylo(): multiple unreferenced parent candidates found for branch_id: ',
        paste(root_ids, collapse=', ')
    )
}

.table2phylo_make_lookup = function(df, columns) {
    node_ids = as.character(df[['branch_id']])
    if (anyDuplicated(node_ids)) {
        stop('Duplicate branch_id values detected in table2phylo().')
    }
    lookup = list()
    for (column_name in columns) {
        values = df[[column_name]]
        if (is.factor(values)) {
            values = as.character(values)
        }
        names(values) = node_ids
        lookup[[column_name]] = values
    }
    lookup
}

.table2phylo_id2value = function(lookup, node_id, column_name) {
    if (!(column_name %in% names(lookup))) {
        return(NA)
    }
    values = lookup[[column_name]][as.character(node_id)]
    if (length(values) == 0) {
        return(NA)
    }
    unname(values[[1]])
}

.table2phylo_add_branch = function(phy, nni, lookup, name_col, dist_col, id2value) {
    out_phy = phy
    nni_name = id2value(lookup, nni, name_col)
    nni_dist = id2value(lookup, nni, dist_col)
    if (is.na(nni_name) || trimws(as.character(nni_name)) == '' || is.na(nni_dist)) {
        stop('Missing branch metadata for node id: ', nni)
    }

    parent_id = id2value(lookup, nni, 'parent')
    if (.table2phylo_is_parent_sentinel(parent_id)) {
        parent_num = integer(0)
    } else {
        parent_name = id2value(lookup, parent_id, name_col)
        parent_num = get_node_num_by_name(out_phy, parent_name)
    }
    sister_id = id2value(lookup, nni, 'sister')
    sister_name = id2value(lookup, sister_id, name_col)
    sister_dist = id2value(lookup, sister_id, dist_col)
    sister_num = get_node_num_by_name(out_phy, sister_name)
    if (length(parent_num) > 1) {
        stop('Ambiguous parent mapping for node id: ', nni)
    }
    if (length(parent_num)==0 && length(sister_num)!=1) {
        stop('Cannot resolve unique sister mapping for node id: ', nni)
    }

    branch = get_single_branch_tree(nni_name, nni_dist)
    if (length(parent_num)==0) {
        if (is.na(sister_dist)) {
            sister_dist = 1e-8
        }
        sister_dist = ifelse(sister_dist<1e-8, 1e-8, sister_dist)
        out_phy = ape::bind.tree(out_phy, branch, where=sister_num, position=sister_dist)
    } else {
        out_phy = ape::bind.tree(out_phy, branch, where=parent_num, position=0)
    }
    return(out_phy)
}

.table2phylo_build_edges = function(df, lookup, phy, root_id, name_col, dist_col, id2value, max_iter) {
    out_phy = phy
    next_node_ids = sort(df[(df$parent==root_id),'branch_id'])
    iter = 0L
    while (length(next_node_ids)>0) {
        iter = iter + 1L
        if (iter > max_iter) {
            stop('Exceeded iteration limit while building edges in table2phylo().')
        }
        for (nni in next_node_ids) {
            out_phy = .table2phylo_add_branch(
                phy=out_phy, nni=nni, lookup=lookup,
                name_col=name_col, dist_col=dist_col, id2value=id2value
            )
        }
        next_node_ids = sort(unique(df[(df$parent %in% next_node_ids),'branch_id']))
    }
    return(out_phy)
}

.table2phylo_assign_node_labels = function(df, lookup, phy, name_col, id2value, max_iter) {
    out_phy = phy
    num_leaf = length(out_phy[['tip.label']])
    num_intnode = as.integer(out_phy[['Nnode']])
    if (length(num_intnode) != 1 || is.na(num_intnode) || num_intnode < 0) {
        num_intnode = max(out_phy[['edge']]) - num_leaf
    }
    out_phy$node.label = rep('placeholder', num_intnode)
    node_ids = seq_len(num_leaf + num_intnode)

    next_node_ids = sort(df[(df[[name_col]] %in% out_phy[['tip.label']]), 'branch_id'])
    iter = 0L
    while (!(length(next_node_ids)==1 && next_node_ids[1] < 0)) {
        iter = iter + 1L
        if (iter > max_iter) {
            stop('Exceeded iteration limit while assigning node labels in table2phylo().')
        }
        tmp_next_node_ids = integer(0)
        for (nni in next_node_ids) {
            if (nni>=0) {
                parent_id = id2value(lookup, nni, 'parent')
                if (.table2phylo_is_parent_sentinel(parent_id)) {
                    next
                }
                nni_name = id2value(lookup, nni, name_col)
                nni_num = node_ids[c(out_phy[['tip.label']], out_phy$node.label)==nni_name]
                if (length(nni_num) != 1) {
                    stop('Ambiguous node label mapping for node id: ', nni)
                }
                parent_num = out_phy[['edge']][(out_phy[['edge']][,2]==nni_num),1]
                if (length(parent_num) != 1) {
                    stop('Ambiguous parent label mapping for node id: ', nni)
                }
                parent_label_index = parent_num - num_leaf
                if (parent_label_index < 1 || parent_label_index > length(out_phy$node.label)) {
                    stop('Parent label index out of range for node id: ', nni)
                }
                parent_name = id2value(lookup, parent_id, name_col)
                out_phy$node.label[parent_label_index] = parent_name
                tmp_next_node_ids = c(tmp_next_node_ids, parent_id)
            }
        }
        next_node_ids = sort(unique(tmp_next_node_ids))
        if (length(next_node_ids)==0) {
            next_node_ids = -999L
        }
    }

    if (sum(out_phy$node.label=='placeholder', na.rm=TRUE)>1) {
        warning('Node label "placeholder" appeared more than once.')
    }
    return(out_phy)
}

table2phylo = function(df, name_col, dist_col) {
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
    df_local[is_root_row, 'sister'] = -999
    df_local[is_root_row, 'parent'] = -999
    lookup = .table2phylo_make_lookup(
        df=df_local,
        columns=unique(c(name_col, dist_col, 'parent', 'sister'))
    )

    root_name = .table2phylo_id2value(lookup, root_id, name_col)
    root_dist = .table2phylo_id2value(lookup, root_id, dist_col)
    if (is.na(root_name) || trimws(as.character(root_name)) == '') {
        stop('Failed to resolve root metadata in table2phylo().')
    }
    if (is.na(root_dist)) {
        root_dist = 0
    }
    phy = get_single_branch_tree(root_name, root_dist)

    max_iter = max(10L, nrow(df_local) * 4L)
    phy = .table2phylo_build_edges(
        df=df_local, lookup=lookup, phy=phy, root_id=root_id,
        name_col=name_col, dist_col=dist_col,
        id2value=.table2phylo_id2value,
        max_iter=max_iter
    )

    phy = remove_redundant_root_edge(phy)
    if (length(phy[['tip.label']]) > 1) {
        phy = ape::ladderize(phy, right=TRUE)
    }
    phy = .table2phylo_assign_node_labels(
        df=df_local, lookup=lookup, phy=phy, name_col=name_col,
        id2value=.table2phylo_id2value,
        max_iter=max_iter
    )
    return(phy)
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
        while (any(lab == out_phy[['node.label']], na.rm=TRUE)) {
            counter = counter + 1L
            lab = paste0('n', counter)
        }
        out_phy[['node.label']][i] = lab
        counter = counter + 1
    }
    return(out_phy)
}
