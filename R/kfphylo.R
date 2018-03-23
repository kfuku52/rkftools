# Title     : TODO
# Objective : TODO
# Created by: kf
# Created on: 3/22/18


get_node_num_by_name = function(phy, node_name) {
    node_names = c(phy$tip.label, phy$node.label)
    node_nums = 1:length(node_names)
    node_num = node_nums[node_names==node_name]
    return(node_num)
}

get_node_name_by_num = function(phy, node_num) {
    node_names = c(phy$tip.label, phy$node.label)
    node_nums = 1:length(node_names)
    node_name = node_names[node_nums==node_num]
    return(node_name)
}

get_root_num = function(phy) {
    root_num = setdiff(phy$edge[,1], phy$edge[,2])
    return(root_num)
}

get_children_num = function(phy, node_num) {
    children_num = phy$edge[(phy$edge[,1]==node_num),2]
    return(children_num)
}

get_tip_labels = function(phy, node_num, out=NULL) {
    if (is.null(out)) {
        out = list()
        out[['nn']] = c(node_num)
        out[['label']] = c()
        node_num = NULL
    }
    if (length(out[['nn']])==0) {
        return(out[['label']])
    }
    if (out[['nn']][1] <= length(phy$tip.label)) {
        out[['label']] = c(out[['label']], phy$tip.label[out[['nn']][1]])
    } else {
        out[['nn']] = c(out[['nn']], get_children_num(phy, out[['nn']][1]))
    }
    out[['nn']] = out[['nn']][out[['nn']]!=out[['nn']][1]]
    Recall(phy=phy, node_num=NULL, out=out)
}

get_duplication_confidence_score = function(phy, node_num) {
    children_num = get_children_num(phy, node_num)
    if (length(children_num)==2) {
        child_leaves = vector(mode='list', length(children_num))
        for (j in 1:length(children_num)) {
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

get_species_overlap_score = function(phy, dc_cutoff=0) {
    # this function assumes that leaf names are: GENUS_SPECIES_GENEID (e.g. Bos_taurus_AF492351.1)
    internal_nodes = (length(phy$tip.label)+1):(max(phy$edge))
    dc_scores = vector(mode='numeric', length(internal_nodes))
    i = 1
    for (int_node in internal_nodes) {
        dc_scores[i] = get_duplication_confidence_score(phy, int_node)
        i = i + 1
    }
    species_overlap_score = sum(dc_scores > dc_cutoff)
    return(species_overlap_score)
}

get_root_position_dependent_species_overlap_scores = function(phy, nslots) {
    species_overlap_scores = vector(mode="numeric", nrow(phy$edge))
    exp_funs = c('reroot', "get_species_overlap_score", "get_duplication_confidence_score",
                "get_children_num", "get_tip_labels")
    cluster = makeCluster(nslots, 'PSOCK', outfile='')
    registerDoParallel(cluster)
    so_score = foreach (i = 1:nrow(phy$edge), .combine=c, .export=exp_funs) %dopar% {
        rt = reroot(tree=phy, node.number=phy$edge[i,2])
        so_score = get_species_overlap_score(phy=rt, dc_cutoff=0)
        names(so_score) = i
        so_score
    }
    stopCluster(cluster)
    species_overlap_scores = so_score[order(as.integer(names(so_score)))]
    return(species_overlap_scores)
}

is_same_root = function(phy1, phy2) {
    if (! is.rooted(phy1)) {
        stop('phy1 is unrooted.')
    }
    if (! is.rooted(phy2)) {
        stop('phy2 is unrooted.')
    }
    if (! identical(sort(phy1$tip.label), sort(phy2$tip.label))) {
        stop('phy1 and phy2 have different sets of leaves.')
    }
    phys = list(phy1, phy2)
    leaf_names = vector(mode='list', length(phys))
    for (i in 1:length(phys)) {
        children = get_children_num(phys[[i]], get_root_num(phys[[i]]))
        leaf_names[[i]] = vector(mode='list', length(children))
        leaf_names[[i]][[1]] = get_tip_labels(phys[[i]], children[1])
        leaf_names[[i]][[2]] = get_tip_labels(phys[[i]], children[2])
    }
    if (identical(sort(leaf_names[[1]][[1]]), sort(leaf_names[[2]][[1]]))) {
        return(TRUE)
    } else if (identical(sort(leaf_names[[1]][[1]]), sort(leaf_names[[2]][[2]]))) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}