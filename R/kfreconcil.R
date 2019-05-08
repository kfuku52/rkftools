# Title     : TODO
# Objective : TODO
# Created by: kef74yk
# Created on: 2019-01-02

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
    num_parallel = min(nslots, nrow(phy$edge))
    cluster = makeCluster(num_parallel, 'PSOCK', outfile='')
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

read_notung_parsable = function(file, mode='D') {
    options(stringsAsFactors=FALSE)
    cols = c('event', 'gn_node', 'lower_sp_node', 'upper_sp_node')
    if (mode=='D') {
        f = file(file,"r")
        event_lines = c()
        repeat {
             str = readLines(con=f,1)
             if (length(str)==0) {
                 break
             }
             event_lines = c(event_lines, str)
        }
        dup_positions = grep("^#D", event_lines)
        if (length(dup_positions)>1) {
            dup_lines = event_lines[dup_positions]
            dup_items = strsplit(dup_lines[2:length(dup_lines)], "\\s")
            event_items = strsplit(dup_lines[2:length(dup_lines)], "\\s")
            df = data.frame(t(data.frame(dup_items)))
            rownames(df) = NULL
            colnames(df) = cols
            df$event = 'D'
        } else {
            df = data.frame(matrix(NA, 0, length(cols)))
            colnames(df) = cols
        }
    } else {
        cat('mode', mode, 'is not supported.')
        df = data.frame(matrix(NA, 0, length(cols)))
        colnames(df) = cols
    }
    close(f)
    return(df)
}