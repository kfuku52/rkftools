# Title     : TODO
# Objective : TODO
# Created by: kf
# Created on: 3/22/18


get_node_num_by_name = function(phy, node_name) {
    node_names = c(phy$tip.label, phy$node.label)
    node_nums = 1:length(node_names)
    node_num = node_nums[node_names %in% node_name]
    return(node_num)
}

get_node_name_by_num = function(phy, node_num) {
    node_names = c(phy$tip.label, phy$node.label)
    node_nums = 1:length(node_names)
    node_name = node_names[node_nums %in% node_num]
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

get_descendent_num = function(phy, node_num) {
    descendent_nums = c()
    children_nums = get_children_num(phy, node_num)
    current_size = length(children_nums)
    next_size = current_size + 1
    while(!all(is.na(children_nums))) {
        children_nums = children_nums[!is.na(children_nums)]
        descendent_nums = c(descendent_nums, children_nums)
        childrens = c()
        for (nn in children_nums) {
            childrens = c(childrens, get_children_num(phy, nn))
        }
        children_nums = childrens
    }
    descendent_nums = sort(descendent_nums)
    return(descendent_nums)
}

get_parent_num = function(phy, node_num) {
    parent_num = phy$edge[(phy$edge[,2]==node_num),1]
    return(parent_num)
}

get_sister_num = function(phy, node_num) {
    parent_num = phy$edge[(phy$edge[,2]==node_num),1]
    sibling_num = phy$edge[(phy$edge[,1]==parent_num),2]
    sister_num = sibling_num[sibling_num!=node_num]
    return(sister_num)
}

# alias
collapse_short_external_edges = function(tree, threshold=1e-6) {
    return(pad_short_edges(tree, threshold=threshold, external_only=TRUE))
}

pad_short_edges = function(tree, threshold=1e-6, external_only=FALSE) {
    stopifnot(ape::is.binary(tree))
    edge_idx = 1:nrow(tree$edge)
    is_target_edge = TRUE
    if (external_only) {
        is_target_edge = is_target_edge & (tree$edge[,2]<=length(tree$tip.label))
    }
    edge_lengths = tree[['edge.length']][is_target_edge]
    min_eel = min(edge_lengths)
    cat('Minimum edge length:', min_eel, '\n')
    is_short_eel = (is_target_edge)&(tree$edge.length<threshold)
    num_short_eel = sum(is_short_eel)
    cat('Number of short edges ( length <', threshold, '):', num_short_eel, '\n')
    if (num_short_eel>0) {
        short_eel_idx = edge_idx[is_short_eel]
        for (i in short_eel_idx) {
            if (tree$edge.length[i]<threshold) {
                shift_value = threshold - tree$edge.length[i]
                sister_node_num = get_sister_num(tree, tree$edge[i,2])
                sister_edge_idx = edge_idx[tree$edge[,2]==sister_node_num]
                root_num = get_root_num(tree)
                flag = TRUE
                flag_root = FALSE
                current_idx = i
                while (flag==TRUE) {
                    parent_node_num = tree$edge[current_idx,1]
                    parent_edge_idx = edge_idx[tree$edge[,2]==parent_node_num]
                    parent_edge_length = tree$edge.length[parent_edge_idx]
                    if (parent_node_num==root_num) {
                        flag = FALSE
                        flag_root = TRUE
                    } else if (parent_edge_length>=threshold+shift_value) {
                        flag = FALSE
                    } else {
                        current_idx = edge_idx[tree$edge[,2]==parent_node_num]
                    }
                }

                tree$edge.length[i] = tree$edge.length[i] +shift_value
                tree$edge.length[sister_edge_idx] = tree$edge.length[sister_edge_idx] + shift_value
                if (flag_root) {
                    cat('Adding branch length to subroot edges,', i, 'and', sister_edge_idx, '\n')
                } else {
                    cat('Transfering branch length from edge', parent_edge_idx, 'to', i, 'and', sister_edge_idx, '\n')
                    tree$edge.length[parent_edge_idx] = tree$edge.length[parent_edge_idx] - shift_value
                }
            }
        }
    }
    return(tree)
}

get_tip_labels = function(phy, node_num, out=NULL) {
    num_leaf = length(phy$tip.label)
    if (node_num > num_leaf) {
        subtree = ape::extract.clade(phy, node_num)
        tip_labels = subtree$tip.label
    } else {
        tip_labels = phy$tip.label[node_num]
    }
    return(tip_labels)
}

get_node_age = function(phy, node_num) {
    stopifnot(is.ultrametric(phy))
    age = 0
    current_node_num = node_num
    while (!is.na(current_node_num)) {
        descendent_node_num = get_children_num(phy, current_node_num)[1]
        if (!is.na(descendent_node_num)) {
            edge_length = phy$edge.length[(phy$edge[,1]==current_node_num)&(phy$edge[,2]==descendent_node_num)]
            age = age + edge_length
        }
        current_node_num = descendent_node_num
    }
    return(age)
}

get_outgroup = function(phy) {
    children_nums = get_children_num(phy, get_root_num(phy))
    outgroup_labels = phy$tip.label
    for (cn in children_nums) {
        og_labels = get_tip_labels(phy, cn)
        if (length(outgroup_labels)>length(og_labels)) {
            outgroup_labels = og_labels
        }
    }
    return(outgroup_labels)
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

get_phy2_root_in_phy1 = function(phy1, phy2, nslots, mode=c("node_num", "index")) {
    if (! is.rooted(phy2)) {
        stop('phy2 is unrooted.')
    }
    if (! identical(sort(phy1$tip.label), sort(phy2$tip.label))) {
        stop('phy1 and phy2 have different sets of leaves.')
    }
    for (i in 1:nrow(phy1$edge)) {
        rphy1 = reroot(tree=phy1, node.number=phy1$edge[i,2])
        if (is_same_root(rphy1, phy2)) {
            if (mode=="node_num") {
                root_pos = phy1$edge[i,2]
            } else if (mode=="index") {
                root_pos = i
            }
            break
        }
    }
    return(root_pos)
}

get_rooted_newick = function(t, madr, rho) {
    notu <- length(t$tip.label)
    dis <- dist.nodes(t)
    pp <- rho[madr]*t$edge.length[madr]
    nn <- t$edge[madr,]
    rt <- reroot(t,nn[2], pos = pp)
    rooted_newick <- write.tree(rt)
    dd <- dis[1:notu,nn]
    sp <- dd[,1]<dd[,2]
    otu2root <- vector(mode="numeric",notu)
    otu2root[sp] <- dd[sp,1] + pp
    otu2root[!sp] <- dd[!sp,1] - pp
    ccv <- 100*sd(otu2root)/mean(otu2root)
    return(list(rooted_newick, rt, ccv))
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
    if (!library('ape',logical.return = TRUE)){
        stop("'ape' package not found, please install it to run MAD")
    }
    if (!library('phytools',logical.return = TRUE)){
        stop("'phytools' package not found, please install it to run MAD")
    }
    t <- NA
    if(class(unrooted_newick)=="phylo"){
        t <- unrooted_newick
    } else {
        t <- read.tree(text=unrooted_newick)
    }
    if(is.rooted(t)){
        t<-unroot(t)
    }
    #t$node.label<-NULL #To allow parsing when identical OTUs are present
    if(!is.binary.tree(t)){
        warning("Input tree is not binary! Internal multifurcations will be converted to branches of length zero and identical OTUs will be collapsed!")
        t<-multi2di(t)
    }
    tf <- t$edge.length<0
    if(any(tf)){
        warning("Input tree contains negative branch lengths. They will be converted to zeros!")
        t$edge.length[tf]<-0
    }

    notu <- length(t$tip.label)
    nbranch <- dim(t$edge)[1]
    npairs <- notu*(notu-1)/2
    nodeids <- 1:(nbranch+1)
    otuids <- 1:notu
    dis <- dist.nodes(t) # phenetic distance. All nodes
    sdis <- dis[1:notu,1:notu] # phenetic distance. otus only

    #### Start recursion to collapse identical OTUs, if present.
    ii<-which(sdis==0,arr.ind=TRUE)
    k<-which(ii[,1]!=ii[,2])
    if(length(k)){
        r<-ii[k[1],1]
        c<-ii[k[1],2]
        vv<-c(paste('@#',t$tip.label[r],'@#',sep=""),paste('(',t$tip.label[r],':0,',t$tip.label[c],':0)',sep=""))
        st<-drop.tip(t,c)
        st$tip.label[st$tip.label==t$tip.label[r]]<-vv[1]
        res<-MAD(st,output_mode)
        if(is.list(res)){
            res[[1]]<-sub(vv[1],vv[2],res[[1]])
        } else{
            res<-sub(vv[1],vv[2],res)
        }
        return(res) #create the list 'res' to return the results
    }
    #### End of recursion

    gc()
    t2 <- t
    t2$edge.length <- rep(1,nbranch)
    disbr <- dist.nodes(t2) # split distance. All nodes
    sdisbr <- disbr[1:notu,1:notu] # split distance. otus only
    rho <- vector(mode = "numeric",length = nbranch) # Store position of the optimized root nodes (branch order as in the input tree)
    bad <- vector(mode = "numeric",length = nbranch) # Store branch ancestor deviations (branch order as in the input tree)
    i2p <- matrix(nrow = nbranch+1, ncol = notu)
    for (br in 1:nbranch){
        #collect the deviations associated with straddling otu pairs
        dij <- t$edge.length[br]
        if(dij==0){
            rho[br]<-NA
            bad[br]<-NA
            next
        }
        rbca <- numeric(npairs)
        i <- t$edge[br,1]
        j <- t$edge[br,2]
        sp <- dis[1:notu,i]<dis[1:notu,j] # otu split for 'br'
        dbc <- matrix(sdis[sp,!sp],nrow=sum(sp),ncol=sum(!sp))
        dbi <- replicate(dim(dbc)[2],dis[(1:notu)[sp],i])

        rho[br] <- sum((dbc-2*dbi)*dbc^-2)/(2*dij*sum(dbc^-2)) # optimized root node relative to 'i' node
        rho[br] <- min(max(0,rho[br]),1)
        dab <- dbi+(dij*rho[br])
        ndab <- length(dab)
        rbca[1:ndab] <- as.vector(2*dab/dbc-1)
        # collect the remaining deviations (non-traversing otus)
        bcsp <- rbind(sp,!sp)
        ij <- c(i,j)
        counter <- ndab
        for (w in c(1,2)){
            if(sum(bcsp[w,])>=2){
            disbrw <- disbr[,ij[w]]
            pairids <- otuids[bcsp[w,]]
            for (z in pairids){
                i2p[,z] <- disbr[z,]+disbrw==disbrw[z]
                }
                for (z in 1:(length(pairids)-1)){
                    p1 <- pairids[z]
                    disp1 <- dis[p1,]
                    pan <- nodeids[i2p[,p1]]
                    for (y in (z+1):length(pairids)){
                        p2 <- pairids[y]
                        pan1 <- pan[i2p[pan,p2]]
                        an <- pan1[which.max(disbrw[pan1])]
                        counter <- counter+1
                        rbca[counter] <- 2*disp1[an]/disp1[p2]-1
                    }
                }
            }
        }
        if(length(rbca)!=npairs){
            stop("Unexpected number of pairs.")
        }
        bad[br] <- sqrt(mean(rbca^2)) # branch ancestor deviation
    }
    # Select the branch with the minum ancestor deviation and calculate the root ambiguity index
    jj <- sort(bad,index.return = TRUE)
    tf<-bad==jj$x[1]
    tf[is.na(tf)]<-FALSE
    nroots <- sum(tf)
    if (nroots>1){
        warning("More than one possible root position. Multiple newick strings printed")
    }
    madr <- which(tf) # Index of the mad root branch(es)
    rai <- jj$x[1]/jj$x[2] # Root ambiguity index
    badr <- bad[tf] # Branch ancestor deviations value for the root(s)
    #Root the tree object, calculate the clock CV and retrieve the newick string
    rt <- vector(mode = "list",nroots) # Rooted tree object
    ccv <- vector(mode = "numeric",nroots) # Clock CV
    rooted_newick <- vector(mode = "character",nroots)

    for (i in 1:length(madr)){
        out = get_rooted_newick(t, madr[i], rho)
        rooted_newick[i] = out[[1]]
        rt[[i]] = out[[2]]
        ccv[i] = out[[3]]
    }
    rooted_newick<-sub(')Root;',');',rooted_newick)
    # Output the result(s)
    if(missing(output_mode)) {
        return(rooted_newick)
    } else {
        if(output_mode=='newick'){
            return(rooted_newick)
        } else if (output_mode=='stats'){ # Rooted newick and stats
            root_stats <- data.frame(ambiguity_index=rai,clock_cv=ccv,ancestor_deviation=badr,n_roots=nroots)
            return(list(rooted_newick,root_stats))
        } else if (output_mode=='full'){ #Rooted newick, stats, unrooted tree object, index of the branch root, ancestor deviations, rooted tree object
            root_stats <- data.frame(ambiguity_index=rai,clock_cv=ccv,ancestor_deviation=badr,n_roots=nroots)
            return(list(rooted_newick,root_stats,t,madr,bad,rt))
        } else if (output_mode=='custom'){ #Rooted newick, stats, unrooted tree object, index of the branch root, ancestor deviations, rooted tree object, rho
            root_stats <- data.frame(ambiguity_index=rai,clock_cv=ccv,ancestor_deviation=badr,n_roots=nroots)
            return(list(rooted_newick,root_stats,t,madr,bad,rt,rho))
        } else{
            return(rooted_newick)
        }
    }
}

transfer_node_labels = function(phy_from, phy_to) {
    for (t in 1:length(phy_to$node.label)) {
        to_node = phy_to$node.label[t]
        to_clade = extract.clade(phy=phy_to, node=to_node, root.edge = 0, interactive = FALSE)
        to_leaves = to_clade$tip.label
        for (f in 1:length(phy_from$node.label)) {
            from_node = phy_from$node.label[f]
            from_clade = extract.clade(phy=phy_from, node=from_node, root.edge = 0, interactive = FALSE)
            from_leaves = from_clade$tip.label
            if (setequal(to_leaves, from_leaves)) {
                phy_to$node.label[t] = from_node
                break
            }
        }
    }
    return(phy_to)
}

get_species_names = function(phy, sep='_') {
    split_names = strsplit(phy$tip.label, sep)
    species_names = c()
    for (sn in split_names) {
        species_names = c(species_names, paste0(sn[1], sep, sn[2]))
    }
    return(species_names)
}

leaf2species = function(leaf_names) {
    split = strsplit(leaf_names, '_')
    species_names = c()
    for (i in 1:length(split)) {
        if (length(split[[i]])>=3) {
            species_names = c(
                species_names,
                paste(split[[i]][[1]], split[[i]][[2]])
            )
        } else {
            warning('leaf name could not be interpreted as genus_species_gene: ', split[[i]], '\n')
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
    stopifnot(all(sort(phy1$tip.label)==sort(phy2$tip.label)))
    phy1_leaves = sort(rkftools::get_tip_labels(phy1, phy1_node))
    phy2_leaves = sort(rkftools::get_tip_labels(phy2, phy2_node))
    is_same_leaves = all(phy1_leaves==phy2_leaves)
    return(is_same_leaves)
}

multi2bi_node_number_transfer = function(multifurcated_tree, bifurcated_tree) {
    mtree = multifurcated_tree
    btree = bifurcated_tree
    stopifnot(all(mtree[['tip.label']]==btree[['tip.label']]))
    stopifnot(rkftools::is_same_root(mtree, btree))
    stopifnot(as.logical(ape::dist.topo(unroot(mtree), unroot(btree), method='PH85')))
    internal_node_counts = table(mtree[['edge']][,1])
    polytomy_parents = as.integer(names(internal_node_counts)[internal_node_counts > 2])
    cat('Polytomy parent nodes:', polytomy_parents, '\n')
    df = data.frame()
    for (mtree_pp in polytomy_parents) {
        mtree_pp_leaves = sort(rkftools::get_tip_labels(mtree, mtree_pp))
        for (btree_internal_node in btree$edge[,1]) {
            btree_in_leaves = sort(rkftools::get_tip_labels(btree, btree_internal_node))
            if (length(mtree_pp_leaves)==length(btree_in_leaves)) {
                if (all(mtree_pp_leaves==btree_in_leaves)) {
                    df = rbind(df, data.frame(mtree_node=mtree_pp, btree_node=btree_internal_node))
                    break
                }
            }
        }
    }
    return(df)
}

collapse_short_branches = function(tree, tol=1e-8) {
    if (any(abs(tree$edge.length) < tol)) {
        cat('Extremely short branches ( n =', sum(abs(tree$edge.length) < tol), ') were collapsed. tol =', tol, '\n')
        tree = ape::di2multi(tree, tol=tol)
    } else {
        cat('No extremely short branch was detected. tol =', tol, '\n')
    }
    return(tree)
}

force_ultrametric = function(tree, stop_if_larger_change=0.01) {
    if (ape::is.ultrametric(tree)) {
        cat('The tree is ultrametric.\n')
    } else {
        cat('The tree is not ultrametric. Adjusting the branch length.\n')
        edge_length_before = tree[['edge.length']]
        tree = ape::chronoMPL(tree)
        edge_length_after = tree[['edge.length']]
        sum_adjustment = sum(abs(edge_length_after-edge_length_before))
        cat('Total branch length difference between before- and after-adjustment:', sum_adjustment, '\n')
        stopifnot(sum_adjustment<(sum(tree[['edge.length']]) * stop_if_larger_change))
    }
    return(tree)
}
