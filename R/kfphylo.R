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
            df$event = 'dup'
        } else {
            df = data.frame(matrix(NA, 0, length(cols)))
            colnames(df) = cols
        }
    } else {
        cat('mode', mode, 'is not supported.')
        df = data.frame(matrix(NA, 0, length(cols)))
        colnames(df) = cols
    }
    return(df)
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
        species_names = c(
            species_names,
            paste(split[[i]][[1]], split[[i]][[2]])
        )
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

