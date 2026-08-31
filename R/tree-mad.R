
#' Root a tree at a scored MAD edge
#'
#' @param t A `phylo` tree.
#' @param madr A single edge index.
#' @param rho Root position proportions, one per edge.
#' @return A list containing Newick text, the rooted tree, and clock CV.
#' @export
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


.calc_mad_branch_stats = function(br, t, dis, sdis, disbr, nodeids, otuids, npairs, notu, nbranch) {
    dij = t$edge.length[br]
    if (dij == 0) {
        return(c(rho=NA_real_, bad=NA_real_))
    }

    rbca = numeric(npairs)
    i = t$edge[br,1]
    j = t$edge[br,2]
    sp = dis[otuids,i] < dis[otuids,j]
    if (!any(sp) || all(sp)) return(c(rho=NA_real_, bad=NA_real_))
    dbc = matrix(sdis[sp,!sp], nrow=sum(sp), ncol=sum(!sp))
    dbi = matrix(dis[otuids[sp],i], nrow=sum(sp), ncol=sum(!sp))
    rho_br = sum((dbc - 2 * dbi) * dbc^-2) / (2 * dij * sum(dbc^-2))
    rho_br = min(max(0, rho_br), 1)

    dab = dbi + (dij * rho_br)
    ndab = length(dab)
    rbca[seq_len(ndab)] = as.vector(2 * dab / dbc - 1)

    bcsp = rbind(sp, !sp)
    ij = c(i, j)
    counter = ndab
    i2p = matrix(FALSE, nrow=nbranch + 1, ncol=length(t$tip.label))
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
    if (counter != npairs) {
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
    .validate_phylo_input(t, context='Input tree')
    if (ape::is.rooted(t)) {
        t = ape::unroot(t)
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
    nbranch = nrow(t$edge)
    dis = ape::dist.nodes(t)
    tip_ids = seq_along(t$tip.label)
    # Zero-distance tips contribute one representative to the MAD objective,
    # as in the historical collapsed calculation. Keep the ORIGINAL topology
    # and edge numbering in every public result instead of using placeholders.
    covered = logical(length(tip_ids))
    representatives = logical(length(tip_ids))
    for (tip in rev(tip_ids)) {
        if (!covered[[tip]]) {
            representatives[[tip]] = TRUE
            covered[dis[tip,tip_ids] == 0] = TRUE
        }
    }
    otuids = which(representatives)
    notu = length(otuids)
    if (notu < 2L) {
        stop('MAD requires at least two tips with positive pairwise distances.')
    }
    sdis = dis[otuids,otuids,drop=FALSE]
    t2 = t
    t2$edge.length = rep(1, nbranch)
    disbr = ape::dist.nodes(t2)
    nodeids = seq_len(nbranch + 1)
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


.run_mad_with_tree = function(t, output_mode, ncpu, use_parallel) {
    scores = .compute_mad_scores(t=t, ncpu=ncpu, use_parallel=use_parallel)
    .format_mad_result(t=t, rho=scores[['rho']], bad=scores[['bad']], output_mode=output_mode)
}


#' Root a tree using minimal ancestor deviation
#'
#' Multifurcations are scored directly and are not resolved into random binary
#' trees. Zero-distance tips share one representative in the scoring objective;
#' all original tips remain in returned trees. In `full` and `custom` results,
#' root indices, deviations, and proportions refer to the returned unrooted
#' tree's edge rows. Clock CV uses all tips of each returned rooted tree.
#'
#' @param unrooted_newick A Newick string or `phylo` tree.
#' @param output_mode One of `"newick"`, `"stats"`, `"full"`, or `"custom"`.
#' @return Newick text or a list whose detail depends on `output_mode`.
#' @export
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
        t=t, output_mode=mode, ncpu=1, use_parallel=FALSE
    ))
}


#' Root a tree using parallel minimal ancestor deviation
#'
#' Multifurcations are scored directly and are not resolved into random binary
#' trees. Zero-distance tips share one representative in the scoring objective;
#' all original tips remain in returned trees. In `full` and `custom` results,
#' root indices, deviations, and proportions refer to the returned unrooted
#' tree's edge rows. Clock CV uses all tips of each returned rooted tree.
#'
#' @param unrooted_newick A Newick string or `phylo` tree.
#' @param output_mode One of `"newick"`, `"stats"`, `"full"`, or `"custom"`.
#' @param ncpu Requested worker count. Automatic parallelism is capped by
#'   `options("rkftools.max_cores")` and skipped for small trees.
#' @return Newick text or a list whose detail depends on `output_mode`.
#' @export
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
        t=t, output_mode=mode, ncpu=num_parallel, use_parallel=TRUE
    ))
}
