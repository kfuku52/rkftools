test_that("MAD, root edges, and root label transfer", {
    withr::local_seed(20260831)
    tr = fixture_gene_tree()
    tr_unlabeled = fixture_trait_tree()
    rerooted = phytools::reroot(tree=tr, node.number=tr$edge[2,2])
    root_idx = get_phy2_root_in_phy1(tr, rerooted, nslots=1, mode="index")
    expect_true(!is.na(root_idx))
    root_idx_vec_slots = get_phy2_root_in_phy1(tr, rerooted, nslots=c(1L, 2L), mode="index")
    expect_true(!is.na(root_idx_vec_slots))
    root_idx_na_slots = get_phy2_root_in_phy1(tr, rerooted, nslots=NA_integer_, mode="index")
    expect_true(!is.na(root_idx_na_slots))

    mad_n1 = MAD_parallel(tr, output_mode="newick", ncpu=1)
    mad_auto = MAD_parallel(tr, output_mode="newick")
    expect_true(length(mad_n1) == length(mad_auto))
    mad_vector_ncpu = MAD_parallel(tr, output_mode="newick", ncpu=c(1L, 2L))
    expect_true(length(mad_vector_ncpu) == length(mad_n1))
    mad_serial = MAD(tr, output_mode="newick")
    expect_true(length(mad_serial) == length(mad_n1))
    tr_no_edge_length = tr
    tr_no_edge_length$edge.length = NULL
    mad_no_edge_length_err = tryCatch(
        {
            MAD(tr_no_edge_length, output_mode="newick")
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(mad_no_edge_length_err))
    expect_true(grepl("no branch lengths", conditionMessage(mad_no_edge_length_err), fixed=TRUE))
    mad_mode_na_err = tryCatch(
        {
            MAD(tr, output_mode=NA_character_)
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(mad_mode_na_err))
    expect_true(grepl("output_mode must be a single non-missing string", conditionMessage(mad_mode_na_err), fixed=TRUE))
    mad_mode_vec_err = tryCatch(
        {
            MAD_parallel(tr, output_mode=c("newick", "stats"))
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(mad_mode_vec_err))
    expect_true(grepl("output_mode must be a single non-missing string", conditionMessage(mad_mode_vec_err), fixed=TRUE))

    old_core_limit_env = Sys.getenv("_R_CHECK_LIMIT_CORES_", unset=NA_character_)
    Sys.setenv("_R_CHECK_LIMIT_CORES_"="TRUE")
    mad_limited = MAD_parallel(tr, output_mode="newick")
    expect_true(length(mad_limited) == length(mad_n1))
    if (is.na(old_core_limit_env)) {
        Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
    } else {
        Sys.setenv("_R_CHECK_LIMIT_CORES_"=old_core_limit_env)
    }

    outgroup_labels = get_outgroup(tr_unlabeled)
    expect_true(identical(as.character(outgroup_labels), "C"))
    outgroup_unrooted_err = tryCatch(
        {
            get_outgroup(ape::unroot(tr_unlabeled))
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(outgroup_unrooted_err))
    expect_true(grepl("requires a rooted tree", conditionMessage(outgroup_unrooted_err), fixed=TRUE))
    single_branch = get_single_branch_tree("A", 0.5)
    expect_true(inherits(single_branch, "phylo"))
    expect_true(identical(as.character(single_branch$tip.label), "A"))
    rooted_newick_ok = get_rooted_newick(tr_unlabeled, madr=1, rho=rep(0.5, nrow(tr_unlabeled$edge)))
    expect_true(length(rooted_newick_ok) == 3)
    rooted_newick_err = tryCatch(
        {
            get_rooted_newick(tr_unlabeled, madr=999, rho=rep(0.5, nrow(tr_unlabeled$edge)))
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(rooted_newick_err))
    expect_true(grepl("must be a single edge index", conditionMessage(rooted_newick_err), fixed=TRUE))
    single_branch_name_err = tryCatch(
        {
            get_single_branch_tree(c("A", "B"), 0.5)
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(single_branch_name_err))
    expect_true(grepl("single non-empty tip label", conditionMessage(single_branch_name_err), fixed=TRUE))

    tr_two_tip = ape::read.tree(text="(A:1,B:1);")
    tr_two_tip_rr = remove_redundant_root_edge(tr_two_tip)
    expect_true(inherits(tr_two_tip_rr, "phylo"))
    expect_true(setequal(tr_two_tip_rr$tip.label, c("A", "B")))
    expect_true(identical(as.integer(tr_two_tip_rr$Nnode), 1L))
    tr_unary_root = ape::read.tree(text="((A:1):1,B:1);")
    tr_unary_root_rr = remove_redundant_root_edge(tr_unary_root)
    expect_true(inherits(tr_unary_root_rr, "phylo"))
    expect_true(setequal(tr_unary_root_rr$tip.label, c("A", "B")))
    expect_true(identical(as.integer(tr_unary_root_rr$Nnode), 1L))
    tr_from_label = ape::read.tree(text="((A:1,B:1):1,C:1);")
    tr_from_label$node.label = c("X", "Y")
    tr_to_unlabeled = ape::read.tree(text="((A:1,B:1):1,C:1);")
    tr_to_unlabeled$node.label = NULL
    tr_transferred = transfer_node_labels(tr_from_label, tr_to_unlabeled)
    expect_true(identical(as.character(tr_transferred$node.label), c("X", "Y")))
    transfer_tip_mismatch_err = tryCatch(
        {
            transfer_node_labels(tr_from_label, ape::read.tree(text="((A:1,D:1):1,C:1);"))
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(transfer_tip_mismatch_err))
    expect_true(grepl("must contain the same tip labels", conditionMessage(transfer_tip_mismatch_err), fixed=TRUE))


    old_max_cores = getOption("rkftools.max_cores")
    options(rkftools.max_cores=1L)
    mad_auto_capped = MAD_parallel(tr, output_mode="newick")
    expect_true(identical(mad_n1, mad_auto_capped))
    options(rkftools.max_cores=old_max_cores)

    # MAD rejects ambiguous modes and trees without any positive distance.
    expect_error_contains(MAD(tr, output_mode="invalid"), "should be one of")
    all_zero_tree = ape::read.tree(text="((A:0,B:0):0,C:0);")
    expect_error_contains(MAD(all_zero_tree, output_mode="newick"), "no positive branch lengths")
})
