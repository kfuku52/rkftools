test_that("short branches and topology transformations", {
    withr::local_seed(20260831)
    tr_unlabeled = fixture_trait_tree()
    tr_na_edge = tr_unlabeled
    tr_na_edge$edge.length[1] = NA_real_
    collapsed_na_edge = collapse_short_branches(tr_na_edge, tol=1e-8)
    expect_true(inherits(collapsed_na_edge, "phylo"))
    expect_true(is.na(collapsed_na_edge$edge.length[1]))
    padded_na_edge = pad_short_edges(tr_na_edge, threshold=1e-6)
    expect_true(inherits(padded_na_edge, "phylo"))
    expect_true(is.na(padded_na_edge$edge.length[1]))
    pad_external_only_na_err = tryCatch(
        {
            pad_short_edges(tr_unlabeled, threshold=1e-6, external_only=NA)
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(pad_external_only_na_err))
    expect_true(grepl("external_only must be a single non-missing logical value", conditionMessage(pad_external_only_na_err), fixed=TRUE))
    force_ultra_na_err = tryCatch(
        {
            force_ultrametric(tr_na_edge)
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(force_ultra_na_err))
    expect_true(grepl("contains NA edge lengths", conditionMessage(force_ultra_na_err), fixed=TRUE))

    short_ext = collapse_short_external_edges(ape::read.tree(text="((A:1e-9,B:1e-9):1,C:1);"), threshold=1e-6)
    expect_true(inherits(short_ext, "phylo"))
    deepest_node = get_deepest_node_num(tr_unlabeled, c(4L, 5L))
    expect_true(identical(as.integer(deepest_node), 4L))
    tr_reorder = ape::read.tree(text="((B:1,A:1):1,C:1);")
    expect_true(isTRUE(has_same_leaves(tr_unlabeled, get_root_num(tr_unlabeled), tr_reorder, get_root_num(tr_reorder))))
    expect_true(isTRUE(is_same_root(tr_unlabeled, tr_reorder)))
    expect_true(isFALSE(is_same_root(tr_unlabeled, ape::read.tree(text="((A:1,C:1):1,B:1);"))))
    tr_poly_root1 = ape::read.tree(text="(A:1,B:1,C:1,D:1);")
    tr_poly_root2 = ape::read.tree(text="(C:1,D:1,A:1,B:1);")
    tr_poly_root1$root.edge = 0
    tr_poly_root2$root.edge = 0
    expect_true(isTRUE(is_same_root(tr_poly_root1, tr_poly_root2)))
    multi2bi_empty = multi2bi_node_number_transfer(tr_unlabeled, tr_unlabeled)
    expect_true(nrow(multi2bi_empty) == 0)
    set.seed(1)
    tr_poly = ape::read.tree(text="((A:1,B:1,C:1):1,D:1);")
    tr_poly_bi = ape::multi2di(tr_poly)
    tr_poly_bi$tip.label = tr_poly$tip.label
    multi2bi_map = multi2bi_node_number_transfer(tr_poly, tr_poly_bi)
    expect_true(nrow(multi2bi_map) >= 1)


    # Only the selected near-zero internal edge is collapsed; a separate negative
    # edge is preserved with its clade and original length.
    tr_mixed_lengths = ape::read.tree(
        text="(((A:1,B:1):-1,(C:1,D:1):0.000000001):1,E:1);"
    )
    tr_mixed_collapsed = collapse_short_branches(tr_mixed_lengths, tol=1e-8)
    expect_true(any(tr_mixed_collapsed$edge.length == -1))
    expect_true(isTRUE(ape::is.monophyletic(tr_mixed_collapsed, c("A", "B"))))
    expect_true(!isTRUE(ape::is.monophyletic(tr_mixed_collapsed, c("C", "D"))))
})
