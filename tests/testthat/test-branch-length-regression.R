test_that("short branches and topology transformations", {
    withr::local_seed(20260831)
    tr_unlabeled = fixture_trait_tree()
    tr_na_edge = tr_unlabeled
    tr_na_edge$edge.length[1] = NA_real_
    collapsed_na_edge = collapse_short_branches(tr_na_edge, tol=1e-8)
    expect_true(is.na(collapsed_na_edge$edge.length[1]))
    expect_error(
        force_ultrametric(tr_na_edge),
        "contains NA edge lengths", fixed=TRUE
    )

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

test_that("short-edge contraction preserves individual unary edges", {
    tree <- ape::read.tree(text="((((A:1)U:2)W:3,B:6)V:0,C:7)R;")
    expected <- ape::read.tree(text="(((A:1)U:2)W:3,B:6,C:7)R;")
    for (order in c("cladewise", "postorder")) {
        input <- ape::reorder.phylo(tree, order)
        result <- collapse_short_branches(input)
        expect_equal(ape::write.tree(result), ape::write.tree(expected))
        expect_equal(ape::cophenetic.phylo(result), ape::cophenetic.phylo(tree))
        expect_identical(input, ape::reorder.phylo(tree, order))
    }
    # Edge identity must not depend on labels or finite, positive lengths.
    tree$node.label[] <- "duplicate"
    tree$edge.length[tree$edge[,2] == 1L] <- NA_real_
    tree$edge.length[tree$edge[,2] == 7L] <- -2
    result <- collapse_short_branches(tree)
    expect_true(is.na(result$edge.length[result$edge[,2] == 1L]))
    expect_equal(sort(result$edge.length, na.last=TRUE), c(-2, 3, 6, 7, NA_real_))
    expect_true(all(result$node.label == "duplicate"))
})
