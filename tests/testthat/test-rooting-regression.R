test_that("root partitions distinguish labels containing signature delimiters", {
    first <- ape::read.tree(text="((a:1,e:1):1,(d:1,b:1,c:1):1);")
    second <- ape::read.tree(text="((a:1,b:1,e:1):1,(d:1,c:1):1);")
    labels <- c(a="a", b="b", c="c", d="a\rb", e="b\rc")
    first$tip.label <- unname(labels[first$tip.label])
    second$tip.label <- unname(labels[second$tip.label])
    expect_false(is_same_root(first, second))
    expect_true(is_same_root(first, ape::rotate(first, get_root_num(first))))
})

test_that("MAD, root edges, and root label transfer", {
    withr::local_seed(20260831)
    tr = fixture_gene_tree()
    tr_unlabeled = fixture_trait_tree()
    tr_no_edge_length = tr
    tr_no_edge_length$edge.length = NULL
    expect_error(
        MAD(tr_no_edge_length, output_mode="newick"),
        "no branch lengths", fixed=TRUE
    )
    expect_error(
        MAD(tr, output_mode=NA_character_),
        "output_mode must be a single non-missing string", fixed=TRUE
    )
    expect_error(
        MAD_parallel(tr, output_mode=c("newick", "stats")),
        "output_mode must be a single non-missing string", fixed=TRUE
    )

    outgroup_labels = get_outgroup(tr_unlabeled)
    expect_true(identical(as.character(outgroup_labels), "C"))
    expect_error(
        get_outgroup(ape::unroot(tr_unlabeled)),
        "requires a rooted tree", fixed=TRUE
    )
    expect_error(
        get_rooted_newick(tr_unlabeled, madr=999, rho=rep(0.5, nrow(tr_unlabeled$edge))),
        "must be a single edge index", fixed=TRUE
    )
    expect_error(
        get_single_branch_tree(c("A", "B"), 0.5),
        "single non-empty tip label", fixed=TRUE
    )

    tr_two_tip = ape::read.tree(text="(A:1,B:1);")
    tr_two_tip_rr = remove_redundant_root_edge(tr_two_tip)
    expect_true(setequal(tr_two_tip_rr$tip.label, c("A", "B")))
    expect_true(identical(as.integer(tr_two_tip_rr$Nnode), 1L))
    tr_unary_root = ape::read.tree(text="((A:1):1,B:1);")
    tr_unary_root_rr = remove_redundant_root_edge(tr_unary_root)
    expect_true(setequal(tr_unary_root_rr$tip.label, c("A", "B")))
    expect_true(identical(as.integer(tr_unary_root_rr$Nnode), 1L))
    tr_from_label = ape::read.tree(text="((A:1,B:1):1,C:1);")
    tr_from_label$node.label = c("X", "Y")
    tr_to_unlabeled = ape::read.tree(text="((A:1,B:1):1,C:1);")
    tr_to_unlabeled$node.label = NULL
    tr_transferred = transfer_node_labels(tr_from_label, tr_to_unlabeled)
    expect_true(identical(as.character(tr_transferred$node.label), c("X", "Y")))
    expect_error(
        transfer_node_labels(tr_from_label, ape::read.tree(text="((A:1,D:1):1,C:1);")),
        "must contain the same tip labels", fixed=TRUE
    )

    # MAD rejects ambiguous modes and trees without any positive distance.
    expect_error(MAD(tr, output_mode="invalid"), "should be one of", fixed=TRUE)
    all_zero_tree = ape::read.tree(text="((A:0,B:0):0,C:0);")
    expect_error(MAD(all_zero_tree, output_mode="newick"), "no positive branch lengths", fixed=TRUE)
})
