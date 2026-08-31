test_that("padding preserves all tip depths through several ancestors", {
    tree <- ape::read.tree(text="(((A:0.1,B:0.1):0.2,C:0.3):2,D:2.3);")
    for (external in c(TRUE, FALSE)) {
        padded <- pad_short_edges(tree, 0.5, external_only=external)
        expect_equal(ape::node.depth.edgelength(padded)[1:4], rep(2.3, 4))
        expect_true(ape::is.ultrametric(padded))
        expect_identical(padded$edge, tree$edge)
        target <- !external | padded$edge[,2] <= 4L
        expect_gte(min(padded$edge.length[target]), 0.5)
        expect_equal(pad_short_edges(padded, 0.5, external), padded)
    }
})

test_that("root movement extends every path equally on all topologies", {
    for (newick in c("(A:0.1,B:0.7,C:0.8);",
                     "(((A:0.1,B:0.2):0.1):0.1,C:0.5,D:0.8);",
                     "((A:0.1,B:0.2,C:0.3):0.1,D:0.8);")) {
        tree <- ape::read.tree(text=newick)
        tips <- seq_along(tree$tip.label)
        before <- ape::node.depth.edgelength(tree)[tips]
        expect_warning(padded <- pad_short_edges(tree, 0.5), NA)
        delta <- ape::node.depth.edgelength(padded)[tips] - before
        expect_equal(delta, rep(delta[[1]], length(tips)))
        expect_gte(min(padded$edge.length), 0.5)
        expect_identical(padded$edge, tree$edge)
    }
})

test_that("short internal edges can borrow from children without moving the root", {
    tree <- ape::read.tree(text="((A:2,B:2):0.1,C:2.1);")
    padded <- pad_short_edges(tree, 0.5)
    expect_equal(ape::node.depth.edgelength(padded)[1:3], rep(2.1, 3))
    expect_equal(padded$edge.length, c(0.5, 1.6, 1.6, 2.1))
    expect_true(ape::is.ultrametric(padded))
})

test_that("padding attains the minimum extension implied by every tip path", {
    withr::local_seed(20260831)
    for (iteration in seq_len(12L)) {
        tree <- ape::rtree(16L)
        tips <- seq_along(tree$tip.label)
        before <- ape::node.depth.edgelength(tree)[tips]
        unit_tree <- tree
        unit_tree$edge.length <- rep(1, nrow(tree$edge))
        path_edges <- ape::node.depth.edgelength(unit_tree)[tips]
        for (external in c(FALSE, TRUE)) {
            # Every path needs one terminal minimum, or one per edge.
            needed <- if (external) rep(0.4, length(tips)) else 0.4 * path_edges
            extension <- max(0, needed - before)
            padded <- pad_short_edges(tree, 0.4, external_only=external)
            expect_equal(ape::node.depth.edgelength(padded)[tips] - before,
                rep(extension, length(tips)), tolerance=1e-12)
            target <- !external | padded$edge[,2] <= length(tips)
            expect_gte(min(padded$edge.length[target]), 0.4)
        }
    }
})

test_that("padding preserves missing lengths without contaminating known edges", {
    tree <- ape::read.tree(text="((A:0.1,B:0.2):1,C:2);")
    tree$edge.length[[1]] <- NA_real_
    padded <- pad_short_edges(tree, 0.5)
    expect_identical(is.na(padded$edge.length), is.na(tree$edge.length))
    expect_gte(min(padded$edge.length, na.rm=TRUE), 0.5)
    expect_equal(padded$edge.length[[4]], 2)
})
