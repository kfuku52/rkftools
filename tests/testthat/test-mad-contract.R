test_that("MAD clamps negative input edges before unrooting combines them", {
    tree <- ape::read.tree(text="((A:1,B:3):-10,C:2);")
    nonnegative <- tree
    nonnegative$edge.length[nonnegative$edge.length < 0] <- 0
    expected <- MAD(nonnegative, "custom")
    expect_warning(actual <- MAD(tree, "custom"), "negative branch lengths")
    expect_equal(actual, expected)
    expect_warning(parallel <- MAD_parallel(tree, "custom", ncpu=2L),
        "negative branch lengths")
    expect_equal(parallel, expected)
    expect_equal(ape::cophenetic.phylo(actual[[3]]),
        ape::cophenetic.phylo(nonnegative))
    expect_equal(tree$edge.length, c(-10, 1, 3, 2))
    tree$edge.length[c(1, 4)] <- .Machine$double.xmax
    expect_error(MAD(tree), "Unrooting produced non-finite branch lengths", fixed=TRUE)
})

test_that("zero-distance duplicates survive every MAD output consistently", {
    for (newick in c("((A:0,B:0):1,C:2,D:3);", "((A:0,B:0):1,(C:0,D:0):2,E:3);")) {
        tree <- ape::read.tree(text=newick)
        for (mode in c("full", "custom")) {
            result <- suppressWarnings(MAD(tree, mode))
            expect_setequal(result[[3]]$tip.label, tree$tip.label)
            expect_length(result[[5]], nrow(result[[3]]$edge))
            expect_equal(result[[5]][result[[4]]],
                rep(min(result[[5]], na.rm=TRUE), length(result[[4]])))
            for (i in seq_along(result[[6]])) {
                restored <- result[[6]][[i]]
                expect_setequal(restored$tip.label, tree$tip.label)
                parsed <- ape::read.tree(text=result[[1]][[i]])
                expect_setequal(parsed$tip.label, tree$tip.label)
                expect_equal(ape::cophenetic.phylo(parsed)[tree$tip.label,tree$tip.label],
                    ape::cophenetic.phylo(restored)[tree$tip.label,tree$tip.label], tolerance=1e-8)
            }
            if (mode == "custom") {
                expect_length(result[[7]], nrow(result[[3]]$edge))
                for (i in seq_along(result[[4]])) {
                    rooted <- get_rooted_newick(result[[3]], result[[4]][[i]], result[[7]])
                    expect_equal(rooted[[2]], result[[6]][[i]])
                }
            }
        }
        expect_equal(suppressWarnings(MAD_parallel(tree, "custom", ncpu=2L)),
            suppressWarnings(MAD(tree, "custom")))
    }
})
