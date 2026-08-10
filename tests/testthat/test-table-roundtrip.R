test_that("unary internal nodes round-trip through branch tables", {
    tab <- data.frame(
        branch_id=c("r", "u", "A", "B"),
        parent=c("-999", "r", "u", "r"),
        sister=c("-999", "B", "-999", "u"),
        label=c("Root", "Unary", "A", "B"),
        dist=c(0, 1, 1, 1),
        stringsAsFactors=FALSE
    )

    phy <- table2phylo(tab, "label", "dist")
    restored <- phylo2table(phy)
    roundtrip <- table2phylo(restored, "label", "dist")

    expect_equal(ape::write.tree(roundtrip), ape::write.tree(phy))
    unary_row <- restored[restored$label == "A",,drop=FALSE]
    expect_identical(as.integer(unary_row$sister), -999L)
})

test_that("root edge lengths and exact zero edges are preserved", {
    tab <- data.frame(
        branch_id=c(2, 0, 1),
        parent=c(-999, 2, 2),
        sister=c(-999, 1, 0),
        label=c("Root", "A", "B"),
        dist=c(0.75, 0, 1),
        stringsAsFactors=FALSE
    )

    phy <- table2phylo(tab, "label", "dist")
    restored <- phylo2table(phy)

    expect_equal(phy$root.edge, 0.75)
    expect_equal(restored$dist[restored$parent == -999], 0.75)
    expect_true(any(restored$dist == 0))
})

test_that("non-finite branch lengths are rejected", {
    phy <- ape::read.tree(text="(A:1,B:1);")
    phy$edge.length[[1]] <- Inf

    expect_error(phylo2table(phy), "finite")
    expect_error(get_single_branch_tree("A", Inf), "finite")
})
