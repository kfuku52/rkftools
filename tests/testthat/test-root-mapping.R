test_that("single-pass root mapping agrees with rerooted trees", {
    set.seed(20260810)
    for (n_tip in c(10L, 25L, 50L)) {
        phy <- ape::rtree(n_tip)
        edge_index <- nrow(phy$edge)
        rerooted <- suppressWarnings(phytools::reroot(
            phy,
            node.number=phy$edge[edge_index, 2]
        ))

        result_index <- get_phy2_root_in_phy1(
            phy, rerooted, nslots=1L, mode="index"
        )
        expect_false(is.na(result_index))
        candidate <- suppressWarnings(phytools::reroot(
            phy,
            node.number=phy$edge[result_index, 2]
        ))
        expect_true(is_same_root(candidate, rerooted))
    }
})

test_that("legacy nslots values do not change root mapping", {
    phy <- ape::rtree(20)
    rerooted <- suppressWarnings(phytools::reroot(
        phy,
        node.number=phy$edge[nrow(phy$edge), 2]
    ))

    expect_identical(
        get_phy2_root_in_phy1(phy, rerooted, nslots=1L),
        get_phy2_root_in_phy1(phy, rerooted, nslots=2L)
    )
})
