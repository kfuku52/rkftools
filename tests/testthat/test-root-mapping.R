test_that("root mapping locates the target split", {
    withr::local_seed(20260810)
    phy <- ape::rtree(25L)
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
})

test_that("legacy nslots values do not change root mapping", {
    phy <- fixture_gene_tree()
    rerooted <- suppressWarnings(phytools::reroot(
        phy,
        node.number=phy$edge[nrow(phy$edge), 2]
    ))

    expect_identical(
        get_phy2_root_in_phy1(phy, rerooted, nslots=1L),
        get_phy2_root_in_phy1(phy, rerooted, nslots=2L)
    )
})
