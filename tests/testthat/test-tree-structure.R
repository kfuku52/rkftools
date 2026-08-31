test_that("malformed phylo graphs are rejected before traversal", {
    bad <- structure(list(
        edge=matrix(c(5,1,5,2,6,7,6,3,7,6,7,4), ncol=2, byrow=TRUE),
        tip.label=c("A", "B", "C", "D"), Nnode=3L, edge.length=rep(1, 6)),
        class="phylo")
    expect_error(phylo2table(bad), "disconnected or cyclic")
    expect_error(get_descendent_num(bad, 6L), "disconnected or cyclic")
    expect_error(MAD(bad), "disconnected or cyclic")
    tree <- ape::read.tree(text="((A:1,B:1):1,C:2);")
    fractional <- tree
    fractional$edge[1,2] <- 4.5
    expect_error(phylo2table(fractional), "integer node numbers")
    duplicated <- tree
    duplicated$edge[4,2] <- duplicated$edge[3,2]
    expect_error(phylo2table(duplicated), "counts|multiply-parented")
    bad_tip <- tree
    bad_tip$edge[2,1] <- 1L
    expect_error(phylo2table(bad_tip), "tip node used as a parent")
    bad_count <- tree
    bad_count$Nnode <- 3L
    expect_error(phylo2table(bad_count), "counts inconsistent")
    bad_count$Nnode <- 1.5
    expect_error(phylo2table(bad_count), "integer Nnode")
})

test_that("deep species-overlap scoring does not depend on the recursion limit", {
    tree <- ape::stree(4096L, type="left")
    tree$tip.label <- paste0("A_sp_g", seq_len(4096L))
    expect_equal(get_species_overlap_score(tree), 4095)
    expect_equal(get_species_overlap_score(tree, dc_cutoff=1), 0)
})
