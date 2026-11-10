test_that("the bundled example plots scores in the input edge order", {
    withr::local_pdf(withr::local_tempfile(fileext=".pdf"), width=12, height=7)
    tree = ape::unroot(species_overlap_example_tree)
    scores = c(2, 3, 3, 2, 3, 3, 1, 2, 3, 3, 2, 3, 3)
    result = withVisible(plot_root_position_species_overlap(tree))
    expect_false(result$visible)
    expect_identical(result$value, scores)

    reordered = ape::reorder.phylo(tree, order="postorder")
    edge_index = match(reordered$edge[,2], tree$edge[,2])
    expect_false(identical(edge_index, seq_along(scores)))
    expect_identical(plot_root_position_species_overlap(reordered), scores[edge_index])
})

test_that("a rooted tree can plot tied zero scores with a custom separator", {
    withr::local_pdf(withr::local_tempfile(fileext=".pdf"), width=12, height=7)
    tree = ape::read.tree(text="(Homo|sapiens|geneA:1,Mus|musculus|geneA:1);")
    old_par = graphics::par(c("mar", "oma", "mfrow"))
    expect_identical(plot_root_position_species_overlap(tree, sep="|"), c(0, 0))
    expect_equal(graphics::par(names(old_par)), old_par)
})
