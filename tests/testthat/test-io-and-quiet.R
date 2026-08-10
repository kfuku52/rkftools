test_that("NOTUNG parser rejects malformed records", {
    path <- tempfile()
    on.exit(unlink(path), add=TRUE)
    writeLines("#D gene species", path)

    expect_error(read_notung_parsable(path), "Malformed NOTUNG")
    expect_error(read_notung_parsable(paste0(path, "-missing")), "does not exist")
})

test_that("transformation helpers are quiet by default", {
    phy <- ape::read.tree(text="((A:1e-9,B:1):1,C:1);")
    traits <- data.frame(x=c(1, 1, 2), row.names=phy$tip.label)

    expect_message(remove_invariant_traits(traits), NA)
    expect_message(pad_short_edges(phy), NA)
    expect_message(fill_node_labels(phy), NA)
})
