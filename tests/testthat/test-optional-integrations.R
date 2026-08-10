test_that("optional PhylogeneticEM API remains available", {
    skip_if_not_installed("PhylogeneticEM")

    expect_true(is.function(getFromNamespace("params_process", "PhylogeneticEM")))
    expect_true(is.function(getFromNamespace("imputed_traits", "PhylogeneticEM")))
})

test_that("Rphylopars integration imputes tip traits", {
    skip_if_not_installed("Rphylopars")

    phy <- ape::read.tree(text="((A:1,B:1):1,C:1);")
    traits <- data.frame(value=c(1, NA, 3), row.names=phy$tip.label)
    result <- phylogenetic_imputation(phy, traits)

    expect_s3_class(result, "data.frame")
    expect_identical(rownames(result), phy$tip.label)
    expect_false(anyNA(result$value))
})
