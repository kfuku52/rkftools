test_that("real PhylogeneticEM fits satisfy every adapter contract", {
    skip_if_not_installed("PhylogeneticEM")
    fixture <- readRDS(test_path("fixtures", "phyloem-fit.rds"))
    expect_equal(fixture$backend_version, "1.8.1")
    for (name in names(fixture$models)) {
        model <- fixture$models[[name]]
        traits <- rownames(model$Y_data)
        tips <- model$phylo$tip.label
        parameters <- PhylogeneticEM::params_process(model)
        tree_table <- get_tree_table(model, "PhylogeneticEM")
        expect_equal(tree_table$num_species, 16)
        expect_equal(tree_table$num_leaf, 16)
        expected_shifts <- if (name == "one_shift") 1L else 0L
        expect_equal(tree_table$num_shift, expected_shifts)
        regimes <- get_regime_table(model, "PhylogeneticEM")
        expect_equal(sum(regimes$param == "shift_value"), expected_shifts)
        expect_true(all(vapply(regimes[traits], is.numeric, logical(1))))
        expect_equal(as.numeric(regimes[regimes$param == "alpha",traits]),
            unname(diag(as.matrix(parameters$selection.strength))))
        if (expected_shifts) {
            nodes <- model$phylo$edge[parameters$shifts$edges,2]
            expected_names <- get_node_name_by_num(fill_node_labels(model$phylo), nodes)
            expect_equal(regimes$node_name[regimes$param == "shift_value"], expected_names)
        }
        leaves <- get_leaf_table(model, "PhylogeneticEM")
        expect_identical(names(leaves), c("regime", "node_name", "param", traits))
        expect_true(all(vapply(leaves[traits], is.numeric, logical(1))))
        expect_true(all(is.finite(as.matrix(leaves[traits]))))
        for (param in c("imputed", "expectations")) {
            rows <- leaves[leaves$param == param,,drop=FALSE]
            expect_identical(rows$node_name, tips)
            expected <- t(PhylogeneticEM::imputed_traits(model,
                trait=seq_along(traits), where="tips", what=param))
            expect_equal(unname(as.matrix(rows[traits])), unname(expected))
        }
        restored <- restore_imputed_leaves(leaves, t(model$Y_data))
        rows <- restored[restored$param == "imputed",traits,drop=FALSE]
        observed <- t(model$Y_data)[tips,,drop=FALSE]
        has_value <- !is.na(observed)
        expect_equal(as.matrix(rows)[has_value], observed[has_value])
        expect_false(anyNA(rows))
    }
})

test_that("Rphylopars fits numeric matrices and data frames equivalently", {
    skip_if_not_installed("Rphylopars")
    phy <- ape::read.tree(text="((A:1,B:1):1,C:1);")
    traits <- data.frame("trait one"=c(1, NA, 3), row.names=phy$tip.label, check.names=FALSE)
    result <- withr::with_seed(123, phylogenetic_imputation(phy, traits))
    matrix_result <- withr::with_seed(123, phylogenetic_imputation(phy, as.matrix(traits)))
    expect_s3_class(result, "data.frame")
    expect_identical(rownames(result), phy$tip.label)
    expect_identical(names(result), "trait one")
    expect_true(is.numeric(result[[1]]))
    expect_true(all(is.finite(result[[1]])))
    expect_equal(matrix_result, result)
})
