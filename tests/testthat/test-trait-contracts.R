test_that("restoring observations never erases imputed missing cells", {
    original <- data.frame(x=c(1, NA), y=c(NA, 4), row.names=c("A", "B"))
    leaves <- data.frame(node_name=c("B", "A", "A"),
        param=c("imputed", "imputed", "observed"), x=c(20, 10, 9), y=c(40, 30, 8))
    out <- restore_imputed_leaves(leaves, original)
    expect_equal(out$x, c(20, 1, 9))
    expect_equal(out$y, c(4, 30, 8))
    expect_identical(out$param, leaves$param)
    expect_equal(restore_imputed_leaves(leaves, as.matrix(original)), out)
})

test_that("replicate means retain dimensions and distinguish NA from NaN", {
    one <- matrix(c(1, 3), nrow=1, dimnames=list("A", c("x_1", "x_2")))
    expect_equal(merge_replicates(one, "_"), data.frame(x=2, row.names="A"))
    missing <- rbind(one, c(NA, NA), c(NA, 5))
    rownames(missing) <- c("A", "B", "C")
    out <- merge_replicates(missing, "_")
    expect_equal(out$x, c(2, NA_real_, 5))
    expect_false(any(is.nan(out$x)))
})

test_that("PhylogeneticEM adapters preserve trait names, types, and species counts", {
    tree <- ape::read.tree(text="((A_sp:1,B_sp:1):1,C_sp:2);")
    traits <- c("trait one", "b")
    pp <- list(shifts=list(edges=1L,
            values=matrix(c(20, 10), 2, dimnames=list(rev(traits), NULL))),
        selection.strength=diag(c(0.2, 0.1)), variance=diag(c(2, 1)),
        optimal.value=setNames(c(4, 3), rev(traits)))
    dimnames(pp$selection.strength) <- dimnames(pp$variance) <- list(rev(traits), rev(traits))
    attr(pp, "log_likelihood") <- -7
    local_mocked_bindings(.phylogeneticem_params_process=function(...) pp)
    model <- list(phylo=tree, Y_data=matrix(1:6, 2, dimnames=list(traits, tree$tip.label)))
    out <- get_regime_table(model, "PhylogeneticEM")
    expect_identical(names(out), c("regime", "node_name", "param", traits))
    expect_true(all(vapply(out[traits], is.numeric, logical(1))))
    expect_equal(out[["trait one"]], c(10, 0.1, 1, 3))
    expect_equal(out$b, c(20, 0.2, 2, 4))
    expect_equal(out$node_name[[1]], fill_node_labels(tree)$node.label[[2]])
    expect_equal(get_tree_table(model, "PhylogeneticEM")$num_species, 3)
    pp$shifts <- list(edges=integer(0), values=matrix(numeric(0), 2, 0))
    pp$selection.strength <- matrix(0.5, 1, 1)
    empty <- get_regime_table(model, "PhylogeneticEM")
    expect_equal(empty$param, c("alpha", "sigma2", "intercept"))
    expect_equal(empty[1,traits], data.frame("trait one"=0.5, b=0.5, check.names=FALSE))
    expect_equal(get_tree_table(model, "PhylogeneticEM")$num_shift, 0)
})

test_that("model and placeholder tables reject ambiguous trait names", {
    tree <- fixture_trait_tree()
    for (trait_names in list("regime", "node_name", "param", c("x","x"), "", NA_character_)) {
        y <- matrix(1, 3, length(trait_names), dimnames=list(tree$tip.label, trait_names))
        fit <- list(tree=tree, Y=y, optima=y, mu=y, residuals=y,
            shift.configuration=integer(), nShifts=0, alpha=1, sigma2=1, intercept=0, logLik=-1)
        expect_error(get_regime_table(fit, "l1ou"), "trait names")
        expect_error(get_leaf_table(fit, "l1ou"), "trait names")
        expect_error(get_placeholder_leaf(tree, y), "trait names")
        expect_error(get_placeholder_regime(tree, y), "trait names")
    }
})
