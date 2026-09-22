test_that("clade collapse, node mapping, and trait similarity", {
    withr::local_seed(20260831)
    tr_unlabeled = fixture_trait_tree()
    collapse_single_trait = collapse_clades(
        tr_unlabeled,
        data.frame(v=c(1, 2, 3), row.names=tr_unlabeled$tip.label),
        collapse_node_nums=5L
    )
    expect_true(length(collapse_single_trait$tree$tip.label) == 2)
    expect_true(ncol(collapse_single_trait$trait) == 1)
    expect_true(all(collapse_single_trait$tree$tip.label %in% rownames(collapse_single_trait$trait)))

    collapse_single_trait_dup = suppressWarnings(collapse_clades(
        tr_unlabeled,
        data.frame(v=c(1, 2, 3), row.names=tr_unlabeled$tip.label),
        collapse_node_nums=rep(as.integer(collapse_single_trait$tree$tip.label[collapse_single_trait$tree$tip.label != "C"]), 2)
    ))
    expect_true(length(collapse_single_trait_dup$tree$tip.label) == length(collapse_single_trait$tree$tip.label))
    expect_error(
        map_node_num(
            tr_unlabeled,
            ape::read.tree(text="((A:1,D:1):1,C:1);")
        ),
        "Tip labels are inconsistent", fixed=TRUE
    )
    expect_error(
        collapse_clades(
            tr_unlabeled,
            data.frame(v=c(1, 2, 3), row.names=tr_unlabeled$tip.label),
            collapse_node_nums="x"
        ),
        "must be integer node numbers", fixed=TRUE
    )
    tr_nested = ape::read.tree(text="((A:1,B:1):1,(C:1,D:1):1);")
    collapse_nested = suppressWarnings(collapse_clades(
        tr_nested,
        data.frame(v=c(1, 2, 3, 4), row.names=tr_nested$tip.label),
        collapse_node_nums=c(5L, 7L)
    ))
    expect_true(length(collapse_nested$collapse_leaf_names) == 1)
    expect_true(identical(names(collapse_nested$collapse_leaf_names), "5"))
    expect_true(setequal(collapse_nested$collapse_leaf_names[[1]], c("A", "B", "C", "D")))
    expect_error(
        get_high_similarity_clades(
            tr_unlabeled,
            data.frame(v=c(1, 2, 3), row.names=tr_unlabeled$tip.label),
            method="pearson",
            threshold=0.5,
            num_test=NA
        ),
        "single non-negative integer", fixed=TRUE
    )
    const_trait_warns = character(0)
    collapse_const_trait = withCallingHandlers(
        collapse_clades(
            tr_nested,
            data.frame(v=rep(5, 4), row.names=tr_nested$tip.label),
            collapse_node_nums=5L
        ),
        warning=function(w) {
            const_trait_warns <<- c(const_trait_warns, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    expect_true(!any(grepl("NaNs produced", const_trait_warns, fixed=TRUE)))
    expect_true(isTRUE(all.equal(as.numeric(collapse_const_trait$trait["5", "v"]), 5, tolerance=1e-10)))

    expect_error(
        collapse_clades(tr_unlabeled, data.frame(v=c(1,2,3), row.names=tr_unlabeled$tip.label), collapse_node_nums=c(1)),
        "collapse_node_nums must contain only internal node numbers", fixed=TRUE
    )
    dup_trait_rows = matrix(c(1, 2, 3, 4), ncol=1)
    rownames(dup_trait_rows) = c("A", "A", "B", "C")
    expect_error(
        collapse_clades(
            tr_unlabeled,
            dup_trait_rows,
            collapse_node_nums=5L
        ),
        "duplicated row name", fixed=TRUE
    )

    # Undefined correlations do not make a clade look perfectly similar.
    undefined_similarity_traits = data.frame(
        t1=c(1, 2, 3),
        t2=c(1, 2, 4),
        row.names=c("A", "B", "C")
    )
    undefined_similarity = get_high_similarity_clades(
        tr_unlabeled,
        undefined_similarity_traits,
        method="pearson",
        threshold=0.9
    )
    expect_true(length(undefined_similarity) == 0)
    expect_error(
        get_high_similarity_clades(
            tr_unlabeled,
            data.frame(t1=c(1, 2, 3), label=c("x", "y", "z"), row.names=c("A", "B", "C")),
            method="pearson",
            threshold=0.9
        ),
        "numeric or logical trait columns"
    )
})
