test_that("clade collapse, node mapping, and trait similarity", {
    withr::local_seed(20260831)
    tr_unlabeled = fixture_trait_tree()
    collapse_single_trait = collapse_clades(
        tr_unlabeled,
        data.frame(v=c(1, 2, 3), row.names=tr_unlabeled$tip.label),
        collapse_node_nums={
            internal_nodes = sort(unique(tr_unlabeled$edge[,1]))
            internal_nodes = internal_nodes[internal_nodes > length(tr_unlabeled$tip.label)]
            size_two_node = internal_nodes[vapply(internal_nodes, function(nn) {
                length(get_tip_labels(tr_unlabeled, nn)) == 2
            }, logical(1))]
            size_two_node[1]
        }
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
    map_identity = map_node_num(tr_unlabeled, tr_unlabeled, collapse_leaf_names=list())
    expect_true(nrow(map_identity) == max(tr_unlabeled$edge[,1]))
    expect_true(identical(as.integer(map_identity$tree_original), as.integer(map_identity$tree_collapsed)))
    map_identity_default = map_node_num(tr_unlabeled, tr_unlabeled)
    expect_true(identical(as.integer(map_identity_default$tree_original), as.integer(map_identity_default$tree_collapsed)))
    map_verbose_na_err = tryCatch(
        {
            map_node_num(tr_unlabeled, tr_unlabeled, verbose=NA)
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(map_verbose_na_err))
    expect_true(grepl("verbose must be a single non-missing logical value", conditionMessage(map_verbose_na_err), fixed=TRUE))
    map_tip_mismatch_err = tryCatch(
        {
            map_node_num(
                tr_unlabeled,
                ape::read.tree(text="((A:1,D:1):1,C:1);")
            )
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(map_tip_mismatch_err))
    expect_true(grepl("Tip labels are inconsistent", conditionMessage(map_tip_mismatch_err), fixed=TRUE))
    collapse_single_trait_chr = collapse_clades(
        tr_unlabeled,
        data.frame(v=c(1, 2, 3), row.names=tr_unlabeled$tip.label),
        collapse_node_nums=as.character({
            internal_nodes = sort(unique(tr_unlabeled$edge[,1]))
            internal_nodes = internal_nodes[internal_nodes > length(tr_unlabeled$tip.label)]
            size_two_node = internal_nodes[vapply(internal_nodes, function(nn) {
                length(get_tip_labels(tr_unlabeled, nn)) == 2
            }, logical(1))]
            size_two_node[1]
        })
    )
    expect_true(length(collapse_single_trait_chr$tree$tip.label) == 2)
    collapse_invalid_type_error = tryCatch(
        {
            collapse_clades(
                tr_unlabeled,
                data.frame(v=c(1, 2, 3), row.names=tr_unlabeled$tip.label),
                collapse_node_nums="x"
            )
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(collapse_invalid_type_error))
    expect_true(grepl("must be integer node numbers", conditionMessage(collapse_invalid_type_error), fixed=TRUE))
    tr_nested = ape::read.tree(text="((A:1,B:1):1,(C:1,D:1):1);")
    collapse_nested = suppressWarnings(collapse_clades(
        tr_nested,
        data.frame(v=c(1, 2, 3, 4), row.names=tr_nested$tip.label),
        collapse_node_nums=c(5L, 7L)
    ))
    expect_true(length(collapse_nested$collapse_leaf_names) == 1)
    expect_true(identical(names(collapse_nested$collapse_leaf_names), "5"))
    expect_true(setequal(collapse_nested$collapse_leaf_names[[1]], c("A", "B", "C", "D")))
    high_sim_num_test_na_err = tryCatch(
        {
            get_high_similarity_clades(
                tr_unlabeled,
                data.frame(v=c(1, 2, 3), row.names=tr_unlabeled$tip.label),
                method="pearson",
                threshold=0.5,
                num_test=NA
            )
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(high_sim_num_test_na_err))
    expect_true(grepl("single non-negative integer", conditionMessage(high_sim_num_test_na_err), fixed=TRUE))
    high_sim_verbose_na_err = tryCatch(
        {
            get_high_similarity_clades(
                tr_unlabeled,
                data.frame(v=c(1, 2, 3), row.names=tr_unlabeled$tip.label),
                method="pearson",
                threshold=0.5,
                verbose=NA
            )
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(high_sim_verbose_na_err))
    expect_true(grepl("verbose must be a single non-missing logical value", conditionMessage(high_sim_verbose_na_err), fixed=TRUE))
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

    collapse_leaf_error = tryCatch(
        {
            collapse_clades(tr_unlabeled, data.frame(v=c(1,2,3), row.names=tr_unlabeled$tip.label), collapse_node_nums=c(1))
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(collapse_leaf_error))
    expect_true(grepl("collapse_node_nums must contain only internal node numbers", conditionMessage(collapse_leaf_error), fixed=TRUE))
    dup_trait_rows = matrix(c(1, 2, 3, 4), ncol=1)
    rownames(dup_trait_rows) = c("A", "A", "B", "C")
    dup_trait_high_sim_err = tryCatch(
        {
            get_high_similarity_clades(
                tr_unlabeled,
                dup_trait_rows,
                method="pearson",
                threshold=0.5
            )
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(dup_trait_high_sim_err))
    expect_true(grepl("duplicated row name", conditionMessage(dup_trait_high_sim_err), fixed=TRUE))
    dup_trait_collapse_err = tryCatch(
        {
            collapse_clades(
                tr_unlabeled,
                dup_trait_rows,
                collapse_node_nums=5L
            )
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(dup_trait_collapse_err))
    expect_true(grepl("duplicated row name", conditionMessage(dup_trait_collapse_err), fixed=TRUE))


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
    expect_error_contains(
        get_high_similarity_clades(
            tr_unlabeled,
            data.frame(t1=c(1, 2, 3), label=c("x", "y", "z"), row.names=c("A", "B", "C")),
            method="pearson",
            threshold=0.9
        ),
        "numeric or logical trait columns"
    )
})
