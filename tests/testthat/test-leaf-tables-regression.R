test_that("leaf tables, collapsed outputs, and placeholders", {
    withr::local_seed(20260831)
    tr = fixture_gene_tree()
    tr_unlabeled = fixture_trait_tree()
    pcm_leaf_no_col = list(
        tree=tr_unlabeled,
        Y=as.data.frame(matrix(c(1, 2, 3), ncol=1, dimnames=list(tr_unlabeled$tip.label, NULL))),
        shift.configuration=structure(integer(0), names=NULL),
        optima=as.data.frame(matrix(c(1, 2, 3), ncol=1, dimnames=list(tr_unlabeled$tip.label, NULL))),
        mu=as.data.frame(matrix(c(1, 2, 3), ncol=1, dimnames=list(tr_unlabeled$tip.label, NULL))),
        residuals=as.data.frame(matrix(c(1, 2, 3), ncol=1, dimnames=list(tr_unlabeled$tip.label, NULL)))
    )
    colnames(pcm_leaf_no_col$Y) = NULL
    colnames(pcm_leaf_no_col$optima) = NULL
    colnames(pcm_leaf_no_col$mu) = NULL
    colnames(pcm_leaf_no_col$residuals) = NULL
    leaf_tbl_no_col = get_leaf_table(pcm_leaf_no_col, mode="l1ou")
    expect_true(!any(is.na(colnames(leaf_tbl_no_col))))
    expect_true("trait1" %in% colnames(leaf_tbl_no_col))
    leaf_tbl_mode_na_err = tryCatch(
        {
            get_leaf_table(pcm_leaf_no_col, mode=NA_character_)
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(leaf_tbl_mode_na_err))
    expect_true(grepl("mode must be a single non-missing string", conditionMessage(leaf_tbl_mode_na_err), fixed=TRUE))
    pcm_leaf_reordered = list(
        tree=tr_unlabeled,
        Y=data.frame(t1=c(30, 20, 10), row.names=rev(tr_unlabeled$tip.label)),
        shift.configuration=structure(integer(0), names=NULL),
        optima=data.frame(t1=c(30, 20, 10), row.names=rev(tr_unlabeled$tip.label)),
        mu=data.frame(t1=c(30, 20, 10), row.names=rev(tr_unlabeled$tip.label)),
        residuals=data.frame(t1=c(30, 20, 10), row.names=rev(tr_unlabeled$tip.label))
    )
    leaf_tbl_reordered = get_leaf_table(pcm_leaf_reordered, mode="l1ou")
    leaf_tbl_reordered_y = leaf_tbl_reordered[as.character(leaf_tbl_reordered$param) == "Y",,drop=FALSE]
    expect_true(identical(as.character(leaf_tbl_reordered_y$node_name), as.character(tr_unlabeled$tip.label)))
    expect_true(identical(as.numeric(leaf_tbl_reordered_y$t1), c(10, 20, 30)))
    pcm_leaf_reordered_cols = list(
        tree=tr_unlabeled,
        Y=data.frame(
            t1=c(1, 2, 3),
            t2=c(4, 5, 6),
            row.names=tr_unlabeled$tip.label
        ),
        shift.configuration=structure(integer(0), names=NULL),
        optima=data.frame(
            t2=c(40, 50, 60),
            t1=c(10, 20, 30),
            row.names=tr_unlabeled$tip.label
        ),
        mu=data.frame(
            t2=c(400, 500, 600),
            t1=c(100, 200, 300),
            row.names=tr_unlabeled$tip.label
        ),
        residuals=data.frame(
            t2=c(4000, 5000, 6000),
            t1=c(1000, 2000, 3000),
            row.names=tr_unlabeled$tip.label
        )
    )
    leaf_tbl_reordered_cols = get_leaf_table(pcm_leaf_reordered_cols, mode="l1ou")
    leaf_tbl_reordered_cols_optima = leaf_tbl_reordered_cols[
        as.character(leaf_tbl_reordered_cols$param) == "optima",,
        drop=FALSE
    ]
    expect_true(identical(as.numeric(leaf_tbl_reordered_cols_optima$t1), c(10, 20, 30)))
    expect_true(identical(as.numeric(leaf_tbl_reordered_cols_optima$t2), c(40, 50, 60)))
    pcm_leaf_bad_col = pcm_leaf_reordered_cols
    colnames(pcm_leaf_bad_col$optima) = c("t2", "wrong_trait")
    leaf_tbl_bad_col_err = tryCatch(
        {
            get_leaf_table(pcm_leaf_bad_col, mode="l1ou")
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(leaf_tbl_bad_col_err))
    expect_true(grepl("must match pcm_out$Y", conditionMessage(leaf_tbl_bad_col_err), fixed=TRUE))
    tree_table_collapsed = data.frame(
        num_shift=0,
        num_regime=1,
        num_conv_regime=0,
        num_uniq_regime=1,
        num_species=0,
        num_leaf=0,
        model_score=0,
        stringsAsFactors=FALSE
    )
    tree_table_restored = tree_table_collapse2original(tree_table_collapsed, tr_unlabeled)
    expect_true(identical(as.integer(tree_table_restored$num_leaf), 3L))
    expect_true(identical(as.integer(tree_table_restored$num_species), 3L))


    tree_collapsed = ape::read.tree(text="(1:1,C:1);")
    node_num_mapping = data.frame(tree_original=c(1L, 2L, 3L), tree_collapsed=c(1L, 1L, 2L))
    regime_input = data.frame(regime="r1", node_name="1", param="shift_value", t1=0.1, stringsAsFactors=FALSE)
    regime_converted = regime_table_collapse2original(
        regime_table=regime_input,
        tree_original=tr_unlabeled,
        tree_collapsed=tree_collapsed,
        node_num_mapping=node_num_mapping
    )
    expect_true(nrow(regime_converted) == 1)
    expect_true(regime_converted$node_name %in% c("A", "B"))

    leaf_input = data.frame(
        regime=c("r1", "r1", NA),
        node_name=c("1", "C", NA),
        param=c("imputed", "imputed", NA),
        t1=c(10, 20, 30),
        stringsAsFactors=FALSE
    )
    leaf_converted = leaf_table_collapse2original(
        leaf_table=leaf_input,
        tree_original=tr_unlabeled,
        tree_collapsed=tree_collapsed,
        node_num_mapping=node_num_mapping
    )
    expect_true(!any(is.na(leaf_converted$param)))
    expect_true(nrow(leaf_converted) == 3)

    restored = restore_imputed_leaves(
        leaf_table=data.frame(regime=0, node_name="A", param="imputed", trait1=NA_real_, stringsAsFactors=FALSE),
        original_trait_table=data.frame(trait1=9, row.names="A")
    )
    expect_true(identical(as.numeric(restored$trait1), 9))

    placeholder_leaf = get_placeholder_leaf(
        tree=tr_unlabeled,
        original_trait_table=data.frame(trait1=c(5, 6, 7), row.names=tr_unlabeled$tip.label)
    )
    expect_true("node_name" %in% colnames(placeholder_leaf))
    expect_true(!("label" %in% colnames(placeholder_leaf)))
    placeholder_regime = get_placeholder_regime(
        tree=tr_unlabeled,
        original_trait_table=data.frame(trait1=c(5, 6, 7), row.names=tr_unlabeled$tip.label)
    )
    expect_true("trait1" %in% colnames(placeholder_regime))
    placeholder_regime_no_col = get_placeholder_regime(
        tree=tr_unlabeled,
        original_trait_table={
            tbl = as.data.frame(matrix(c(5, 6, 7), ncol=1, dimnames=list(tr_unlabeled$tip.label, NULL)))
            colnames(tbl) = NULL
            tbl
        }
    )
    expect_true("trait1" %in% colnames(placeholder_regime_no_col))
    placeholder_tree = get_placeholder_tree(
        tree=tr_unlabeled,
        original_trait_table=data.frame(trait1=c(5, 6, 7), row.names=tr_unlabeled$tip.label)
    )
    expect_true("num_leaf" %in% colnames(placeholder_tree))
})
