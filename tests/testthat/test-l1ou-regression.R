test_that("l1ou model summaries, regimes, and parameter alignment", {
    withr::local_seed(20260831)
    tr = fixture_gene_tree()
    tr_unlabeled = fixture_trait_tree()
    pcm_mock = list(
        tree=tr,
        shift.configuration=structure(c(1L, 2L, 3L), names=c("1", "1", "2")),
        nShifts=3L,
        score=1.23
    )
    tree_table_mock = get_tree_table(pcm_mock, mode="l1ou")
    expect_true(identical(as.integer(tree_table_mock$num_regime), 3L))
    expect_true(identical(as.integer(tree_table_mock$num_conv_regime), 1L))
    expect_true(identical(as.integer(tree_table_mock$num_uniq_regime), 2L))
    tree_table_mode_na_err = tryCatch(
        {
            get_tree_table(pcm_mock, mode=NA_character_)
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(tree_table_mode_na_err))
    expect_true(grepl("mode must be a single non-missing string", conditionMessage(tree_table_mode_na_err), fixed=TRUE))


    pcm_reg = list(
        Y=data.frame(t1=c(1, 2, 3), row.names=tr_unlabeled$tip.label),
        tree=tr_unlabeled,
        shift.configuration=structure(c(1L, 2L), names=c("regA", "regB")),
        nShifts=2L,
        shift.values=matrix(c(0.1, 0.2), ncol=1),
        shift.means=matrix(c(1.1, 1.2), ncol=1),
        alpha=1,
        sigma2=1,
        intercept=0,
        logLik=-1
    )
    regime_tbl = get_regime_table(pcm_reg, mode="l1ou")
    expect_true(any(as.character(regime_tbl$regime) == "regA"))
    expect_true(any(as.character(regime_tbl$regime) == "regB"))
    regime_tbl_shift_values = regime_tbl[as.character(regime_tbl$param) == "shift_value",,drop=FALSE]
    expect_true(identical(as.character(regime_tbl_shift_values$regime), c("regA", "regB")))
    expect_true(identical(as.numeric(regime_tbl_shift_values$t1), c(0.1, 0.2)))
    pcm_reg_reversed = pcm_reg
    pcm_reg_reversed$shift.configuration = structure(c(2L, 1L), names=c("regB", "regA"))
    pcm_reg_reversed$shift.values = matrix(c(0.2, 0.1), ncol=1)
    pcm_reg_reversed$shift.means = matrix(c(1.2, 1.1), ncol=1)
    regime_tbl_reversed = get_regime_table(pcm_reg_reversed, mode="l1ou")
    regime_tbl_reversed_shift_values = regime_tbl_reversed[
        as.character(regime_tbl_reversed$param) == "shift_value",,
        drop=FALSE
    ]
    expect_true(identical(as.character(regime_tbl_reversed_shift_values$regime), c("regB", "regA")))
    expect_true(identical(as.numeric(regime_tbl_reversed_shift_values$t1), c(0.2, 0.1)))
    regime_tbl_mode_na_err = tryCatch(
        {
            get_regime_table(pcm_reg, mode=NA_character_)
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(regime_tbl_mode_na_err))
    expect_true(grepl("mode must be a single non-missing string", conditionMessage(regime_tbl_mode_na_err), fixed=TRUE))
    pcm_reg_multi = pcm_reg
    pcm_reg_multi$Y = data.frame(
        t1=c(1, 2, 3),
        t2=c(4, 5, 6),
        row.names=tr_unlabeled$tip.label
    )
    pcm_reg_multi$shift.values = matrix(c(0.1, 0.2, 0.3, 0.4), ncol=2, byrow=TRUE)
    pcm_reg_multi$shift.means = matrix(c(1.1, 1.2, 1.3, 1.4), ncol=2, byrow=TRUE)
    pcm_reg_multi$alpha = 2
    pcm_reg_multi$sigma2 = 3
    pcm_reg_multi$intercept = 4
    pcm_reg_multi$logLik = -5
    regime_tbl_multi = get_regime_table(pcm_reg_multi, mode="l1ou")
    alpha_row_multi = regime_tbl_multi[as.character(regime_tbl_multi$param) == "alpha", c("t1", "t2"), drop=FALSE]
    expect_true(nrow(alpha_row_multi) == 1)
    expect_true(identical(as.numeric(alpha_row_multi[1,]), c(2, 2)))
    pcm_reg_bad_alpha = pcm_reg_multi
    pcm_reg_bad_alpha$alpha = c(1, 2, 3)
    regime_tbl_bad_alpha_err = tryCatch(
        {
            get_regime_table(pcm_reg_bad_alpha, mode="l1ou")
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(regime_tbl_bad_alpha_err))
    expect_true(grepl("must have length 1 or 2", conditionMessage(regime_tbl_bad_alpha_err), fixed=TRUE))
    pcm_reg_bad_shift_type = pcm_reg
    pcm_reg_bad_shift_type$shift.configuration = structure(c("x", "2"), names=c("regA", "regB"))
    regime_tbl_shift_type_err = tryCatch(
        {
            get_regime_table(pcm_reg_bad_shift_type, mode="l1ou")
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(regime_tbl_shift_type_err))
    expect_true(grepl("must contain integer edge indices", conditionMessage(regime_tbl_shift_type_err), fixed=TRUE))
    pcm_reg_na_shift = pcm_reg
    pcm_reg_na_shift$shift.configuration = structure(c(NA_integer_, 2L), names=c("regA", "regB"))
    regime_tbl_na_shift_err = tryCatch(
        {
            get_regime_table(pcm_reg_na_shift, mode="l1ou")
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(regime_tbl_na_shift_err))
    expect_true(grepl("must contain integer edge indices", conditionMessage(regime_tbl_na_shift_err), fixed=TRUE))
    leaf_regimes_ok = get_leaf_regimes(pcm_reg, mode="l1ou")
    expect_true(identical(as.character(leaf_regimes_ok$label), as.character(tr_unlabeled$tip.label)))
    leaf_regime_by_label = setNames(as.character(leaf_regimes_ok$regime), as.character(leaf_regimes_ok$label))
    expect_true(identical(leaf_regime_by_label[c("A", "B", "C")], c(A="regB", B="regA", C="0")))
    leaf_regimes_reversed = get_leaf_regimes(pcm_reg_reversed, mode="l1ou")
    leaf_regime_reversed_by_label = setNames(
        as.character(leaf_regimes_reversed$regime),
        as.character(leaf_regimes_reversed$label)
    )
    expect_true(identical(leaf_regime_reversed_by_label[c("A", "B", "C")], c(A="regB", B="regA", C="0")))
    leaf_regimes_mode_na_err = tryCatch(
        {
            get_leaf_regimes(pcm_reg, mode=NA_character_)
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(leaf_regimes_mode_na_err))
    expect_true(grepl("mode must be a single non-missing string", conditionMessage(leaf_regimes_mode_na_err), fixed=TRUE))
    pcm_reg_bad_labels = pcm_reg
    rownames(pcm_reg_bad_labels$Y) = c("x", "y", "z")
    leaf_regimes_err = tryCatch(
        {
            get_leaf_regimes(pcm_reg_bad_labels, mode="l1ou")
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(leaf_regimes_err))
    expect_true(grepl("must match tree tip labels", conditionMessage(leaf_regimes_err), fixed=TRUE))
    pcm_reg_bad_shift = pcm_reg
    pcm_reg_bad_shift$shift.configuration = structure(999L, names="regA")
    leaf_regimes_shift_err = tryCatch(
        {
            get_leaf_regimes(pcm_reg_bad_shift, mode="l1ou")
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(leaf_regimes_shift_err))
    expect_true(grepl("invalid edge index", conditionMessage(leaf_regimes_shift_err), fixed=TRUE))
    leaf_regimes_na_shift_err = tryCatch(
        {
            get_leaf_regimes(pcm_reg_na_shift, mode="l1ou")
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(leaf_regimes_na_shift_err))
    expect_true(grepl("must contain integer edge indices", conditionMessage(leaf_regimes_na_shift_err), fixed=TRUE))

    # L1OU parameter values are aligned by trait name and remain numeric.
    pcm_reg_named = pcm_reg_multi
    pcm_reg_named$shift.values = matrix(
        c(20, 10, 40, 30),
        nrow=2,
        byrow=TRUE,
        dimnames=list(NULL, c("t2", "t1"))
    )
    pcm_reg_named$shift.means = matrix(
        c(120, 110, 140, 130),
        nrow=2,
        byrow=TRUE,
        dimnames=list(NULL, c("t2", "t1"))
    )
    pcm_reg_named$alpha = c(t2=2, t1=1)
    pcm_reg_named$sigma2 = c(t2=4, t1=3)
    pcm_reg_named$intercept = c(t2=6, t1=5)
    named_regimes = get_regime_table(pcm_reg_named, mode="l1ou")
    named_shift = named_regimes[named_regimes$param == "shift_value",,drop=FALSE]
    named_alpha = named_regimes[named_regimes$param == "alpha",,drop=FALSE]
    expect_true(identical(as.numeric(named_shift[1, c("t1", "t2")]), c(10, 20)))
    expect_true(identical(as.numeric(named_alpha[1, c("t1", "t2")]), c(1, 2)))
    expect_true(is.numeric(named_regimes$t1))
    expect_true(is.numeric(named_regimes$t2))
    expect_error_contains(get_tree_table(pcm_mock, mode="unsupported"), "mode must be one of")
})
