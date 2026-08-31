test_that("bootstrap table shape, names, and types", {
    withr::local_seed(20260831)
    tr_unlabeled = fixture_trait_tree()
    pcm_bt = list(tree=tr_unlabeled)
    bt_tbl = get_bootstrap_table(pcm_bt, list(detection.rate=c(1, 2, 3, 4)), mode="l1ou")
    expect_true(nrow(bt_tbl) == (length(tr_unlabeled$tip.label) + tr_unlabeled$Nnode))
    expect_true(sum(is.na(bt_tbl$bootstrap_support)) == 1)
    bt_mode_na_err = tryCatch(
        {
            get_bootstrap_table(pcm_bt, list(detection.rate=c(1, 2, 3, 4)), mode=NA_character_)
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(bt_mode_na_err))
    expect_true(grepl("mode must be a single non-missing string", conditionMessage(bt_mode_na_err), fixed=TRUE))
    bt_tbl_factor = get_bootstrap_table(
        pcm_bt,
        list(detection.rate=factor(c("0.1", "0.2", "0.3", "0.4"))),
        mode="l1ou"
    )
    expect_true(isTRUE(all.equal(
        as.numeric(bt_tbl_factor$bootstrap_support),
        c(0.1, 0.2, 0.3, NA, 0.4),
        tolerance=1e-10
    )))
    bt_non_numeric_err = tryCatch(
        {
            get_bootstrap_table(
                pcm_bt,
                list(detection.rate=c("0.1", "bad", "0.3", "0.4")),
                mode="l1ou"
            )
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(bt_non_numeric_err))
    expect_true(grepl("must be numeric or coercible to numeric", conditionMessage(bt_non_numeric_err), fixed=TRUE))
    bt_err = tryCatch(
        {
            get_bootstrap_table(pcm_bt, list(detection.rate=1:5), mode="l1ou")
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(bt_err))
    expect_true(grepl("Length mismatch", conditionMessage(bt_err), fixed=TRUE))
})
