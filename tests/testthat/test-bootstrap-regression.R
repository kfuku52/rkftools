test_that("bootstrap support is placed by node and numeric factors are decoded", {
    withr::local_seed(20260831)
    tr_unlabeled = fixture_trait_tree()
    pcm_bt = list(tree=tr_unlabeled)
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
    expect_error(
        get_bootstrap_table(
            pcm_bt,
            list(detection.rate=c("0.1", "bad", "0.3", "0.4")),
            mode="l1ou"
        ),
        "must be numeric or coercible to numeric", fixed=TRUE
    )
    expect_error(
        get_bootstrap_table(pcm_bt, list(detection.rate=1:5), mode="l1ou"),
        "Length mismatch", fixed=TRUE
    )
})
