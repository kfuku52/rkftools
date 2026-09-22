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

test_that("bootstrap modes are validated before reading model inputs", {
    expect_identical(
        get_bootstrap_table(list(tree=fixture_trait_tree()), list(detection.rate=1:4)),
        get_bootstrap_table(list(tree=fixture_trait_tree()), list(detection.rate=1:4), "l1ou")
    )
    expect_error(get_bootstrap_table(NULL, NULL, "other"),
        'mode must be one of: "l1ou".', fixed=TRUE)
    for (mode in list(NULL, NA_character_, c("l1ou", "l1ou"), 1)) {
        expect_error(get_bootstrap_table(NULL, NULL, mode),
            "mode must be a single non-missing string.", fixed=TRUE)
    }
    expect_error(get_bootstrap_table(NULL, NULL, " "),
        "mode must be a single non-empty string.", fixed=TRUE)
})

test_that("bootstrap tables preserve missing support on a star tree", {
    tree = ape::read.tree(text="(A:1,B:1,C:1)root;")
    result = get_bootstrap_table(list(tree=tree), list(detection.rate=c(0, NA, 1)))
    expect_identical(result, data.frame(
        node_name=c("A", "B", "C", "root"),
        bootstrap_support=c(0, NA_real_, 1, NA_real_)
    ))
})
