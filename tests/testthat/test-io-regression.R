test_that("NOTUNG records and command-line parsing", {
    withr::local_seed(20260831)
    notung_tmp = tempfile()
    on.exit(unlink(notung_tmp), add=TRUE)
    writeLines(c("#D g1 s1 s2", "  #D g2 s3 s4"), notung_tmp)
    records = read_notung_parsable(notung_tmp, mode="D")
    expect_identical(as.character(records$gn_node), c("g1", "g2"))
    expect_error(
        read_notung_parsable(notung_tmp, mode=NA_character_),
        "mode must be a single non-missing string", fixed=TRUE
    )

    expect_true(isTRUE(is.blank(c("", NA_character_))))
    expect_true(isFALSE(is.blank(c("", "x"))))
    parsed_args = get_parsed_args(c("--threads=2", "--dry-run"), print=FALSE)
    expect_true(identical(as.numeric(parsed_args[["threads"]]), 2))
    expect_true(isTRUE(parsed_args[["dry-run"]]))
    expect_error(
        get_parsed_args(c("--=1"), print=FALSE),
        "parameter name is empty", fixed=TRUE
    )
    expect_error(
        get_parsed_args(c(NA_character_), print=FALSE),
        "single non-missing string", fixed=TRUE
    )

    # Command-line parsing is silent by default and never prints credential values.
    parsed_default_output = capture.output(parsed_default <- get_parsed_args(c("--threads=3")))
    expect_true(identical(parsed_default_output, character(0)))
    expect_true(identical(parsed_default$threads, 3))
    expect_true(identical(get_parsed_args("--sample=0012")$sample, "0012"))
    parsed_secret_output = capture.output(
        parsed_secret <- get_parsed_args(c("--api-key=do-not-print", "--threads=2"), print=TRUE)
    )
    expect_true(identical(parsed_secret[["api-key"]], "do-not-print"))
    expect_true(any(grepl("<redacted>", parsed_secret_output, fixed=TRUE)))
    expect_true(!any(grepl("do-not-print", parsed_secret_output, fixed=TRUE)))
    expect_error(get_parsed_args("threads=2"), 'expected a "--name"', fixed=TRUE)
    expect_error(get_parsed_args(c("--threads=1", "--threads=2")), "Duplicate long argument", fixed=TRUE)
})
