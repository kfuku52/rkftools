test_that("NOTUNG records and command-line parsing", {
    withr::local_seed(20260831)
    notung_tmp = tempfile()
    writeLines(c("#D g1 s1 s2"), notung_tmp)
    notung_df1 = read_notung_parsable(notung_tmp, mode="D")
    expect_true(nrow(notung_df1) == 1)
    expect_true(identical(as.character(notung_df1$gn_node), "g1"))
    writeLines(c("#D g1 s1 s2", "#D g2 s3 s4"), notung_tmp)
    notung_df2 = read_notung_parsable(notung_tmp, mode="D")
    expect_true(nrow(notung_df2) == 2)
    expect_true(identical(as.character(notung_df2$gn_node), c("g1", "g2")))
    writeLines(c("  #D g3 s5 s6"), notung_tmp)
    notung_df3 = read_notung_parsable(notung_tmp, mode="D")
    expect_true(nrow(notung_df3) == 1)
    expect_true(identical(as.character(notung_df3$gn_node), "g3"))
    notung_mode_na_err = tryCatch(
        {
            read_notung_parsable(notung_tmp, mode=NA_character_)
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(notung_mode_na_err))
    expect_true(grepl("mode must be a single non-missing string", conditionMessage(notung_mode_na_err), fixed=TRUE))
    notung_mode_vec_err = tryCatch(
        {
            read_notung_parsable(notung_tmp, mode=c("D", "X"))
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(notung_mode_vec_err))
    expect_true(grepl("mode must be a single non-missing string", conditionMessage(notung_mode_vec_err), fixed=TRUE))
    unlink(notung_tmp)


    expect_true(isTRUE(is.blank(c("", NA_character_))))
    expect_true(isFALSE(is.blank(c("", "x"))))
    is_blank_false_trigger_err = tryCatch(
        {
            is.blank(c("", NA_character_), false.triggers=NA)
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(is_blank_false_trigger_err))
    expect_true(grepl("false.triggers must be a single non-missing logical value", conditionMessage(is_blank_false_trigger_err), fixed=TRUE))
    parsed_args = get_parsed_args(c("--threads=2", "--dry-run"), print=FALSE)
    expect_true(identical(as.numeric(parsed_args[["threads"]]), 2))
    expect_true(isTRUE(parsed_args[["dry-run"]]))
    parsed_arg_err = tryCatch(
        {
            get_parsed_args(c("--=1"), print=FALSE)
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(parsed_arg_err))
    expect_true(grepl("parameter name is empty", conditionMessage(parsed_arg_err), fixed=TRUE))
    parsed_arg_na_err = tryCatch(
        {
            get_parsed_args(c(NA_character_), print=FALSE)
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(parsed_arg_na_err))
    expect_true(grepl("single non-missing string", conditionMessage(parsed_arg_na_err), fixed=TRUE))
    parsed_print_na_err = tryCatch(
        {
            get_parsed_args(c("--threads=2"), print=NA)
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(parsed_print_na_err))
    expect_true(grepl("print must be a single non-missing logical value", conditionMessage(parsed_print_na_err), fixed=TRUE))

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
    expect_error_contains(get_parsed_args("threads=2"), 'expected a "--name"')
    expect_error_contains(get_parsed_args(c("--threads=1", "--threads=2")), "Duplicate long argument")
})
