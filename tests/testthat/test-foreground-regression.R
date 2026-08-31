test_that("independent foreground lineage counts", {
    withr::local_seed(20260831)
    tr_fg = ape::read.tree(text="((A:1,B:1):1,C:1);")
    trait_numeric = data.frame(species=c("A", "B", "C"), fg=c(1, 1, 0), stringsAsFactors=FALSE)
    fg_count = count_foreground_lineage(tr_fg, trait_numeric)
    expect_true(identical(as.integer(fg_count$fg), 1L))
    fg_unknown_err = tryCatch(
        {
            count_foreground_lineage(
                tr_fg,
                data.frame(species=c("A", "X", "C"), fg=c(1, 1, 0), stringsAsFactors=FALSE)
            )
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(fg_unknown_err))
    expect_true(grepl("not present in tree tip labels", conditionMessage(fg_unknown_err), fixed=TRUE))


    tr_fg_qualified = ape::read.tree(text="((Genus_species:1,Genus_cf_species:1):1,Other_species:1);")
    trait_qualified_numeric = data.frame(
        species=c("Genus_species", "Genus_cf_species", "Other_species"),
        fg=c(0, 1, 0),
        stringsAsFactors=FALSE
    )
    fg_count_qualified = count_foreground_lineage(tr_fg_qualified, trait_qualified_numeric)
    expect_true(identical(as.integer(fg_count_qualified$fg), 1L))


    fg_na = suppressWarnings(count_foreground_lineage(
        tr_fg,
        data.frame(species=c("A", "B", "C"), fg=c(1, NA, 0), stringsAsFactors=FALSE)
    ))
    expect_true(is.na(fg_na$fg))
    expect_error_contains(
        remove_invariant_traits(data.frame(a=c("x", "y"))),
        "numeric or logical trait columns"
    )
    expect_error_contains(
        remove_invariant_traits(data.frame(a=c(1, Inf))),
        "non-finite trait values"
    )
})
