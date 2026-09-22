source("tools/load-local.R")
invisible(load_local_package())
filter <- Sys.getenv("RKFTOOLS_TEST_FILTER", "")
testthat::test_local(".", filter=if (nzchar(filter)) filter else NULL,
                     stop_on_failure=TRUE, stop_on_warning=TRUE)
