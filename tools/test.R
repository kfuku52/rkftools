source("tools/load-local.R")
invisible(load_local_package())
testthat::test_local(".", stop_on_warning=TRUE)
