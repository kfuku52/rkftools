args <- commandArgs(trailingOnly=TRUE)
if (length(setdiff(args, "--minimal"))) stop("Usage: Rscript tools/setup.R [--minimal]")
dev_library <- Sys.getenv("RKFTOOLS_DEV_LIBRARY", ".local/R-library")
dir.create(dev_library, recursive=TRUE, showWarnings=FALSE)
.libPaths(c(normalizePath(dev_library), .libPaths()))
description <- read.dcf("DESCRIPTION")
packages <- trimws(unlist(strsplit(paste(description[1,c("Imports", "Suggests")], collapse=","), ",")))
packages <- sub("[[:space:]]*[(].*$", "", packages)
if ("--minimal" %in% args) packages <- setdiff(packages, c("PhylogeneticEM", "Rphylopars"))
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
if (requireNamespace("testthat", quietly=TRUE) && utils::packageVersion("testthat") < "3.1.8") {
    missing <- union(missing, "testthat")
}
if (length(missing)) {
    repos <- getOption("repos")
    if (!length(repos)) repos <- c(CRAN="https://cloud.r-project.org")
    repos[repos == "@CRAN@"] <- "https://cloud.r-project.org"
    utils::install.packages(missing, lib=dev_library, repos=repos)
}
unavailable <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
if (length(unavailable)) stop("Dependency installation failed: ", paste(unavailable, collapse=", "))
message("Development library: ", normalizePath(dev_library))
message("Ready for ", if ("--minimal" %in% args) "minimal" else "full", " checks.")
