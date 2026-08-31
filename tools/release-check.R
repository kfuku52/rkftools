args <- commandArgs(trailingOnly=TRUE)
expected_tag <- if (length(args)) args[[1]] else if (Sys.getenv("GITHUB_REF_TYPE") == "tag") Sys.getenv("GITHUB_REF_NAME") else ""

description <- read.dcf("DESCRIPTION")
version <- unname(description[[1, "Version"]])
expected_news <- paste0("# rkftools ", version)
news <- readLines("NEWS.md", warn=FALSE)
readme <- readLines("README.md", warn=FALSE)

errors <- character()
if (!length(news) || !identical(news[[1]], expected_news)) {
    errors <- c(errors, paste0("NEWS.md must start with: ", expected_news))
}
version_badge <- paste0("version-", version, "-informational")
if (!any(grepl(version_badge, readme, fixed=TRUE))) {
    errors <- c(errors, paste0("README version badge must contain: ", version_badge))
}
if (nzchar(expected_tag) && !identical(expected_tag, paste0("v", version))) {
    errors <- c(errors, paste0(
        "Release tag ", expected_tag, " does not match DESCRIPTION version v", version
    ))
}

if (length(errors)) {
    stop(paste(errors, collapse="\n"), call.=FALSE)
}
message("Release metadata is consistent for rkftools ", version, ".")
