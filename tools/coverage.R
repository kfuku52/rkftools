args <- commandArgs(trailingOnly=TRUE)
threshold <- if (length(args)) as.numeric(args[[1]]) else 80

if (length(threshold) != 1L || is.na(threshold) || !is.finite(threshold) ||
    threshold < 0 || threshold > 100) {
    stop("Coverage threshold must be a finite number between 0 and 100.",
         call.=FALSE)
}

coverage <- covr::package_coverage(type="tests", quiet=FALSE)
percentage <- as.numeric(covr::percent_coverage(coverage))

message(sprintf("Package coverage: %.2f%% (required: %.2f%%)",
                percentage, threshold))
if (percentage < threshold) {
    stop(sprintf("Coverage %.2f%% is below the %.2f%% threshold.",
                 percentage, threshold), call.=FALSE)
}
