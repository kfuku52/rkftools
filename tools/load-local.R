load_local_package <- function(path=".") {
    dev_library <- Sys.getenv("RKFTOOLS_DEV_LIBRARY", file.path(path, ".local", "R-library"))
    if (dir.exists(dev_library)) .libPaths(c(normalizePath(dev_library), .libPaths()))
    if (!requireNamespace("pkgload", quietly=TRUE)) {
        stop("Install development dependencies with make setup.", call.=FALSE)
    }
    path <- normalizePath(path, mustWork=TRUE)
    pkgload::load_all(path, quiet=TRUE, export_all=FALSE, helpers=FALSE)
    namespace_path <- normalizePath(getNamespaceInfo(asNamespace("rkftools"), "path"))
    if (!identical(namespace_path, path)) stop("Loaded package is not the requested source tree.")
    commit <- Sys.getenv("RKFTOOLS_BENCHMARK_COMMIT", "unknown")
    dirty <- NA
    if (file.exists(file.path(path, ".git")) && nzchar(Sys.which("git"))) {
        commit <- system2("git", c("-C", shQuote(path), "rev-parse", "HEAD"), stdout=TRUE)
        dirty <- length(system2("git", c("-C", shQuote(path), "status", "--porcelain"), stdout=TRUE)) > 0L
    }
    metadata <- list(version=as.character(utils::packageVersion("rkftools")),
        commit=commit, dirty=dirty, source=namespace_path, R=R.version.string,
        platform=R.version$platform)
    message("rkftools ", metadata$version, " | ", commit,
        if (isTRUE(dirty)) " (working tree changes)" else "", " | ", namespace_path)
    metadata
}
