library(rkftools)

benchmark_root_mapping <- function(n_tip) {
    set.seed(20260810 + n_tip)
    tree <- ape::rtree(n_tip)
    rerooted <- suppressWarnings(phytools::reroot(
        tree,
        node.number=tree$edge[nrow(tree$edge), 2]
    ))
    elapsed <- system.time({
        result <- get_phy2_root_in_phy1(tree, rerooted, nslots=1L)
    })[["elapsed"]]
    stopifnot(!is.na(result))
    data.frame(
        benchmark="get_phy2_root_in_phy1",
        tips=n_tip,
        edges=nrow(tree$edge),
        elapsed_seconds=unname(elapsed)
    )
}

result <- do.call(rbind, lapply(c(200L, 400L, 800L), benchmark_root_mapping))
print(result, row.names=FALSE)

summary_path <- Sys.getenv("GITHUB_STEP_SUMMARY")
if (nzchar(summary_path)) {
    write("## rkftools benchmark\n", file=summary_path, append=TRUE)
    write("```\n", file=summary_path, append=TRUE)
    suppressWarnings(write.table(
        result,
        file=summary_path,
        append=TRUE,
        row.names=FALSE,
        quote=FALSE
    ))
    write("```\n", file=summary_path, append=TRUE)
}
