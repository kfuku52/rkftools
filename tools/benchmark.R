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

benchmark_root_overlap <- function(n_tip) {
    tree <- ape::stree(n_tip, type="star")
    tree$edge.length <- rep(1, nrow(tree$edge))
    tree$root.edge <- 0
    tree$tip.label <- paste0("S", seq_len(n_tip), "_sp_g")
    elapsed <- system.time({
        result <- get_root_position_dependent_species_overlap_scores(tree)
    })[["elapsed"]]
    stopifnot(length(result) == nrow(tree$edge))
    data.frame(
        benchmark="get_root_position_dependent_species_overlap_scores",
        tips=n_tip,
        edges=nrow(tree$edge),
        elapsed_seconds=unname(elapsed)
    )
}

result <- rbind(
    do.call(rbind, lapply(c(200L, 400L, 800L), benchmark_root_mapping)),
    do.call(rbind, lapply(c(50L, 100L, 200L), benchmark_root_overlap))
)
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
