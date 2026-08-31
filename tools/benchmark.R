# Measure the requested source checkout, never a globally installed rkftools.
arguments <- commandArgs(trailingOnly=TRUE)
options <- list(source=".", output="benchmark/results.csv", sizes="200,400,800",
    repetitions="3", baseline="", `max-ratio`="")
for (argument in arguments) {
    if (!grepl("^--[^=]+=", argument)) stop("Expected --name=value: ", argument)
    name <- sub("=.*$", "", sub("^--", "", argument))
    if (!name %in% names(options)) stop("Unknown benchmark option: ", name)
    options[[name]] <- sub("^--[^=]+=", "", argument)
}
sizes <- as.integer(strsplit(options$sizes, ",", fixed=TRUE)[[1]])
repetitions <- as.integer(options$repetitions)
if (anyNA(sizes) || any(sizes < 4L) || length(repetitions) != 1L ||
    is.na(repetitions) || repetitions < 2L) stop("Use sizes >= 4 and at least two repetitions.")
if (!grepl("[.]csv$", options$output)) stop("--output must end in .csv")
driver <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][[1]])
source(file.path(dirname(driver), "load-local.R"))
metadata <- load_local_package(options$source)
cases <- list()
add_case <- function(id, fun, check, tips, shape) {
    cases[[id]] <<- list(fun=fun, check=check, tips=tips, shape=shape)
}
tip_side <- function(tree, node) {
    if (node <= length(tree$tip.label)) tree$tip.label[[node]] else ape::extract.clade(tree, node)$tip.label
}
same_distances <- function(tree1, tree2) {
    tips <- sort(tree1$tip.label)
    stopifnot(setequal(tips, tree2$tip.label), isTRUE(all.equal(
        ape::cophenetic.phylo(tree1)[tips,tips], ape::cophenetic.phylo(tree2)[tips,tips], tolerance=1e-8)))
}
make_mapping_and_table_cases <- function(n, shape) {
    set.seed(20260831L + n)
    tree <- if (shape == "random") ape::rtree(n) else ape::stree(n, type="left")
    if (is.null(tree$edge.length)) tree$edge.length <- rep(1, nrow(tree$edge))
    target <- suppressWarnings(phytools::reroot(tree, tree$edge[nrow(tree$edge),2],
        pos=tree$edge.length[[nrow(tree$edge)]] / 2))
    root <- setdiff(target$edge[,1], target$edge[,2])
    reference <- tip_side(target, target$edge[which(target$edge[,1] == root)[[1]],2])
    add_case(paste("root_mapping", shape, n, sep="/"),
        function() get_phy2_root_in_phy1(tree, target, mode="index"),
        function(result) {
            stopifnot(length(result) == 1L, !is.na(result), result %in% seq_len(nrow(tree$edge)))
            side <- tip_side(tree, tree$edge[result,2])
            stopifnot(setequal(side, reference) || setequal(setdiff(tree$tip.label, side), reference))
        }, n, shape)
    add_case(paste("phylo2table", shape, n, sep="/"), function() phylo2table(tree),
        function(result) same_distances(tree, table2phylo(result, "label", "dist")), n, shape)
}
for (n in sizes) for (shape in c("random", "left")) make_mapping_and_table_cases(n, shape)
make_padding_case <- function(shape) {
    tree <- ape::stree(128L, type=shape)
    tree$edge.length <- ifelse(tree$edge[,2] <= 128L, 0.1, 1)
    before <- ape::node.depth.edgelength(tree)[1:128]
    add_case(paste("padding", shape, 128, sep="/"),
        function() pad_short_edges(tree, 0.2, external_only=TRUE),
        function(result) {
            delta <- ape::node.depth.edgelength(result)[1:128] - before
            stopifnot(identical(result$edge, tree$edge), min(result$edge.length[result$edge[,2] <= 128]) >= 0.2,
                isTRUE(all.equal(delta, rep(delta[[1]], 128))))
        }, 128L, shape)
}
for (shape in c("balanced", "left")) make_padding_case(shape)
make_mad_case <- function(n) {
    set.seed(20260831L + n)
    tree <- ape::rtree(n)
    add_case(paste("MAD", "random", n, sep="/"), function() MAD(tree, "custom"),
        function(result) {
            stopifnot(length(result[[5]]) == nrow(result[[3]]$edge),
                all(result[[5]][result[[4]]] == min(result[[5]], na.rm=TRUE)))
            for (rooted in result[[6]]) same_distances(tree, rooted)
        }, n, "random")
}
for (n in c(12L, 24L)) make_mad_case(n)
make_overlap_case <- function(n) {
    tree <- ape::stree(n, type="star")
    tree$edge.length <- rep(1, nrow(tree$edge))
    tree$root.edge <- 0
    tree$tip.label <- paste0("S", seq_len(n), "_sp_g")
    add_case(paste("root_overlap", "star", n, sep="/"),
        function() get_root_position_dependent_species_overlap_scores(tree),
        function(result) stopifnot(length(result) == nrow(tree$edge), all(result == 0)), n, "star")
}
for (n in c(50L, 100L, 200L)) make_overlap_case(n)

results <- outputs <- samples <- list()
for (id in names(cases)) {
    case <- cases[[id]]
    case$check(case$fun()) # Untimed warm-up and independent correctness check.
    elapsed <- heap <- numeric(repetitions)
    for (iteration in seq_len(repetitions)) {
        gc(reset=TRUE)
        elapsed[[iteration]] <- system.time(value <- case$fun(), gcFirst=FALSE)[["elapsed"]]
        memory <- gc()
        heap[[iteration]] <- sum(memory[,ncol(memory)])
    }
    case$check(value)
    outputs[[id]] <- value
    results[[id]] <- data.frame(case=id, tips=case$tips, shape=case$shape,
        repetitions=repetitions, median_seconds=median(elapsed),
        min_seconds=min(elapsed), max_seconds=max(elapsed), peak_r_heap_mib=max(heap))
    samples[[id]] <- data.frame(case=id, repetition=seq_len(repetitions),
        elapsed_seconds=elapsed, peak_r_heap_mib=heap)
    message(sprintf("%-32s %.3fs median; %.1f MiB peak R heap", id, median(elapsed), max(heap)))
}
result <- do.call(rbind, results)
rownames(result) <- NULL
result$version <- metadata$version
result$commit <- metadata$commit
result$platform <- metadata$platform
result$r_version <- as.character(getRversion())
dir.create(dirname(options$output), recursive=TRUE, showWarnings=FALSE)
utils::write.csv(result, options$output, row.names=FALSE)
utils::write.csv(do.call(rbind, samples), sub("[.]csv$", ".samples.csv", options$output), row.names=FALSE)
saveRDS(outputs, sub("[.]csv$", ".outputs.rds", options$output))
saveRDS(metadata, sub("[.]csv$", ".metadata.rds", options$output))

if (nzchar(options$baseline)) {
    baseline <- utils::read.csv(options$baseline, stringsAsFactors=FALSE)
    expected <- readRDS(sub("[.]csv$", ".outputs.rds", options$baseline))
    if (!setequal(names(outputs), names(expected))) stop("Baseline cases differ from this run.")
    for (id in names(outputs)) {
        difference <- all.equal(outputs[[id]], expected[[id]], tolerance=1e-8)
        if (!isTRUE(difference)) stop("Output changed for ", id, ": ", paste(difference, collapse="; "))
    }
    comparison <- merge(result, baseline[c("case", "median_seconds")], by="case", suffixes=c("", "_before"))
    comparison$ratio <- comparison$median_seconds / pmax(comparison$median_seconds_before, 0.001)
    utils::write.csv(comparison, sub("[.]csv$", ".comparison.csv", options$output), row.names=FALSE)
    message("Every benchmark output matches the baseline.")
    if (nzchar(options[["max-ratio"]])) {
        maximum <- as.numeric(options[["max-ratio"]])
        if (!is.finite(maximum) || maximum <= 1) stop("--max-ratio must be greater than 1.")
        slow <- comparison$ratio > maximum & comparison$median_seconds > 0.05
        if (any(slow)) stop("Performance regression: ", paste(comparison$case[slow], collapse=", "))
    }
}
print(result, row.names=FALSE)
summary <- Sys.getenv("GITHUB_STEP_SUMMARY")
if (nzchar(summary)) {
    cat("## Benchmark\n\n```\n", paste(capture.output(print(result, row.names=FALSE)), collapse="\n"),
        "\n```\n\nHeap columns report R's allocation high-water marks, not process RSS.\n",
        file=summary, append=TRUE)
}
