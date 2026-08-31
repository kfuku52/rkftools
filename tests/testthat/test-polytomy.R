make_internal_polytomy <- function() {
    ape::read.tree(text="((A_a_g1:1,A_a_g2:1,B_b_g3:1)P:2,C_c_g4:3)R;")
}

make_root_polytomy <- function() {
    tree <- ape::read.tree(text="(A_a_g1:1,A_a_g2:1,B_b_g3:1,C_c_g4:1)R;")
    tree$root.edge <- 0
    tree
}

expect_same_phylo <- function(actual, expected) {
    expect_true(isTRUE(ape::all.equal.phylo(
        actual,
        expected,
        use.edge.length=TRUE,
        use.tip.label=TRUE
    )))
}

test_that("multifurcations round-trip through branch tables", {
    for (tree in list(make_internal_polytomy(), make_root_polytomy())) {
        branch_table <- phylo2table(tree)
        restored <- table2phylo(branch_table, "label", "dist")

        expect_same_phylo(restored, tree)
        expect_true(ape::is.rooted(restored))
        expect_no_error(phylo2table(restored))
        child_counts <- table(branch_table$parent[branch_table$parent != -999])
        multifurcating_parent <- names(child_counts)[child_counts > 2L]
        expect_true(length(multifurcating_parent) >= 1L)
        expect_true(all(
            branch_table$sister[
                as.character(branch_table$parent) %in% multifurcating_parent
            ] == -999L
        ))
    }
})

test_that("species overlap uses maximum pairwise child overlap", {
    tree <- make_internal_polytomy()
    poly_node <- as.integer(names(which(table(tree$edge[,1]) > 2L))[[1]])

    expect_equal(get_duplication_confidence_score(tree, poly_node), 1)
    expect_equal(get_species_overlap_score(tree, dc_cutoff=0.5), 1)

    reordered <- ape::rotate(tree, poly_node)
    expect_equal(
        get_species_overlap_score(reordered, dc_cutoff=0.5),
        get_species_overlap_score(tree, dc_cutoff=0.5)
    )
    root_scores <- get_root_position_dependent_species_overlap_scores(tree)
    expect_length(root_scores, nrow(tree$edge))
    expect_true(all(is.finite(root_scores)))

    star <- ape::stree(60, type="star")
    star$edge.length <- rep(1, nrow(star$edge))
    star$root.edge <- 0
    star$tip.label <- c(
        "A_sp_g1",
        "A_sp_g2",
        paste0("S", 3:60, "_sp_g", 3:60)
    )
    expect_equal(get_species_overlap_score(star), 1)
    expect_equal(
        get_root_position_dependent_species_overlap_scores(star),
        rep(1, nrow(star$edge))
    )
})

test_that("root and node mapping support multifurcations", {
    root_poly <- make_root_polytomy()
    reordered <- ape::rotate(root_poly, get_root_num(root_poly))

    expect_identical(
        get_phy2_root_in_phy1(reordered, root_poly),
        get_root_num(reordered)
    )
    expect_true(is.na(get_phy2_root_in_phy1(
        reordered, root_poly, mode="index"
    )))

    search_tree <- ape::read.tree(
        text="(A_a_g1:1,(A_a_g2:1,(B_b_g3:1,C_c_g4:1):1)X:1)R;"
    )
    target_tree <- ape::read.tree(
        text="(A_a_g1:1,A_a_g2:1,(B_b_g3:1,C_c_g4:1):1)T;"
    )
    target_tree$root.edge <- 0
    matched_node <- get_node_num_by_name(search_tree, "X")
    expect_identical(
        get_phy2_root_in_phy1(search_tree, target_tree),
        matched_node
    )
    expect_identical(
        get_phy2_root_in_phy1(search_tree, target_tree, mode="index"),
        which(search_tree$edge[,2] == matched_node)[[1]]
    )

    resolved <- ape::multi2di(root_poly, random=FALSE)
    mapping <- multi2bi_node_number_transfer(root_poly, resolved)
    expect_true(any(
        mapping$mtree_node == get_root_num(root_poly) &
            mapping$btree_node == get_root_num(resolved)
    ))
})

test_that("branch length transformations preserve multifurcating topology", {
    tree <- ape::read.tree(text="((A:0.1,B:0.7,C:0.8)P:2,D:2.8)R;")
    before_depth <- ape::node.depth.edgelength(tree)[seq_along(tree$tip.label)]
    padded <- pad_short_edges(tree, threshold=0.5, external_only=TRUE)
    after_depth <- ape::node.depth.edgelength(padded)[seq_along(tree$tip.label)]

    expect_false(ape::is.binary(padded))
    expect_equal(after_depth, before_depth)
    expect_gte(min(padded$edge.length[padded$edge[,2] <= length(tree$tip.label)]), 0.5)

    root_three <- ape::read.tree(text="(A:0.1,B:0.6,C:0.8)R;")
    root_three$root.edge <- 0
    padded_root <- pad_short_edges(
        root_three,
        threshold=0.5,
        external_only=TRUE
    )
    expect_gte(min(padded_root$edge.length), 0.5)
    expect_equal(
        padded_root$edge.length - root_three$edge.length,
        rep(0.4, 3)
    )

    non_ultrametric <- ape::read.tree(text="((A:1,B:1,C:1)P:1,D:3)R;")
    ultrametric <- force_ultrametric(
        non_ultrametric,
        stop_if_larger_change=10
    )
    expect_true(ape::is.ultrametric(ultrametric))
    expect_false(ape::is.binary(ultrametric))
})

test_that("unrooted binary trees retain binary transformation behavior", {
    tree <- ape::read.tree(text="((A:0.1,B:0.7):0.2,C:0.8,D:0.9);")
    expect_true(ape::is.binary(tree))
    expect_false(ape::is.rooted(tree))
    expect_false(contains_polytomy(tree))

    padded <- pad_short_edges(tree, threshold=0.5, external_only=TRUE)
    expect_equal(padded$edge.length, c(0, 0.5, 1.1, 1, 1.1))
    expect_equal(ape::node.depth.edgelength(padded)[1:4] -
        ape::node.depth.edgelength(tree)[1:4], rep(0.2, 4))

    expected_ultrametric <- ape::chronoMPL(tree)
    actual_ultrametric <- force_ultrametric(
        tree,
        stop_if_larger_change=100
    )
    expect_same_phylo(actual_ultrametric, expected_ultrametric)

    species_tree <- tree
    species_tree$tip.label <- c(
        "A_sp_g1", "A_sp_g2", "B_sp_g3", "C_sp_g4"
    )
    expect_true(is.na(get_duplication_confidence_score(
        species_tree,
        get_root_num(species_tree)
    )))
    expect_equal(get_species_overlap_score(species_tree), 1)

    unrooted_star <- ape::stree(4, type="star")
    unrooted_star$edge.length <- rep(1, nrow(unrooted_star$edge))
    expect_true(contains_polytomy(unrooted_star))
})

test_that("collapsed multifurcations use a phylogenetic root estimate", {
    tree <- ape::read.tree(text="((A:1,B:2,C:4)P:1,D:5)R;")
    poly_node <- as.integer(names(which(table(tree$edge[,1]) > 2L))[[1]])
    traits <- data.frame(value=c(1, 2, 8, 0), row.names=tree$tip.label)

    expect_warning(
        collapsed <- collapse_clades(tree, traits, poly_node),
        NA
    )
    expected_gls <- sum(c(1, 2, 8) / c(1, 2, 4)) / sum(1 / c(1, 2, 4))
    expect_equal(collapsed$trait[as.character(poly_node), "value"], expected_gls)
})

test_that("multifurcating GLS retains short and zero-length information", {
    for (short_length in c(1e-10, 0)) {
        tree <- ape::read.tree(text=paste0(
            "((A:", format(short_length, scientific=TRUE),
            ",B:1,C:1)P:1,D:2)R;"
        ))
        poly_node <- get_node_num_by_name(tree, "P")
        traits <- data.frame(
            value=c(10, 0, 0, 0),
            row.names=tree$tip.label
        )
        expect_warning(
            collapsed <- collapse_clades(tree, traits, poly_node),
            NA
        )
        expected <- if (short_length == 0) {
            10
        } else {
            sum(c(10, 0, 0) / c(short_length, 1, 1)) /
                sum(1 / c(short_length, 1, 1))
        }
        expect_equal(
            collapsed$trait[as.character(poly_node), "value"],
            expected,
            tolerance=1e-7
        )
    }
})

test_that("MAD scores multifurcations without random resolution", {
    tree <- make_internal_polytomy()

    set.seed(1)
    result1 <- MAD(tree, "full")
    set.seed(999)
    result2 <- MAD(tree, "full")
    parallel_result <- MAD_parallel(tree, "full", ncpu=1)

    expect_equal(result1[[1]], result2[[1]])
    expect_equal(result1[[2]], result2[[2]])
    expect_equal(result1[[1]], parallel_result[[1]])
    expect_equal(result1[[2]], parallel_result[[2]])
    expect_true(contains_polytomy(result1[[3]]))
    expect_true(all(vapply(
        result1[[6]],
        contains_polytomy,
        logical(1)
    )))
})
