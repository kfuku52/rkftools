test_that("mapping a fully collapsed tree ignores its synthetic root", {
    tree <- ape::read.tree(text="((A:1,B:1):1,C:2);")
    traits <- data.frame(value=c(1, 2, 3), row.names=tree$tip.label)
    collapsed <- collapse_clades(tree, traits, get_root_num(tree))
    result <- map_node_num(tree, collapsed$tree, collapsed$collapse_leaf_names)
    expect_identical(result$tree_original, seq_len(5L))
    expect_identical(result$tree_collapsed, rep(1L, 5L))
})

test_that("unary nodes map in ancestor order without conflating tips", {
    tree <- ape::read.tree(text="(((A:1):1,B:2):1,C:3);")
    result <- map_node_num(tree, tree)
    expect_identical(result$tree_original, result$tree_collapsed)
    reordered <- ape::reorder.phylo(tree, "postorder")
    expect_identical(map_node_num(tree, reordered), result)
    clade <- get_parent_num(tree, match("A", tree$tip.label))
    # A named collapse map also distinguishes an internal unary node from A.
    collapsed <- tree
    collapsed$tip.label[[1]] <- as.character(clade)
    collapsed <- ape::collapse.singles(collapsed)
    result <- map_node_num(tree, collapsed, setNames(list("A"), as.character(clade)))
    expect_equal(result$tree_collapsed[result$tree_original %in% c(1L, clade)], c(1L, 1L))
})

test_that("regime restoration uses original names for every replacement", {
    original <- ape::read.tree(text="((A:1,B:1)7:1,(C:1,D:1)Z:1)R;")
    collapsed <- collapse_clades(original,
        data.frame(value=rep(1,4), row.names=original$tip.label), c(6L,7L))
    mapping <- map_node_num(original, collapsed$tree, collapsed$collapse_leaf_names)
    tab <- data.frame(node_name=c("6", "7", NA), param="shift_value", value=c(10,20,30))
    result <- regime_table_collapse2original(tab, original, collapsed$tree, mapping)
    expect_identical(result$node_name, c("7", "Z", NA_character_))
    expect_identical(result$value, tab$value)
    # Swapping two labels must also be simultaneous.
    swapped <- original
    swapped$node.label <- c("7", "R", "Z")
    tab <- data.frame(node_name=c("R", "7"), value=c(1,2))
    mapping <- data.frame(tree_original=1:7, tree_collapsed=1:7)
    expect_identical(regime_table_collapse2original(tab, original, swapped, mapping)$node_name,
        c("7", "R"))
})
