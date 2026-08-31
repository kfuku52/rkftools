test_that("node lookup, traversal, and ages", {
    withr::local_seed(20260831)
    tr_unlabeled = ape::read.tree(text="((A:1,B:1):1,C:1);")
    expect_true(is.null(tr_unlabeled$node.label))
    tr_filled = fill_node_labels(tr_unlabeled)
    expect_true(length(tr_filled$node.label) == tr_filled$Nnode)
    expect_true(!any(is.na(tr_filled$node.label) | tr_filled$node.label == ""))
    root_num = get_root_num(tr_unlabeled)
    expect_true(identical(as.integer(root_num), 4L))
    expect_true(isTRUE(is_root(tr_unlabeled, root_num)))
    expect_true(isTRUE(is_leaf(tr_unlabeled, 1L)))
    expect_true(isFALSE(is_leaf(tr_unlabeled, root_num)))
    all_node_names = c(tr_filled$tip.label, tr_filled$node.label)
    mapped_nums = get_node_num_by_name(tr_filled, c("C", "A", "C"))
    expect_true(identical(as.integer(mapped_nums), c(3L, 1L, 3L)))
    mapped_nums_na = get_node_num_by_name(tr_filled, c("C", NA_character_, "X", "A"))
    expect_true(identical(as.integer(mapped_nums_na), c(3L, 1L)))
    mapped_names = get_node_name_by_num(tr_filled, c(5L, 4L, 5L))
    expect_true(identical(as.character(mapped_names), as.character(all_node_names[c(5L, 4L, 5L)])))
    mapped_names_na = get_node_name_by_num(tr_filled, c(5L, NA_integer_, 4L, 999L))
    expect_true(identical(as.character(mapped_names_na), as.character(all_node_names[c(5L, 4L)])))
    internal_nodes = sort(unique(tr_unlabeled$edge[,1]))
    internal_nodes = internal_nodes[internal_nodes > length(tr_unlabeled$tip.label)]
    size_two_node = internal_nodes[vapply(internal_nodes, function(nn) {
        length(get_tip_labels(tr_unlabeled, nn)) == 2
    }, logical(1))][1]
    tip_labels_vec = get_tip_labels(tr_unlabeled, c(3L, size_two_node))
    expect_true(identical(as.character(tip_labels_vec), c("C", "A", "B")))
    tip_labels_invalid_err = tryCatch(
        {
            get_tip_labels(tr_unlabeled, 999L)
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(tip_labels_invalid_err))
    expect_true(grepl("outside valid node range", conditionMessage(tip_labels_invalid_err), fixed=TRUE))
    children_vec = get_children_num(tr_unlabeled, c(4L, 5L))
    expect_true(identical(as.integer(children_vec), c(5L, 3L, 1L, 2L)))
    children_na = get_children_num(tr_unlabeled, NA_integer_)
    expect_true(identical(as.integer(children_na), integer(0)))
    parent_vec = get_parent_num(tr_unlabeled, c(1L, 2L, 3L))
    expect_true(identical(as.integer(parent_vec), c(5L, 5L, 4L)))
    parent_na = get_parent_num(tr_unlabeled, NA_integer_)
    expect_true(identical(as.integer(parent_na), integer(0)))
    sister_vec = get_sister_num(tr_unlabeled, c(1L, 2L))
    expect_true(identical(as.integer(sister_vec), c(2L, 1L)))
    sister_na = get_sister_num(tr_unlabeled, NA_integer_)
    expect_true(identical(as.integer(sister_na), integer(0)))
    ancestor_a = get_ancestor_num(tr_unlabeled, 1L)
    expect_true(identical(as.integer(ancestor_a), c(5L, 4L)))
    ancestor_invalid = get_ancestor_num(tr_unlabeled, 99L)
    expect_true(length(ancestor_invalid) == 0)
    desc_leaf = get_descendent_num(tr_unlabeled, c(4L, 5L), leaf_only=TRUE)
    expect_true(identical(as.integer(desc_leaf), c(1L, 2L, 3L)))
    desc_leaf_only_na_err = tryCatch(
        {
            get_descendent_num(tr_unlabeled, c(4L, 5L), leaf_only=NA)
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(desc_leaf_only_na_err))
    expect_true(grepl("leaf_only must be a single non-missing logical value", conditionMessage(desc_leaf_only_na_err), fixed=TRUE))
    mrca_tr_unlabeled = ape::mrca(tr_unlabeled)
    nearest_out = get_nearest_tips(
        tr_unlabeled,
        query="A",
        subjects=c("B", "C"),
        mrca_matrix=mrca_tr_unlabeled
    )
    expect_true(identical(as.character(nearest_out$nearests), "B"))
    nearest_query_err = tryCatch(
        {
            get_nearest_tips(tr_unlabeled, query="X", subjects=c("B", "C"), mrca_matrix=mrca_tr_unlabeled)
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(nearest_query_err))
    expect_true(grepl("query must be present in mrca_matrix row names", conditionMessage(nearest_query_err), fixed=TRUE))
    nearest_subject_err = tryCatch(
        {
            get_nearest_tips(tr_unlabeled, query="A", subjects=c("B", "X"), mrca_matrix=mrca_tr_unlabeled)
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(nearest_subject_err))
    expect_true(grepl("subjects are missing in mrca_matrix column names", conditionMessage(nearest_subject_err), fixed=TRUE))
    nearest_query_vec_err = tryCatch(
        {
            get_nearest_tips(tr_unlabeled, query=c("A", "B"), subjects=c("B", "C"), mrca_matrix=mrca_tr_unlabeled)
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(nearest_query_vec_err))
    expect_true(grepl("query must be a single non-missing string", conditionMessage(nearest_query_vec_err), fixed=TRUE))
    nearest_subject_empty_err = tryCatch(
        {
            get_nearest_tips(tr_unlabeled, query="A", subjects=character(0), mrca_matrix=mrca_tr_unlabeled)
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(nearest_subject_empty_err))
    expect_true(grepl("subjects must contain at least one non-missing tip label", conditionMessage(nearest_subject_empty_err), fixed=TRUE))
    ultra_tree = ape::compute.brlen(ape::stree(4), 1)
    ultra_tree = force_ultrametric(ultra_tree, stop_if_larger_change=1)
    tip_age = get_node_age(ultra_tree, 1)
    expect_true(isTRUE(all.equal(as.numeric(tip_age), 0, tolerance=1e-10)))
    node_age_err = tryCatch(
        {
            get_node_age(ultra_tree, 999)
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(node_age_err))
    expect_true(grepl("must be a single integer", conditionMessage(node_age_err), fixed=TRUE))
})
