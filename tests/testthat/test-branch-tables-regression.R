test_that("branch table schema, distances, and graph validation", {
    withr::local_seed(20260831)
    tbl = data.frame(
        branch_id=c(1L, 2L, 3L),
        parent=c(3L, 3L, -999L),
        sister=c(2L, 1L, -999L),
        label=c("A", "B", "Root"),
        dist=c(0.1, 0.2, 0.0),
        stringsAsFactors=FALSE
    )
    tbl_phy = table2phylo(tbl, name_col="label", dist_col="dist")
    expect_true(inherits(tbl_phy, "phylo"))
    expect_true(setequal(tbl_phy$tip.label, c("A", "B")))
    tbl_phy_dist = cophenetic(tbl_phy)[c("A", "B"), c("A", "B")]
    tbl_phy_expected = matrix(c(0, 0.3, 0.3, 0), nrow=2, byrow=TRUE, dimnames=list(c("A", "B"), c("A", "B")))
    expect_true(isTRUE(all.equal(tbl_phy_dist, tbl_phy_expected, tolerance=1e-10)))

    tbl_phy2table = phylo2table(tbl_phy, name_col="label", dist_col="dist")
    expect_true(identical(as.character(colnames(tbl_phy2table)), c("branch_id", "parent", "sister", "label", "dist")))
    expect_true(!any(is.na(suppressWarnings(as.integer(tbl_phy2table$branch_id)))))
    expect_true(!any(is.na(suppressWarnings(as.integer(tbl_phy2table$parent)))))
    expect_true(!any(is.na(suppressWarnings(as.integer(tbl_phy2table$sister)))))
    expect_true(setequal(as.character(tbl_phy2table$label), c("A", "B", "Root")))
    tbl_phy_roundtrip = table2phylo(tbl_phy2table, name_col="label", dist_col="dist")
    expect_true(inherits(tbl_phy_roundtrip, "phylo"))
    expect_true(setequal(tbl_phy_roundtrip$tip.label, tbl_phy$tip.label))
    roundtrip_dist = cophenetic(tbl_phy_roundtrip)[c("A", "B"), c("A", "B")]
    expect_true(isTRUE(all.equal(roundtrip_dist, tbl_phy_expected, tolerance=1e-10)))

    tbl_readme = data.frame(
        branch_id=c(6L, 2L, 0L, 1L, 5L, 3L, 4L),
        parent=c(-999L, 6L, 2L, 2L, 6L, 5L, 5L),
        sister=c(-999L, 5L, 1L, 0L, 2L, 4L, 3L),
        node_name=c("Root", "Clade_AB", "A", "B", "Clade_CD", "C", "D"),
        dist=c(0, 0.42, 0.16, 0.18, 0.50, 0.22, 0.25),
        stringsAsFactors=FALSE
    )
    tbl_readme_phy = table2phylo(tbl_readme, name_col="node_name", dist_col="dist")
    tbl_readme_newick = ape::write.tree(tbl_readme_phy)
    tbl_readme_phy_from_newick = ape::read.tree(text=tbl_readme_newick)
    tbl_readme_dist = cophenetic(tbl_readme_phy)[c("A", "B", "C", "D"), c("A", "B", "C", "D")]
    tbl_readme_expected = matrix(
        c(
            0, 0.34, 1.30, 1.33,
            0.34, 0, 1.32, 1.35,
            1.30, 1.32, 0, 0.47,
            1.33, 1.35, 0.47, 0
        ),
        nrow=4,
        byrow=TRUE,
        dimnames=list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))
    )
    tbl_readme_newick_dist = cophenetic(tbl_readme_phy_from_newick)[c("A", "B", "C", "D"), c("A", "B", "C", "D")]
    expect_true(isTRUE(all.equal(tbl_readme_dist, tbl_readme_expected, tolerance=1e-10)))
    expect_true(isTRUE(all.equal(tbl_readme_newick_dist, tbl_readme_expected, tolerance=1e-10)))
    tbl_readme_roundtrip = phylo2table(tbl_readme_phy_from_newick, name_col="node_name", dist_col="dist")
    expect_true(identical(as.character(colnames(tbl_readme_roundtrip)), c("branch_id", "parent", "sister", "node_name", "dist")))
    expect_true(!any(is.na(suppressWarnings(as.integer(tbl_readme_roundtrip$branch_id)))))
    expect_true(!any(is.na(suppressWarnings(as.integer(tbl_readme_roundtrip$parent)))))
    expect_true(!any(is.na(suppressWarnings(as.integer(tbl_readme_roundtrip$sister)))))
    expect_true(isTRUE(all.equal(tbl_readme_roundtrip, tbl_readme, check.attributes=FALSE, tolerance=1e-10)))
    tbl_readme_phy_roundtrip = table2phylo(tbl_readme_roundtrip, name_col="node_name", dist_col="dist")
    tbl_readme_roundtrip_dist = cophenetic(tbl_readme_phy_roundtrip)[c("A", "B", "C", "D"), c("A", "B", "C", "D")]
    expect_true(isTRUE(all.equal(tbl_readme_roundtrip_dist, tbl_readme_expected, tolerance=1e-10)))

    tbl_single_phylo2table = phylo2table(get_single_branch_tree("A", 0.1), name_col="label", dist_col="dist")
    tbl_single_phylo2table_expected = data.frame(
        branch_id=c(1L, 0L),
        parent=c(-999L, 1L),
        sister=c(-999L, -999L),
        label=c("n0", "A"),
        dist=c(0, 0.1),
        stringsAsFactors=FALSE
    )
    expect_true(isTRUE(all.equal(tbl_single_phylo2table, tbl_single_phylo2table_expected, check.attributes=FALSE, tolerance=1e-10)))
    tbl_single_phylo2table_phy = table2phylo(tbl_single_phylo2table, name_col="label", dist_col="dist")
    expect_true(identical(as.character(tbl_single_phylo2table_phy$tip.label), "A"))

    tbl_zero_sister_dist = tbl
    tbl_zero_sister_dist$dist[tbl_zero_sister_dist$branch_id == 1L] = 0
    tbl_zero_sister_dist_phy = table2phylo(tbl_zero_sister_dist, name_col="label", dist_col="dist")
    expect_true(inherits(tbl_zero_sister_dist_phy, "phylo"))
    expect_true(setequal(tbl_zero_sister_dist_phy$tip.label, c("A", "B")))
    tbl_zero_sister_dist_dist = cophenetic(tbl_zero_sister_dist_phy)[c("A", "B"), c("A", "B")]
    tbl_zero_sister_dist_expected = matrix(c(0, 0.2, 0.2, 0), nrow=2, byrow=TRUE, dimnames=list(c("A", "B"), c("A", "B")))
    expect_true(isTRUE(all.equal(tbl_zero_sister_dist_dist, tbl_zero_sister_dist_expected, tolerance=1e-10)))

    tbl_short_sister_dist = tbl
    tbl_short_sister_dist$dist[tbl_short_sister_dist$branch_id == 1L] = 5e-9
    tbl_short_sister_dist_phy = table2phylo(tbl_short_sister_dist, name_col="label", dist_col="dist")
    expect_true(inherits(tbl_short_sister_dist_phy, "phylo"))
    expect_true(setequal(tbl_short_sister_dist_phy$tip.label, c("A", "B")))
    tbl_short_sister_dist_dist = cophenetic(tbl_short_sister_dist_phy)[c("A", "B"), c("A", "B")]
    tbl_short_sister_dist_expected = matrix(c(0, 0.200000005, 0.200000005, 0), nrow=2, byrow=TRUE, dimnames=list(c("A", "B"), c("A", "B")))
    expect_true(isTRUE(all.equal(tbl_short_sister_dist_dist, tbl_short_sister_dist_expected, tolerance=1e-10)))

    tbl_root_dist_na = tbl
    tbl_root_dist_na$dist[tbl_root_dist_na$branch_id == 3L] = NA_real_
    tbl_root_dist_na_phy = table2phylo(tbl_root_dist_na, name_col="label", dist_col="dist")
    expect_true(inherits(tbl_root_dist_na_phy, "phylo"))
    expect_true(setequal(tbl_root_dist_na_phy$tip.label, c("A", "B")))

    tbl_parent_na = data.frame(
        branch_id=c(10L, 20L, 5L),
        parent=c(5L, 5L, NA),
        sister=c(20L, 10L, -999L),
        label=c("A", "B", "Root"),
        dist=c(0.1, 0.2, NA_real_),
        stringsAsFactors=FALSE
    )
    tbl_parent_na_phy = table2phylo(tbl_parent_na, name_col="label", dist_col="dist")
    expect_true(inherits(tbl_parent_na_phy, "phylo"))
    expect_true(setequal(tbl_parent_na_phy$tip.label, c("A", "B")))
    tbl_factor_ids = data.frame(
        branch_id=factor(c(1L, 2L, 3L)),
        parent=factor(c(3L, 3L, -999L)),
        sister=factor(c(2L, 1L, -999L)),
        label=c("A", "B", "Root"),
        dist=c(0.1, 0.2, 0.0),
        stringsAsFactors=TRUE
    )
    tbl_factor_ids_phy = table2phylo(tbl_factor_ids, name_col="label", dist_col="dist")
    expect_true(inherits(tbl_factor_ids_phy, "phylo"))
    expect_true(setequal(tbl_factor_ids_phy$tip.label, c("A", "B")))
    tbl_missing_branch_id_err = tryCatch(
        {
            table2phylo(
                data.frame(
                    branch_id=c(1L, NA_integer_, 3L),
                    parent=c(3L, 3L, -999L),
                    sister=c(2L, 1L, -999L),
                    label=c("A", "B", "Root"),
                    dist=c(0.1, 0.2, 0.0),
                    stringsAsFactors=FALSE
                ),
                name_col="label",
                dist_col="dist"
            )
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(tbl_missing_branch_id_err))
    expect_true(grepl("branch_id contains missing/blank", conditionMessage(tbl_missing_branch_id_err), fixed=TRUE))
    tbl_single_tip = data.frame(
        branch_id=1L,
        parent=-999L,
        sister=-999L,
        label="A",
        dist=0.1,
        stringsAsFactors=FALSE
    )
    tbl_single_tip_phy = table2phylo(tbl_single_tip, name_col="label", dist_col="dist")
    expect_true(inherits(tbl_single_tip_phy, "phylo"))
    expect_true(identical(as.character(tbl_single_tip_phy$tip.label), "A"))
    tbl_single_tip_root = data.frame(
        branch_id=c(1L, 2L),
        parent=c(2L, -999L),
        sister=c(-999L, -999L),
        label=c("A", "Root"),
        dist=c(0.1, 0.0),
        stringsAsFactors=FALSE
    )
    tbl_single_tip_root_phy = table2phylo(tbl_single_tip_root, name_col="label", dist_col="dist")
    expect_true(inherits(tbl_single_tip_root_phy, "phylo"))
    expect_true(identical(as.character(tbl_single_tip_root_phy$tip.label), "A"))
    tbl_partial_parent = data.frame(
        branch_id=c("11", "d", "1", "3"),
        parent=c("12", "12", "", "a"),
        sister=c("3", "1", "2", "4"),
        label=c("J", "F", "O", "O"),
        dist=c("NaN", "Inf", "x", "0"),
        stringsAsFactors=FALSE
    )
    tbl_partial_parent_err = tryCatch(
        {
            table2phylo(tbl_partial_parent, name_col="label", dist_col="dist")
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(tbl_partial_parent_err))
    expect_true(grepl("contains non-numeric value", conditionMessage(tbl_partial_parent_err), fixed=TRUE))

    tbl_multi_root = tbl_parent_na
    tbl_multi_root$parent[1] = -999L
    multi_root_error = tryCatch(
        {
            table2phylo(tbl_multi_root, name_col="label", dist_col="dist")
            NULL
        },
        error=function(e) e
    )
    expect_true(!is.null(multi_root_error))
    expect_true(grepl("Ambiguous root candidate", conditionMessage(multi_root_error), fixed=TRUE))


    # Graph validation rejects orphaned/cyclic tables and asymmetric sisters.
    tbl_cycle = data.frame(
        branch_id=c(1L, 2L, 3L),
        parent=c(2L, 1L, -999L),
        sister=c(-999L, -999L, -999L),
        label=c("A", "B", "Root"),
        dist=c(0.1, 0.2, 0),
        stringsAsFactors=FALSE
    )
    expect_error_contains(table2phylo(tbl_cycle, "label", "dist"), "Disconnected or cyclic")
    tbl_bad_sister = tbl
    tbl_bad_sister$sister[tbl_bad_sister$branch_id == 2L] = 2L
    expect_error_contains(table2phylo(tbl_bad_sister, "label", "dist"), "Non-reciprocal sister")
    tbl_bad_root_sister = tbl
    tbl_bad_root_sister$sister[tbl_bad_root_sister$branch_id == 3L] = 1L
    expect_error_contains(table2phylo(tbl_bad_root_sister, "label", "dist"), "root row must have a sentinel sister")
    tbl_unknown_parent = tbl
    tbl_unknown_parent$parent[tbl_unknown_parent$branch_id == 1L] = 99L
    expect_error_contains(table2phylo(tbl_unknown_parent, "label", "dist"), "Unknown parent")

    tr_tip_name_collision = ape::read.tree(text="(n0:0,B:0.2);")
    tr_tip_name_collision$node.label = NULL
    collision_table = phylo2table(tr_tip_name_collision)
    expect_true(!anyDuplicated(collision_table$label))
    expect_true("n0" %in% collision_table$label)
    expect_true(isTRUE(all.equal(
        cophenetic(table2phylo(collision_table, "label", "dist")),
        cophenetic(tr_tip_name_collision),
        tolerance=1e-12
    )))
    tr_repeated_internal = ape::read.tree(text="((A:1,B:1)100:1,(C:1,D:1)100:1)100;")
    repeated_internal_table = phylo2table(tr_repeated_internal)
    repeated_internal_restored = table2phylo(repeated_internal_table, "label", "dist")
    expect_true(identical(
        as.character(repeated_internal_restored$node.label),
        as.character(tr_repeated_internal$node.label)
    ))
})
