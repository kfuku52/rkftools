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
    tbl_readme = data.frame(
        branch_id=c(6L, 2L, 0L, 1L, 5L, 3L, 4L),
        parent=c(-999L, 6L, 2L, 2L, 6L, 5L, 5L),
        sister=c(-999L, 5L, 1L, 0L, 2L, 4L, 3L),
        node_name=c("Root", "Clade_AB", "A", "B", "Clade_CD", "C", "D"),
        dist=c(0, 0.42, 0.16, 0.18, 0.50, 0.22, 0.25),
        stringsAsFactors=FALSE
    )
    tbl_readme_phy = table2phylo(tbl_readme, name_col="node_name", dist_col="dist")
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
    expect_true(isTRUE(all.equal(tbl_readme_dist, tbl_readme_expected, tolerance=1e-10)))
    tbl_readme_roundtrip = phylo2table(tbl_readme_phy, name_col="node_name", dist_col="dist")
    expect_true(isTRUE(all.equal(tbl_readme_roundtrip, tbl_readme, check.attributes=FALSE, tolerance=1e-10)))
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

    tbl_short_sister_dist = tbl
    tbl_short_sister_dist$dist[tbl_short_sister_dist$branch_id == 1L] = 5e-9
    tbl_short_sister_dist_phy = table2phylo(tbl_short_sister_dist, name_col="label", dist_col="dist")
    expect_true(setequal(tbl_short_sister_dist_phy$tip.label, c("A", "B")))
    tbl_short_sister_dist_dist = cophenetic(tbl_short_sister_dist_phy)[c("A", "B"), c("A", "B")]
    tbl_short_sister_dist_expected = matrix(c(0, 0.200000005, 0.200000005, 0), nrow=2, byrow=TRUE, dimnames=list(c("A", "B"), c("A", "B")))
    expect_true(isTRUE(all.equal(tbl_short_sister_dist_dist, tbl_short_sister_dist_expected, tolerance=1e-10)))

    tbl_parent_na = data.frame(
        branch_id=c(10L, 20L, 5L),
        parent=c(5L, 5L, NA),
        sister=c(20L, 10L, -999L),
        label=c("A", "B", "Root"),
        dist=c(0.1, 0.2, NA_real_),
        stringsAsFactors=FALSE
    )
    tbl_parent_na_phy = table2phylo(tbl_parent_na, name_col="label", dist_col="dist")
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
    expect_true(setequal(tbl_factor_ids_phy$tip.label, c("A", "B")))
    expect_error(
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
        ),
        "branch_id contains missing/blank", fixed=TRUE
    )
    tbl_single_tip = data.frame(
        branch_id=1L,
        parent=-999L,
        sister=-999L,
        label="A",
        dist=0.1,
        stringsAsFactors=FALSE
    )
    tbl_single_tip_phy = table2phylo(tbl_single_tip, name_col="label", dist_col="dist")
    expect_true(identical(as.character(tbl_single_tip_phy$tip.label), "A"))
    tbl_partial_parent = data.frame(
        branch_id=c("11", "d", "1", "3"),
        parent=c("12", "12", "", "a"),
        sister=c("3", "1", "2", "4"),
        label=c("J", "F", "O", "O"),
        dist=c("NaN", "Inf", "x", "0"),
        stringsAsFactors=FALSE
    )
    expect_error(
        table2phylo(tbl_partial_parent, name_col="label", dist_col="dist"),
        "contains non-numeric value", fixed=TRUE
    )

    tbl_multi_root = tbl_parent_na
    tbl_multi_root$parent[1] = -999L
    expect_error(
        table2phylo(tbl_multi_root, name_col="label", dist_col="dist"),
        "Ambiguous root candidate", fixed=TRUE
    )

    # Graph validation rejects orphaned/cyclic tables and asymmetric sisters.
    tbl_cycle = data.frame(
        branch_id=c(1L, 2L, 3L),
        parent=c(2L, 1L, -999L),
        sister=c(-999L, -999L, -999L),
        label=c("A", "B", "Root"),
        dist=c(0.1, 0.2, 0),
        stringsAsFactors=FALSE
    )
    expect_error(table2phylo(tbl_cycle, "label", "dist"), "Disconnected or cyclic", fixed=TRUE)
    tbl_bad_sister = tbl
    tbl_bad_sister$sister[tbl_bad_sister$branch_id == 2L] = 2L
    expect_error(table2phylo(tbl_bad_sister, "label", "dist"), "Non-reciprocal sister", fixed=TRUE)
    tbl_bad_root_sister = tbl
    tbl_bad_root_sister$sister[tbl_bad_root_sister$branch_id == 3L] = 1L
    expect_error(table2phylo(tbl_bad_root_sister, "label", "dist"), "root row must have a sentinel sister", fixed=TRUE)
    tbl_unknown_parent = tbl
    tbl_unknown_parent$parent[tbl_unknown_parent$branch_id == 1L] = 99L
    expect_error(table2phylo(tbl_unknown_parent, "label", "dist"), "Unknown parent", fixed=TRUE)

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
