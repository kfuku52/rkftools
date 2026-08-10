test_that("node numbers must be integer-valued", {
    phy <- ape::read.tree(text="((A:1,B:1):1,C:1);")

    expect_error(get_tip_labels(phy, 1.9), "integer")
    expect_error(get_children_num(phy, 4.2), "integer")
    expect_error(
        collapse_clades(
            phy,
            data.frame(value=1:3, row.names=phy$tip.label),
            collapse_node_nums=4.5
        ),
        "integer node numbers"
    )
})

test_that("species overlap validates topology and cutoff", {
    binary <- ape::read.tree(text="((A_a_g1:1,A_a_g2:1):1,B_b_g3:1);")
    polytomy <- ape::read.tree(text="((A_a_g1:1,B_b_g2:1,C_c_g3:1):1,D_d_g4:1);")

    expect_error(get_species_overlap_score(binary, dc_cutoff=NA_real_), "dc_cutoff")
    expect_error(get_species_overlap_score(binary, dc_cutoff=1.1), "dc_cutoff")
    expect_error(get_species_overlap_score(polytomy), "binary")
})

test_that("logical controls are strict", {
    expect_error(get_parsed_args(character(), print=1), "logical")
    expect_error(is.blank(NULL, false.triggers="yes"), "logical")
})
