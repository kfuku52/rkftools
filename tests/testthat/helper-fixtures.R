# Small immutable starting points shared by the former smoke regressions.
fixture_gene_tree <- function() {
    ape::read.tree(text="((A_a_x:1,B_b_y:1):1,(A_a_z:1,C_c_w:1):1);")
}

fixture_trait_tree <- function() {
    ape::read.tree(text="((A:1,B:1):1,C:1);")
}

expect_error_contains <- function(expr, text) {
    expect_error(force(expr), regexp=text, fixed=TRUE)
}
