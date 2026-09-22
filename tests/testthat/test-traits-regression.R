test_that("trait filtering, replicates, and alignment", {
    withr::local_seed(20260831)
    comp = calc_complementarity(c(1, 2, 3), c(1, 2, 0), method="weighted")
    expect_true(isTRUE(all.equal(comp, 1 / 3, tolerance=1e-10)))

    ri = remove_invariant_traits(data.frame(a=c(1, 1), b=c(1, 2)))
    expect_true(identical(as.character(colnames(ri$trait_table)), "b"))

    rep_tbl = data.frame("a.1"=c(1, 2), "a.2"=c(3, 4), check.names=FALSE)
    rownames(rep_tbl) = c("x", "y")
    rep_merged = merge_replicates(rep_tbl, replicate_sep=".")
    expect_true(identical(colnames(rep_merged), "a"))
    expect_true(isTRUE(all.equal(as.numeric(rep_merged[,"a"]), c(2, 3), tolerance=1e-10)))
    rep_bases = get_expression_bases(rep_tbl, replicate_sep=".")
    expect_true(identical(as.character(rep_bases), "a"))

    rep_tbl_mixed = data.frame("gene_1"=c(1, 2), "gene_2"=c(3, 4), "control"=c(5, 6), check.names=FALSE)
    rownames(rep_tbl_mixed) = c("x", "y")
    rep_merged_mixed = merge_replicates(rep_tbl_mixed, replicate_sep="_")
    expect_true(identical(colnames(rep_merged_mixed), c("gene", "control")))
    expect_true(isTRUE(all.equal(as.numeric(rep_merged_mixed[,"gene"]), c(2, 3), tolerance=1e-10)))
    expect_true(identical(as.numeric(rep_merged_mixed[,"control"]), c(5, 6)))

    rep_tbl_qualified = data.frame(
        "Genus_species.1"=c(1, 2),
        "Genus_species.2"=c(3, 4),
        "Genus_cf_species.1"=c(5, 6),
        "Genus_cf_species.2"=c(7, 8),
        check.names=FALSE
    )
    rownames(rep_tbl_qualified) = c("x", "y")
    rep_merged_qualified = merge_replicates(rep_tbl_qualified, replicate_sep=".")
    expect_true(identical(colnames(rep_merged_qualified), c("Genus_species", "Genus_cf_species")))
    expect_true(isTRUE(all.equal(as.numeric(rep_merged_qualified[, "Genus_cf_species"]), c(6, 7), tolerance=1e-10)))
    expect_error(
        merge_replicates(rep_tbl, replicate_sep=NA_character_),
        "replicate_sep must be a single non-missing string", fixed=TRUE
    )

    sorted_single_col = sort_exp(
        exp=data.frame(gene_id=c("B", "A"), stringsAsFactors=FALSE),
        tree=ape::read.tree(text="(A:1,B:1);"),
        col="gene_id"
    )
    expect_true(identical(as.character(sorted_single_col$gene_id), c("A", "B")))
    expect_error(
        sort_exp(
            exp=data.frame(gene_id="A", value=1, stringsAsFactors=FALSE),
            tree=ape::read.tree(text="(A:1,B:1);"),
            col="gene_id"
        ),
        "missing rows for tree tip label", fixed=TRUE
    )
    expect_error(
        sort_exp(
            exp=data.frame(gene_id=c("A", "A"), value=c(1, 2), stringsAsFactors=FALSE),
            tree=ape::read.tree(text="(A:1,B:1);"),
            col="gene_id"
        ),
        "duplicated values in key column", fixed=TRUE
    )
    # Replicate suffix removal retains every part of a multi-separator base name.
    rep_tbl_multi_sep = data.frame(
        "alpha_beta_1"=c(1, 3),
        "alpha_beta_2"=c(3, 5),
        check.names=FALSE
    )
    rep_multi_sep = merge_replicates(rep_tbl_multi_sep, "_")
    expect_true(identical(colnames(rep_multi_sep), "alpha_beta"))
    expect_true(identical(get_expression_bases(rep_tbl_multi_sep, "_"), "alpha_beta"))
})

test_that("complementarity weights remain finite at extreme magnitudes", {
    x <- c(1, 0.5, 0)
    y <- c(0, 1, 0.25)
    for (scale in c(1e-300, 1, 1e308)) {
        expect_equal(calc_complementarity(x*scale, y*scale), calc_complementarity(x,y))
    }
    expect_equal(calc_complementarity(c(1e308,1e308), c(0,0)), 1)
    expect_equal(calc_complementarity(c(1e308,1e308), c(1e308,0)), 1/3)
    expect_equal(calc_complementarity(c(1e308,1e308), c(1e308,1e308)), 0)
    expect_equal(calc_complementarity(x*1e308, y*1e308, "independent"),
        calc_complementarity(x,y,"independent"))
})
