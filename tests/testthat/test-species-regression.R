test_that("species-label parsing and validation", {
    withr::local_seed(20260831)
    expect_true(identical(as.character(get_species_name("A_B_gene1", species_parser="legacy")), "A B"))
    expect_true(identical(as.character(get_species_name("A_cf_B_gene1", species_parser="taxonomic")), "A cf B"))
    expect_true(identical(as.character(get_species_name("Amoeba_sp_JDSRuffled_gene1", species_parser="taxonomic")), "Amoeba sp JDSRuffled"))
    expect_true(identical(as.character(get_species_name("Bacillus_subtilis_subsp_168_gene1", species_parser="taxonomic")), "Bacillus subtilis subsp 168"))
    expect_true(identical(as.character(get_species_name("Solanum_lycopersicum_cultivar_Heinz1706_gene1", species_parser="taxonomic")), "Solanum lycopersicum cultivar Heinz1706"))
    expect_true(identical(as.character(get_species_name("Escherichia_coli_serovar_O157_gene1", species_parser="taxonomic")), "Escherichia coli serovar O157"))
    tr_dot = ape::read.tree(text="(A.B.g1:1,C.D.g2:1);")
    dot_species = get_species_names(tr_dot, sep=".")
    expect_true(identical(as.character(dot_species), c("A.B", "C.D")))
    tr_taxonomic = ape::read.tree(text="(A_B_gene1:1,A_cf_B_gene2:1,Amoeba_sp_JDSRuffled_gene3:1);")
    taxonomic_species = get_species_names(tr_taxonomic, species_parser="taxonomic")
    expect_true(identical(as.character(taxonomic_species), c("A_B", "A_cf_B", "Amoeba_sp_JDSRuffled")))
    tr_bad_species = ape::read.tree(text="(A:1,B_c_d:1);")
    bad_species = suppressWarnings(get_species_names(tr_bad_species, sep="_"))
    expect_true(is.na(bad_species[1]))
    expect_true(identical(as.character(bad_species[2]), "B_c"))
    expect_error(
        get_species_names(tr_dot, sep=NA_character_),
        "sep must be a single non-missing string", fixed=TRUE
    )

    leaf2species_bad = suppressWarnings(leaf2species(c("A_B_g1", "bad")))
    expect_true(length(leaf2species_bad) == 2)
    expect_true(identical(as.character(leaf2species_bad[1]), "A B"))
    expect_true(is.na(leaf2species_bad[2]))
    leaf2species_taxonomic = suppressWarnings(leaf2species(
        c("A_B_gene1", "A_cf_B_gene2", "Amoeba_sp_JDSRuffled_gene3", "Bacillus_subtilis_subsp_168_gene4"),
        species_parser="taxonomic"
    ))
    expect_true(identical(as.character(leaf2species_taxonomic), c("A B", "A cf B", "Amoeba sp JDSRuffled", "Bacillus subtilis subsp 168")))
    leaf2species_taxonomic_underbar = suppressWarnings(leaf2species(
        c("A_B_gene1", "A_cf_B_gene2", "Amoeba_sp_JDSRuffled_gene3"),
        use_underbar=TRUE,
        species_parser="taxonomic"
    ))
    expect_true(identical(as.character(leaf2species_taxonomic_underbar), c("A_B", "A_cf_B", "Amoeba_sp_JDSRuffled")))
})

test_that("restored model species counts match species-only and gene labels", {
    for (parser in c("legacy", "taxonomic")) {
        labels <- if (parser == "legacy") {
            c("Homo_sapiens", "Homo_sapiens_g1", "Mus_musculus_g1", "unparsed")
        } else c("A_cf_B", "A_cf_B_g1", "C_D_g1", "unparsed")
        for (sep in c("_", ".")) {
            tree <- ape::read.tree(text="((A:1,B:1):1,(C:1,D:1):1);")
            tree$tip.label <- gsub("_", sep, labels, fixed=TRUE)
            fit <- list(tree=tree, shift.configuration=integer(), nShifts=0, score=0)
            summary <- get_tree_table(fit, "l1ou", parser, sep)
            restored <- tree_table_collapse2original(summary, tree, parser, sep)
            expect_equal(summary$num_species, 3L)
            expect_identical(restored, summary)
        }
    }
})
