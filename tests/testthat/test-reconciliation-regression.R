test_that("species-overlap reference scores and parser modes", {
    withr::local_seed(20260831)
    tr_unlabeled = fixture_trait_tree()
    tr_fg = fixture_trait_tree()
    tr = ape::read.tree(text="((A_a_x:1,B_b_y:1):1,(A_a_z:1,C_c_w:1):1);")
    so = get_species_overlap_score(tr, dc_cutoff=0)
    expect_true(identical(as.numeric(so), 1))

    so_by_root = get_root_position_dependent_species_overlap_scores(tr, nslots=1)
    expect_true(length(so_by_root) == nrow(tr$edge))
    so_by_root_vec_slots = get_root_position_dependent_species_overlap_scores(tr, nslots=c(1L, 2L))
    expect_true(length(so_by_root_vec_slots) == nrow(tr$edge))
    dup_score_root = get_duplication_confidence_score(tr, get_root_num(tr))
    expect_true(isTRUE(all.equal(as.numeric(dup_score_root), 1 / 3, tolerance=1e-10)))

    tr_vertebrate_duplication = ape::read.tree(text=paste0(
        "(((Homo_sapiens_geneA:0.12,Mus_musculus_geneA:0.12):0.18,",
        "(Danio_rerio_geneA:0.16,Oryzias_latipes_geneA:0.16):0.14):0.25,",
        "((Homo_sapiens_geneB:0.11,Mus_musculus_geneB:0.11):0.19,",
        "(Danio_rerio_geneB:0.15,Oryzias_latipes_geneB:0.15):0.15):0.25);"
    ))
    vertebrate_species = get_species_names(tr_vertebrate_duplication)
    expect_true(setequal(
        unique(as.character(vertebrate_species)),
        c("Homo_sapiens", "Mus_musculus", "Danio_rerio", "Oryzias_latipes")
    ))
    vertebrate_root_scores = get_root_position_dependent_species_overlap_scores(
        tr_vertebrate_duplication,
        nslots=1
    )
    vertebrate_root_child_edges = which(
        tr_vertebrate_duplication$edge[,1] == get_root_num(tr_vertebrate_duplication)
    )
    expect_true(identical(as.integer(vertebrate_root_child_edges), c(1L, 8L)))
    vertebrate_unrooted_scores = c(
        vertebrate_root_scores[2:7],
        vertebrate_root_scores[vertebrate_root_child_edges[1]],
        vertebrate_root_scores[9:14]
    )
    expect_true(identical(
        as.numeric(vertebrate_root_scores),
        c(1, 2, 3, 3, 2, 3, 3, 1, 2, 3, 3, 2, 3, 3)
    ))
    expect_true(identical(
        as.numeric(vertebrate_unrooted_scores),
        c(2, 3, 3, 2, 3, 3, 1, 2, 3, 3, 2, 3, 3)
    ))
    expect_true(identical(which(vertebrate_unrooted_scores == min(vertebrate_unrooted_scores)), 7L))

    tr_taxonomic_overlap = ape::read.tree(text="((Genus_species_gene1:1,Genus_cf_species_gene2:1):1,(Genus_species_gene3:1,Other_species_gene4:1):1);")
    so_taxonomic_legacy = get_species_overlap_score(tr_taxonomic_overlap, dc_cutoff=0, species_parser="legacy")
    so_taxonomic = get_species_overlap_score(tr_taxonomic_overlap, dc_cutoff=0, species_parser="taxonomic")
    expect_true(identical(as.numeric(so_taxonomic_legacy), 1))
    expect_true(identical(as.numeric(so_taxonomic), 1))
    taxonomic_dup_node = get_parent_num(tr_taxonomic_overlap, get_node_num_by_name(tr_taxonomic_overlap, "Genus_species_gene1"))
    dup_score_taxonomic_legacy = get_duplication_confidence_score(
        tr_taxonomic_overlap,
        taxonomic_dup_node,
        species_parser="legacy"
    )
    dup_score_taxonomic = get_duplication_confidence_score(
        tr_taxonomic_overlap,
        taxonomic_dup_node,
        species_parser="taxonomic"
    )
    expect_true(identical(as.numeric(dup_score_taxonomic_legacy), 0))
    expect_true(identical(as.numeric(dup_score_taxonomic), 0))
    so_by_root_taxonomic = get_root_position_dependent_species_overlap_scores(
        tr_taxonomic_overlap,
        nslots=1,
        species_parser="taxonomic"
    )
    expect_true(length(so_by_root_taxonomic) == nrow(tr_taxonomic_overlap$edge))
    expect_true(isTRUE(contains_polytomy(ape::read.tree(text="((A:1,B:1,C:1):1,D:1);"))))
    expect_true(isFALSE(contains_polytomy(tr_unlabeled)))

    # Reconciliation requires parseable species labels rather than silently using
    # malformed gene labels as distinct species.
    tr_malformed_species = ape::read.tree(text="(A:1,B_c_gene:1);")
    expect_error_contains(get_species_overlap_score(tr_malformed_species), "Unable to parse species")
    expect_error_contains(
        get_duplication_confidence_score(tr_malformed_species, get_root_num(tr_malformed_species)),
        "Unable to parse species"
    )
    tr_blank_species = ape::read.tree(text="(A__gene:1,B_c_gene:1);")
    expect_true(is.na(suppressWarnings(get_species_names(tr_blank_species))[1]))

    expect_error_contains(
        count_foreground_lineage(
            tr_fg,
            data.frame(species=c("A", "B"), fg=c(1, 0), stringsAsFactors=FALSE)
        ),
        "missing rows for tree tip label"
    )

    # The dynamic root-score implementation agrees with the previous rerooting
    # definition over randomized trees.
    set.seed(20260809)
    for (iteration in seq_len(5L)) {
        random_tree = ape::rtree(10)
        random_tree$tip.label = paste0(
            "Genus", rep(seq_len(4), length.out=10),
            "_species", rep(seq_len(4), length.out=10),
            "_gene", seq_len(10)
        )
        fast_scores = get_root_position_dependent_species_overlap_scores(random_tree)
        slow_scores = vapply(seq_len(nrow(random_tree$edge)), function(edge_index) {
            rerooted_tree = suppressWarnings(phytools::reroot(
                tree=random_tree,
                node.number=random_tree$edge[edge_index, 2],
                position=random_tree$edge.length[edge_index] / 2
            ))
            get_species_overlap_score(rerooted_tree)
        }, numeric(1))
        expect_true(identical(as.numeric(fast_scores), as.numeric(slow_scores)))
    }
})
