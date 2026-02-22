library(rkftools)

comp = calc_complementarity(c(1, 2, 3), c(1, 2, 0), method="weighted")
stopifnot(isTRUE(all.equal(comp, 1 / 3, tolerance=1e-10)))

tr = ape::read.tree(text="((A_a_x:1,B_b_y:1):1,(A_a_z:1,C_c_w:1):1);")
so = get_species_overlap_score(tr, dc_cutoff=0)
stopifnot(identical(as.numeric(so), 1))

so_by_root = get_root_position_dependent_species_overlap_scores(tr, nslots=1)
stopifnot(length(so_by_root) == nrow(tr$edge))
so_by_root_vec_slots = get_root_position_dependent_species_overlap_scores(tr, nslots=c(1L, 2L))
stopifnot(length(so_by_root_vec_slots) == nrow(tr$edge))
dup_score_root = get_duplication_confidence_score(tr, get_root_num(tr))
stopifnot(isTRUE(all.equal(as.numeric(dup_score_root), 1 / 3, tolerance=1e-10)))
tr_dot = ape::read.tree(text="(A.B.g1:1,C.D.g2:1);")
dot_species = get_species_names(tr_dot, sep=".")
stopifnot(identical(as.character(dot_species), c("A.B", "C.D")))
tr_bad_species = ape::read.tree(text="(A:1,B_c_d:1);")
bad_species = suppressWarnings(get_species_names(tr_bad_species, sep="_"))
stopifnot(is.na(bad_species[1]))
stopifnot(identical(as.character(bad_species[2]), "B_c"))
species_sep_na_err = tryCatch(
    {
        get_species_names(tr_dot, sep=NA_character_)
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(species_sep_na_err))
stopifnot(grepl("sep must be a single non-missing string", conditionMessage(species_sep_na_err), fixed=TRUE))
species_sep_vec_err = tryCatch(
    {
        get_species_names(tr_dot, sep=c("_", "."))
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(species_sep_vec_err))
stopifnot(grepl("sep must be a single non-missing string", conditionMessage(species_sep_vec_err), fixed=TRUE))

notung_tmp = tempfile()
writeLines(c("#D g1 s1 s2"), notung_tmp)
notung_df1 = read_notung_parsable(notung_tmp, mode="D")
stopifnot(nrow(notung_df1) == 1)
stopifnot(identical(as.character(notung_df1$gn_node), "g1"))
writeLines(c("#D g1 s1 s2", "#D g2 s3 s4"), notung_tmp)
notung_df2 = read_notung_parsable(notung_tmp, mode="D")
stopifnot(nrow(notung_df2) == 2)
stopifnot(identical(as.character(notung_df2$gn_node), c("g1", "g2")))
writeLines(c("  #D g3 s5 s6"), notung_tmp)
notung_df3 = read_notung_parsable(notung_tmp, mode="D")
stopifnot(nrow(notung_df3) == 1)
stopifnot(identical(as.character(notung_df3$gn_node), "g3"))
notung_mode_na_err = tryCatch(
    {
        read_notung_parsable(notung_tmp, mode=NA_character_)
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(notung_mode_na_err))
stopifnot(grepl("mode must be a single non-missing string", conditionMessage(notung_mode_na_err), fixed=TRUE))
notung_mode_vec_err = tryCatch(
    {
        read_notung_parsable(notung_tmp, mode=c("D", "X"))
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(notung_mode_vec_err))
stopifnot(grepl("mode must be a single non-missing string", conditionMessage(notung_mode_vec_err), fixed=TRUE))
unlink(notung_tmp)

rerooted = phytools::reroot(tree=tr, node.number=tr$edge[2,2])
root_idx = get_phy2_root_in_phy1(tr, rerooted, nslots=1, mode="index")
stopifnot(!is.na(root_idx))
root_idx_vec_slots = get_phy2_root_in_phy1(tr, rerooted, nslots=c(1L, 2L), mode="index")
stopifnot(!is.na(root_idx_vec_slots))
root_idx_na_slots = get_phy2_root_in_phy1(tr, rerooted, nslots=NA_integer_, mode="index")
stopifnot(!is.na(root_idx_na_slots))

mad_n1 = MAD_parallel(tr, output_mode="newick", ncpu=1)
mad_auto = MAD_parallel(tr, output_mode="newick")
stopifnot(length(mad_n1) == length(mad_auto))
mad_vector_ncpu = MAD_parallel(tr, output_mode="newick", ncpu=c(1L, 2L))
stopifnot(length(mad_vector_ncpu) == length(mad_n1))
mad_serial = MAD(tr, output_mode="newick")
stopifnot(length(mad_serial) == length(mad_n1))
mad_mode_na_err = tryCatch(
    {
        MAD(tr, output_mode=NA_character_)
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(mad_mode_na_err))
stopifnot(grepl("output_mode must be a single non-missing string", conditionMessage(mad_mode_na_err), fixed=TRUE))
mad_mode_vec_err = tryCatch(
    {
        MAD_parallel(tr, output_mode=c("newick", "stats"))
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(mad_mode_vec_err))
stopifnot(grepl("output_mode must be a single non-missing string", conditionMessage(mad_mode_vec_err), fixed=TRUE))

old_core_limit_env = Sys.getenv("_R_CHECK_LIMIT_CORES_", unset=NA_character_)
Sys.setenv("_R_CHECK_LIMIT_CORES_"="TRUE")
mad_limited = MAD_parallel(tr, output_mode="newick")
stopifnot(length(mad_limited) == length(mad_n1))
if (is.na(old_core_limit_env)) {
    Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
} else {
    Sys.setenv("_R_CHECK_LIMIT_CORES_"=old_core_limit_env)
}
stopifnot(isTRUE(is.blank(c("", NA_character_))))
stopifnot(isFALSE(is.blank(c("", "x"))))
is_blank_false_trigger_err = tryCatch(
    {
        is.blank(c("", NA_character_), false.triggers=NA)
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(is_blank_false_trigger_err))
stopifnot(grepl("false.triggers must be a single non-missing logical value", conditionMessage(is_blank_false_trigger_err), fixed=TRUE))
parsed_args = get_parsed_args(c("--threads=2", "--dry-run"), print=FALSE)
stopifnot(identical(as.numeric(parsed_args[["threads"]]), 2))
stopifnot(isTRUE(parsed_args[["dry-run"]]))
parsed_arg_err = tryCatch(
    {
        get_parsed_args(c("--=1"), print=FALSE)
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(parsed_arg_err))
stopifnot(grepl("parameter name is empty", conditionMessage(parsed_arg_err), fixed=TRUE))
parsed_arg_na_err = tryCatch(
    {
        get_parsed_args(c(NA_character_), print=FALSE)
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(parsed_arg_na_err))
stopifnot(grepl("single non-missing string", conditionMessage(parsed_arg_na_err), fixed=TRUE))
parsed_print_na_err = tryCatch(
    {
        get_parsed_args(c("--threads=2"), print=NA)
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(parsed_print_na_err))
stopifnot(grepl("print must be a single non-missing logical value", conditionMessage(parsed_print_na_err), fixed=TRUE))
ri = remove_invariant_traits(data.frame(a=c(1, 1), b=c(1, 2)))
stopifnot(identical(as.character(colnames(ri$trait_table)), "b"))

pcm_mock = list(
    tree=tr,
    shift.configuration=structure(c(1L, 2L, 3L), names=c("1", "1", "2")),
    nShifts=3L,
    score=1.23
)
tree_table_mock = get_tree_table(pcm_mock, mode="l1ou")
stopifnot(identical(as.integer(tree_table_mock$num_regime), 3L))
stopifnot(identical(as.integer(tree_table_mock$num_conv_regime), 1L))
stopifnot(identical(as.integer(tree_table_mock$num_uniq_regime), 2L))
tree_table_mode_na_err = tryCatch(
    {
        get_tree_table(pcm_mock, mode=NA_character_)
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(tree_table_mode_na_err))
stopifnot(grepl("mode must be a single non-missing string", conditionMessage(tree_table_mode_na_err), fixed=TRUE))

tr_fg = ape::read.tree(text="((A:1,B:1):1,C:1);")
trait_numeric = data.frame(species=c("A", "B", "C"), fg=c(1, 1, 0), stringsAsFactors=FALSE)
fg_count = count_foreground_lineage(tr_fg, trait_numeric)
stopifnot(identical(as.integer(fg_count$fg), 1L))
fg_unknown_err = tryCatch(
    {
        count_foreground_lineage(
            tr_fg,
            data.frame(species=c("A", "X", "C"), fg=c(1, 1, 0), stringsAsFactors=FALSE)
        )
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(fg_unknown_err))
stopifnot(grepl("not present in tree tip labels", conditionMessage(fg_unknown_err), fixed=TRUE))

tr_unlabeled = ape::read.tree(text="((A:1,B:1):1,C:1);")
stopifnot(is.null(tr_unlabeled$node.label))
tr_filled = fill_node_labels(tr_unlabeled)
stopifnot(length(tr_filled$node.label) == tr_filled$Nnode)
stopifnot(!any(is.na(tr_filled$node.label) | tr_filled$node.label == ""))
root_num = get_root_num(tr_unlabeled)
stopifnot(identical(as.integer(root_num), 4L))
stopifnot(isTRUE(is_root(tr_unlabeled, root_num)))
stopifnot(isTRUE(is_leaf(tr_unlabeled, 1L)))
stopifnot(isFALSE(is_leaf(tr_unlabeled, root_num)))
all_node_names = c(tr_filled$tip.label, tr_filled$node.label)
mapped_nums = get_node_num_by_name(tr_filled, c("C", "A", "C"))
stopifnot(identical(as.integer(mapped_nums), c(3L, 1L, 3L)))
mapped_nums_na = get_node_num_by_name(tr_filled, c("C", NA_character_, "X", "A"))
stopifnot(identical(as.integer(mapped_nums_na), c(3L, 1L)))
mapped_names = get_node_name_by_num(tr_filled, c(5L, 4L, 5L))
stopifnot(identical(as.character(mapped_names), as.character(all_node_names[c(5L, 4L, 5L)])))
mapped_names_na = get_node_name_by_num(tr_filled, c(5L, NA_integer_, 4L, 999L))
stopifnot(identical(as.character(mapped_names_na), as.character(all_node_names[c(5L, 4L)])))
internal_nodes = sort(unique(tr_unlabeled$edge[,1]))
internal_nodes = internal_nodes[internal_nodes > length(tr_unlabeled$tip.label)]
size_two_node = internal_nodes[vapply(internal_nodes, function(nn) {
    length(get_tip_labels(tr_unlabeled, nn)) == 2
}, logical(1))][1]
tip_labels_vec = get_tip_labels(tr_unlabeled, c(3L, size_two_node))
stopifnot(identical(as.character(tip_labels_vec), c("C", "A", "B")))
tip_labels_invalid_err = tryCatch(
    {
        get_tip_labels(tr_unlabeled, 999L)
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(tip_labels_invalid_err))
stopifnot(grepl("outside valid node range", conditionMessage(tip_labels_invalid_err), fixed=TRUE))
children_vec = get_children_num(tr_unlabeled, c(4L, 5L))
stopifnot(identical(as.integer(children_vec), c(5L, 3L, 1L, 2L)))
children_na = get_children_num(tr_unlabeled, NA_integer_)
stopifnot(identical(as.integer(children_na), integer(0)))
parent_vec = get_parent_num(tr_unlabeled, c(1L, 2L, 3L))
stopifnot(identical(as.integer(parent_vec), c(5L, 5L, 4L)))
parent_na = get_parent_num(tr_unlabeled, NA_integer_)
stopifnot(identical(as.integer(parent_na), integer(0)))
sister_vec = get_sister_num(tr_unlabeled, c(1L, 2L))
stopifnot(identical(as.integer(sister_vec), c(2L, 1L)))
sister_na = get_sister_num(tr_unlabeled, NA_integer_)
stopifnot(identical(as.integer(sister_na), integer(0)))
ancestor_a = get_ancestor_num(tr_unlabeled, 1L)
stopifnot(identical(as.integer(ancestor_a), c(5L, 4L)))
ancestor_invalid = get_ancestor_num(tr_unlabeled, 99L)
stopifnot(length(ancestor_invalid) == 0)
desc_leaf = get_descendent_num(tr_unlabeled, c(4L, 5L), leaf_only=TRUE)
stopifnot(identical(as.integer(desc_leaf), c(1L, 2L, 3L)))
desc_leaf_only_na_err = tryCatch(
    {
        get_descendent_num(tr_unlabeled, c(4L, 5L), leaf_only=NA)
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(desc_leaf_only_na_err))
stopifnot(grepl("leaf_only must be a single non-missing logical value", conditionMessage(desc_leaf_only_na_err), fixed=TRUE))
mrca_tr_unlabeled = ape::mrca(tr_unlabeled)
nearest_out = get_nearest_tips(
    tr_unlabeled,
    query="A",
    subjects=c("B", "C"),
    mrca_matrix=mrca_tr_unlabeled
)
stopifnot(identical(as.character(nearest_out$nearests), "B"))
nearest_query_err = tryCatch(
    {
        get_nearest_tips(tr_unlabeled, query="X", subjects=c("B", "C"), mrca_matrix=mrca_tr_unlabeled)
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(nearest_query_err))
stopifnot(grepl("query must be present in mrca_matrix row names", conditionMessage(nearest_query_err), fixed=TRUE))
nearest_subject_err = tryCatch(
    {
        get_nearest_tips(tr_unlabeled, query="A", subjects=c("B", "X"), mrca_matrix=mrca_tr_unlabeled)
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(nearest_subject_err))
stopifnot(grepl("subjects are missing in mrca_matrix column names", conditionMessage(nearest_subject_err), fixed=TRUE))
nearest_query_vec_err = tryCatch(
    {
        get_nearest_tips(tr_unlabeled, query=c("A", "B"), subjects=c("B", "C"), mrca_matrix=mrca_tr_unlabeled)
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(nearest_query_vec_err))
stopifnot(grepl("query must be a single non-missing string", conditionMessage(nearest_query_vec_err), fixed=TRUE))
nearest_subject_empty_err = tryCatch(
    {
        get_nearest_tips(tr_unlabeled, query="A", subjects=character(0), mrca_matrix=mrca_tr_unlabeled)
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(nearest_subject_empty_err))
stopifnot(grepl("subjects must contain at least one non-missing tip label", conditionMessage(nearest_subject_empty_err), fixed=TRUE))
ultra_tree = ape::compute.brlen(ape::stree(4), 1)
ultra_tree = force_ultrametric(ultra_tree, stop_if_larger_change=1)
tip_age = get_node_age(ultra_tree, 1)
stopifnot(isTRUE(all.equal(as.numeric(tip_age), 0, tolerance=1e-10)))
node_age_err = tryCatch(
    {
        get_node_age(ultra_tree, 999)
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(node_age_err))
stopifnot(grepl("must be a single integer", conditionMessage(node_age_err), fixed=TRUE))
outgroup_labels = get_outgroup(tr_unlabeled)
stopifnot(identical(as.character(outgroup_labels), "C"))
outgroup_unrooted_err = tryCatch(
    {
        get_outgroup(ape::unroot(tr_unlabeled))
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(outgroup_unrooted_err))
stopifnot(grepl("requires a rooted tree", conditionMessage(outgroup_unrooted_err), fixed=TRUE))
single_branch = get_single_branch_tree("A", 0.5)
stopifnot(inherits(single_branch, "phylo"))
stopifnot(identical(as.character(single_branch$tip.label), "A"))
rooted_newick_ok = get_rooted_newick(tr_unlabeled, madr=1, rho=rep(0.5, nrow(tr_unlabeled$edge)))
stopifnot(length(rooted_newick_ok) == 3)
rooted_newick_err = tryCatch(
    {
        get_rooted_newick(tr_unlabeled, madr=999, rho=rep(0.5, nrow(tr_unlabeled$edge)))
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(rooted_newick_err))
stopifnot(grepl("must be a single edge index", conditionMessage(rooted_newick_err), fixed=TRUE))
single_branch_name_err = tryCatch(
    {
        get_single_branch_tree(c("A", "B"), 0.5)
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(single_branch_name_err))
stopifnot(grepl("single non-empty tip label", conditionMessage(single_branch_name_err), fixed=TRUE))
tr_na_edge = tr_unlabeled
tr_na_edge$edge.length[1] = NA_real_
collapsed_na_edge = collapse_short_branches(tr_na_edge, tol=1e-8)
stopifnot(inherits(collapsed_na_edge, "phylo"))
stopifnot(is.na(collapsed_na_edge$edge.length[1]))
padded_na_edge = pad_short_edges(tr_na_edge, threshold=1e-6)
stopifnot(inherits(padded_na_edge, "phylo"))
stopifnot(is.na(padded_na_edge$edge.length[1]))
pad_external_only_na_err = tryCatch(
    {
        pad_short_edges(tr_unlabeled, threshold=1e-6, external_only=NA)
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(pad_external_only_na_err))
stopifnot(grepl("external_only must be a single non-missing logical value", conditionMessage(pad_external_only_na_err), fixed=TRUE))
force_ultra_na_err = tryCatch(
    {
        force_ultrametric(tr_na_edge)
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(force_ultra_na_err))
stopifnot(grepl("contains NA edge lengths", conditionMessage(force_ultra_na_err), fixed=TRUE))
leaf2species_bad = suppressWarnings(leaf2species(c("A_B_g1", "bad")))
stopifnot(length(leaf2species_bad) == 2)
stopifnot(identical(as.character(leaf2species_bad[1]), "A B"))
stopifnot(is.na(leaf2species_bad[2]))
leaf2species_use_underbar_na_err = tryCatch(
    {
        leaf2species(c("A_B_g1"), use_underbar=NA)
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(leaf2species_use_underbar_na_err))
stopifnot(grepl("use_underbar must be a single non-missing logical value", conditionMessage(leaf2species_use_underbar_na_err), fixed=TRUE))
stopifnot(identical(as.character(get_species_name("A_B_gene1")), "A B"))
stopifnot(isTRUE(contains_polytomy(ape::read.tree(text="((A:1,B:1,C:1):1,D:1);"))))
stopifnot(isFALSE(contains_polytomy(tr_unlabeled)))
short_ext = collapse_short_external_edges(ape::read.tree(text="((A:1e-9,B:1e-9):1,C:1);"), threshold=1e-6)
stopifnot(inherits(short_ext, "phylo"))
deepest_node = get_deepest_node_num(tr_unlabeled, c(4L, 5L))
stopifnot(identical(as.integer(deepest_node), 4L))
tr_reorder = ape::read.tree(text="((B:1,A:1):1,C:1);")
stopifnot(isTRUE(has_same_leaves(tr_unlabeled, get_root_num(tr_unlabeled), tr_reorder, get_root_num(tr_reorder))))
stopifnot(isTRUE(is_same_root(tr_unlabeled, tr_reorder)))
stopifnot(isFALSE(is_same_root(tr_unlabeled, ape::read.tree(text="((A:1,C:1):1,B:1);"))))
tr_poly_root1 = ape::read.tree(text="(A:1,B:1,C:1,D:1);")
tr_poly_root2 = ape::read.tree(text="(C:1,D:1,A:1,B:1);")
tr_poly_root1$root.edge = 0
tr_poly_root2$root.edge = 0
stopifnot(isTRUE(is_same_root(tr_poly_root1, tr_poly_root2)))
multi2bi_empty = multi2bi_node_number_transfer(tr_unlabeled, tr_unlabeled)
stopifnot(nrow(multi2bi_empty) == 0)
set.seed(1)
tr_poly = ape::read.tree(text="((A:1,B:1,C:1):1,D:1);")
tr_poly_bi = ape::multi2di(tr_poly)
tr_poly_bi$tip.label = tr_poly$tip.label
multi2bi_map = multi2bi_node_number_transfer(tr_poly, tr_poly_bi)
stopifnot(nrow(multi2bi_map) >= 1)

rep_tbl = data.frame("a.1"=c(1, 2), "a.2"=c(3, 4), check.names=FALSE)
rownames(rep_tbl) = c("x", "y")
rep_merged = merge_replicates(rep_tbl, replicate_sep=".")
stopifnot(identical(colnames(rep_merged), "a"))
stopifnot(isTRUE(all.equal(as.numeric(rep_merged[,"a"]), c(2, 3), tolerance=1e-10)))
rep_bases = get_expression_bases(rep_tbl, replicate_sep=".")
stopifnot(identical(as.character(rep_bases), "a"))

rep_tbl_mixed = data.frame("gene_1"=c(1, 2), "gene_2"=c(3, 4), "control"=c(5, 6), check.names=FALSE)
rownames(rep_tbl_mixed) = c("x", "y")
rep_merged_mixed = merge_replicates(rep_tbl_mixed, replicate_sep="_")
stopifnot(identical(colnames(rep_merged_mixed), c("gene", "control")))
stopifnot(isTRUE(all.equal(as.numeric(rep_merged_mixed[,"gene"]), c(2, 3), tolerance=1e-10)))
stopifnot(identical(as.numeric(rep_merged_mixed[,"control"]), c(5, 6)))

rep_tbl_na = data.frame("g_1"=c(NA, 1), "g_2"=c(NA, 3), check.names=FALSE)
rownames(rep_tbl_na) = c("x", "y")
rep_merged_na = merge_replicates(rep_tbl_na, replicate_sep="_")
stopifnot(is.na(rep_merged_na["x", "g"]))
stopifnot(!is.nan(rep_merged_na["x", "g"]))
rep_sep_na_err = tryCatch(
    {
        merge_replicates(rep_tbl, replicate_sep=NA_character_)
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(rep_sep_na_err))
stopifnot(grepl("replicate_sep must be a single non-missing string", conditionMessage(rep_sep_na_err), fixed=TRUE))
rep_bases_sep_na_err = tryCatch(
    {
        get_expression_bases(rep_tbl, replicate_sep=NA_character_)
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(rep_bases_sep_na_err))
stopifnot(grepl("replicate_sep must be a single non-missing string", conditionMessage(rep_bases_sep_na_err), fixed=TRUE))

sorted_single_col = sort_exp(
    exp=data.frame(gene_id=c("B", "A"), stringsAsFactors=FALSE),
    tree=ape::read.tree(text="(A:1,B:1);"),
    col="gene_id"
)
stopifnot(identical(as.character(sorted_single_col$gene_id), c("A", "B")))
sort_missing_err = tryCatch(
    {
        sort_exp(
            exp=data.frame(gene_id="A", value=1, stringsAsFactors=FALSE),
            tree=ape::read.tree(text="(A:1,B:1);"),
            col="gene_id"
        )
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(sort_missing_err))
stopifnot(grepl("missing rows for tree tip label", conditionMessage(sort_missing_err), fixed=TRUE))
sort_dup_err = tryCatch(
    {
        sort_exp(
            exp=data.frame(gene_id=c("A", "A"), value=c(1, 2), stringsAsFactors=FALSE),
            tree=ape::read.tree(text="(A:1,B:1);"),
            col="gene_id"
        )
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(sort_dup_err))
stopifnot(grepl("duplicated values in key column", conditionMessage(sort_dup_err), fixed=TRUE))
imp_in = data.frame(t1=c(1, NA, 3), row.names=tr_unlabeled$tip.label)
imp_out = tryCatch(
    {
        phylogenetic_imputation(tr_unlabeled, imp_in)
    },
    error=function(e) e
)
if (inherits(imp_out, "error")) {
    stopifnot(grepl("Rphylopars", conditionMessage(imp_out), fixed=TRUE))
} else {
    stopifnot(is.data.frame(imp_out))
    stopifnot(nrow(imp_out) == length(tr_unlabeled$tip.label))
    stopifnot(ncol(imp_out) == 1)
    stopifnot(identical(rownames(imp_out), tr_unlabeled$tip.label))
}

pcm_reg = list(
    Y=data.frame(t1=c(1, 2, 3), row.names=tr_unlabeled$tip.label),
    tree=tr_unlabeled,
    shift.configuration=structure(c(1L, 2L), names=c("regA", "regB")),
    nShifts=2L,
    shift.values=matrix(c(0.1, 0.2), ncol=1),
    shift.means=matrix(c(1.1, 1.2), ncol=1),
    alpha=1,
    sigma2=1,
    intercept=0,
    logLik=-1
)
regime_tbl = get_regime_table(pcm_reg, mode="l1ou")
stopifnot(any(as.character(regime_tbl$regime) == "regA"))
stopifnot(any(as.character(regime_tbl$regime) == "regB"))
regime_tbl_mode_na_err = tryCatch(
    {
        get_regime_table(pcm_reg, mode=NA_character_)
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(regime_tbl_mode_na_err))
stopifnot(grepl("mode must be a single non-missing string", conditionMessage(regime_tbl_mode_na_err), fixed=TRUE))
pcm_reg_multi = pcm_reg
pcm_reg_multi$Y = data.frame(
    t1=c(1, 2, 3),
    t2=c(4, 5, 6),
    row.names=tr_unlabeled$tip.label
)
pcm_reg_multi$shift.values = matrix(c(0.1, 0.2, 0.3, 0.4), ncol=2, byrow=TRUE)
pcm_reg_multi$shift.means = matrix(c(1.1, 1.2, 1.3, 1.4), ncol=2, byrow=TRUE)
pcm_reg_multi$alpha = 2
pcm_reg_multi$sigma2 = 3
pcm_reg_multi$intercept = 4
pcm_reg_multi$logLik = -5
regime_tbl_multi = get_regime_table(pcm_reg_multi, mode="l1ou")
alpha_row_multi = regime_tbl_multi[as.character(regime_tbl_multi$param) == "alpha", c("t1", "t2"), drop=FALSE]
stopifnot(nrow(alpha_row_multi) == 1)
stopifnot(identical(as.numeric(alpha_row_multi[1,]), c(2, 2)))
pcm_reg_bad_alpha = pcm_reg_multi
pcm_reg_bad_alpha$alpha = c(1, 2, 3)
regime_tbl_bad_alpha_err = tryCatch(
    {
        get_regime_table(pcm_reg_bad_alpha, mode="l1ou")
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(regime_tbl_bad_alpha_err))
stopifnot(grepl("must have length 1 or 2", conditionMessage(regime_tbl_bad_alpha_err), fixed=TRUE))
pcm_reg_bad_shift_type = pcm_reg
pcm_reg_bad_shift_type$shift.configuration = structure(c("x", "2"), names=c("regA", "regB"))
regime_tbl_shift_type_err = tryCatch(
    {
        get_regime_table(pcm_reg_bad_shift_type, mode="l1ou")
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(regime_tbl_shift_type_err))
stopifnot(grepl("must contain integer edge indices", conditionMessage(regime_tbl_shift_type_err), fixed=TRUE))
pcm_reg_na_shift = pcm_reg
pcm_reg_na_shift$shift.configuration = structure(c(NA_integer_, 2L), names=c("regA", "regB"))
regime_tbl_na_shift_err = tryCatch(
    {
        get_regime_table(pcm_reg_na_shift, mode="l1ou")
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(regime_tbl_na_shift_err))
stopifnot(grepl("must contain integer edge indices", conditionMessage(regime_tbl_na_shift_err), fixed=TRUE))
leaf_regimes_ok = get_leaf_regimes(pcm_reg, mode="l1ou")
stopifnot(identical(as.character(leaf_regimes_ok$label), as.character(tr_unlabeled$tip.label)))
leaf_regimes_mode_na_err = tryCatch(
    {
        get_leaf_regimes(pcm_reg, mode=NA_character_)
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(leaf_regimes_mode_na_err))
stopifnot(grepl("mode must be a single non-missing string", conditionMessage(leaf_regimes_mode_na_err), fixed=TRUE))
pcm_reg_bad_labels = pcm_reg
rownames(pcm_reg_bad_labels$Y) = c("x", "y", "z")
leaf_regimes_err = tryCatch(
    {
        get_leaf_regimes(pcm_reg_bad_labels, mode="l1ou")
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(leaf_regimes_err))
stopifnot(grepl("must match tree tip labels", conditionMessage(leaf_regimes_err), fixed=TRUE))
pcm_reg_bad_shift = pcm_reg
pcm_reg_bad_shift$shift.configuration = structure(999L, names="regA")
leaf_regimes_shift_err = tryCatch(
    {
        get_leaf_regimes(pcm_reg_bad_shift, mode="l1ou")
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(leaf_regimes_shift_err))
stopifnot(grepl("invalid edge index", conditionMessage(leaf_regimes_shift_err), fixed=TRUE))
leaf_regimes_na_shift_err = tryCatch(
    {
        get_leaf_regimes(pcm_reg_na_shift, mode="l1ou")
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(leaf_regimes_na_shift_err))
stopifnot(grepl("must contain integer edge indices", conditionMessage(leaf_regimes_na_shift_err), fixed=TRUE))
pcm_leaf_no_col = list(
    tree=tr_unlabeled,
    Y=as.data.frame(matrix(c(1, 2, 3), ncol=1, dimnames=list(tr_unlabeled$tip.label, NULL))),
    shift.configuration=structure(integer(0), names=NULL),
    optima=as.data.frame(matrix(c(1, 2, 3), ncol=1, dimnames=list(tr_unlabeled$tip.label, NULL))),
    mu=as.data.frame(matrix(c(1, 2, 3), ncol=1, dimnames=list(tr_unlabeled$tip.label, NULL))),
    residuals=as.data.frame(matrix(c(1, 2, 3), ncol=1, dimnames=list(tr_unlabeled$tip.label, NULL)))
)
colnames(pcm_leaf_no_col$Y) = NULL
colnames(pcm_leaf_no_col$optima) = NULL
colnames(pcm_leaf_no_col$mu) = NULL
colnames(pcm_leaf_no_col$residuals) = NULL
leaf_tbl_no_col = get_leaf_table(pcm_leaf_no_col, mode="l1ou")
stopifnot(!any(is.na(colnames(leaf_tbl_no_col))))
stopifnot("trait1" %in% colnames(leaf_tbl_no_col))
leaf_tbl_mode_na_err = tryCatch(
    {
        get_leaf_table(pcm_leaf_no_col, mode=NA_character_)
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(leaf_tbl_mode_na_err))
stopifnot(grepl("mode must be a single non-missing string", conditionMessage(leaf_tbl_mode_na_err), fixed=TRUE))
pcm_leaf_reordered = list(
    tree=tr_unlabeled,
    Y=data.frame(t1=c(30, 20, 10), row.names=rev(tr_unlabeled$tip.label)),
    shift.configuration=structure(integer(0), names=NULL),
    optima=data.frame(t1=c(30, 20, 10), row.names=rev(tr_unlabeled$tip.label)),
    mu=data.frame(t1=c(30, 20, 10), row.names=rev(tr_unlabeled$tip.label)),
    residuals=data.frame(t1=c(30, 20, 10), row.names=rev(tr_unlabeled$tip.label))
)
leaf_tbl_reordered = get_leaf_table(pcm_leaf_reordered, mode="l1ou")
leaf_tbl_reordered_y = leaf_tbl_reordered[as.character(leaf_tbl_reordered$param) == "Y",,drop=FALSE]
stopifnot(identical(as.character(leaf_tbl_reordered_y$node_name), as.character(tr_unlabeled$tip.label)))
stopifnot(identical(as.numeric(leaf_tbl_reordered_y$t1), c(10, 20, 30)))
tree_table_collapsed = data.frame(
    num_shift=0,
    num_regime=1,
    num_conv_regime=0,
    num_uniq_regime=1,
    num_species=0,
    num_leaf=0,
    model_score=0,
    stringsAsFactors=FALSE
)
tree_table_restored = tree_table_collapse2original(tree_table_collapsed, tr_unlabeled)
stopifnot(identical(as.integer(tree_table_restored$num_leaf), 3L))
stopifnot(identical(as.integer(tree_table_restored$num_species), 3L))

pcm_bt = list(tree=tr_unlabeled)
bt_tbl = get_bootstrap_table(pcm_bt, list(detection.rate=c(1, 2, 3, 4)), mode="l1ou")
stopifnot(nrow(bt_tbl) == (length(tr_unlabeled$tip.label) + tr_unlabeled$Nnode))
stopifnot(sum(is.na(bt_tbl$bootstrap_support)) == 1)
bt_mode_na_err = tryCatch(
    {
        get_bootstrap_table(pcm_bt, list(detection.rate=c(1, 2, 3, 4)), mode=NA_character_)
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(bt_mode_na_err))
stopifnot(grepl("mode must be a single non-missing string", conditionMessage(bt_mode_na_err), fixed=TRUE))
bt_tbl_factor = get_bootstrap_table(
    pcm_bt,
    list(detection.rate=factor(c("0.1", "0.2", "0.3", "0.4"))),
    mode="l1ou"
)
stopifnot(isTRUE(all.equal(
    as.numeric(bt_tbl_factor$bootstrap_support),
    c(0.1, 0.2, 0.3, NA, 0.4),
    tolerance=1e-10
)))
bt_non_numeric_err = tryCatch(
    {
        get_bootstrap_table(
            pcm_bt,
            list(detection.rate=c("0.1", "bad", "0.3", "0.4")),
            mode="l1ou"
        )
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(bt_non_numeric_err))
stopifnot(grepl("must be numeric or coercible to numeric", conditionMessage(bt_non_numeric_err), fixed=TRUE))
bt_err = tryCatch(
    {
        get_bootstrap_table(pcm_bt, list(detection.rate=1:5), mode="l1ou")
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(bt_err))
stopifnot(grepl("Length mismatch", conditionMessage(bt_err), fixed=TRUE))

tree_collapsed = ape::read.tree(text="(1:1,C:1);")
node_num_mapping = data.frame(tree_original=c(1L, 2L, 3L), tree_collapsed=c(1L, 1L, 2L))
regime_input = data.frame(regime="r1", node_name="1", param="shift_value", t1=0.1, stringsAsFactors=FALSE)
regime_converted = regime_table_collapse2original(
    regime_table=regime_input,
    tree_original=tr_unlabeled,
    tree_collapsed=tree_collapsed,
    node_num_mapping=node_num_mapping
)
stopifnot(nrow(regime_converted) == 1)
stopifnot(regime_converted$node_name %in% c("A", "B"))

leaf_input = data.frame(
    regime=c("r1", "r1", NA),
    node_name=c("1", "C", NA),
    param=c("imputed", "imputed", NA),
    t1=c(10, 20, 30),
    stringsAsFactors=FALSE
)
leaf_converted = leaf_table_collapse2original(
    leaf_table=leaf_input,
    tree_original=tr_unlabeled,
    tree_collapsed=tree_collapsed,
    node_num_mapping=node_num_mapping
)
stopifnot(!any(is.na(leaf_converted$param)))
stopifnot(nrow(leaf_converted) == 3)

restored = restore_imputed_leaves(
    leaf_table=data.frame(regime=0, node_name="A", param="imputed", trait1=NA_real_, stringsAsFactors=FALSE),
    original_trait_table=data.frame(trait1=9, row.names="A")
)
stopifnot(identical(as.numeric(restored$trait1), 9))

placeholder_leaf = get_placeholder_leaf(
    tree=tr_unlabeled,
    original_trait_table=data.frame(trait1=c(5, 6, 7), row.names=tr_unlabeled$tip.label)
)
stopifnot("node_name" %in% colnames(placeholder_leaf))
stopifnot(!("label" %in% colnames(placeholder_leaf)))
placeholder_regime = get_placeholder_regime(
    tree=tr_unlabeled,
    original_trait_table=data.frame(trait1=c(5, 6, 7), row.names=tr_unlabeled$tip.label)
)
stopifnot("trait1" %in% colnames(placeholder_regime))
placeholder_regime_no_col = get_placeholder_regime(
    tree=tr_unlabeled,
    original_trait_table={
        tbl = as.data.frame(matrix(c(5, 6, 7), ncol=1, dimnames=list(tr_unlabeled$tip.label, NULL)))
        colnames(tbl) = NULL
        tbl
    }
)
stopifnot("trait1" %in% colnames(placeholder_regime_no_col))
placeholder_tree = get_placeholder_tree(
    tree=tr_unlabeled,
    original_trait_table=data.frame(trait1=c(5, 6, 7), row.names=tr_unlabeled$tip.label)
)
stopifnot("num_leaf" %in% colnames(placeholder_tree))

collapse_single_trait = collapse_clades(
    tr_unlabeled,
    data.frame(v=c(1, 2, 3), row.names=tr_unlabeled$tip.label),
    collapse_node_nums={
        internal_nodes = sort(unique(tr_unlabeled$edge[,1]))
        internal_nodes = internal_nodes[internal_nodes > length(tr_unlabeled$tip.label)]
        size_two_node = internal_nodes[vapply(internal_nodes, function(nn) {
            length(get_tip_labels(tr_unlabeled, nn)) == 2
        }, logical(1))]
        size_two_node[1]
    }
)
stopifnot(length(collapse_single_trait$tree$tip.label) == 2)
stopifnot(ncol(collapse_single_trait$trait) == 1)
stopifnot(all(collapse_single_trait$tree$tip.label %in% rownames(collapse_single_trait$trait)))

collapse_single_trait_dup = suppressWarnings(collapse_clades(
    tr_unlabeled,
    data.frame(v=c(1, 2, 3), row.names=tr_unlabeled$tip.label),
    collapse_node_nums=rep(as.integer(collapse_single_trait$tree$tip.label[collapse_single_trait$tree$tip.label != "C"]), 2)
))
stopifnot(length(collapse_single_trait_dup$tree$tip.label) == length(collapse_single_trait$tree$tip.label))
map_identity = map_node_num(tr_unlabeled, tr_unlabeled, collapse_leaf_names=list())
stopifnot(nrow(map_identity) == max(tr_unlabeled$edge[,1]))
stopifnot(identical(as.integer(map_identity$tree_original), as.integer(map_identity$tree_collapsed)))
map_identity_default = map_node_num(tr_unlabeled, tr_unlabeled)
stopifnot(identical(as.integer(map_identity_default$tree_original), as.integer(map_identity_default$tree_collapsed)))
map_verbose_na_err = tryCatch(
    {
        map_node_num(tr_unlabeled, tr_unlabeled, verbose=NA)
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(map_verbose_na_err))
stopifnot(grepl("verbose must be a single non-missing logical value", conditionMessage(map_verbose_na_err), fixed=TRUE))
map_tip_mismatch_err = tryCatch(
    {
        map_node_num(
            tr_unlabeled,
            ape::read.tree(text="((A:1,D:1):1,C:1);")
        )
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(map_tip_mismatch_err))
stopifnot(grepl("Tip labels are inconsistent", conditionMessage(map_tip_mismatch_err), fixed=TRUE))
collapse_single_trait_chr = collapse_clades(
    tr_unlabeled,
    data.frame(v=c(1, 2, 3), row.names=tr_unlabeled$tip.label),
    collapse_node_nums=as.character({
        internal_nodes = sort(unique(tr_unlabeled$edge[,1]))
        internal_nodes = internal_nodes[internal_nodes > length(tr_unlabeled$tip.label)]
        size_two_node = internal_nodes[vapply(internal_nodes, function(nn) {
            length(get_tip_labels(tr_unlabeled, nn)) == 2
        }, logical(1))]
        size_two_node[1]
    })
)
stopifnot(length(collapse_single_trait_chr$tree$tip.label) == 2)
collapse_invalid_type_error = tryCatch(
    {
        collapse_clades(
            tr_unlabeled,
            data.frame(v=c(1, 2, 3), row.names=tr_unlabeled$tip.label),
            collapse_node_nums="x"
        )
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(collapse_invalid_type_error))
stopifnot(grepl("must be integer node numbers", conditionMessage(collapse_invalid_type_error), fixed=TRUE))
tr_nested = ape::read.tree(text="((A:1,B:1):1,(C:1,D:1):1);")
collapse_nested = suppressWarnings(collapse_clades(
    tr_nested,
    data.frame(v=c(1, 2, 3, 4), row.names=tr_nested$tip.label),
    collapse_node_nums=c(5L, 7L)
))
stopifnot(length(collapse_nested$collapse_leaf_names) == 1)
stopifnot(identical(names(collapse_nested$collapse_leaf_names), "5"))
stopifnot(setequal(collapse_nested$collapse_leaf_names[[1]], c("A", "B", "C", "D")))
high_sim_num_test_na_err = tryCatch(
    {
        get_high_similarity_clades(
            tr_unlabeled,
            data.frame(v=c(1, 2, 3), row.names=tr_unlabeled$tip.label),
            method="pearson",
            threshold=0.5,
            num_test=NA
        )
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(high_sim_num_test_na_err))
stopifnot(grepl("single non-negative integer", conditionMessage(high_sim_num_test_na_err), fixed=TRUE))
high_sim_verbose_na_err = tryCatch(
    {
        get_high_similarity_clades(
            tr_unlabeled,
            data.frame(v=c(1, 2, 3), row.names=tr_unlabeled$tip.label),
            method="pearson",
            threshold=0.5,
            verbose=NA
        )
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(high_sim_verbose_na_err))
stopifnot(grepl("verbose must be a single non-missing logical value", conditionMessage(high_sim_verbose_na_err), fixed=TRUE))
const_trait_warns = character(0)
collapse_const_trait = withCallingHandlers(
    collapse_clades(
        tr_nested,
        data.frame(v=rep(5, 4), row.names=tr_nested$tip.label),
        collapse_node_nums=5L
    ),
    warning=function(w) {
        const_trait_warns <<- c(const_trait_warns, conditionMessage(w))
        invokeRestart("muffleWarning")
    }
)
stopifnot(!any(grepl("NaNs produced", const_trait_warns, fixed=TRUE)))
stopifnot(isTRUE(all.equal(as.numeric(collapse_const_trait$trait["5", "v"]), 5, tolerance=1e-10)))

collapse_leaf_error = tryCatch(
    {
        collapse_clades(tr_unlabeled, data.frame(v=c(1,2,3), row.names=tr_unlabeled$tip.label), collapse_node_nums=c(1))
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(collapse_leaf_error))
stopifnot(grepl("collapse_node_nums must contain only internal node numbers", conditionMessage(collapse_leaf_error), fixed=TRUE))
dup_trait_rows = matrix(c(1, 2, 3, 4), ncol=1)
rownames(dup_trait_rows) = c("A", "A", "B", "C")
dup_trait_high_sim_err = tryCatch(
    {
        get_high_similarity_clades(
            tr_unlabeled,
            dup_trait_rows,
            method="pearson",
            threshold=0.5
        )
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(dup_trait_high_sim_err))
stopifnot(grepl("duplicated row name", conditionMessage(dup_trait_high_sim_err), fixed=TRUE))
dup_trait_collapse_err = tryCatch(
    {
        collapse_clades(
            tr_unlabeled,
            dup_trait_rows,
            collapse_node_nums=5L
        )
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(dup_trait_collapse_err))
stopifnot(grepl("duplicated row name", conditionMessage(dup_trait_collapse_err), fixed=TRUE))

tbl = data.frame(
    branch_id=c(1L, 2L, 3L),
    parent=c(3L, 3L, -999L),
    sister=c(2L, 1L, -999L),
    label=c("A", "B", "Root"),
    dist=c(0.1, 0.2, 0.0),
    stringsAsFactors=FALSE
)
tbl_phy = table2phylo(tbl, name_col="label", dist_col="dist")
stopifnot(inherits(tbl_phy, "phylo"))
stopifnot(setequal(tbl_phy$tip.label, c("A", "B")))
tbl_phy_dist = cophenetic(tbl_phy)[c("A", "B"), c("A", "B")]
tbl_phy_expected = matrix(c(0, 0.3, 0.3, 0), nrow=2, byrow=TRUE, dimnames=list(c("A", "B"), c("A", "B")))
stopifnot(isTRUE(all.equal(tbl_phy_dist, tbl_phy_expected, tolerance=1e-10)))

tbl_root_dist_na = tbl
tbl_root_dist_na$dist[tbl_root_dist_na$branch_id == 3L] = NA_real_
tbl_root_dist_na_phy = table2phylo(tbl_root_dist_na, name_col="label", dist_col="dist")
stopifnot(inherits(tbl_root_dist_na_phy, "phylo"))
stopifnot(setequal(tbl_root_dist_na_phy$tip.label, c("A", "B")))

tbl_parent_na = data.frame(
    branch_id=c(10L, 20L, 5L),
    parent=c(5L, 5L, NA),
    sister=c(20L, 10L, -999L),
    label=c("A", "B", "Root"),
    dist=c(0.1, 0.2, NA_real_),
    stringsAsFactors=FALSE
)
tbl_parent_na_phy = table2phylo(tbl_parent_na, name_col="label", dist_col="dist")
stopifnot(inherits(tbl_parent_na_phy, "phylo"))
stopifnot(setequal(tbl_parent_na_phy$tip.label, c("A", "B")))
tbl_factor_ids = data.frame(
    branch_id=factor(c(1L, 2L, 3L)),
    parent=factor(c(3L, 3L, -999L)),
    sister=factor(c(2L, 1L, -999L)),
    label=c("A", "B", "Root"),
    dist=c(0.1, 0.2, 0.0),
    stringsAsFactors=TRUE
)
tbl_factor_ids_phy = table2phylo(tbl_factor_ids, name_col="label", dist_col="dist")
stopifnot(inherits(tbl_factor_ids_phy, "phylo"))
stopifnot(setequal(tbl_factor_ids_phy$tip.label, c("A", "B")))
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
stopifnot(!is.null(tbl_missing_branch_id_err))
stopifnot(grepl("branch_id contains missing/blank", conditionMessage(tbl_missing_branch_id_err), fixed=TRUE))
tbl_single_tip = data.frame(
    branch_id=1L,
    parent=-999L,
    sister=-999L,
    label="A",
    dist=0.1,
    stringsAsFactors=FALSE
)
tbl_single_tip_phy = table2phylo(tbl_single_tip, name_col="label", dist_col="dist")
stopifnot(inherits(tbl_single_tip_phy, "phylo"))
stopifnot(identical(as.character(tbl_single_tip_phy$tip.label), "A"))
tbl_single_tip_root = data.frame(
    branch_id=c(1L, 2L),
    parent=c(2L, -999L),
    sister=c(-999L, -999L),
    label=c("A", "Root"),
    dist=c(0.1, 0.0),
    stringsAsFactors=FALSE
)
tbl_single_tip_root_phy = table2phylo(tbl_single_tip_root, name_col="label", dist_col="dist")
stopifnot(inherits(tbl_single_tip_root_phy, "phylo"))
stopifnot(identical(as.character(tbl_single_tip_root_phy$tip.label), "A"))
tbl_partial_parent = data.frame(
    branch_id=c("11", "d", "1", "3"),
    parent=c("12", "12", "", "a"),
    sister=c("3", "1", "2", "4"),
    label=c("J", "F", "O", "O"),
    dist=c("NaN", "Inf", "x", "0"),
    stringsAsFactors=FALSE
)
tbl_partial_parent_phy = table2phylo(tbl_partial_parent, name_col="label", dist_col="dist")
stopifnot(inherits(tbl_partial_parent_phy, "phylo"))
tr_two_tip = ape::read.tree(text="(A:1,B:1);")
tr_two_tip_rr = remove_redundant_root_edge(tr_two_tip)
stopifnot(inherits(tr_two_tip_rr, "phylo"))
stopifnot(setequal(tr_two_tip_rr$tip.label, c("A", "B")))
stopifnot(identical(as.integer(tr_two_tip_rr$Nnode), 1L))
tr_unary_root = ape::read.tree(text="((A:1):1,B:1);")
tr_unary_root_rr = remove_redundant_root_edge(tr_unary_root)
stopifnot(inherits(tr_unary_root_rr, "phylo"))
stopifnot(setequal(tr_unary_root_rr$tip.label, c("A", "B")))
stopifnot(identical(as.integer(tr_unary_root_rr$Nnode), 1L))
tr_from_label = ape::read.tree(text="((A:1,B:1):1,C:1);")
tr_from_label$node.label = c("X", "Y")
tr_to_unlabeled = ape::read.tree(text="((A:1,B:1):1,C:1);")
tr_to_unlabeled$node.label = NULL
tr_transferred = transfer_node_labels(tr_from_label, tr_to_unlabeled)
stopifnot(identical(as.character(tr_transferred$node.label), c("X", "Y")))
transfer_tip_mismatch_err = tryCatch(
    {
        transfer_node_labels(tr_from_label, ape::read.tree(text="((A:1,D:1):1,C:1);"))
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(transfer_tip_mismatch_err))
stopifnot(grepl("must contain the same tip labels", conditionMessage(transfer_tip_mismatch_err), fixed=TRUE))

tbl_multi_root = tbl_parent_na
tbl_multi_root$parent[1] = -999L
multi_root_error = tryCatch(
    {
        table2phylo(tbl_multi_root, name_col="label", dist_col="dist")
        NULL
    },
    error=function(e) e
)
stopifnot(!is.null(multi_root_error))
stopifnot(grepl("Ambiguous root candidate", conditionMessage(multi_root_error), fixed=TRUE))

old_max_cores = getOption("rkftools.max_cores")
options(rkftools.max_cores=1L)
mad_auto_capped = MAD_parallel(tr, output_mode="newick")
stopifnot(identical(mad_n1, mad_auto_capped))
options(rkftools.max_cores=old_max_cores)
