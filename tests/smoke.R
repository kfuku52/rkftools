library(rkftools)

comp = calc_complementarity(c(1, 2, 3), c(1, 2, 0), method="weighted")
stopifnot(isTRUE(all.equal(comp, 1 / 3, tolerance=1e-10)))

tr = ape::read.tree(text="((A_a_x:1,B_b_y:1):1,(A_a_z:1,C_c_w:1):1);")
so = get_species_overlap_score(tr, dc_cutoff=0)
stopifnot(identical(as.numeric(so), 1))

so_by_root = get_root_position_dependent_species_overlap_scores(tr, nslots=1)
stopifnot(length(so_by_root) == nrow(tr$edge))

rerooted = phytools::reroot(tree=tr, node.number=tr$edge[2,2])
root_idx = get_phy2_root_in_phy1(tr, rerooted, nslots=1, mode="index")
stopifnot(!is.na(root_idx))

mad_n1 = MAD_parallel(tr, output_mode="newick", ncpu=1)
mad_auto = MAD_parallel(tr, output_mode="newick")
stopifnot(length(mad_n1) == length(mad_auto))

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
