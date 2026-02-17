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

old_max_cores = getOption("rkftools.max_cores")
options(rkftools.max_cores=1L)
mad_auto_capped = MAD_parallel(tr, output_mode="newick")
stopifnot(identical(mad_n1, mad_auto_capped))
options(rkftools.max_cores=old_max_cores)
