\name{rkftools-functions}
\alias{rkftools-functions}
\alias{calc_complementarity}
\alias{collapse_clades}
\alias{collapse_short_branches}
\alias{collapse_short_external_edges}
\alias{contains_polytomy}
\alias{count_foreground_lineage}
\alias{fill_node_labels}
\alias{force_ultrametric}
\alias{get_ancestor_num}
\alias{get_bootstrap_table}
\alias{get_children_num}
\alias{get_deepest_node_num}
\alias{get_descendent_num}
\alias{get_duplication_confidence_score}
\alias{get_expression_bases}
\alias{get_high_similarity_clades}
\alias{get_leaf_regimes}
\alias{get_leaf_table}
\alias{get_nearest_tips}
\alias{get_node_age}
\alias{get_node_name_by_num}
\alias{get_node_num_by_name}
\alias{get_outgroup}
\alias{get_parent_num}
\alias{get_parsed_args}
\alias{get_phy2_root_in_phy1}
\alias{get_placeholder_leaf}
\alias{get_placeholder_regime}
\alias{get_placeholder_tree}
\alias{get_regime_table}
\alias{get_root_num}
\alias{get_root_position_dependent_species_overlap_scores}
\alias{get_rooted_newick}
\alias{get_single_branch_tree}
\alias{get_sister_num}
\alias{get_species_name}
\alias{get_species_names}
\alias{get_species_overlap_score}
\alias{get_tip_labels}
\alias{get_tree_table}
\alias{has_same_leaves}
\alias{is_leaf}
\alias{is_root}
\alias{is_same_root}
\alias{is.blank}
\alias{leaf_table_collapse2original}
\alias{leaf2species}
\alias{MAD}
\alias{MAD_parallel}
\alias{map_node_num}
\alias{merge_replicates}
\alias{multi2bi_node_number_transfer}
\alias{pad_short_edges}
\alias{phylogenetic_imputation}
\alias{phylo2table}
\alias{read_notung_parsable}
\alias{regime_table_collapse2original}
\alias{remove_invariant_traits}
\alias{remove_redundant_root_edge}
\alias{restore_imputed_leaves}
\alias{sort_exp}
\alias{table2phylo}
\alias{transfer_node_labels}
\alias{tree_table_collapse2original}
\title{Evolutionary Biology and Phylogenetic Utilities}
\description{
Functions for phylogenetic tree traversal and conversion, gene-tree
reconciliation, comparative-method output processing, and trait-table
preparation.
}
\usage{
calc_complementarity(array1, array2, method = "weighted")
remove_invariant_traits(trait_table, small_dif = 0.001)
merge_replicates(trait_table, replicate_sep)
get_high_similarity_clades(tree, trait_table, method, threshold,
  verbose = FALSE, num_test = 0)
collapse_clades(tree, trait_table, collapse_node_nums)
map_node_num(tree_original, tree_collapsed, collapse_leaf_names = list(),
  verbose = FALSE)
get_tree_table(pcm_out, mode, species_parser = "legacy", sep = "_")
get_regime_table(pcm_out, mode)
get_leaf_regimes(pcm_out, mode)
get_leaf_table(pcm_out, mode)
get_bootstrap_table(pcm_out, bootstrap_result, mode = "l1ou")
tree_table_collapse2original(tree_table, tree_original,
  species_parser = "legacy", sep = "_")
get_deepest_node_num(tree, node_nums)
regime_table_collapse2original(regime_table, tree_original, tree_collapsed,
  node_num_mapping)
leaf_table_collapse2original(leaf_table, tree_original, tree_collapsed,
  node_num_mapping)
restore_imputed_leaves(leaf_table, original_trait_table)
get_placeholder_leaf(tree, original_trait_table)
get_placeholder_regime(tree, original_trait_table)
get_placeholder_tree(tree, original_trait_table)
sort_exp(exp, tree, col = "gene_id")
phylogenetic_imputation(tree, trait_table)
get_expression_bases(trait_table, replicate_sep)
count_foreground_lineage(tree, trait)

get_node_num_by_name(phy, node_name)
get_node_name_by_num(phy, node_num)
get_root_num(phy)
get_children_num(phy, node_num)
is_root(phy, node_num)
is_leaf(phy, node_num)
get_descendent_num(phy, node_num, leaf_only = FALSE)
get_parent_num(phy, node_num)
get_sister_num(phy, node_num)
get_ancestor_num(phy, node_num)
collapse_short_external_edges(tree, threshold = 1e-6)
pad_short_edges(tree, threshold = 1e-6, external_only = FALSE)
get_tip_labels(phy, node_num, out = NULL)
get_nearest_tips(phy, query, subjects, mrca_matrix)
get_node_age(phy, node_num)
get_outgroup(phy)
is_same_root(phy1, phy2)
get_phy2_root_in_phy1(phy1, phy2, nslots = NULL,
  mode = c("node_num", "index"))
get_rooted_newick(t, madr, rho)
MAD(unrooted_newick, output_mode)
MAD_parallel(unrooted_newick, output_mode, ncpu = NULL)
transfer_node_labels(phy_from, phy_to)
get_species_name(a, species_parser = "legacy", sep = "_")
get_species_names(phy, sep = "_", species_parser = "legacy")
leaf2species(leaf_names, use_underbar = FALSE,
  species_parser = "legacy", sep = "_")
contains_polytomy(phy)
has_same_leaves(phy1, phy1_node, phy2, phy2_node)
multi2bi_node_number_transfer(multifurcated_tree, bifurcated_tree)
collapse_short_branches(tree, tol = 1e-8)
force_ultrametric(tree, stop_if_larger_change = 0.01)
get_single_branch_tree(name, dist)
remove_redundant_root_edge(phy)
table2phylo(df, name_col, dist_col)
phylo2table(phy, name_col = "label", dist_col = "dist")
fill_node_labels(phy)

get_duplication_confidence_score(phy, node_num,
  species_parser = "legacy", sep = "_")
get_species_overlap_score(phy, dc_cutoff = 0,
  species_parser = "legacy", sep = "_")
get_root_position_dependent_species_overlap_scores(phy, nslots = NULL,
  species_parser = "legacy", sep = "_")
read_notung_parsable(file, mode = "D")

get_parsed_args(args, print = FALSE)
is.blank(x, false.triggers = FALSE)
}
\arguments{
\item{array1, array2}{Equal-length, non-negative numeric vectors.}
\item{method}{A supported complementarity or correlation method.}
\item{trait_table, original_trait_table}{A numeric trait matrix or data frame,
usually with tree tip labels as row names.}
\item{small_dif}{Non-negative range below which a trait is invariant.}
\item{replicate_sep}{Literal separator before replicate suffixes.}
\item{tree, phy, phy1, phy2, phy_from, phy_to, multifurcated_tree,
bifurcated_tree}{An object of class \code{"phylo"}.}
\item{threshold, tol}{Finite branch-length or similarity threshold.}
\item{external_only}{Whether short-edge padding is restricted to terminal
branches.}
\item{verbose}{Whether to print progress information.}
\item{num_test}{Optional maximum number of similarity comparisons; zero means
no limit.}
\item{collapse_node_nums}{Internal nodes whose clades should be collapsed.}
\item{tree_original, tree_collapsed}{Trees before and after clade collapse.}
\item{collapse_leaf_names}{Named descriptions of collapsed leaf groups.}
\item{pcm_out}{An output object from a supported comparative-method package.}
\item{mode}{A supported processing, lookup, or NOTUNG event mode.}
\item{species_parser}{Species-label convention: \code{"legacy"} or
\code{"taxonomic"}.}
\item{sep}{Literal separator used in gene labels.}
\item{bootstrap_result}{Bootstrap result containing detection rates.}
\item{node_nums, node_num, phy1_node, phy2_node}{One or more phylogenetic node
numbers.}
\item{regime_table, leaf_table, tree_table}{Long-format comparative-method
result tables.}
\item{node_num_mapping}{Mapping between original and collapsed node numbers.}
\item{exp}{Expression data frame to order by tree tips.}
\item{col}{Name of the expression-table key column.}
\item{trait}{Data frame with a \code{species} column and foreground traits.}
\item{node_name}{One or more tip or internal-node labels.}
\item{leaf_only}{Whether descendant lookup should return only tips.}
\item{out}{Internal accumulator; normally leave as \code{NULL}.}
\item{query, subjects}{Tip labels used for nearest-tip lookup.}
\item{mrca_matrix}{Named matrix returned by \code{ape::mrca()}.}
\item{nslots, ncpu}{Requested worker count. Values are capped by task count and
\code{options("rkftools.max_cores")}. For root-position species-overlap scoring,
\code{nslots} is retained for compatibility and no worker pool is required.}
\item{t}{A phylogenetic tree.}
\item{madr}{Edge index selected by MAD rooting.}
\item{rho}{Root position proportions, one per edge.}
\item{unrooted_newick}{A Newick string, path, or \code{"phylo"} tree.}
\item{output_mode}{MAD output detail: \code{"newick"}, \code{"stats"},
\code{"full"}, or \code{"custom"}.}
\item{a}{A single gene label from which to parse a species.}
\item{leaf_names}{Gene labels from which to parse species names.}
\item{use_underbar}{Whether parsed species names retain underscores.}
\item{dc_cutoff}{Minimum duplication-confidence score to count.}
\item{stop_if_larger_change}{Maximum allowed ultrametric adjustment.}
\item{name, dist}{Label and branch length for a single-tip tree.}
\item{df}{Branch table with \code{branch_id}, \code{parent}, \code{sister},
label, and distance columns.}
\item{name_col, dist_col}{Names of the label and distance columns.}
\item{file}{Path to a NOTUNG parsable-output file.}
\item{args}{Character vector of \code{--name} or \code{--name=value}
arguments.}
\item{print}{Whether parsed arguments should be printed; sensitive values are
redacted.}
\item{x}{Object tested for blank values.}
\item{false.triggers}{If true, logical \code{FALSE} is also blank.}
}
\details{
Branch-table conversion validates a connected rooted graph with at most two
children per node and reciprocal sister relationships. Species parsing requires
at least genus and species components for reconciliation functions. Optional
integrations raise an informative error when their package is not installed.
}
\value{
The return type depends on the function: a modified \code{"phylo"} object,
numeric or logical vector, mapping or result data frame, parsed argument list,
or a list containing transformed trait data. See each function's name and
examples for the expected family of result.
}
\examples{
tree <- ape::read.tree(text = "((A:1,B:1):1,C:1);")
get_root_num(tree)
get_tip_labels(tree, get_root_num(tree))

tab <- phylo2table(tree)
restored <- table2phylo(tab, "label", "dist")
stopifnot(setequal(restored$tip.label, tree$tip.label))

get_parsed_args(c("--threads=2", "--dry-run"))
}
\keyword{utilities}
