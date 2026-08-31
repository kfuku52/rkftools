# Usage

## Convert between Newick, ape::phylo, and a branch table

`table2phylo()` and `phylo2table()` convert between a branch table and an
`ape::phylo` object. Newick strings or files can be converted to and from
`ape::phylo` with `ape::read.tree()` and `ape::write.tree()`. `branch_id`,
`parent`, and `sister` are numerical labels. `parent` and `sister` refer to
`branch_id`, while `node_name` stores the tip or internal-node label.
`phylo2table()` returns the same schema with numerical labels generated from
clade signatures, matching genegalleon's `numerical_label` convention.
`table2phylo()` validates that the table describes one connected rooted tree.
Binary children use reciprocal sister references; unary and multifurcating
children use sister sentinels. Exact zero-length branches are preserved.

```r
library(rkftools)

branch_table = data.frame(
    branch_id=c(6, 2, 0, 1, 5, 3, 4),
    parent=c(-999, 6, 2, 2, 6, 5, 5),
    sister=c(-999, 5, 1, 0, 2, 4, 3),
    node_name=c("Root", "Clade_AB", "A", "B", "Clade_CD", "C", "D"),
    dist=c(0, 0.42, 0.16, 0.18, 0.50, 0.22, 0.25),
    stringsAsFactors=FALSE
)

tree = table2phylo(branch_table, name_col="node_name", dist_col="dist")
newick = ape::write.tree(tree)
tree_from_newick = ape::read.tree(text=newick)
roundtrip_table = phylo2table(tree_from_newick, name_col="node_name", dist_col="dist")
```

![Newick, ape::phylo, and branch table conversion example](../man/figures/table2phylo_roundtrip.png)

## Score candidate root positions by species overlap

```r
library(rkftools)

gene_tree = ape::read.tree(text=paste0(
    "(((Homo_sapiens_geneA:0.12,Mus_musculus_geneA:0.12):0.18,",
    "(Danio_rerio_geneA:0.16,Oryzias_latipes_geneA:0.16):0.14):0.25,",
    "((Homo_sapiens_geneB:0.11,Mus_musculus_geneB:0.11):0.19,",
    "(Danio_rerio_geneB:0.15,Oryzias_latipes_geneB:0.15):0.15):0.25);"
))

root_scores = get_root_position_dependent_species_overlap_scores(
    gene_tree,
    nslots=1
)
```

The figure uses unique branch IDs after unrooting the tree; the minimum-score
candidate root branch is branch 7.

![root-position species-overlap score example](../man/figures/root_position_species_overlap.png)

# Parallel tuning
`MAD_parallel()` automatically uses available CPU cores when `ncpu` is omitted,
caps automatic parallelism at eight cores, and avoids parallel overhead for
small trees. Root-position species-overlap scoring now uses one bidirectional
tree traversal, so its legacy `nslots` argument is accepted for compatibility
but no worker pool is needed. `get_phy2_root_in_phy1()` likewise matches edge
bipartitions in one traversal and retains `nslots` only for compatibility.

To cap cores globally:
```r
options(rkftools.max_cores = 8)
```

# Input validation and optional integrations

Tree and trait helpers reject disconnected trees, duplicated tip labels,
malformed species identifiers, duplicated trait rows, and non-finite values
where these would make a result unreliable. Repeated internal-node labels such
as bootstrap values remain supported. `get_parsed_args()` does not print by
default; when printing is requested, credential-like values are redacted.

Functions backed by optional packages produce an installation error naming the
required package. Install `PhylogeneticEM` or `Rphylopars` when using their
corresponding modes or imputation helper.

# Transformation contracts

`pad_short_edges()` transfers length through every ancestor while keeping all
root-to-tip distances unchanged when possible. If padding requires moving the
root, every finite root-to-tip distance increases by the same minimum amount.
`external_only=TRUE` applies the threshold to tips; internal branches still
transfer length and remain non-negative. Missing lengths remain missing and do
not constrain transfers across those edges.

`map_node_num()` distinguishes tips from unary internal nodes and matches unary
chains in ancestor order. `collapse_clades()` records original node numbers as
collapse-map keys; preserve those keys when passing its result to
`map_node_num()`. A completely collapsed tree maps to its single tip, rather
than the artificial root needed by the phylo representation.

MAD treats zero-distance tip groups as one representative in its scoring
objective, while preserving all tips in every returned tree. In `full` and
`custom`, root indices, per-edge deviations, and root proportions correspond to
the edge rows of the returned unrooted tree. Clock CV uses all returned tips.

Model summaries accept both gene labels and species-only labels. As in the
existing l1ou adapter, labels without a recognizable genus/species pair count
verbatim; missing labels are rejected instead of being counted as a species.
Trait imputation accepts numeric matrices and data frames. Restoring observed
leaves only replaces cells that had non-missing observations.
