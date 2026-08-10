# rkftools 0.1.10

- Preserve zero-length root edges so root polytomies remain rooted across
  repeated branch-table conversions.
- Stabilize Brownian GLS ancestral-state estimates for very short and zero
  branches without discarding informative covariance directions.
- Restore v0.1.8 behavior for unrooted binary trees in branch-length
  transformations and reconciliation scores, and distinguish trivalent binary
  roots from true polytomies.
- Reduce root-position species-overlap scoring on high-degree nodes from cubic
  to quadratic time and add a star-tree benchmark.
- Align the README branch-table contract with multifurcation support.

# rkftools 0.1.9

- Preserve unary, binary, and multifurcating topologies in branch-table
  conversion, short-edge padding, and ultrametric adjustment.
- Generalize duplication-confidence and root-position species-overlap scoring
  to any number of child clades using maximum pairwise Jaccard overlap.
- Match multifurcating root partitions and map polytomy nodes into rooted binary
  refinements, including root polytomies.
- Estimate collapsed multifurcating clade traits with Brownian-motion GLS.
- Score multifurcating trees directly in MAD without random binary resolution.
- Add regression coverage for internal and root polytomies across conversions,
  reconciliation, branch lengths, ancestral traits, mapping, and rooting.

# rkftools 0.1.8

- Establish an 80% coverage floor from the initial 81.3% baseline.
- Make branch-table conversion symmetric for unary nodes, preserve root-edge
  lengths, and reject non-finite branch lengths.
- Add strict reusable validators for tree topology, integer node identifiers,
  logical controls, finite thresholds, and optional integration inputs.
- Replace repeated rerooting in `get_phy2_root_in_phy1()` with single-pass edge
  bipartition matching and cache pairwise trait similarities across clades.
- Centralize phylogenetic traversal indices and preallocate collapsed-table and
  placeholder-table construction.
- Make transformation helpers quiet by default with opt-in `verbose` messages,
  aggregate parsing warnings, and reject malformed NOTUNG records.
- Add focused testthat suites, optional-package integration checks, roxygen2
  function documentation, scheduled CRAN checks, R-devel CI, benchmarks, and
  release metadata automation.

# rkftools 0.1.7

- Validate branch tables before conversion and preserve zero-length branches.
- Fix trait-name alignment, replicate parsing, missing-value similarity, and
  collapsed-clade reconstruction.
- Replace repeated rerooting and clade extraction with linear-time tree
  traversals for root scoring and node-label transfer.
- Make short-branch collapse topology-aware so unrelated negative branches are
  not removed.
- Harden species parsing, command-line parsing, MAD edge cases, and input
  validation throughout the public API.
- Add optional dependency metadata, multi-platform CI, and expanded regression
  tests and documentation.
