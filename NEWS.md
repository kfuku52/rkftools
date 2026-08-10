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
