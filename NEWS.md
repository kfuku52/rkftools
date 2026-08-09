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
