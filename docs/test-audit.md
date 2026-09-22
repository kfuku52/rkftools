# Test audit (0.1.13)

Reviewed every `tests/testthat/test-*.R` file, the shared fixtures, the test
entrypoint, and the implementations behind overlapping checks. Keep a test when
its removal would lose a concrete failure signal; assertion count and coverage
are not retention criteria. Production code and CI thresholds are unchanged.

| Area / test files | Decision and remaining failure signal |
| --- | --- |
| `bootstrap-regression` | Remove the plain numeric shape check and repeated mode validation. The factor-input result checks actual support placement, including the root's missing value; invalid values and length mismatches remain. |
| `branch-tables-regression`, `table-roundtrip` | Use the nontrivial four-tip table for binary distances and table values; custom-column contracts retain schema checks. Remove the two-tip duplicate, repeated shape/coercion checks, and ape's own Newick round-trip. Retain unary nodes, root/zero/tiny lengths, factor IDs, malformed graphs, custom columns, and label collisions. |
| `branch-length-regression`, `edge-padding`, `polytomy` | Remove class-only smoke checks, repeated padding cases, and ultrametric checks already implied by explicit equal tip depths. One seeded path-bound example replaces twelve iterations; keep hand-built borrowing, missing-length, unary, and multifurcating cases. Retain negative-edge collapse, GLS, and root-extension behavior. |
| `clade-collapse-regression`, `collapse-mapping` | Remove ordinary identity mappings covered by the harder unary identity/reordering case, and repeated coercion/validator checks. Retain full collapse, synthetic roots, duplicate/nested collapse requests, constant traits, and undefined correlations. |
| `rooting-regression`, `root-mapping`, `mad-contract` | Remove length-only MAD/worker smoke tests and three sizes of the same mapping scenario. Retain one nontrivial root match, legacy argument equivalence, actual serial/parallel MAD output equality, zero-distance duplicates, output edge consistency, and root-label transfer. |
| `reconciliation-regression`, `tree-structure` | Remove assertions about fixture indexing and arrays computed solely by the test. One seeded rerooting comparison retains an independent reference for dynamic root scoring; fixed biological scores and the deep-tree recursion regression remain. Test cyclic-graph rejection through traversal rather than repeating it through three entrypoints. |
| `species-regression`, `foreground-regression` | Remove repeated scalar/vector examples, retaining distinct qualifier grammar, separators, malformed labels, lineage merging, and missing/unknown species. |
| `io-regression`, `io-and-quiet` | Combine ordinary and whitespace-prefixed NOTUNG records. Retain malformed/missing-file checks, quiet defaults, leading-zero identifiers, duplicate arguments, and credential redaction. |
| `l1ou-regression`, `leaf-tables-regression` | Remove repeated mode checks, redundant regime membership assertions, and simple observation restoration covered by mixed missing/observed cells. Retain shift ordering, trait-name alignment, scalar parameter expansion, collapsed outputs, and placeholders. |
| `traits-regression`, `trait-contracts`, `optional-integrations` | Remove duplicated all-missing replicate and permissive imputation smoke tests. Replace the Rphylopars internal-argument mock with reordered matrix input in the real integration. Remove the fixture-version assertion. Keep the PhylogeneticEM mock for reversed trait names and scalar expansion, which the fitted fixture does not exercise; its expected transformed values are independent of the mock return layout. |
| `traversal-regression`, `validation` | Remove fixture-only checks and repeated shared logical/string validator permutations. Retain ordered traversal, nearest-tip lookup, node ages, integer-node validation, cutoff boundaries, and representative strict logical controls. |

Error cases that remain use `expect_error()` directly, replacing a captured
error plus separate existence/message checks. The one-line custom error helper
is removed. Real optional-backend tests still skip explicitly in a minimal
installation and are required by the existing full delivery check.
