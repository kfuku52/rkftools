# Development

Run all commands below from the repository root (the directory containing
`DESCRIPTION` and `Makefile`). First check `R --version` and `Rscript --version`;
activate an existing compatible R environment or put its bin directory on PATH
if either is missing or has the wrong CPU architecture. `make setup-minimal`
installs R packages, not R itself. No CONTRIBUTING, lockfile, container, separate
lint target, or static type-check configuration is maintained here.

Use a working R installation with matching C/C++ and Fortran toolchains when
installing the optional backends from source. Required runtime packages remain
`ape` and `phytools`; development tools and backend packages are Suggests.

```sh
make setup-minimal  # Runtime and development packages, excluding optional backends
make test           # Every regression; absent optional backends explicitly skip
make check          # Build/check; permits absent Suggests
```

Setup adds missing packages to `.local/R-library`, which is ignored by Git and
excluded from package builds. Existing libraries are not overwritten. Use
`make DEV_LIB=/absolute/path/to/library ...` to choose another library. The
Makefile exports that library for every R process while preserving existing
library search paths. The setup, test, benchmark, fixture-generation, and
documentation-figure scripts also accept `RKFTOOLS_DEV_LIBRARY` when run
directly. The coverage script uses R's standard library paths instead: run
`make coverage`, or set `R_LIBS` to the development library for a standalone
`Rscript tools/coverage.R` invocation.

Setup requires network access only for missing packages; review missing
dependencies before starting an installation. `make setup` also installs
PhylogeneticEM and Rphylopars. Optional backend tests explicitly skip when those packages are absent; adapter contract
tests still run without them. `make check-full` and `make check-as-cran` require
all Suggests, matching the strict CI configuration. Missing dependencies are
errors in those checks, not silent skips.

`make document` regenerates the Rd files with roxygen2, with Markdown enabled
in DESCRIPTION so inline code and function references render as R help.
`make figures` executes the R examples in `docs/usage.md` and regenerates both
figures in `man/figures/`. It uses the actual tree labels, edge rows, branch
table, and computed scores, so figure values are not maintained separately.
Regenerate and visually inspect the figures after changing the examples or
tree conversion/scoring code. No plotting packages beyond base R and the
existing runtime dependencies are needed.

`make release-check` checks DESCRIPTION, NEWS, and the README badge. An explicit
tag argument to `tools/release-check.R` also verifies the tag version. Bump the
package version before pushing, following the repository policy.

## Choose verification for the change

`make test` loads this checkout, uses small synthetic inputs, and does not
fetch data or regenerate fitted fixtures. Installed optional backends are
exercised too (Rphylopars performs small fits). For a focused iteration use
`make test TEST_FILTER='table-roundtrip|branch-tables'`: the filter is a regex
on test filenames without `test-` and `.R`. An empty filter runs everything;
no matching files, failed assertions, and unexpected warnings fail the command.
For direct execution, use `RKFTOOLS_TEST_FILTER` with `Rscript tools/test.R`.
Filtering does not alter `R CMD check` or CI coverage.

| Changed area | Focused `TEST_FILTER` starting point |
| --- | --- |
| Branch-table conversion | `table-roundtrip|branch-tables` |
| Shared tree indexing, traversal, labels | `tree-structure|traversal|validation|root-mapping|collapse-mapping|table-roundtrip` |
| Branch transforms, collapse | `branch-length|edge-padding|polytomy|clade-collapse|collapse-mapping` |
| Rooting, MAD, reconciliation | `rooting|root-mapping|mad-contract|reconciliation` |
| Species/foreground, traits | `species|foreground|traits|trait-contracts` |
| Model adapters or imputation | `l1ou|leaf-tables|bootstrap|trait-contracts|optional-integrations` |
| Argument parsing or NOTUNG input | `io-` |

These are starting points, not coverage guarantees: read affected callers and
contracts, add a regression for changed behavior, then run unfiltered
`make test` for R/test/runner changes. Backend changes require their real
`optional-integrations` tests without missing-package skips.

For delivery of R code, tests, test tooling, package metadata, or generated
help, run `make check-full` and `make coverage` (80% floor) when all Suggests
are available. If dependencies are missing, run `make check`, report the
missing checks/skips, and leave full validation to the existing CI; do not
relax thresholds or silently treat minimal validation as full validation.
Plain prose/Skill edits need link/command review and execution of newly
introduced procedures, not fixture regeneration. Before a push, also run
`make release-check` after synchronizing DESCRIPTION, NEWS, and the README
version badge. `make check-as-cran` is the extended scheduled/manual check,
not the routine local iteration command.

Successful tests have no failures or unexpected warnings. Package checks
should have no errors or warnings; inspect and report any NOTEs. Coverage
must meet the printed floor and release-check must report consistent metadata.
There is no separate lint/type command; package checking covers R code
analysis, help consistency, and examples. Tests use temporary files for I/O.
Existing `make check*` targets write ignored tarballs and `rkftools.Rcheck/`
inside the working directory; `.agents/` is excluded from the built package.
Use a temporary source copy if those existing
artifacts must be preserved; never edit generated check copies as source.

Benchmarks, `--as-cran` network checks, optional dependency installation, and
`tools/generate-phyloem-fixture.R` are separate, deliberate operations. See
[performance](performance.md) for benchmarks. No large research dataset is
required for routine regression checks.

## Tests and fixtures

All tests live in `tests/testthat/`. The previous `tests/smoke.R` assertions are
organized by feature and run through both `make test` and `R CMD check`.
`helper-fixtures.R` provides small, independent starting trees. Backend adapters
also have focused contracts for trait alignment, missing cells, and zero shifts.
See the [test audit](test-audit.md) for the 0.1.13 consolidation decisions.

`fixtures/phyloem-fit.rds` contains real PhylogeneticEM fits: two traits with no
selected shift, a selected shift, and a univariate model. The data are synthetic
and deterministic. The fixture records the backend version and seed. Regenerate
it deliberately, with all optional dependencies installed:

```sh
Rscript tools/generate-phyloem-fixture.R
make test
```

The generator uses PhylogeneticEM's public fitting API. Backend deprecation
warnings during fitting are visible; they are not suppressed in the generator.
Normal tests read the small fitted objects and exercise the real backend
extraction/reconstruction functions without repeating fitting in every run.
Rphylopars integration tests actually fit numeric matrix and data-frame inputs
and compare their reconstructions.

The 0.1.11 local verification used R 4.4.3, PhylogeneticEM 1.8.1, Rphylopars
0.3.10, and roxygen2 8.1.0. CRAN Deriv 4.3.0 failed to compile on that R version
because `R_ClosureFormals` was unavailable. The isolated verification library
used the archived Deriv 4.2.0 instead; no runtime dependency bound was added to
rkftools. Prefer current R for full optional-backend setup. If reproducing that
older environment, the workaround is only needed until Deriv restores that R
compatibility or declares an appropriate R requirement. Local compiler overrides
were confined to a temporary Makevars file; global R libraries were not changed.

## CI coverage

| Trigger | Coverage |
| --- | --- |
| Relevant push or PR | Linux R 4.1 and release; release includes all optional integrations and coverage in the same job |
| Weekly or manual full run | Also Linux R 4.2, 4.3, devel, macOS release, Windows release; Linux release uses `--as-cran` |
| Manual `compatibility_only=true` | Only the five extended configurations, for a commit whose base checks already passed |
| Version tag | Build one tarball, check that exact tarball with all Suggests, then publish it |
| Weekly or manual benchmark | Repeated timings, correctness checks, heap/RSS measurements, CSV and typed-output artifacts |

Push checks are limited to default branches; PR checks cover incoming changes.
Superseded checks are cancelled. Routine checks use two environment setups,
instead of the previous seven check setups plus a separate coverage setup.
Platform-specific coverage remains on GitHub-hosted runners. Failed check logs
are retained for five days; benchmark artifacts for fourteen days.

## Source layout

Tree indexing and validation are in `R/tree-index.R`. Traversal, transformations,
rooting, MAD, conversion, labels, and collapse have separate modules. Trait
tables, similarity, foreground states, and imputation likewise have focused
files. Public comparative-model entry points are in `R/model-tables.R`, with
backend-specific conversion in `model-l1ou.R` and `model-phylogeneticem.R`.
These internal boundaries do not change the exported API.
