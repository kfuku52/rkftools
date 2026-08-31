# Development

Use a working R installation with matching C/C++ and Fortran toolchains when
installing the optional backends from source. Required runtime packages remain
`ape` and `phytools`; development tools and backend packages are Suggests.

```sh
make setup          # All direct runtime, development, and optional dependencies
make test           # Every regression, including the former smoke tests
make check-full     # Build and check with all Suggests required
make coverage       # Enforce the 80% coverage floor
```

Setup adds missing packages to `.local/R-library`, which is ignored by Git and
excluded from package builds. Existing libraries are not overwritten. Use
`make DEV_LIB=/absolute/path/to/library ...` to choose another library. The
Makefile exports that library for every R process while preserving existing
library search paths. Standalone scripts also accept `RKFTOOLS_DEV_LIBRARY`.

For a smaller environment, run `make setup-minimal` and `make check`. Optional
backend tests explicitly skip when those packages are absent; adapter contract
tests still run without them. `make check-full` and `make check-as-cran` require
all Suggests, matching the strict CI configuration. Missing dependencies are
errors in those checks, not silent skips.

`make document` regenerates the Rd files with roxygen2. `make release-check`
checks DESCRIPTION, NEWS, and the README badge; an explicit tag argument to
`tools/release-check.R` also verifies the tag version. Bump the package version
before pushing, following the repository policy.

## Tests and fixtures

All tests live in `tests/testthat/`. The previous `tests/smoke.R` assertions are
organized by feature and run through both `make test` and `R CMD check`.
`helper-fixtures.R` provides small, independent starting trees. Backend adapters
also have lightweight tests for types, names, missing cells, and zero shifts.

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
