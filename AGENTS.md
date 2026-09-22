<!-- BEGIN KF AGENT POLICY: source=https://github.com/kfuku52/kf-agent-policy; version=10; sha256=82e3c0eb467582a414d9a6b2feaaaf6f5c8ae330d30f2e3efbf8c303155d0e2e -->
# Common agent policy

Repository-specific instructions override these defaults.

- Follow the user's task scope within higher-priority instructions and execution
  permissions. Complete implementation through affected verification and a result
  report; a plan or investigation ends with its requested deliverable. Continue
  authorized work without repeated approval; identify actual blocking boundaries.
- Inspect the worktree and preserve unrelated changes. Refresh remote information
  when needed; do not merge, rebase, or switch branches merely to inspect it.
- Prefer the default branch when starting work without an established branch.
  Preserve an existing task branch; follow explicit user branch instructions.
  Never create or switch branches solely for a commit, push, release, or PR.
- Change or recommend branch protection only when explicitly asked. Honor explicit
  repository-specific direct-push exceptions; otherwise report a rejected push
  without bypassing protection or inventing a branch or PR.
- Unpublished implementation details may be redesigned; preserve existing public
  APIs, file formats, and saved-data compatibility unless a breaking change is
  authorized. Update affected producers, consumers, tests, examples, and docs.
- Fix verified root causes; do not hide failures with fallbacks or weaker checks.
  Document unavoidable workarounds and their removal conditions.
- Read relevant docs and run the repository's check entrypoint for the change and
  phase. Verify affected behavior; report checks run and omitted. Repeat or broaden
  successful checks only for new changes, failures, or unresolved concerns.
- For library metadata, require demonstrated incompatibility for exact pins or
  upper bounds; keep reproducibility locks separate.
- When editing READMEs, keep them concise with useful visuals inline; put extended
  guides in linked documentation.
- For GitHub push/release work, use `prepare-github-push` in `.agents/skills/`.
  Local-only commits need no version bump; GitHub pushes require one.
- For software performance work, use `benchmark-performance` in `.agents/skills/`.
  Performance claims require comparable measurements and equivalent output.
- For GitHub Actions edits, use `optimize-github-actions` in `.agents/skills/`.
  Preserve required coverage; never run untrusted PR code on self-hosted runners.
<!-- END KF AGENT POLICY -->

# Working on rkftools

- Start with [README](README.md) and [development](docs/development.md).
  The latter owns setup, test selection, and delivery commands. Read
  [usage contracts](docs/usage.md#transformation-contracts) before changing
  scientific behavior or output structure; consult [test audit](docs/test-audit.md)
  before removing or consolidating regression cases.
- This is an R package, not a standalone CLI. Start shared tree changes at
  `R/tree-index.R`, conversions at `R/tree-table.R`, and model adapters at
  `R/model-tables.R`. Use the development guide's test-selection table for
  related modules. Public exports are in `NAMESPACE`; help originates in R
  roxygen comments.
- Run commands from the repository root with R and Rscript on PATH:
  `make setup-minimal`, `make test`, `make check`. Focus an iteration with
  `make test TEST_FILTER='table-roundtrip|branch-tables'`. Before delivery use
  the guide's checks, including `make release-check` for a push. There is no
  separate configured lint or static type checker; `R CMD check` is the
  existing code/help/example check. Do not claim it proves static type safety.
- Use [verify-rkftools-change](.agents/skills/verify-rkftools-change/SKILL.md)
  to choose and report checks after changes. Setup may download packages;
  benchmark runs and fitted-fixture regeneration are deliberate workflows,
  not prerequisites for every edit.
- Preserve numerical-label/edge-row correspondence, zero/missing branch
  lengths, trait names and missing cells, species parsing, and MAD mode/list
  contracts. Do not change seeds, fitting models, thresholds, tolerances, or
  legacy arguments merely to make verification pass. See usage/help and the
  affected regression for the precise contract.
- Keep `.local/`, `rkftools.Rcheck/`, package tarballs, `benchmark/`, IDE
  checkpoints, and personal R settings out of commits. Do not hand-edit
  generated `man/*.Rd` or `NAMESPACE`: use `make document` when needed.
  Regenerate figures or `tests/testthat/fixtures/phyloem-fit.rds` only for an
  intentional corresponding change, then review them.
- Finish with the intended diff, checks and skips, environment limitations,
  and (when requested) pushed commit/version. Preserve unrelated files and
  distinguish commands actually run from static review or unavailable checks.
