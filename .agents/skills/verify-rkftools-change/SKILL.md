---
name: verify-rkftools-change
description: Select and run rkftools regression and package checks after source, test, tooling, or documentation changes. Use for change verification, not performance benchmarking or release publication itself.
---

# Verify an rkftools change

Input: the working-tree diff (or an explicit revision range), the intended
behavior, and any environment or execution limits. Work from the repository
root. The output is a concise verification report, not a new report file.

1. Inspect `git status --short` and the relevant diff. Read
   [development](../../../docs/development.md#choose-verification-for-the-change)
   for commands and the area-to-test map. For changed semantics, read the
   affected R implementation, callers, regression, and
   [usage contracts](../../../docs/usage.md#transformation-contracts).
2. Check that R and Rscript run, and use the existing development library.
   If unavailable, first look for an already installed compatible environment;
   do not install another runtime or modify global configuration implicitly.
   Setup instructions and the network boundary are in the development guide.
3. Select a filename regex from the affected contracts and run the existing
   test entrypoint. For example, a branch-table change starts with
   `make test TEST_FILTER='table-roundtrip|branch-tables'` and should execute
   both contexts. Inspect the context list: a successful unrelated selection
   is not evidence for the change. A no-match selection must fail.
4. Apply the guide's delivery checks for the changed area, including unfiltered
   tests after R/test/runner changes. Real backend changes need their installed
   integration tests; mocks alone do not establish compatibility. Do not
   regenerate fitted objects, alter seeds/tolerances, or run benchmarks just
   to obtain a green result. New diagnostic I/O belongs in a temporary directory;
   preserve pre-existing check artifacts by using a temporary source copy.
5. Review `git diff --check` and the final diff for accidental generated files.
   Report the source revision/dirty state, R version, commands, selected
   contexts, failures/warnings/skips, package-check status and coverage when
   applicable. Distinguish executed checks, static inspection, and checks
   blocked by missing dependencies. Link existing logs when useful.

If a command fails, preserve its nonzero result and diagnose the first relevant
failure. Correct an invalid filter rather than accepting an empty run. If the
environment prevents a required check, report the missing prerequisite and
exact remaining command; do not lower the check or call minimal coverage full.
Publication, when requested, follows the existing prepare-github-push skill.
