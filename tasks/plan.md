# Implementation Plan: Close Remaining GRHayL M1 Review Findings

## Contract

Resolve the findings that remain in the current checkout: declare and dispatch
the M1/flux-source errors expected by the implementation and existing test,
format the affected changed C/header code to the repository style, and add the
missing Radiation M1 ownership route to `wiki/source-map.md`.

The GRHayLib Radiation symlink and Bash 3.2 runner compatibility were present
in the reviewed tree, so this plan does not change them.

## Current State and Prerequisite

The initial planning read found the M1/flux-source declarations and fatal
dispatch missing, while the existing `unit_test_m1_error_handling.c` referenced
all 24 statuses. During implementation, the shared checkout was concurrently
updated; the current source now contains all 24 identifiers and matching
`GHL_CASE_ERROR` entries. The previously missing `ghl_error_m1_con2prim_failure`
is in the order shown by the existing test and archived PR source, with its
exact test diagnostic. The assigned agent's patch did not apply because that
shared update had already occurred. The checkout's Git worktree reference
remains unavailable, so a base-to-head diff cannot be audited here.

## Scoped Build

Compile and run only the current M1 runner targets and their GRHayL dependencies
in two separate disposable configurations: default HDF5 and `--disable-hdf5`.
The runner currently selects `unit_test_m1_closure_fallback`,
`unit_test_m1_diffusion_flux`, `unit_test_m1_error_handling`,
`unit_test_m1_fd_jacobian`, `unit_test_m1_neutrino_rusanov_flux`,
`unit_test_m1_neutrino_seeded_invariants`,
`unit_test_m1_neutrino_source_update`, `unit_test_m1_rate_provider`,
`unit_test_m1_thcm1_blended_rusanov`, and `unit_test_rusanov_flux`.

Do not build the whole toolkit or a Cactus/Einstein Toolkit configuration for
these findings.

## Ordered Tasks

### Task 1: Restore the public error codes and fatal mappings

**Description:** Recover the authoritative enum mapping, restore the 24
M1/flux-source values expected by the current source and error test, and add an
identifiable `GHL_CASE_ERROR` entry for each in `abort_if_error.c`.

**Acceptance criteria:**

- Preserve the authoritative numeric values and enum order; `ghl_success` is
  the only status that returns without calling `ghl_Error`.
- Every M1/flux-source entry emits the exact diagnostic already specified in
  `m1_error_cases`; do not use a generic fallback for these codes.
- The existing error test verifies each exit status and diagnostic, and the
  dispatcher compiles without an unhandled-enumerator warning.

**Files:** `GRHayL/include/ghl.h`,
`GRHayL/GRHayL_Core/abort_if_error.c`; reuse
`Unit_Tests/unit_test_m1_error_handling.c` without adding cases unless the
authoritative mapping proves the existing inventory incomplete.

**Dependency:** The enum mapping prerequisite above.

### Task 2: Format the affected changed C and header paths

**Description:** Use the checked-in `.clang-format` on the C/header paths in the
Issue 4 inventory that are part of the PR diff. Limit the existing
`NRPyLeakage_nucleon_blocking.h` cleanup to its changed hunk. Review the diff to
confirm formatting-only edits preserve behavior.

**Acceptance criteria:**

- `clang-format --dry-run --Werror` reports no diagnostics on changed C/header
  paths under the same formatter version used for the formatting pass.
- `git diff --check` passes after a valid PR diff is available.
- Existing M1 runner targets pass in both scoped HDF5 configurations after the
  formatting pass.

**Files:** The affected C/header paths enumerated in the Issue 4 review,
including the currently reported diagnostics in
`GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_implicit_solve.c`,
`GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_source_update.c`,
`Unit_Tests/m1_test_utils.h`,
`Unit_Tests/unit_test_m1_diffusion_flux.c`,
`Unit_Tests/unit_test_m1_neutrino_seeded_invariants.c`, and
`Unit_Tests/unit_test_m1_neutrino_source_update.c`.

**Dependency:** A valid PR diff is needed to distinguish changed lines from
pre-existing formatting in modified files.

### Task 3: Add Radiation M1 to the source map

**Description:** Add one concise `GRHayL/Radiation/` row to
`wiki/source-map.md`, linking the M1 owner hub and the authoritative public
headers, Doxygen contract, scoped runner, and fixture/test route.

**Acceptance criteria:**

- The row routes source ownership and build/test impact without duplicating the
  M1 owner pages.
- Every repo-relative link in the new row resolves within the checkout.
- The repository-local broken Markdown link check documented in
  `wiki/lint/CHECKS.md` reports no missing or escaping target.

**Files:** `wiki/source-map.md`.

## Final Checkpoint

- Run `bash Unit_Tests/run_m1_tests.sh --build` in separate default-HDF5 and
  `--disable-hdf5` disposable configurations.
- Confirm `unit_test_m1_error_handling` observes all expected fatal exit codes
  and exact diagnostics; confirm no enum warning from `abort_if_error.c`.
- Run the formatter check on the 25 C/header paths enumerated by Issue 4 and
  the repository-local Markdown link audit; both pass. `git diff --check` is
  blocked because the checkout's `.git` worktree reference points to a missing
  directory.
- Reconfirm the existing `src/Radiation -> ../../../GRHayL/Radiation` link and
  that the runner uses Bash 3.2-compatible indexed arrays; make no changes to
  these already-resolved findings unless the final diff shows a regression.

## Risks and Open Dependency

The current checkout lacks both a usable Git worktree reference and the
reported enum declarations. The exact public numeric mapping must come from an
authoritative PR revision or artifact before implementation. The formatter
version is also not pinned by the reviewed files; keep formatter verification
and formatting on the same selected version, and do not treat one environment's
diagnostic count as a project-wide threshold.
