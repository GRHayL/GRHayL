# Task List: Remaining GRHayL M1 Review Findings

- [x] Confirm the canonical `ghl_error_m1_con2prim_failure` order and exact
      diagnostic against the existing test and archived PR source; the shared
      checkout already contained the declaration and fatal mapping when the
      implementation agent's patch ran.
- [x] Build/run the existing M1 inventory in separate default-HDF5 and
      `--disable-hdf5` disposable configurations; both passed.
- [x] Format the 25 C/header paths from the Issue 4 inventory with
      clang-format 23.1.1; the dry-run check passes.
- [x] Add the Radiation M1 ownership/build/test route to `wiki/source-map.md`;
      the repository-local Markdown link audit passes.
- [x] Reconfirm `implementations/GRHayLib/src/Radiation` resolves to
      `../../../GRHayL/Radiation`, and the runner uses Bash-3.2-compatible
      indexed arrays rather than `mapfile` or associative arrays.
- [ ] `git diff --check` remains unavailable because `.git` points to a
      missing worktree metadata directory.
