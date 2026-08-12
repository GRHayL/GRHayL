# Build And CI

Repo ground truth: `README.md`, `configure`, `generate_makefile.sh`,
`Doxyfile`, `.github/run_tests.sh`, `.github/actions/*/action.yml`, and
`.github/workflows/*.yml`.

## Configure Flow

`configure` is the primary build entry point. It probes host OS, required
commands, compiler, math library linking, OpenMP support for unit-test links,
and optional HDF5 support. It then parses `GRHayL/make.code.defn` recursively
through `scripts/parser`, creates the build tree, symlinks public headers into
the configured build include area, and writes the top-level `Makefile`.

Common flags visible in `configure`:

| Flag | Effect |
| --- | --- |
| `-h`, `--help` | Print flag help. |
| `-u`, `--usage` | Print basic usage, or named examples for `ubuntu`, `mac`, `falcon`, `sawtooth`, or `lemhi`. |
| `-l`, `--license` | Print `LICENSE`. |
| `-r`, `--reconfigure` | Overwrite a previous configured build directory. |
| `-s`, `--silent` | Generate less verbose Makefile command output. |
| `--noomp` | Disable OpenMP flags when linking unit tests. |
| `--disable-hdf5` | Disable HDF5 and omit tabulated sources/tests selected by the script. |
| `--prefix=<dir>` | Installation prefix. |
| `--builddir=<dir>` | Build directory; default is `build`. |
| `--buildtype=<type>` | Compiler flag preset. Current help and parser disagree: help advertises `nocflags`, but the parser rejects it; the parser accepts undocumented `plain`, which supplies no preset flags. |
| `--cflags="<flags>"` | Extra compiler flags. |
| `--clibs="<libs>"` | Extra linker flags. |
| `--hdf5dir=<dir>` | HDF5 base directory containing include and lib subdirectories. |
| `--hdf5inc=<dir>` | HDF5 include directory. Must be paired with `--hdf5lib` for custom paths. |
| `--hdf5lib=<dir>` | HDF5 lib directory. Must be paired with `--hdf5inc` for custom paths. |

`configure` writes Makefile targets for `all`, `grhayl`, `tests`, `datagen`,
`clean`, `realclean`, `install`, and `uninstall`. `tests` and `datagen` are
populated from `Unit_Tests/unit_test_*.c` and
`Unit_Tests/data_gen/unit_test_data_*.c`, filtered by the HDF5 setting.

Build-type flags have another live help/parser mismatch. Help advertises
`-fno-finite-math-only` for `production`; the parser currently generates
`-Wall -std=c99 -march=native -O3` without that flag. Treat `nocflags` and the
documented production flag string as broken documentation routes. `plain` and
the parser's emitted production flags describe current behavior; maintainer
intent remains unknown.

Configuration, compilation, installation, consumer linking, and execution are
separate evidence classes:

- successful `configure` proves generated target selection only;
- `make grhayl`, `make tests`, and `make datagen` prove compilation/linking of
  selected targets, not test execution;
- `make install` proves only files copied into the chosen prefix;
- a consumer compile/link against that prefix proves a referenced symbol is
  linkable in that host/mode; and
- only running the consumer or test proves the exercised runtime path.

Do not infer a later class from an earlier one. Workflow YAML and a configured
target likewise prove intent/selection, not a successful historical run.

## Legacy Makefile Generator

`generate_makefile.sh` is a separate, currently **broken** Makefile generator.
It finds every path containing `make.code.defn`, excludes paths matching
`ET/`, scans `SRCS` blocks, and writes a top-level `Makefile` intended to build
`lib/libgrhayl.a`.

This script is not the same path as `configure`. A clean disposable-tree
reproduction currently:

- turns comment words in `GRHayL/Induction/make.code.defn` into source/object
  targets;
- includes `implementations/GRHayLib/src/initialize_and_shutdown.c` because it
  scans manifests outside `GRHayL/`;
- emits `-I./include`, although public headers live under `GRHayL/include/` in
  this checkout;
- has no shebang while using Bash `[[ ... ]]` syntax, so interpreter selection
  depends on caller-shell fallback; and
- has no no-HDF5 source filtering.

Generation itself exits successfully, so exit status alone does not validate
the Makefile. Do not run `make` from this output. Exact repair-versus-retirement
intent is unknown and needs maintainer confirmation.

## `make.code.defn` Inclusion

`make.code.defn` files are the source registry for the configured build.
`configure` starts at `GRHayL/make.code.defn`, follows `subdirs`, and appends
each child `sources`, `headers`, and `install_headers` entry. The generated
Makefile compiles listed C sources under `GRHayL/` into the configured build
directory and links `build/lib/libghl.so` or the host shared-library extension.

During configuration, every `GRHayL/include/*.h` is symlinked into
`<builddir>/include/ghl`. `make install` instead copies the headers parsed from
`GRHayL/include/make.code.defn` into `<prefix>/include/ghl`. Both mechanisms
currently select the same 16 headers, including `ghl_unit_tests.h`, but they
are separate lists and can drift. Installation then copies the versioned
shared library and symlink into `<prefix>/lib`. Installed presence does not by
itself classify a header as production versus test-only API.

`ghl_unit_tests.h` is a concrete installed-surface caveat. Its `static inline`
helpers are caller-compiled, but its non-inline test-helper declarations are
not part of any GRHayL library source manifest. Several are defined only by
helper files under `Unit_Tests/`; no repo definition is visible for the
declared binary read/write family or `ghl_initial_random_data`. Therefore an
installed header compile is not proof those declarations link against
`libghl`. Treat them as test-only or unresolved, not production library API,
until install intent and definitions are reconciled.

`generate_makefile.sh` instead scans all repo `make.code.defn` files found by
`find`, except `ET/` paths, and turns listed sources into static-library object
rules under `build/`.

## HDF5 Contract

Default `configure` builds use HDF5. When HDF5 is enabled, `configure` first
tries `pkg-config hdf5`; custom include/lib paths can be provided with
`--hdf5dir` or `--hdf5inc` plus `--hdf5lib`. It also compiles a small HDF5
probe program.

With `--disable-hdf5`, `configure`:

- adds `-DGHL_DISABLE_HDF5` to `CFLAGS`;
- filters implementation sources with the exact path/name predicate in
  `configure`; despite their paths, it explicitly retains
  `Con2Prim/Tabulated/tabulated_primitive_guess_helpers.c` and sources under
  `Con2Prim/Tabulated/neural_network_guess/`;
- excludes `unit_test_*tabulated*.c`, `unit_test_con2prim_debug.c`, and the
  NRPyLeakage unit tests from the generated unit-test list;
- excludes tabulated data generators from the generated data-generator list.

For Neutrinos-specific HDF5/EOS details, route public table-backed API behavior
through [Neutrinos API and data](gems/neutrinos/api-and-data.md), and fixture or
SLy4 table setup through [Neutrinos tests and fixtures](gems/neutrinos/tests-and-fixtures.md).

Manual or downstream no-HDF5 builds must mirror current script behavior: define
`GHL_DISABLE_HDF5` and reproduce its source-selection predicate. The broader
README wording that all tabulated implementation sources are omitted is not an
exact description of the current generated source list.

GRHayLib is separate implementation-specific build routing. Its Cactus
`configuration.ccl` hard-codes `requires HDF5`; that thorn requirement is not
the same contract as core GRHayL configured `--disable-hdf5` support. Route
Cactus CCL, aggregate header, and thorn source-registry questions through
[GRHayLib Cactus build boundary](implementations/grhaylib/cactus-build-boundary.md).

## GitHub Actions Matrix

Workflows live in `.github/workflows/`:

| Workflow | Compiler | OS matrix | Coverage step status |
| --- | --- | --- | --- |
| `github-actions-Ubuntu-gcc.yml` | `gcc` | `ubuntu-22.04`, `ubuntu-24.04` | 13 existing job groups invoke the shared coverage action; the focused CompOSE job uploads only its Python XML |
| `github-actions-Ubuntu-clang.yml` | `clang` | `ubuntu-22.04`, `ubuntu-24.04` | all 13 jobs invoke coverage action |
| `github-actions-Ubuntu-intel.yml` | `intel` / `icx` | `ubuntu-22.04`, `ubuntu-24.04` | 2 of 13 jobs invoke coverage action |
| `github-actions-MacOS-gcc.yml` | Homebrew GCC | `macos-15`, `macos-26` | all 13 jobs invoke coverage action; local collection body is commented |
| `github-actions-MacOS-clang.yml` | Homebrew LLVM clang | `macos-15`, `macos-26` | no jobs invoke coverage action |

Each workflow ignores pushes and pull requests when **all** changed paths match
its `paths-ignore` list, including `docs/**`, `wiki/**`, Markdown/reStructuredText
patterns, and `implementations/**`. A mixed change with any non-ignored path can
trigger the workflow; path filters apply to `push`/`pull_request`, while the
separately declared schedule remains eligible independently. These semantics
come from the
[GitHub Actions workflow syntax](https://docs.github.com/en/actions/reference/workflows-and-actions/workflow-syntax#onpushpull_requestpull_request_targetpathspaths-ignore),
not merely from local YAML key names. All five workflows use cron
`33 15 1,15 * *`. Their `push` event is restricted to branch `main`; their
`pull_request` event has no branch filter in local YAML. Do not infer project
support beyond the OS/compiler
pairs encoded in these workflow matrices and the usage examples in `configure`.
Because `implementations/**` is path-ignored, repository CI does not validate
GRHayLib implementation-only changes unless another touched path or a scheduled
run causes jobs to execute. Even then, the listed jobs are core
`configure`/unit-test jobs, not a GRHayLib Cactus thorn build.

Common job groups across workflows:

| Job | Test scope |
| --- | --- |
| `ET-Legacy` | `conservs`, `primitives`, `induction_gauge_rhs`, `HLL_flux`, `reconstruction`, `flux_source` |
| `c2p-routines` | `apply_conservative_limits`, `con2prim_multi_method_hybrid`, `enforce_primitive_limits_and_compute_u0`, `compute_conservs_and_Tmunu` |
| `c2p-failure` | `hybrid_failure`, `c2p_nn_guess` |
| `tabulated-eos` | tabulated EOS table read/interpolation |
| `piecewise-polytrope-eos` | piecewise-polytrope EOS |
| `grhayl-core` | core struct/metric/stress-energy suite |
| `flux` | `hybrid_flux`, `tabulated_flux` |
| `reconstruction` | `PLM_reconstruction`, `WENOZ_reconstruction` |
| `neutrinos` | NRPyLeakage optically thin gas, constant density sphere, luminosities; see [Neutrinos tests and fixtures](gems/neutrinos/tests-and-fixtures.md) |
| `con-to-prim-tabulated` | tabulated C2P routines |
| `code-failure` | expected error-code failures |
| `induction-interpolators` | cell/vertex interpolation variants; see [Induction verification workflows](gems/induction/verification-workflows.md) |
| `induction-flux` | vector-potential HLL flux variants; see [Induction verification workflows](gems/induction/verification-workflows.md) |
| `compose-regularized-eos` | 100% Python line/branch coverage, synthetic fixed-profile conversion, and unchanged StellarCollapse C integration |

Composite actions:

- `.github/actions/OS_setup/action.yml` installs compiler/HDF5 dependencies.
- `.github/actions/compile_GRHayL/action.yml` runs `configure`, `make tests
  datagen`, and `make install`.
- `.github/actions/code-coverage/action.yml` selects compiler/OS-specific
  collection steps, several of which contain only comments, then invokes
  `codecov/codecov-action@v5` unconditionally when the composite action itself
  is called. Action invocation is not proof that a usable coverage artifact was
  collected or uploaded.
- `codecov.yml` at repo root configures Codecov report behavior for those
  uploads, including per-gem coverage components and ignored test paths. Its
  header comment requires validating any change with
  `curl -X POST --data-binary @codecov.yml https://codecov.io/validate`.
  The CompOSE flag and component require project and patch coverage of 100%
  with zero threshold; global patch coverage also targets 100%.

## `.github/run_tests.sh`

`.github/run_tests.sh` is a broad repository replay driver, not a complete
workflow matrix and not a fixture generator:

1. Runs `./configure -r`.
2. Runs `make tests`.
3. Exports `LD_LIBRARY_PATH` with `build/lib`.
4. Downloads binary fixtures from the repo-visible `GRHayL/TestData` raw URL
   base.
5. Downloads EOS tables from the repo-visible `stellarcollapse.org/EOS` URLs
   where needed, decompressing `*.bz2` files.
6. Runs the compiled tests under `test/`, including the direct
   `unit_test_c2p_nn_guess` route.
7. Runs `unit_test_code_error` over error-code keys `0` through `85`, expecting
   each invocation to fail at process level.
8. Runs `pyghl append SLy4_3335_rho391_temp163_ye66.h5` before
   `./test/unit_test_con2prim_tabulated SLy4_3335_rho391_temp163_ye66.h5 1`;
   this records only the visible runner/workflow setup command for NN-enabled
   tabulated replay.
9. Continues selected compiled-test runs.
10. Removes root-level `*.bin`, `*.h5`, and `*.bz2` files if execution reaches
    the final cleanup command. Because `set -Eeuxo pipefail` exits on an earlier
    failure, cleanup is not guaranteed.

That cleanup is indiscriminate: an already-present matching file is skipped by
the downloader but still removed by the final glob. Run this driver only in a
disposable checkout without user-owned root-level `.bin`, `.h5`, or `.bz2`
files.

The runner directly invokes 27 of the 29 configured default test binaries. It
does not invoke `unit_test_WENOZ_reconstruction` (workflow matrices do) or
`unit_test_con2prim_debug` (no runner/workflow invocation is visible). The
composite-action YAML configures `tests` and `datagen` compilation, but neither
that action nor the local runner executes data-generator binaries. Tracked YAML
therefore establishes a workflow-configured compile route only. After an
observed successful action or local `make datagen`, those binaries are
`compiled-unrun` until a separate command executes them.

## Coverage Caveats

Coverage is configured through workflow flags and `.github/actions/code-coverage/action.yml`.
Repo evidence shows these caveats:

- Linux GCC uses `gcovr`; Ubuntu image handling differs for `ubuntu22`.
- Linux clang uses `llvm-profdata` and `llvm-cov`, with a comment that expected
  coverage files are still not generated.
- Linux Intel action body is commented, with a note questioning compatibility.
- macOS GCC and clang coverage commands are commented, with notes about needing
  installed tool versions or `llvm-cov gcov`.
- Some workflow coverage steps are commented out, especially macOS clang and
  most Ubuntu Intel jobs.
- The focused Ubuntu GCC CompOSE job bypasses coverage-file discovery: it
  uploads only `compose-coverage.xml` under the `compose` flag, disables
  search, and fails the job on an upload error.
- Workflows ignore docs-only and implementation-only pull-request changes, so
  CI coverage does not prove those paths are exercised.

## Ground Truth References

- GitHub Actions path-filter semantics:
  https://docs.github.com/en/actions/reference/workflows-and-actions/workflow-syntax#onpushpull_requestpull_request_targetpathspaths-ignore
