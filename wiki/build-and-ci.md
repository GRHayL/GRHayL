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
| `--enable-m1-debug` | Define `GRHAYL_M1_DEBUG`, compiling the additional expensive M1 runtime-parameter validation. Off by default, and orthogonal to `--buildtype`. |
| `--prefix=<dir>` | Installation prefix. |
| `--builddir=<dir>` | Build directory; default is `build`. |
| `--buildtype=<type>` | Compiler flag preset. Current help and parser disagree: help advertises `nocflags`, but the parser rejects it; the parser accepts undocumented `plain`, which supplies no preset flags. |
| `--cflags="<flags>"` | Extra safe compiler flags. Unsafe floating-point modes are rejected. |
| `--clibs="<libs>"` | Extra linker flags. |
| `--hdf5dir=<dir>` | HDF5 base directory containing include and lib subdirectories. |
| `--hdf5inc=<dir>` | HDF5 include directory. Must be paired with `--hdf5lib` for custom paths. |
| `--hdf5lib=<dir>` | HDF5 lib directory. Must be paired with `--hdf5inc` for custom paths. |

`configure` writes Makefile targets for `all`, `grhayl`, `tests`, `datagen`,
`clean`, `realclean`, `install`, and `uninstall`. `tests` and `datagen` are
populated from `Unit_Tests/unit_test_*.c` and
`Unit_Tests/data_gen/unit_test_data_*.c`, filtered by the HDF5 setting.

The remaining build-type mismatch concerns the no-flags name: help advertises
`nocflags`, while the parser rejects it and accepts undocumented `plain`.
Production emits its documented
`-Wall -std=c99 -march=native -fno-finite-math-only -O3` flags so runtime
finite-value checks retain their required semantics.
All build types reject explicit unsafe floating-point flag tokens listed by
`./configure --help` when supplied through `CC` or `--cflags`, including
fast-math and finite-only umbrella modes and the individual transformations on
which those modes rely. A compile probe also rejects final flag sets that define
fast-math or positive finite-only macros, and an execution probe requires
gradual underflow. For Intel LLVM, `configure` appends
`-fp-model=precise -no-ftz` after user flags because optimized ICX builds
otherwise use fast arithmetic and flush subnormals. Build routes that bypass
`configure`, including Cactus and downstream executable links, must enforce the
same semantics; `-no-ftz` must reach the program containing `main`.

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
`<builddir>/include/ghl`. `make install` instead copies the public headers
parsed from `GRHayL/include/make.code.defn` into `<prefix>/include/ghl`; it
excludes `ghl_unit_tests.h` and the internal storage-definition companion
`ghl_eos_functions_declaration.h`. Installation then copies the versioned
shared library and symlink into `<prefix>/lib`. Installed presence does not by
itself classify a header as production versus test-only API.
`make install` also copies the FDINT notice, stored as `THIRD_PARTY_NOTICES`, into
`<prefix>/share/doc/grhayl` so binary installations retain the FDINT BSD-3
notice.

`ghl_unit_tests.h` remains available to source-tree tests only; its helpers are
not production library API.

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
  `Con2Prim/Tabulated/tabulated_primitive_guess_helpers.c`, sources under
  `Con2Prim/Tabulated/neural_network_guess/`, and the direct tabulated HLLE
  flux implementations. The disabled direct-tabulated-solver stubs also remain
  because their source path does not match the exclusion tokens;
- excludes `unit_test_*tabulated*.c`, `unit_test_con2prim_debug.c`, and the
  table-backed NRPyLeakage tests `constant_density_sphere`, `luminosities`, and
  `optically_thin_gas` from the generated unit-test list. The NRPyLeakage
  `physics` and `classifier_fallback` tests remain available;
- excludes tabulated data generators from the generated data-generator list.

For Neutrinos-specific HDF5/EOS details, route public table-backed API behavior
through [Neutrinos API and data](gems/neutrinos/api-and-data.md), and fixture or
SLy4 table setup through [Neutrinos tests and fixtures](gems/neutrinos/tests-and-fixtures.md).

Manual or downstream no-HDF5 builds must mirror current script behavior: define
`GHL_DISABLE_HDF5` and reproduce its source-selection predicate. The README
lists the retained Con2Prim helpers and the exact exclusion patterns. Loader
entry points become disabled-feature stubs; pure NN inference from an
independently valid in-memory model does not inherently require HDF5. The five
public direct tabulated solver entry points are non-mutating disabled-feature
stubs, while the six public direct tabulated HLLE flux variants retain their
real implementations.

GRHayLib is separate implementation-specific build routing. Its Cactus
`configuration.ccl` hard-codes `requires HDF5`; that thorn requirement is not
the same contract as core GRHayL configured `--disable-hdf5` support. Route
Cactus CCL, aggregate header, and thorn source-registry questions through
[GRHayLib Cactus build boundary](implementations/grhaylib/cactus-build-boundary.md).

## GitHub Actions Matrix

Workflows live in `.github/workflows/`:

| Workflow | Compiler | OS matrix | Coverage step status |
| --- | --- | --- | --- |
| `github-actions-Ubuntu-gcc.yml` | `gcc` | `ubuntu-22.04`, `ubuntu-24.04` | the Radiation M1 job and 13 other job groups invoke the shared coverage action; the focused CompOSE job uploads only its Python XML |
| `github-actions-Ubuntu-clang.yml` | `clang` | `ubuntu-22.04`, `ubuntu-24.04` | all 13 jobs invoke coverage action |
| `github-actions-Ubuntu-intel.yml` | `intel` / `icx` | `ubuntu-22.04`, `ubuntu-24.04` | 2 of 13 jobs invoke coverage action |
| `github-actions-MacOS-gcc.yml` | Homebrew GCC | `macos-15`, `macos-26` | all 13 jobs invoke coverage action; local collection body is commented |
| `github-actions-MacOS-clang.yml` | Homebrew LLVM clang | `macos-15`, `macos-26` | no jobs invoke coverage action |

The Ubuntu-Clang `c2p-failure` matrix configures its Ubuntu 24.04
`c2p_nn_guess` variant without HDF5. That existing job variant exercises and
uploads coverage for the disabled-feature path; the other variants retain
HDF5-enabled builds.

Each workflow ignores pushes and pull requests when **all** changed paths match
its `paths-ignore` list, including `docs/**`, `wiki/**`, and
Markdown/reStructuredText patterns. The Ubuntu-GCC workflow no longer ignores
`implementations/**`; the other four compiler workflows still do. A mixed
change with any non-ignored path can trigger the workflow; path filters apply
to `push`/`pull_request`, while the
separately declared schedule remains eligible independently. These semantics
come from the
[GitHub Actions workflow syntax](https://docs.github.com/en/actions/reference/workflows-and-actions/workflow-syntax#onpushpull_requestpull_request_targetpathspaths-ignore),
not merely from local YAML key names. Every compiler workflow uses cron
`33 15 1,15 * *`. Their `push` event is restricted to branch `main`; their
`pull_request` event has no branch filter in local YAML. Do not infer project
support beyond the OS/compiler
pairs encoded in these workflow matrices and the usage examples in `configure`.
An implementation-only change therefore triggers Ubuntu-GCC, but its listed
test jobs are core `configure`/unit-test jobs, not a GRHayLib Cactus thorn build
or direct GRHayLib validation. The other compiler workflows remain skipped for
such a change.

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
| `neutrinos` | Table-free NRPyLeakage physics checks (blocking, Fermi edges, detailed balance, invalid inputs, asymmetric optical-depth stencil), optically thin gas, constant density sphere, and luminosities; see [Neutrinos tests and fixtures](gems/neutrinos/tests-and-fixtures.md) |
| `con-to-prim-tabulated` | tabulated C2P routines |
| `code-failure` | expected error-code failures |
| `induction-interpolators` | cell/vertex interpolation variants; see [Induction verification workflows](gems/induction/verification-workflows.md) |
| `induction-flux` | vector-potential HLL flux variants; see [Induction verification workflows](gems/induction/verification-workflows.md) |
| `compose-regularized-eos` | 100% Python line/branch coverage, synthetic fixed-profile conversion, and unchanged StellarCollapse C integration |

Flux fixture downloads use the immutable TestData reference recorded in
`.github/et-legacy-testdata-ref`; only the `tabulated_flux` matrix leg downloads
the LS220 EOS table.

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
- `codecov.yml` defers Codecov notifications until an explicit final trigger.
  The `codecov-finalize` job in the Ubuntu GCC workflow waits for every local
  job and for the macOS GCC, Ubuntu Clang, and Ubuntu Intel coverage workflows
  for the same source head recorded by the current workflow run, then sends the
  single final Codecov notification.
  macOS Clang is excluded because its coverage collection steps are disabled.

## `.github/run_tests.sh`

`.github/run_tests.sh` is a broad repository replay driver, not a complete
workflow matrix and not a fixture generator:

1. Runs `./configure -r`.
2. Runs `make tests datagen` (data generators are compiled, not executed).
3. Exports `LD_LIBRARY_PATH` with `build/lib`.
4. Downloads binary fixtures from the repo-visible `GRHayL/TestData` raw URL
   base.
5. Downloads EOS tables from the repo-visible `stellarcollapse.org/EOS` URLs
   where needed, decompressing `*.bz2` files.
6. Runs the compiled tests under `test/`, including the direct
   `unit_test_c2p_nn_guess` route.
7. Runs `unit_test_code_error` over error-code keys `0` through `88`, expecting
   each invocation to fail at process level.
8. Runs `pyghl append SLy4_3335_rho391_temp163_ye66.h5` before
   `./test/unit_test_con2prim_tabulated SLy4_3335_rho391_temp163_ye66.h5 1`;
   this records only the visible runner/workflow setup command for NN-enabled
   tabulated replay.
9. Continues selected compiled-test runs.
10. An `EXIT` trap removes only paths created by that run, including downloaded
    or decompressed files and its private expected-error work directory;
    preexisting paths are preserved, including on early failure.

The broad runner does not invoke the scoped Radiation M1 suite,
`unit_test_WENOZ_reconstruction`, `unit_test_tabulated_eos_compose`, or
`unit_test_con2prim_debug`. Radiation M1 runs through its dedicated action,
WENOZ through the reconstruction workflow matrix, and CompOSE through its
focused workflow; no runner/workflow invocation is visible for the debug
binary. The composite-action YAML configures `tests` and `datagen` compilation,
but neither that action nor the local runner executes data-generator binaries. Tracked YAML
therefore establishes a workflow-configured compile route only. After an
observed successful action or local `make datagen`, those binaries are
`compiled-unrun` until a separate command executes them.

## Coverage Caveats

Coverage is configured through workflow flags and `.github/actions/code-coverage/action.yml`.
Repo evidence shows these caveats:

- Linux GCC uses `gcovr`; Ubuntu image handling differs for `ubuntu22`.
- The Ubuntu GCC Radiation M1 jobs upload a gcovr Cobertura report filtered to
  `GRHayL/Radiation/`, with automatic file search disabled for that upload. The
  `radiation_m1` Codecov component has a 100% project coverage target.
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
- All workflows ignore docs-only pull-request changes. Only Ubuntu-GCC runs its
  normal GRHayL jobs for implementation-only changes; it does not inspect
  GRHayLib symlinks or build/test the Cactus thorn.

## Ground Truth References

- GitHub Actions path-filter semantics:
  https://docs.github.com/en/actions/reference/workflows-and-actions/workflow-syntax#onpushpull_requestpull_request_targetpathspaths-ignore
