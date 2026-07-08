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
| `--buildtype=<type>` | Compiler flag preset. Supported values in help are `nocflags`, `debug`, `debug-opt`, `opt`, and `production`; implementation also accepts `plain` while help names `nocflags`. |
| `--cflags="<flags>"` | Extra compiler flags. |
| `--clibs="<libs>"` | Extra linker flags. |
| `--hdf5dir=<dir>` | HDF5 base directory containing include and lib subdirectories. |
| `--hdf5inc=<dir>` | HDF5 include directory. Must be paired with `--hdf5lib` for custom paths. |
| `--hdf5lib=<dir>` | HDF5 lib directory. Must be paired with `--hdf5inc` for custom paths. |

`configure` writes Makefile targets for `all`, `grhayl`, `tests`, `datagen`,
`clean`, `realclean`, `install`, and `uninstall`. `tests` and `datagen` are
populated from `Unit_Tests/unit_test_*.c` and
`Unit_Tests/data_gen/unit_test_data_*.c`, filtered by the HDF5 setting.

## Legacy Makefile Generator

`generate_makefile.sh` is a separate Makefile generator. It finds
`make.code.defn` files, excludes paths matching `ET/`, scans `SRCS` blocks, and
writes a top-level `Makefile` that builds `lib/libgrhayl.a`.

This script is not the same path as `configure`: it creates a static archive,
expects HDF5 include/lib locations in the generated Makefile, and warns the
user to add HDF5 paths when no HDF5 directory argument is given. It does not
implement the no-HDF5 source filtering done by `configure`.

## `make.code.defn` Inclusion

`make.code.defn` files are the source registry for the configured build.
`configure` starts at `GRHayL/make.code.defn`, follows `subdirs`, and appends
each child `sources`, `headers`, and `install_headers` entry. The generated
Makefile compiles listed C sources under `GRHayL/` into the configured build
directory and links `build/lib/libghl.so` or the host shared-library extension.

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
- removes tabulated implementation sources by path/name pattern;
- excludes `unit_test_*tabulated*.c`, `unit_test_con2prim_debug.c`, and the
  NRPyLeakage unit tests from the generated unit-test list;
- excludes tabulated data generators from the generated data-generator list.

For Neutrinos-specific HDF5/EOS details, route public table-backed API behavior
through [Neutrinos API and data](gems/neutrinos/api-and-data.md), and fixture or
SLy4 table setup through [Neutrinos tests and fixtures](gems/neutrinos/tests-and-fixtures.md).

Manual or downstream no-HDF5 builds must mirror both parts visible in repo
docs and script behavior: define `GHL_DISABLE_HDF5` and exclude HDF5/tabulated
sources that require HDF5 runtime paths.

## GitHub Actions Matrix

Workflows live in `.github/workflows/`:

| Workflow | Compiler | OS matrix | Coverage step status |
| --- | --- | --- | --- |
| `github-actions-Ubuntu-gcc.yml` | `gcc` | `ubuntu-22.04`, `ubuntu-24.04` | enabled in jobs |
| `github-actions-Ubuntu-clang.yml` | `clang` | `ubuntu-22.04`, `ubuntu-24.04` | enabled in jobs |
| `github-actions-Ubuntu-intel.yml` | `intel` / `icx` | `ubuntu-22.04`, `ubuntu-24.04` | mixed; most coverage steps are commented, but `c2p-failure` and `reconstruction` are enabled |
| `github-actions-MacOS-gcc.yml` | Homebrew GCC | `macos-15`, `macos-26` | enabled in jobs |
| `github-actions-MacOS-clang.yml` | Homebrew LLVM clang | `macos-15`, `macos-26` | commented out in jobs |

Each workflow ignores pushes and pull requests that touch only paths listed in
`paths-ignore`, including `docs/**`, `*.md`, and `implementations/**`. Scheduled
runs are also configured. Do not infer project support beyond the OS/compiler
pairs encoded in these workflow matrices and the usage examples in `configure`.

Common job groups across workflows:

| Job | Test scope |
| --- | --- |
| `ET-Legacy` | `conservs`, `primitives`, `induction_gauge_rhs`, `HLL_flux`, `reconstruction`, `flux_source` |
| `c2p-routines` | `apply_conservative_limits`, `con2prim_multi_method_hybrid`, `enforce_primitive_limits_and_compute_u0`, `compute_conservs_and_Tmunu` |
| `c2p-failure` | `hybrid_failure` |
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

Composite actions:

- `.github/actions/OS_setup/action.yml` installs compiler/HDF5 dependencies.
- `.github/actions/compile_GRHayL/action.yml` runs `configure`, `make tests
  datagen`, and `make install`.
- `.github/actions/code-coverage/action.yml` runs coverage tooling and then
  `codecov/codecov-action@v5`.

## `.github/run_tests.sh`

`.github/run_tests.sh` is a full local-style CI driver:

1. Runs `./configure -r`.
2. Runs `make tests`.
3. Exports `LD_LIBRARY_PATH` with `build/lib`.
4. Downloads binary fixtures from the repo-visible `GRHayL/TestData` raw URL
   base.
5. Downloads EOS tables from the repo-visible `stellarcollapse.org/EOS` URLs
   where needed, decompressing `*.bz2` files.
6. Runs the compiled tests under `test/`.
7. Runs `unit_test_code_error` over error-code keys `0` through `82`, expecting
   each invocation to fail at process level.
8. Removes downloaded `*.bin`, `*.h5`, and `*.bz2` files from the repo root.

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
- Workflows ignore docs-only and implementation-only pull-request changes, so
  CI coverage does not prove those paths are exercised.
