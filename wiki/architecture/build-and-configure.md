# Build And Configure

> Supported `./configure` build system, recursive `make.code.defn` parsing, `--disable-hdf5` filtering, floating-point checks, style settings, and the broken legacy `generate_makefile.sh` script. · Status: contested
> Up: [Architecture](index.md)

## Summary

`configure` is the supported GRHayL build generator: it probes compiler, math, OpenMP, and HDF5 toolchains, parses recursive `GRHayL/**/make.code.defn` manifests through `scripts/parser`, symlinks public headers into `<builddir>/include/ghl`, and emits a top-level `Makefile` with `all`, `grhayl`, `tests`, `datagen`, `clean`, `realclean`, `install`, and `uninstall` targets. Repository formatting rules live in `.clang-format` and `.editorconfig`. Legacy `generate_makefile.sh` is currently broken and must not be built with `make`.

## Detail

### `configure` Options And Floating-Point Guards

`configure` probes host OS, commands, C compilers, `-lm`, OpenMP (for test and data-generator compilation/linking), and optional HDF5 support, then populates `tests` and `datagen` targets from `Unit_Tests/unit_test_*.c` and `Unit_Tests/data_gen/unit_test_data_*.c` filtered by HDF5 mode:

| Flag | Effect |
| --- | --- |
| `-h`, `--help`, `-u`, `--usage`, `-l`, `--license` | Print flag help, basic/host usage (`ubuntu`, `mac`, `falcon`, `sawtooth`, `lemhi`), or `LICENSE`. (Ev: `configure::--help`) |
| `-r`, `--reconfigure`, `-s`, `--silent` | Allow reconfiguration without restoring prior flags or cleaning old objects; reduce command echo. (Ev: `configure::--reconfigure`) |
| `--noomp`, `--disable-hdf5` | Disable OpenMP for tests/generators; disable HDF5 and filter tabulated sources/tests. (Ev: `configure::--disable-hdf5`) |
| `--prefix=<dir>`, `--builddir=<dir>` | Set install prefix and build directory (default `build`). (Ev: `configure::--builddir`) |
| `--buildtype=<type>`, `--cflags="<flags>"`, `--clibs="<libs>"` | Select compiler flag preset and append extra safe compile/link flags. (Ev: `configure::--buildtype`) |
| `--hdf5dir=<dir>`, `--hdf5inc=<dir>`, `--hdf5lib=<dir>` | Set custom HDF5 root (`--hdf5dir` takes precedence) or paired include/lib paths. (Ev: `configure::--hdf5dir`) |

Claim status: contested; contradiction: CONTR-0001. See [CONTR-0001](../contradictions.md#contr-0001).

The `production` build type emits `-Wall -std=c99 -march=native -fno-finite-math-only -O3` so runtime finite-value checks keep IEEE semantics. Every build type rejects explicit unsafe floating-point flags listed by `./configure --help` (fast-math, finite-math-only, and component flags) supplied via `CC` or `--cflags`, runs a compile probe rejecting fast-math/finite-only predefined macros, and runs an execution probe requiring gradual underflow. For Intel LLVM (`icx`), `configure` appends `-fp-model=precise -no-ftz` after user flags so optimized builds preserve IEEE arithmetic and subnormals; build routes bypassing `configure` (including Cactus and downstream executables linking `main`) must enforce equivalent semantics.

Evidence (descriptive behavior): `configure::--buildtype`; corroboration: `configure::-fno-finite-math-only`.

### Recursive `make.code.defn` Parsing And Installation

`configure` starts at `GRHayL/make.code.defn`, parses `subdirs`, `sources`, `headers`, and `install_headers` via `scripts/parser` and `scripts/match_all.awk`, and writes `Makefile` rules; `make` (`make grhayl`) then compiles the listed C files into `<builddir>` and links `<builddir>/lib/libghl.so` (or host shared-library extension).

- During `configure`, every `GRHayL/include/*.h` is symlinked into `<builddir>/include/ghl`.
- During `make install`, only headers listed in `GRHayL/include/make.code.defn` (`install_headers`) are copied into `<prefix>/include/ghl`, excluding `ghl_unit_tests.h` and internal storage companion `ghl_eos_functions_declaration.h`.
- `make install` also copies the versioned shared library and symlink into `<prefix>/lib` and installs `THIRD_PARTY_NOTICES` (FDINT BSD-3 notice) into `<prefix>/share/doc/grhayl`.

Evidence (descriptive behavior): `GRHayL/include/make.code.defn::install_headers`; corroboration: `scripts/parser` and `configure::THIRD_PARTY_NOTICES`.

### HDF5 Detection And `--disable-hdf5` Filtering

Default builds enable HDF5: without custom paths, `configure` tries `pkg-config hdf5` when available, then compiles and links (without executing) a small HDF5 probe. With `--disable-hdf5`, `configure`:

- adds `-DGHL_DISABLE_HDF5` to `CFLAGS`;
- filters implementation sources with its path/name predicate while explicitly retaining `Con2Prim/Tabulated/tabulated_primitive_guess_helpers.c`, `Con2Prim/Tabulated/neural_network_guess/`, direct tabulated HLLE flux sources, and the non-mutating disabled-feature direct tabulated Con2Prim stubs;
- excludes `unit_test_*tabulated*.c`, `unit_test_con2prim_debug.c`, and table-backed NRPyLeakage tests (`constant_density_sphere`, `luminosities`, `optically_thin_gas`) while keeping `physics` and `classifier_fallback`;
- excludes tabulated data generators from `DGEXES`.

Manual or downstream no-HDF5 builds must define `GHL_DISABLE_HDF5` and mirror this exact source-selection predicate (`README.md` documents the retained helpers and exclusion tokens). Direct tabulated HLLE flux routines remain linkable in no-HDF5 mode but require caller-supplied EOS callbacks because built-in table dispatch is disabled; see [HLLE Flux Variant Matrix](../gems/flux-source/hlle-flux-variant-matrix.md). By contrast, [GRHayLib Cactus Build](../implementations/grhaylib-cactus-build.md) hard-codes `requires HDF5`.

Evidence (descriptive behavior): `configure::GHL_DISABLE_HDF5`; corroboration: `README.md::GHL_DISABLE_HDF5`.

### Legacy `generate_makefile.sh`, Style Configs, And Read-Only Checks

Claim status: contested; contradiction: CONTR-0002. See [CONTR-0002](../contradictions.md#contr-0002).

`generate_makefile.sh` scans every repository `make.code.defn` outside `ET/`, pulls `implementations/GRHayLib/src/initialize_and_shutdown.c` into library sources, emits checkout-invalid `-I./include` (public headers live under `GRHayL/include/`), lacks a shebang despite Bash `[[ ... ]]` syntax, and performs no `--disable-hdf5` filtering. Inspect it only inside a disposable tree and never invoke `make` on its output.

Repository code-formatting and editor settings live in `.clang-format` and `.editorconfig` (plus `.gitattributes` and `.gitignore`). Side-effect-free build-script checks safe in the shared tree:

```bash
./configure --help
./scripts/parser awk GRHayL/make.code.defn subdirs
./scripts/parser awk GRHayL/include/make.code.defn install_headers
for s in configure scripts/parser; do sh -n "$s" || exit; done
for s in generate_makefile.sh .github/run_tests.sh; do bash -n "$s" || exit; done
```

Evidence (descriptive behavior): `generate_makefile.sh::libgrhayl.a`; corroboration: `.clang-format` and `.editorconfig`.

## Sources

- [configure](../../configure) - options, `--buildtype`, floating-point probes, HDF5 filtering, and `Makefile` generation
- [generate_makefile.sh](../../generate_makefile.sh) - broken legacy static `Makefile` generator
- [scripts/parser](../../scripts/parser) - manifest parser wrapper invoking `scripts/match_all.awk`
- [GRHayL/make.code.defn](../../GRHayL/make.code.defn) - top-level library `subdirs` manifest
- [GRHayL/include/make.code.defn](../../GRHayL/include/make.code.defn) - `install_headers` manifest
- [README.md](../../README.md) - build, install, and `--disable-hdf5` documentation
- [THIRD_PARTY_NOTICES](../../THIRD_PARTY_NOTICES) - FDINT BSD-3 notice installed under `<prefix>/share/doc/grhayl`
- [.clang-format](../../.clang-format) - C/C++ formatting configuration
- [.editorconfig](../../.editorconfig) - editor whitespace and indentation settings

## See Also

- Depends on: [Public API Map](public-api-map.md)
- Depends on: [Generated Boundaries](generated-boundaries.md)
- Contrasts with: [GRHayLib Cactus Build](../implementations/grhaylib-cactus-build.md)
- See also: [CI Workflows](ci-workflows.md)
- See also: [Contradictions](../contradictions.md#contr-0001)
