# Flux Source Generated NRPy Boundary

## Generation

[generate_flux_source.sh](../../../GRHayL/Flux_Source/generate_flux_source.sh)
clones the pinned NRPy 2 commit into a temporary checkout and runs the single
[generate_flux_source.py](../../../GRHayL/Flux_Source/generate_flux_source.py)
driver. Use a Python environment with NRPy 2's runtime dependencies available.
From the repository root, generate into a new or empty directory:

```sh
GRHayL/Flux_Source/generate_flux_source.sh /tmp/grhayl-flux-stage
```

The driver rejects output directories inside `GRHayL/Flux_Source` and nonempty
destinations. It checks its C output names against the root and variant
`make.code.defn` files. Review staged C and compile it before copying it into
the checked-in source tree.

The Python driver uses NRPy 2's ADM metric conversion, GRMHD stress-energy
tensor, magnetic fast-wave estimate, GRHD characteristic roots, and C code
generator. GRHayL's magnetic fields are already scaled by `1/sqrt(4 pi)` as
required by the generated equations. The driver retains the checked C
interfaces: EOS callback errors propagate, legacy wrappers abort on errors,
and HLLE wave speeds are checked before any output is written.

## Checked-In C Kernels

The [root build list](../../../GRHayL/Flux_Source/make.code.defn) compiles the
source-term and three characteristic-speed kernels. The variant build lists
compile three directions each of hybrid, hybrid entropy, tabulated, and
tabulated entropy HLLE fluxes:

| Generated C | Python path | Build list |
| --- | --- | --- |
| `ghl_calculate_source_terms.c` | `generate_source()` | [root](../../../GRHayL/Flux_Source/make.code.defn) |
| `ghl_calculate_characteristic_speed_dirn*.c` | `generate_speeds()` | [root](../../../GRHayL/Flux_Source/make.code.defn) |
| `hybrid/ghl_calculate_HLLE_fluxes_dirn*_hybrid.c` | `generate_fluxes()` | [hybrid](../../../GRHayL/Flux_Source/hybrid/make.code.defn) |
| `hybrid_entropy/ghl_calculate_HLLE_fluxes_dirn*_hybrid_entropy.c` | `generate_fluxes()` | [hybrid entropy](../../../GRHayL/Flux_Source/hybrid_entropy/make.code.defn) |
| `tabulated/ghl_calculate_HLLE_fluxes_dirn*_tabulated.c` | `generate_fluxes()` | [tabulated](../../../GRHayL/Flux_Source/tabulated/make.code.defn) |
| `tabulated_entropy/ghl_calculate_HLLE_fluxes_dirn*_tabulated_entropy.c` | `generate_fluxes()` | [tabulated entropy](../../../GRHayL/Flux_Source/tabulated_entropy/make.code.defn) |

`configure --disable-hdf5` still compiles the tabulated HLLE kernels; it
excludes the table-dependent tests and data generators. The Ubuntu-Clang
`c2p-failure` no-HDF5 job link-checks their public symbols.

## Review Changes

Review the Python equations and affected C together when changing
characteristic-speed roots, HLLE output fields, source-term derivatives,
extrinsic-curvature contractions, C signatures, or build lists. Compare
generated values with the current kernels and measure representative kernel
timing before replacing the checked-in C. Follow
[characteristic speeds](characteristic-speeds-contract.md),
[HLLE variants](hlle-flux-variant-matrix.md), and
[source terms](source-terms-contract.md) for caller and test behavior.

Keep derivations in [GRMHD equations](../../../docs/raw/derivation.md) and
[Flux Source documentation](../../../docs/raw/Flux_Source.dox). The checked-in
C and active build configuration define executable behavior; the Python file
defines the equations used to regenerate it.
