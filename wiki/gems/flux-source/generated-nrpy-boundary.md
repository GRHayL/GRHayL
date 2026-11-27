# Generated NRPy Boundary

> Staging regeneration workflow (`generate_flux_source.sh`, `generate_flux_source.py`), pinned NRPy 2 usage, and manifest parity for `Flux_Source` C kernels. · Status: confirmed
> Up: [Flux_Source](index.md)

## Summary

`GRHayL/Flux_Source/generate_flux_source.py` is the symbolic NRPy 2 generator for GRHayL's source-term, characteristic-speed, and all four HLLE flux kernel families. `generate_flux_source.sh <empty-staging-dir>` clones a pinned NRPy 2 revision into a temporary checkout and runs the driver into an empty directory outside `GRHayL/Flux_Source`, verifying output filenames against the `make.code.defn` manifests. Checked-in C files govern compiled behavior.

## Detail

### Staging Command And Symbolic Structure

From repository root (with NRPy 2 Python dependencies installed):

```sh
GRHayL/Flux_Source/generate_flux_source.sh /tmp/grhayl-flux-stage
```

The Python driver rejects output paths inside `GRHayL/Flux_Source` and nonempty directories, and checks its emitted file list against root and variant `make.code.defn` manifests. It uses NRPy 2's ADM metric conversions, comoving magnetic norm, magnetic fast-wave estimate, GRHD characteristic roots, and C code generator, grouping `rho_b * h + b^2` before multiplying by `u^\mu u^\nu` to preserve precision in small momentum fluxes under GRHayL's `1 / sqrt(4*pi)` magnetic normalization, while emitting both `_checked` and legacy aborting wrappers.

| Generated C file(s) | Python entry point | Build manifest |
| --- | --- | --- |
| `ghl_calculate_source_terms.c` | `generate_source()` | `GRHayL/Flux_Source/make.code.defn` (Ev: `GRHayL/Flux_Source/generate_flux_source.py::generate_source`) |
| `ghl_calculate_characteristic_speed_dirn{0,1,2}.c` | `generate_speeds()` | `GRHayL/Flux_Source/make.code.defn` (Ev: `GRHayL/Flux_Source/generate_flux_source.py::generate_speeds`) |
| `{hybrid,hybrid_entropy,tabulated,tabulated_entropy}/ghl_calculate_HLLE_fluxes_dirn*_*.c` | `generate_fluxes()` | Matching variant `make.code.defn` (Ev: `GRHayL/Flux_Source/generate_flux_source.py::generate_fluxes`) |

Review staged C, numerical differences, and kernel timing before copying generated output over checked-in C files, and keep Python equations and C files synchronized on any formula, signature, or manifest edit.

Evidence (generated evidence): `GRHayL/Flux_Source/generate_flux_source.sh`; corroboration: `GRHayL/Flux_Source/generate_flux_source.py`.

## Sources

- [GRHayL/Flux_Source/generate_flux_source.sh](../../../GRHayL/Flux_Source/generate_flux_source.sh) - pinned NRPy 2 checkout and invocation wrapper
- [GRHayL/Flux_Source/generate_flux_source.py](../../../GRHayL/Flux_Source/generate_flux_source.py) - `generate_source`, `generate_speeds`, `generate_fluxes`, manifest parity check
- [GRHayL/Flux_Source/make.code.defn](../../../GRHayL/Flux_Source/make.code.defn) - root Flux_Source manifest
- [docs/raw/Flux_Source.dox](../../../docs/raw/Flux_Source.dox) - Flux_Source equations and docs

## See Also

- Depends on: [Generated Boundaries](../../architecture/generated-boundaries.md)
- See also: [Characteristic Speeds](characteristic-speeds.md)
- See also: [HLLE Flux Variant Matrix](hlle-flux-variant-matrix.md)
- See also: [Source Terms](source-terms.md)
