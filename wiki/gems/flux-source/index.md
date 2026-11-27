# Flux_Source

> Router for GRMHD characteristic speeds, directional HLLE flux variants, geometric source terms, NRPy 2 generation, and Flux_Source tests.

| Page | Go here when... |
| --- | --- |
| [Characteristic Speeds](characteristic-speeds.md) | You need `ghl_calculate_characteristic_speed_dirn{0,1,2}` (`_checked`/legacy), sign conventions, EOS mutation, or Induction coupling. |
| [HLLE Flux Variant Matrix](hlle-flux-variant-matrix.md) | You need directional hybrid/tabulated/entropy HLLE functions, wave-speed validation, deprecated generic globals, or no-HDF5 behavior. |
| [Source Terms](source-terms.md) | You need `ghl_calculate_source_terms` (`_checked`/legacy), caller-supplied metric derivatives, output fields, or EOS callback mutation. |
| [Generated NRPy Boundary](generated-nrpy-boundary.md) | You need `generate_flux_source.sh` / `generate_flux_source.py`, pinned NRPy 2 staging, or manifest parity across generated C kernels. |
| [Flux_Source Tests And Reference Data](tests-and-reference-data.md) | You need hybrid/tabulated/ET_Legacy/CompOSE flux tests, perturbation scales, synthetic `test_compute_h_and_cs2.c`, or Induction HLL exclusion. |

Back to [GRHayL Gems](../index.md).
