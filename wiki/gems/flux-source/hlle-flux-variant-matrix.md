# Flux Source HLLE Flux Variant Matrix

## Routing Purpose

Use this page when choosing or auditing Flux_Source HLLE flux variants. It maps
the public direct functions, variant directories, build lists, and EOS pointer
surface without copying generated formulas.

## Public Matrix

All direct public functions are declared in
[GRHayL/include/ghl_flux_source.h](../../../GRHayL/include/ghl_flux_source.h).
Each takes reconstructed right/left face primitives, EOS parameters, a face
metric, direction-specific `cmin`/`cmax`, and a conservative flux output.

| Variant | `dirn0` | `dirn1` | `dirn2` | Directory |
| --- | --- | --- | --- | --- |
| hybrid | `ghl_calculate_HLLE_fluxes_dirn0_hybrid` | `ghl_calculate_HLLE_fluxes_dirn1_hybrid` | `ghl_calculate_HLLE_fluxes_dirn2_hybrid` | [GRHayL/Flux_Source/hybrid/](../../../GRHayL/Flux_Source/hybrid/) |
| hybrid entropy | `ghl_calculate_HLLE_fluxes_dirn0_hybrid_entropy` | `ghl_calculate_HLLE_fluxes_dirn1_hybrid_entropy` | `ghl_calculate_HLLE_fluxes_dirn2_hybrid_entropy` | [GRHayL/Flux_Source/hybrid_entropy/](../../../GRHayL/Flux_Source/hybrid_entropy/) |
| tabulated | `ghl_calculate_HLLE_fluxes_dirn0_tabulated` | `ghl_calculate_HLLE_fluxes_dirn1_tabulated` | `ghl_calculate_HLLE_fluxes_dirn2_tabulated` | [GRHayL/Flux_Source/tabulated/](../../../GRHayL/Flux_Source/tabulated/) |
| tabulated entropy | `ghl_calculate_HLLE_fluxes_dirn0_tabulated_entropy` | `ghl_calculate_HLLE_fluxes_dirn1_tabulated_entropy` | `ghl_calculate_HLLE_fluxes_dirn2_tabulated_entropy` | [GRHayL/Flux_Source/tabulated_entropy/](../../../GRHayL/Flux_Source/tabulated_entropy/) |

Variant build lists:

- [GRHayL/Flux_Source/hybrid/make.code.defn](../../../GRHayL/Flux_Source/hybrid/make.code.defn)
- [GRHayL/Flux_Source/hybrid_entropy/make.code.defn](../../../GRHayL/Flux_Source/hybrid_entropy/make.code.defn)
- [GRHayL/Flux_Source/tabulated/make.code.defn](../../../GRHayL/Flux_Source/tabulated/make.code.defn)
- [GRHayL/Flux_Source/tabulated_entropy/make.code.defn](../../../GRHayL/Flux_Source/tabulated_entropy/make.code.defn)

## Direct Functions And Pointer Surface

The direct variant functions above are source-backed public calls. Tests select
them directly by EOS family, entropy mode, and flux direction.

[GRHayL/include/ghl_eos_functions.h](../../../GRHayL/include/ghl_eos_functions.h)
and
[GRHayL/include/ghl_eos_functions_declaration.h](../../../GRHayL/include/ghl_eos_functions_declaration.h)
also declare generic function pointers named
`ghl_calculate_HLLE_fluxes_dirn0`, `ghl_calculate_HLLE_fluxes_dirn1`, and
`ghl_calculate_HLLE_fluxes_dirn2`. Source review of
[GRHayL/GRHayL_Core/initialize_eos.c](../../../GRHayL/GRHayL_Core/initialize_eos.c)
shows EOS initialization assigning `ghl_compute_h_and_cs2` and
`ghl_con2prim_multi_method`; no GRHayL-local assignment to the generic HLLE
flux pointers is present in the checked-in `GRHayL/` tree. Treat direct variant
symbols as the authoritative route unless downstream integration wires the
generic pointers.

## Caller Contract

HLLE flux callers must supply:

- `prims_r` and `prims_l`: reconstructed face primitive states with `u0`,
  velocity, magnetic field, density, pressure, and EOS-specific fields already
  valid.
- `eos`: parameters compatible with the active `ghl_compute_h_and_cs2`
  function pointer.
- `metric_face`: face-centered ADM metric.
- `cmin_dirn*` and `cmax_dirn*`: characteristic speeds computed for the same
  direction and face.
- `cons`: caller-owned conservative output receiving flux components.

Entropy variants read `prims_r->entropy` and `prims_l->entropy` and write
`cons->entropy`. Calling an entropy variant with unset or stale primitive
entropy gives the conservative entropy flux a bad input even when density,
momentum, and energy fields look valid.

## Tabulated And HDF5 Boundary

Tabulated variants live in checked-in source directories regardless of build
mode, but tabulated EOS initialization is guarded by `GHL_DISABLE_HDF5` in
[GRHayL/GRHayL_Core/initialize_eos.c](../../../GRHayL/GRHayL_Core/initialize_eos.c).
Use this as route-level evidence only: HDF5-disabled behavior belongs to build
and EOS pages, not to duplicated flux formulas here.

## Generated-Source Boundary

[GRHayL/Flux_Source/IGM_All_fluxes.py](../../../GRHayL/Flux_Source/IGM_All_fluxes.py)
contains the local generator naming path for all 12
`ghl_calculate_HLLE_fluxes_dirn*_<variant>.c` files. Review generated C and the
Python source together when formulas, variables, or output fields change.

## Evidence Links

- [GRHayL/include/ghl_flux_source.h](../../../GRHayL/include/ghl_flux_source.h)
- [GRHayL/include/ghl_eos_functions.h](../../../GRHayL/include/ghl_eos_functions.h)
- [GRHayL/include/ghl_eos_functions_declaration.h](../../../GRHayL/include/ghl_eos_functions_declaration.h)
- [GRHayL/GRHayL_Core/initialize_eos.c](../../../GRHayL/GRHayL_Core/initialize_eos.c)
- [GRHayL/Flux_Source/IGM_All_fluxes.py](../../../GRHayL/Flux_Source/IGM_All_fluxes.py)
- [Unit_Tests/unit_test_hybrid_flux.c](../../../Unit_Tests/unit_test_hybrid_flux.c)
- [Unit_Tests/unit_test_tabulated_flux.c](../../../Unit_Tests/unit_test_tabulated_flux.c)
- [Unit_Tests/unit_test_ET_Legacy_flux_source.c](../../../Unit_Tests/unit_test_ET_Legacy_flux_source.c)
- [docs/raw/Flux_Source.dox](../../../docs/raw/Flux_Source.dox)
