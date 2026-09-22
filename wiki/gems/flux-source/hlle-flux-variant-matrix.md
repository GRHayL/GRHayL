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

| Variant | Direct functions, `dirn0/1/2` | Fields written | Build modes | Direct replay |
| --- | --- | --- | --- | --- |
| hybrid | `ghl_calculate_HLLE_fluxes_dirn0_hybrid`<br>`ghl_calculate_HLLE_fluxes_dirn1_hybrid`<br>`ghl_calculate_HLLE_fluxes_dirn2_hybrid` | `rho`, `tau`, `SD[0..2]` | default and no-HDF5 | all three in `unit_test_hybrid_flux`; all three in ET Legacy |
| hybrid entropy | `ghl_calculate_HLLE_fluxes_dirn0_hybrid_entropy`<br>`ghl_calculate_HLLE_fluxes_dirn1_hybrid_entropy`<br>`ghl_calculate_HLLE_fluxes_dirn2_hybrid_entropy` | hybrid fields plus `entropy` | default and no-HDF5 | all three in `unit_test_hybrid_flux` |
| tabulated | `ghl_calculate_HLLE_fluxes_dirn0_tabulated`<br>`ghl_calculate_HLLE_fluxes_dirn1_tabulated`<br>`ghl_calculate_HLLE_fluxes_dirn2_tabulated` | hybrid fields plus `Y_e` | default and no-HDF5 | all three in `unit_test_tabulated_flux` |
| tabulated entropy | `ghl_calculate_HLLE_fluxes_dirn0_tabulated_entropy`<br>`ghl_calculate_HLLE_fluxes_dirn1_tabulated_entropy`<br>`ghl_calculate_HLLE_fluxes_dirn2_tabulated_entropy` | hybrid fields plus `Y_e` and `entropy` | default and no-HDF5 | all three in `unit_test_tabulated_flux` |

Every direct symbol has a public declaration, one checked-in definition, and a
variant-manifest entry. [configure](../../../configure) retains the real
tabulated flux definitions under `--disable-hdf5`, and the Ubuntu-Clang
`c2p-failure` no-HDF5 variant checks that they remain linkable. The tabulated
flux test and data generator still require HDF5/table support and are excluded
in that mode.

Variant build lists:

- [GRHayL/Flux_Source/hybrid/make.code.defn](../../../GRHayL/Flux_Source/hybrid/make.code.defn)
- [GRHayL/Flux_Source/hybrid_entropy/make.code.defn](../../../GRHayL/Flux_Source/hybrid_entropy/make.code.defn)
- [GRHayL/Flux_Source/tabulated/make.code.defn](../../../GRHayL/Flux_Source/tabulated/make.code.defn)
- [GRHayL/Flux_Source/tabulated_entropy/make.code.defn](../../../GRHayL/Flux_Source/tabulated_entropy/make.code.defn)

## Direct Functions And Pointer Surface

The direct variant functions above are the public calls. Tests may select them
through test-local function pointers by EOS family, entropy mode, and flux
direction. The former unsuffixed generic globals
`ghl_calculate_HLLE_fluxes_dirn0/1/2` were removed because they were never
initialized and had incompatible `const` primitive signatures. Source clients
that referenced or manually assigned those globals must select a direct
variant and rebuild.

For direction `d`, simple/hybrid callers choose
`ghl_calculate_HLLE_fluxes_dirn<d>_hybrid` or its `_entropy` form; tabulated
callers choose the corresponding `_tabulated` or `_tabulated_entropy` form.
Primitive inputs remain mutable because tabulated thermodynamic callbacks can
clamp them. Supply copies when original face states must be preserved, and
initialize the EOS global callbacks before a direct kernel call.

## Caller Contract

HLLE flux callers must supply:

- `prims_r` and `prims_l`: reconstructed face primitive states with `u0`,
  velocity, magnetic field, density, pressure, and EOS-specific fields already
  valid.
- `eos`: parameters compatible with the active `ghl_compute_h_and_cs2`
  function pointer.
- `metric_face`: face-centered ADM metric.
- `cmin_dirn*` and `cmax_dirn*`: characteristic speeds computed for the same
  direction and face. They are non-negative magnitudes; each kernel divides by
  `cmin + cmax` without a zero-denominator check.
- `cons`: caller-owned conservative output receiving flux components.

Entropy variants read `prims_r->entropy` and `prims_l->entropy` and write
`cons->entropy`. Calling an entropy variant with unset or stale primitive
entropy gives the conservative entropy flux a bad input even when density,
momentum, and energy fields look valid.

Only fields listed in the matrix are written. Other members of `cons` retain
their prior value; these routines do not initialize the whole struct.

All 12 routines call `ghl_compute_h_and_cs2` twice and discard its returned
error code. Primitive arguments are mutable: production tabulated dispatch
clamps `rho`, `Y_e`, and `temperature`, then overwrites `press` and `eps`.
Callers needing immutable reconstructed states must pass copies and validate
EOS/table inputs before calling.

## Tabulated And HDF5 Boundary

Tabulated variants remain built in no-HDF5 targets; only their table-dependent
tests and generators are filtered. Tabulated EOS initialization is separately
guarded by `GHL_DISABLE_HDF5` in
[GRHayL/GRHayL_Core/initialize_eos.c](../../../GRHayL/GRHayL_Core/initialize_eos.c).
Link visibility is not runtime support: these real kernels call the global
`ghl_compute_h_and_cs2` dispatch and discard its error code. A no-HDF5 build
cannot initialize compatible tabulated EOS dispatch, so callers must not invoke
the tabulated variants there; a hybrid or unset global dispatch can otherwise
produce the wrong EOS calculation or a null call.

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

## Evidence Status

- **Build-configured:** every listed direct symbol is in default and no-HDF5
  builds. The Ubuntu-Clang `c2p-failure` no-HDF5 variant link-checks the
  tabulated symbols.
- **Direct replay:** every direction in each row is called by the named test;
  runner and compiler workflows configure those executions.
- **Fixture-generation:** matching data generators call every row/direction,
  but generated outputs use the same implementation and are not an independent
  oracle.
- **Coverage gaps:** ignored EOS error returns; primitive mutation; and zero
  `cmin + cmax` have no focused tests. The
  no-HDF5 matrix variant link-checks the retained algebraic tabulated symbols.
