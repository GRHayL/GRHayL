# Flux Source Characteristic Speeds Contract

## Routing Purpose

Use this page when checking Flux_Source characteristic-speed calls,
direction-specific generated kernels, or Induction HLL speed coupling. This is
a caller contract and route map; formulas stay in source, Doxygen, and the
derivation notes.

## Public API

The public declarations live in
[GRHayL/include/ghl_flux_source.h](../../../GRHayL/include/ghl_flux_source.h):

| Direction | Function | Built source |
| --- | --- | --- |
| `dirn0` | `ghl_calculate_characteristic_speed_dirn0` | [GRHayL/Flux_Source/ghl_calculate_characteristic_speed_dirn0.c](../../../GRHayL/Flux_Source/ghl_calculate_characteristic_speed_dirn0.c) |
| `dirn1` | `ghl_calculate_characteristic_speed_dirn1` | [GRHayL/Flux_Source/ghl_calculate_characteristic_speed_dirn1.c](../../../GRHayL/Flux_Source/ghl_calculate_characteristic_speed_dirn1.c) |
| `dirn2` | `ghl_calculate_characteristic_speed_dirn2` | [GRHayL/Flux_Source/ghl_calculate_characteristic_speed_dirn2.c](../../../GRHayL/Flux_Source/ghl_calculate_characteristic_speed_dirn2.c) |

The parent build list
[GRHayL/Flux_Source/make.code.defn](../../../GRHayL/Flux_Source/make.code.defn)
compiles all three direction files.

## Caller Contract

Each direction function takes:

- `prims_r` and `prims_l`: right and left primitive states already
  reconstructed to the face for that direction.
- `eos`: EOS parameters compatible with the active EOS function table.
- `metric_face`: ADM metric quantities evaluated at the same face.
- `cmin_dirn*` and `cmax_dirn*`: caller-owned output pointers receiving the
  minimum and maximum characteristic speeds.

The kernels call `ghl_compute_h_and_cs2` for both reconstructed states. That
function pointer is declared in
[GRHayL/include/ghl_eos_functions.h](../../../GRHayL/include/ghl_eos_functions.h)
and assigned to hybrid or tabulated enthalpy/sound-speed implementations in
[GRHayL/GRHayL_Core/initialize_eos.c](../../../GRHayL/GRHayL_Core/initialize_eos.c).
Keep EOS behavior routed through EOS pages; this page only records that the
speed kernels depend on `h` and `cs2`.

The outputs are direction-specific `cmin` and `cmax` speeds. HLLE flux kernels
consume these same values through their `cmin_dirn*` and `cmax_dirn*`
arguments.

## Generated-Source Boundary

[GRHayL/Flux_Source/IGM_Characteristic_Speeds.py](../../../GRHayL/Flux_Source/IGM_Characteristic_Speeds.py)
contains the local generator naming path for
`ghl_calculate_characteristic_speed_dirn0`, `dirn1`, and `dirn2`. Treat the
checked-in C files and the Python source as coupled evidence. Do not copy or
hand-expand generated formulas into KB pages.

## Induction Coupling

Induction HLL packing uses characteristic speeds produced under the Flux_Source
direction convention. The coupling is documented in
[GRHayL/include/ghl_induction.h](../../../GRHayL/include/ghl_induction.h),
where `ghl_HLL_vars` stores `c1_min`, `c1_max`, `c2_min`, and `c2_max` and
points readers back to the three `ghl_calculate_characteristic_speed_dirn*`
functions. Direction or sign changes here can therefore affect hydrodynamic
HLLE fluxes and vector-potential HLL flux setup.

## Evidence Links

- [GRHayL/include/ghl_flux_source.h](../../../GRHayL/include/ghl_flux_source.h)
- [GRHayL/include/ghl_eos_functions.h](../../../GRHayL/include/ghl_eos_functions.h)
- [GRHayL/include/ghl_induction.h](../../../GRHayL/include/ghl_induction.h)
- [GRHayL/GRHayL_Core/initialize_eos.c](../../../GRHayL/GRHayL_Core/initialize_eos.c)
- [GRHayL/Flux_Source/make.code.defn](../../../GRHayL/Flux_Source/make.code.defn)
- [GRHayL/Flux_Source/IGM_Characteristic_Speeds.py](../../../GRHayL/Flux_Source/IGM_Characteristic_Speeds.py)
- [Unit_Tests/data_gen/unit_test_data_hybrid_flux.c](../../../Unit_Tests/data_gen/unit_test_data_hybrid_flux.c)
  and
  [Unit_Tests/data_gen/unit_test_data_tabulated_flux.c](../../../Unit_Tests/data_gen/unit_test_data_tabulated_flux.c)
  call the characteristic-speed kernels while generating flux fixtures.
- [Unit_Tests/unit_test_ET_Legacy_flux_source.c](../../../Unit_Tests/unit_test_ET_Legacy_flux_source.c)
  calls the characteristic-speed kernels directly before hybrid HLLE fluxes.
- [Unit_Tests/unit_test_hybrid_flux.c](../../../Unit_Tests/unit_test_hybrid_flux.c)
  and
  [Unit_Tests/unit_test_tabulated_flux.c](../../../Unit_Tests/unit_test_tabulated_flux.c)
  replay precomputed `cmin`/`cmax` fixtures.
- [Unit_Tests/unit_test_HLL_flux.c](../../../Unit_Tests/unit_test_HLL_flux.c)
  is Induction HLL evidence only.
- [docs/raw/Flux_Source.dox](../../../docs/raw/Flux_Source.dox)
- [docs/raw/derivation.md](../../../docs/raw/derivation.md)
