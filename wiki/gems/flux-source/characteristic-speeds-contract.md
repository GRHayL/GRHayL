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
compiles the `dirn0`, `dirn1`, and `dirn2` files.

Every row has a legacy `void` entry point and a matching `_checked` entry point
returning `ghl_error_codes_t`; the legacy wrapper aborts on a checked error.
Both are declared publicly and share a checked-in definition file. The files
have normal and `--disable-hdf5` build membership. Data generators call checked
entry points, while the ET Legacy replay exercises legacy wrappers. The
ordinary hybrid and tabulated flux tests do **not** call these routines; they
read precomputed speed arrays from their fixtures.

## Caller Contract

Each direction function takes:

- `prims_r` and `prims_l`: right and left primitive states already
  reconstructed to the face for that direction.
- `eos`: EOS parameters compatible with the active EOS function table.
- `metric_face`: ADM metric quantities evaluated at the same face.
- `cmin_dirn*` and `cmax_dirn*`: caller-owned output pointers receiving the
  non-negative left-going and right-going speed magnitudes used by the HLL
  formulas. Despite the names, `cmin` is not a signed minimum eigenvalue: the
  generator defines it as the negated minimum of zero and the two left-going
  speeds. Branch-free min/max evaluation can leave negative roundoff residues;
  the HLLE routines floor those residues at zero within their tolerance.

The kernels call `ghl_compute_h_and_cs2` for both reconstructed states. That
function pointer is declared in
[GRHayL/include/ghl_eos_functions.h](../../../GRHayL/include/ghl_eos_functions.h)
and assigned to hybrid or tabulated enthalpy/sound-speed implementations in
[GRHayL/GRHayL_Core/initialize_eos.c](../../../GRHayL/GRHayL_Core/initialize_eos.c).
Keep EOS behavior routed through EOS pages; this page only records that the
speed kernels depend on `h` and `cs2`.

The primitive pointers are intentionally non-`const`. The production tabulated
enthalpy/sound-speed implementation limits `rho`, `Y_e`, and `temperature` to table bounds and
recomputes `press` and `eps` in place. Checked speed routines return the exact
`ghl_error_codes_t` from either callback and leave both speed outputs unchanged
on failure; legacy wrappers abort instead. A successful first callback may mutate its primitive before a
failing second callback; rollback is not promised. Callers needing unchanged
face states must pass copies.

The outputs are direction-specific `cmin` and `cmax` speeds. HLLE flux kernels
consume these same values through their `cmin_dirn*` and `cmax_dirn*`
arguments.

## Generated-Source Boundary

[GRHayL/Flux_Source/generate_flux_source.py](../../../GRHayL/Flux_Source/generate_flux_source.py)
contains the NRPy 2 generator path for
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

This is a contract coupling, not end-to-end test evidence. The Induction HLL
fixture generators fill `cmin` and `cmax` with random values (including negative
values) and never call these characteristic-speed routines. No visible test
passes production characteristic-speed outputs into an Induction HLL routine.

## Evidence Status

- **Built:** each listed declaration, definition, and manifest entry exists in
  both HDF5 modes.
- **Direct execution route:** `unit_test_ET_Legacy_flux_source` calls every
  listed function, and the ordinary runner plus every compiler workflow configures that
  replay. That test installs its own `compute_h_and_cs2` callback rather than
  production EOS dispatch.
- **Fixture-generation only:** hybrid and tabulated data generators call every
  listed function; normal test jobs download rather than regenerate those fixtures.
- **Replay only:** `unit_test_hybrid_flux` and `unit_test_tabulated_flux` consume
  stored speed arrays.
- **Production-EOS evidence:** `unit_test_tabulated_eos_compose` checks an
  analytic magnetized speed bound in every direction and the tabulated
  table-bound/mutation behavior. The oracle is implemented separately but uses the
  same magnetosonic model as the kernel.
- **Focused errors:** `unit_test_hybrid_flux` injects first- and second-callback
  failures into every checked direction and verifies unchanged outputs.
- **Coverage gap:** no Flux_Source-to-Induction end-to-end test.

## Evidence Links

- [GRHayL/include/ghl_flux_source.h](../../../GRHayL/include/ghl_flux_source.h)
- [GRHayL/include/ghl_eos_functions.h](../../../GRHayL/include/ghl_eos_functions.h)
- [GRHayL/include/ghl_induction.h](../../../GRHayL/include/ghl_induction.h)
- [GRHayL/GRHayL_Core/initialize_eos.c](../../../GRHayL/GRHayL_Core/initialize_eos.c)
- [GRHayL/Flux_Source/make.code.defn](../../../GRHayL/Flux_Source/make.code.defn)
- [GRHayL/Flux_Source/generate_flux_source.py](../../../GRHayL/Flux_Source/generate_flux_source.py)
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
