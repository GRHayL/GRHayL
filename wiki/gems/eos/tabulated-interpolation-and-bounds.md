# Tabulated EOS Interpolation And Bounds

## Purpose

This page routes tabulated EOS interpolation, temperature recovery, table-bound
handling, beta-equilibrium rho-map helpers, and enthalpy/sound-speed entry
points. It is not a full API reference; exact signatures stay in
`GRHayL/include/ghl_nrpyeos_tabulated.h` and `GRHayL/include/ghl_eos_functions.h`.
For table layout, table keys, HDF5 adapter, energy shift, and table lifecycle
context, read [tabulated table contract](tabulated-table-contract.md)
(`wiki/gems/eos/tabulated-table-contract.md`).

## Wrapper Families

Public tabulated wrappers are declared in `GRHayL/include/ghl_eos_functions.h`,
implemented through NRPyEOS functions declared in
`GRHayL/include/ghl_nrpyeos_tabulated.h`, and assigned in
`GRHayL/EOS/Tabulated/NRPyEOS_initialize_tabulated_functions.c`.

- Direct from `(rho,Y_e,T)`: routes pressure, energy, sound speed, entropy,
  `depsdT`, chemical potentials, and nucleon fractions through
  `GRHayL/EOS/Tabulated/interpolators/NRPyEOS_from_rho_Ye_T_interpolate_n_quantities.c`
  plus the specific wrapper files in
  `GRHayL/EOS/Tabulated/interpolators/`.
- From `eps`: routes temperature recovery from `(rho,Y_e,eps)`, then
  interpolates requested outputs through
  `GRHayL/EOS/Tabulated/interpolators/NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities.c`.
  Public names include `ghl_tabulated_compute_T_from_eps`,
  `ghl_tabulated_compute_P_T_from_eps`,
  `ghl_tabulated_compute_P_S_T_from_eps`,
  `ghl_tabulated_compute_P_cs2_T_from_eps`, and
  `ghl_tabulated_compute_P_S_depsdT_T_from_eps`.
- Wrappers from `P`: routes temperature recovery from `(rho,Y_e,P)` through the same
  auxiliary temperature finder. Public names include
  `ghl_tabulated_compute_eps_T_from_P`,
  `ghl_tabulated_compute_eps_cs2_T_from_P`, and
  `ghl_tabulated_compute_eps_S_T_from_P`.
- Wrappers from `S`: routes temperature recovery from `(rho,Y_e,S)`. Public names
  include `ghl_tabulated_compute_P_T_from_S` and
  `ghl_tabulated_compute_P_eps_T_from_S`.
- Wrappers from `h`: routes temperature recovery from `(rho,Y_e,h)`. Public names
  include `ghl_tabulated_compute_P_eps_S_T_from_h` and
  `ghl_tabulated_compute_P_eps_dPdrho_dPdeps_T_from_h`.
- Beta-equilibrium rho-map helpers: `ghl_tabulated_compute_Ye_of_rho_beq_constant_T`
  and `ghl_tabulated_compute_Ye_P_eps_of_rho_beq_constant_T` build cached
  `rho` maps in
  `GRHayL/EOS/Tabulated/NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T.c`.
  Consumers include `ghl_tabulated_compute_Ye_from_rho`,
  `ghl_tabulated_compute_P_from_rho`, `ghl_tabulated_compute_rho_from_P`,
  `ghl_tabulated_compute_eps_from_rho`, `ghl_tabulated_compute_dP_drho_from_rho`,
  and `ghl_tabulated_compute_deps_dP_from_rho`.

`GRHayL/EOS/Tabulated/interpolators/make.code.defn` owns the built
interpolator source list for these wrapper families.

## Current StellarCollapse Grid Contract

The live StellarCollapse path stores coordinate arrays as `logrho`, `logtemp`,
and `ye`. For current interpolation to be correct, each stored coordinate
array must be uniformly linearly spaced, strictly increasing, and contain at
least two points. The physical density and temperature values must be
positive, so uniform spacing in `logrho` and `logtemp` means logarithmically
spaced physical `rho` and `T`; `ye` is uniformly spaced directly.

`NRPyEOS_read_table_set_EOS_params.c` takes coordinate bounds from the first
and last array entries and derives each inverse spacing from the first
interval. `NRPyEOS_tabulated_helpers.h` then uses those values to select an
adjacent cell and perform eight-point trilinear interpolation. The current
path does not recover a nonuniform grid from all coordinate entries.

For CompOSE conversion-specific consequences, including the density and
charge-coordinate maps, read the
[CompOSE EOS adapter how-to](../neutrinos/compose-eos-adapter-how-to.md#grid-and-interpolation-contract).

## Temperature Recovery

Auxiliary-input wrappers recover temperature before interpolation. They pass
`eos->root_finding_precision` from `ghl_eos_parameters` into
`NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities`, which then calls
the helper routines in
`GRHayL/EOS/Tabulated/interpolators/NRPyEOS_tabulated_helpers.h`. The helper
path uses `NRPyEOS_findtemp_from_any` and can fall back to `NRPyEOS_bisection`
when interpolation/Newton-style updates do not converge cleanly.

The direct `(rho,Y_e,T)` path does not recover temperature; it checks table
bounds, finds interpolation spots, and interpolates requested table quantities.
For pressure, enthalpy, and internal energy, the direct helper applies the
source-owned log/energy-shift conversions before returning values.

## Bounds Paths

There are two different bounds contracts:

- Hard-error checks: `NRPyEOS_checkbounds` and
  `NRPyEOS_checkbounds_kt0_noTcheck` live in
  `GRHayL/EOS/Tabulated/interpolators/NRPyEOS_tabulated_helpers.h`. Direct
  `(rho,Y_e,T)` interpolation checks `rho`, `Y_e`, and `T`; auxiliary-input
  interpolation checks `rho` and `Y_e` before recovering `T`. Failures return
  table min/max error codes rather than mutating caller inputs.
- Clamp/enforce helpers: `NRPyEOS_enforce_table_bounds_rho_Ye_T`,
  `NRPyEOS_enforce_table_bounds_rho_Ye_eps`,
  `NRPyEOS_enforce_table_bounds_rho_Ye_S`, and
  `NRPyEOS_enforce_table_bounds_rho_Ye_P` live in
  `GRHayL/EOS/Tabulated/NRPyEOS_enforce_table_bounds.c`. These mutate passed
  variables into EOS bounds for `(rho,Y_e,T)`, `(rho,Y_e,eps)`,
  `(rho,Y_e,S)`, and `(rho,Y_e,P)`.

Use hard-error checks when out-of-table input should fail fast. Use enforce
helpers where the caller contract is to clamp to configured `ghl_eos_parameters`
limits before continuing.

## Enthalpy And cs2

`NRPyEOS_tabulated_compute_enthalpy_and_cs2` in
`GRHayL/EOS/Tabulated/NRPyEOS_tabulated_compute_enthalpy_and_cs2.c` can mutate
`ghl_primitive_quantities`: it enforces `rho`, `Y_e`, and `temperature` bounds
through `ghl_tabulated_enforce_bounds_rho_Ye_T`, then computes pressure, energy,
and `cs2` through the tabulated `(rho,Y_e,T)` wrapper before forming enthalpy.

## Coverage Routes

- `Unit_Tests/unit_test_tabulated_eos.c` validates direct `(rho,Y_e,T)`
  wrappers, auxiliary-input wrappers from `eps`, `P`, and `S`, the enthalpy and
  `cs2` path through `ghl_compute_h_and_cs2`, selected enforce behavior, and
  beta-equilibrium rho-map helpers.
- `Unit_Tests/unit_test_code_error.c` links error coverage for out-of-table
  direct interpolation, out-of-table auxiliary temperature-recovery wrappers,
  too-many-variable calls to both shared interpolation helpers, and
  beta-equilibrium rho-map out-of-range cases.

## Source-Of-Truth Paths

- `GRHayL/include/ghl_nrpyeos_tabulated.h`
- `GRHayL/include/ghl_eos_functions.h`
- `GRHayL/EOS/Tabulated/NRPyEOS_initialize_tabulated_functions.c`
- `GRHayL/EOS/Tabulated/interpolators/`
- `GRHayL/EOS/Tabulated/interpolators/NRPyEOS_tabulated_helpers.h`
- `GRHayL/EOS/Tabulated/interpolators/NRPyEOS_from_rho_Ye_T_interpolate_n_quantities.c`
- `GRHayL/EOS/Tabulated/interpolators/NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities.c`
- `GRHayL/EOS/Tabulated/NRPyEOS_enforce_table_bounds.c`
- `GRHayL/EOS/Tabulated/NRPyEOS_tabulated_compute_enthalpy_and_cs2.c`
- `GRHayL/EOS/Tabulated/NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T.c`
- `GRHayL/EOS/Tabulated/NRPyEOS_tabulated_get_index.c`
- `GRHayL/EOS/Tabulated/interpolators/make.code.defn`
- `Unit_Tests/unit_test_tabulated_eos.c`
- `Unit_Tests/unit_test_code_error.c`
