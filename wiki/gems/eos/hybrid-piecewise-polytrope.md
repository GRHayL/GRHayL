# Hybrid Piecewise-Polytrope EOS

## Purpose

This page routes the simple ideal-fluid and hybrid piecewise-polytrope EOS
internals. For wrapper and dispatch context, start with
[EOS initialization and dispatch](initialization-and-dispatch.md), then use this
page for the cold-piece helpers, thermal pressure coupling, bounds, and direct
test routes.

Repo-local source and Doxygen own the equations and API signatures. This page
keeps only routing notes and short behavior summaries.

## Simple Ideal-Fluid Setup

The simple EOS uses the hybrid implementation path with one ideal-fluid piece.
`ghl_initialize_simple_eos_functions_and_params` first initializes EOS function
pointers for `ghl_eos_simple`, then calls `ghl_initialize_simple_eos`
(`GRHayL/GRHayL_Core/initialize_eos.c`).

`ghl_initialize_simple_eos` sets:

- `eos_type = ghl_eos_simple`
- `neos = 1`
- `Gamma_th = Gamma_ppoly[0] = Gamma`
- `K_ppoly[0] = 1`
- `rho_ppoly[0] = 0.0`
- `p_ppoly[0] = 0.0`
- `eps_integ_const[0] = 0.0`

It also fills density and pressure atmosphere, floor, and ceiling fields, then
computes `eps_min`, `eps_max`, `eps_atm`, `entropy_min`, `entropy_max`,
`entropy_atm`, and `tau_atm` through the hybrid entropy helper and ideal-fluid
pressure relation (`GRHayL/GRHayL_Core/initialize_eos.c`).

## Hybrid Piecewise-Polytrope Setup

`ghl_initialize_hybrid_eos_functions_and_params` initializes EOS function
pointers for `ghl_eos_hybrid`, then calls `ghl_initialize_hybrid_eos`
(`GRHayL/GRHayL_Core/initialize_eos.c`). The low-level initializer owns these
piecewise-polytrope fields:

- `neos`: number of polytropic pieces, bounded by `MAX_EOS_PARAMS` in
  `GRHayL/include/ghl.h`.
- `rho_ppoly`: density breakpoints copied from initializer input.
- `Gamma_ppoly`: per-piece adiabatic indices copied from initializer input.
- `K_ppoly`: `K_ppoly[0]` comes from initializer input; later entries are
  filled by `ghl_hybrid_set_K_ppoly_and_eps_integ_consts`.
- `eps_integ_const`: initialized with `eps_integ_const[0] = 0.0` for one-piece
  setup and filled for all pieces by
  `ghl_hybrid_set_K_ppoly_and_eps_integ_consts`.
- `p_ppoly`: pressure breakpoints computed after `K_ppoly` and
  `eps_integ_const` are ready, using
  `ghl_hybrid_compute_P_cold_and_eps_cold` on each positive `rho_ppoly`
  breakpoint.
- `Gamma_th`: thermal adiabatic index copied from initializer input.

After piece setup, the hybrid initializer computes pressure, internal-energy,
entropy, and atmosphere/floor/ceiling values from the cold EOS helper and
entropy helper (`GRHayL/GRHayL_Core/initialize_eos.c`). The field definitions
for `neos`, `rho_ppoly`, `Gamma_ppoly`, `K_ppoly`, `eps_integ_const`,
`p_ppoly`, and `Gamma_th` live in `ghl_eos_parameters`
(`GRHayL/include/ghl.h`).

## Helper Map

- Density index lookup: `NRPyEOS_find_polytropic_index` chooses the piece from
  `rho_ppoly` and returns `0` for one-piece EOSs
  (`GRHayL/EOS/Hybrid/NRPyEOS_find_polytropic_index.c`).
- Pressure index lookup: `NRPyEOS_find_polytropic_index_from_P` chooses the
  piece from `p_ppoly` and returns `0` for one-piece EOSs
  (`GRHayL/EOS/Hybrid/NRPyEOS_find_polytropic_index_from_P.c`).
- Piece constants lookup: `NRPyEOS_get_K_and_Gamma` uses the density index to
  return the active `K_ppoly` and `Gamma_ppoly`
  (`GRHayL/EOS/Hybrid/NRPyEOS_get_K_and_gamma.c`).
- Continuity constants: `NRPyEOS_set_K_ppoly_and_eps_integ_consts` fills
  `K_ppoly` and `eps_integ_const` entries after the initializer supplies
  `K_ppoly[0]`, `rho_ppoly`, and `Gamma_ppoly`
  (`GRHayL/EOS/Hybrid/NRPyEOS_set_K_ppoly_and_eps_integ_consts.c`).
- Cold pressure: `NRPyEOS_compute_P_cold` uses the active `K_ppoly` and
  `Gamma_ppoly` for a density
  (`GRHayL/EOS/Hybrid/NRPyEOS_compute_P_cold.c`).
- Cold pressure and cold energy: `NRPyEOS_compute_P_cold_and_eps_cold` first
  clamps density with `NRPyEOS_hybrid_enforce_bounds__rho`, then uses the
  active piece and `eps_integ_const`
  (`GRHayL/EOS/Hybrid/NRPyEOS_compute_P_cold_and_eps_cold.c`,
  `GRHayL/EOS/Hybrid/NRPyEOS_enforce_bounds.c`).
- Thermal pressure role: `NRPyEOS_hybrid_compute_epsilon` computes cold
  pressure/energy first, then uses total pressure minus `P_cold` with
  `Gamma_th` to add the thermal contribution
  (`GRHayL/EOS/Hybrid/NRPyEOS_hybrid_compute_epsilon.c`).
- Entropy helper: `NRPyEOS_compute_entropy_function` selects `Gamma_ppoly` by
  density and returns the hybrid entropy function used by initialization
  (`GRHayL/EOS/Hybrid/NRPyEOS_compute_entropy_function.c`,
  `GRHayL/GRHayL_Core/initialize_eos.c`).
- Rho from pressure: `NRPyEOS_hybrid_compute_rho_cold_from_P_cold` selects the
  pressure piece from `p_ppoly`, computes the corresponding cold density, then
  enforces rho bounds (`GRHayL/EOS/Hybrid/NRPyEOS_hybrid_compute_rho_cold_from_P_cold.c`).
- Rho bounds: `NRPyEOS_hybrid_enforce_bounds__rho` clamps to `rho_min` and
  `rho_max` and reports whether input was already in range
  (`GRHayL/EOS/Hybrid/NRPyEOS_enforce_bounds.c`).
- Enthalpy and sound speed: `NRPyEOS_hybrid_compute_enthalpy_and_cs2` is the
  hybrid entry point behind `ghl_compute_h_and_cs2`; it combines cold pressure,
  cold energy, thermal energy, `Gamma_th`, and the active `Gamma_ppoly`
  (`GRHayL/EOS/Hybrid/NRPyEOS_hybrid_compute_enthalpy_and_cs2.c`,
  `GRHayL/GRHayL_Core/initialize_eos.c`).

`NRPyEOS_initialize_hybrid_functions` wires the public hybrid function pointers
to these helpers, and `GRHayL/EOS/Hybrid/make.code.defn` lists the compiled
hybrid EOS sources.

## Test Coverage

`Unit_Tests/unit_test_piecewise_polytrope.c` directly initializes a four-piece
hybrid EOS and checks computed `K_ppoly` and `eps_integ_const` entries against
reference values from the NRPyEOS Python path. This is the direct regression
route for `ghl_hybrid_set_K_ppoly_and_eps_integ_consts`.

`Unit_Tests/test_compute_h_and_cs2.c` provides a small test helper for
`ghl_compute_h_and_cs2` call sites; it is not the direct piecewise-polytrope
constant test.

## Source-of-truth Paths

- `GRHayL/include/ghl.h`
- `GRHayL/include/ghl_nrpyeos_hybrid.h`
- `GRHayL/GRHayL_Core/initialize_eos.c`
- `GRHayL/EOS/Hybrid/NRPyEOS_initialize_hybrid_functions.c`
- `GRHayL/EOS/Hybrid/NRPyEOS_find_polytropic_index.c`
- `GRHayL/EOS/Hybrid/NRPyEOS_find_polytropic_index_from_P.c`
- `GRHayL/EOS/Hybrid/NRPyEOS_get_K_and_gamma.c`
- `GRHayL/EOS/Hybrid/NRPyEOS_set_K_ppoly_and_eps_integ_consts.c`
- `GRHayL/EOS/Hybrid/NRPyEOS_compute_P_cold.c`
- `GRHayL/EOS/Hybrid/NRPyEOS_compute_P_cold_and_eps_cold.c`
- `GRHayL/EOS/Hybrid/NRPyEOS_hybrid_compute_epsilon.c`
- `GRHayL/EOS/Hybrid/NRPyEOS_compute_entropy_function.c`
- `GRHayL/EOS/Hybrid/NRPyEOS_hybrid_compute_enthalpy_and_cs2.c`
- `GRHayL/EOS/Hybrid/NRPyEOS_enforce_bounds.c`
- `GRHayL/EOS/Hybrid/NRPyEOS_hybrid_compute_rho_cold_from_P_cold.c`
- `GRHayL/EOS/Hybrid/make.code.defn`
- `docs/raw/EOS.dox`
- `Unit_Tests/unit_test_piecewise_polytrope.c`
- `Unit_Tests/test_compute_h_and_cs2.c`
