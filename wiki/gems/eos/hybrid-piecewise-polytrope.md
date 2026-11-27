# Hybrid Piecewise-Polytrope EOS

> Simple ideal-fluid and hybrid piecewise-polytrope constants, cold/thermal EOS helpers, breakpoint extents, and direct unit tests. · Status: confirmed
> Up: [Equation Of State (EOS)](index.md)

## Summary

`GRHayL/EOS/Hybrid/` implements GRHayL's code-unit simple ideal-fluid (`neos = 1`) and piecewise-polytrope hybrid EOS routines, splitting cold polytropic pressure/energy from a `Gamma_th` thermal component across `neos - 1` density and pressure transitions. Unlike tabulated EOS, hybrid enthalpy/sound-speed evaluation performs no unit conversion, finiteness validation, or sound-speed cleaning.

## Detail

### Piecewise-Polytrope Setup And Helper Map

Simple EOS (`ghl_eos_simple`) sets `neos = 1`, `Gamma_th = Gamma_ppoly[0] = Gamma`, `K_ppoly[0] = 1`, and `rho_ppoly[0] = p_ppoly[0] = eps_integ_const[0] = 0.0`. Hybrid EOS (`ghl_eos_hybrid`) copies `neos` gammas and `neos - 1` density breakpoints (`rho_ppoly[0..neos - 2]`, never reading a one-piece or terminal breakpoint), computes continuous `K_ppoly` and `eps_integ_const` (`eps_integ_const[0] = 0.0`), and computes `neos - 1` pressure transitions in `p_ppoly` (leaving `p_ppoly[neos - 1] = 0.0`). All values stay in caller-supplied GRHayL code units.

| Helper | File | Behavior & preconditions |
| --- | --- | --- |
| `NRPyEOS_find_polytropic_index` / `_from_P` | `NRPyEOS_find_polytropic_index.c`, `NRPyEOS_find_polytropic_index_from_P.c` | Selects piece index from `neos - 1` `rho_ppoly` or `p_ppoly` transitions; returns `0` for one-piece EOS. (Ev: `GRHayL/EOS/Hybrid/NRPyEOS_find_polytropic_index.c`) |
| `NRPyEOS_get_K_and_Gamma` / `NRPyEOS_set_K_ppoly_and_eps_integ_consts` | `NRPyEOS_get_K_and_gamma.c`, `NRPyEOS_set_K_ppoly_and_eps_integ_consts.c` | Looks up active `(K, Gamma)` or enforces pressure/energy continuity across transitions. (Ev: `GRHayL/EOS/Hybrid/NRPyEOS_set_K_ppoly_and_eps_integ_consts.c`) |
| `NRPyEOS_compute_P_cold` / `_and_eps_cold` | `NRPyEOS_compute_P_cold.c`, `NRPyEOS_compute_P_cold_and_eps_cold.c` | Evaluates cold polytropic `P_cold` and `eps_cold`; `_and_eps_cold` clamps density first via `NRPyEOS_hybrid_enforce_bounds__rho`. (Ev: `GRHayL/EOS/Hybrid/NRPyEOS_compute_P_cold_and_eps_cold.c`) |
| `NRPyEOS_hybrid_compute_epsilon` / `NRPyEOS_compute_entropy_function` | `NRPyEOS_hybrid_compute_epsilon.c`, `NRPyEOS_compute_entropy_function.c` | Adds thermal energy `(P - P_cold) / ((Gamma_th - 1) * rho)` using bounded cold values and original `rho`; computes `P / rho^(Gamma_ppoly - 1)`. (Ev: `GRHayL/EOS/Hybrid/NRPyEOS_hybrid_compute_epsilon.c`) |
| `NRPyEOS_hybrid_compute_rho_cold_from_P_cold` / `NRPyEOS_hybrid_enforce_bounds__rho` | `NRPyEOS_hybrid_compute_rho_cold_from_P_cold.c`, `NRPyEOS_enforce_bounds.c` | Inverts `P_cold` (requires nonzero `K_ppoly`) and clamps `rho` to `[rho_min, rho_max]`. (Ev: `GRHayL/EOS/Hybrid/NRPyEOS_enforce_bounds.c`) |
| `NRPyEOS_hybrid_compute_enthalpy_and_cs2` | `NRPyEOS_hybrid_compute_enthalpy_and_cs2.c` | Backs `ghl_compute_h_and_cs2`; returns `ghl_success` unconditionally without clamping `prims->rho`, checking finiteness, or cleaning negative/superluminal `cs2`; mixes bounded cold terms with original `rho` if caller supplies out-of-domain density. (Ev: `GRHayL/EOS/Hybrid/NRPyEOS_hybrid_compute_enthalpy_and_cs2.c`) |

### Preconditions And Direct Test Coverage

Every `ghl_hybrid_*` declaration in `ghl_eos_functions.h` has matching storage and assignment in `NRPyEOS_initialize_hybrid_functions.c`. Successful initialization does not exclude solver-specific singularities (such as Noble recovery with `Gamma_th = 0`).

`Unit_Tests/unit_test_piecewise_polytrope.c` tests four-piece `K_ppoly`, `eps_integ_const`, and all three `p_ppoly` transitions against independent adjacent-piece formulas, plus a one-piece sentinel and an exact-extent two-piece breakpoint case; invalid `neos`, breakpoint ordering, and gamma-domain errors remain test gaps.

Evidence (descriptive behavior): `GRHayL/EOS/Hybrid/make.code.defn`; corroboration: `Unit_Tests/unit_test_piecewise_polytrope.c`.

## Sources

- [GRHayL/include/ghl_nrpyeos_hybrid.h](../../../GRHayL/include/ghl_nrpyeos_hybrid.h) - hybrid NRPyEOS declarations
- [GRHayL/EOS/Hybrid/make.code.defn](../../../GRHayL/EOS/Hybrid/make.code.defn) - compiled hybrid sources
- [GRHayL/EOS/Hybrid/NRPyEOS_find_polytropic_index.c](../../../GRHayL/EOS/Hybrid/NRPyEOS_find_polytropic_index.c) - density piece index
- [GRHayL/EOS/Hybrid/NRPyEOS_find_polytropic_index_from_P.c](../../../GRHayL/EOS/Hybrid/NRPyEOS_find_polytropic_index_from_P.c) - pressure piece index
- [GRHayL/EOS/Hybrid/NRPyEOS_set_K_ppoly_and_eps_integ_consts.c](../../../GRHayL/EOS/Hybrid/NRPyEOS_set_K_ppoly_and_eps_integ_consts.c) - continuity constants
- [GRHayL/EOS/Hybrid/NRPyEOS_compute_P_cold_and_eps_cold.c](../../../GRHayL/EOS/Hybrid/NRPyEOS_compute_P_cold_and_eps_cold.c) - cold `P` and `eps`
- [GRHayL/EOS/Hybrid/NRPyEOS_hybrid_compute_epsilon.c](../../../GRHayL/EOS/Hybrid/NRPyEOS_hybrid_compute_epsilon.c) - thermal energy helper
- [GRHayL/EOS/Hybrid/NRPyEOS_enforce_bounds.c](../../../GRHayL/EOS/Hybrid/NRPyEOS_enforce_bounds.c) - `rho` clamping
- [GRHayL/EOS/Hybrid/NRPyEOS_hybrid_compute_enthalpy_and_cs2.c](../../../GRHayL/EOS/Hybrid/NRPyEOS_hybrid_compute_enthalpy_and_cs2.c) - hybrid enthalpy and sound speed
- [Unit_Tests/unit_test_piecewise_polytrope.c](../../../Unit_Tests/unit_test_piecewise_polytrope.c) - piecewise-polytrope unit test

## See Also

- Depends on: [EOS Initialization And Dispatch](initialization-and-dispatch.md)
- Validated by: [EOS Tests And Reference Data](tests-and-reference-data.md)
- Contrasts with: [Tabulated Table Layout](tabulated-table-layout.md)
