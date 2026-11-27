# Tabulated Interpolator Catalog

> Catalog of built tabulated EOS interpolation wrapper groups, registry audit, beta-equilibrium `rho`-map caches, enthalpy/sound-speed routines, and direct-test coverage distinctions. · Status: confirmed
> Up: [Equation Of State (EOS)](index.md)

## Summary

`GRHayL/EOS/Tabulated/interpolators/make.code.defn` builds tabulated wrapper families keyed by `(rho, Y_e, T)`, `(rho, Y_e, eps)`, `(rho, Y_e, P)`, `(rho, Y_e, S)`, and `(rho, Y_e, h)`, all assigned in `NRPyEOS_initialize_tabulated_functions.c` alongside `ghl_tabulated_get_index_T`, beta-equilibrium `rho`-map caches, and `ghl_tabulated_free_beq_quantities`.

## Detail

### Registry Audit And Wrapper Groups By Input Tuple

All public `ghl_tabulated_*` pointers in `ghl_eos_functions.h` have storage in `ghl_eos_functions_declaration.h` and assignment in `NRPyEOS_initialize_tabulated_functions.c` when HDF5 is enabled, including `ghl_tabulated_free_beq_quantities` (`NRPyEOS_tabulated_free_beq_quantities`, also called by `NRPyEOS_free_memory`) and `ghl_tabulated_get_index_T` (`NRPyEOS_tabulated_get_index.c`, which rejects nonfinite/nonpositive `T` and out-of-range log-temperature endpoints with `-1`). Former prototype-only `rho`/`Y_e` index declarations and the duplicate cleanup prototype were removed.

| Input tuple | Shared engine | Built public wrapper outputs |
| --- | --- | --- |
| `(rho, Y_e, T)` | `NRPyEOS_from_rho_Ye_T_interpolate_n_quantities` | `P`, `eps`, `cs2`, `P+eps`, `P+eps+S`, `P+eps+cs2`, `P+eps+S+cs2`, `P+eps+depsdT`, chemical-potential groups, and nucleon-fraction groups. (Ev: `GRHayL/EOS/Tabulated/interpolators/make.code.defn`) |
| `(rho, Y_e, eps)` | `NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities` (`NRPyEOS_eps_key`) | `T`, `P+T`, `P+S+T`, `P+cs2+T`, `P+S+depsdT+T`. (Ev: `GRHayL/EOS/Tabulated/NRPyEOS_initialize_tabulated_functions.c`) |
| `(rho, Y_e, P)` | Auxiliary engine (`NRPyEOS_press_key`) | `eps+T`, `eps+S+T`, `eps+cs2+T`. |
| `(rho, Y_e, S)` | Auxiliary engine (`NRPyEOS_entropy_key`) | `P+T`, `P+eps+T`. |
| `(rho, Y_e, h)` | Auxiliary engine (`NRPyEOS_enthalpy_key`) | `P+eps+S+T`, `P+eps+dPdrho+dPdeps+T`. |

### Beta-Equilibrium `rho`-Map Caches, Derived Enthalpy, And Test Coverage

- **Beta-equilibrium caches (`NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T.c`):**
  - `NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T` builds `Ye_of_lr` from the temperature-grid bin chosen by `get_index_T` (it does **not** interpolate `munu` to an off-grid `T`), recording the first strict adjacent sign crossing per density (storing minimum `Y_e` and succeeding when no strict crossing exists; an exact zero node is not a crossing).
  - `NRPyEOS_tabulated_compute_Ye_P_eps_of_rho_beq_constant_T` extends the cache with `lp_of_lr`, `le_of_lr`, and `lh_of_lr` evaluated at the requested temperature (so off-grid calls do not establish exact beta equilibrium at `T`).
  - Downstream readers (`Ye_from_rho`, `P_from_rho`, `rho_from_P`, `eps_from_rho`, `dP_drho_from_rho`, `deps_dP_from_rho`) do **not** allocate on demand and assume finite positive arguments, valid caches, monotone pressure for `rho_from_P`, and nonzero `dP_drho` for `deps_dP_from_rho`. A failed combined rebuild can update `Ye_of_lr` while retaining older pressure/energy caches, and a `Ye`-only rebuild leaves older thermodynamic caches stale; rebuild the combined cache whenever all maps must match.
- **Derived enthalpy (`NRPyEOS_tabulate_enthalpy.c`):** `NRPyEOS_tabulate_enthalpy_checked` returns `ghl_error_invalid_eos_table` before `log(h)` on invalid points; legacy `void NRPyEOS_tabulate_enthalpy` discards status and stops at the first invalid entry.
- **Direct-test distinctions (`Unit_Tests/unit_test_tabulated_eos.c`):** does not directly invoke `ghl_tabulated_compute_P_S_T_from_eps`, either `from_h` wrapper, `ghl_tabulated_compute_Ye_of_rho_beq_constant_T`, the `T`/`eps`/`S` enforce pointers, `ghl_tabulated_get_index_T`, `ghl_tabulated_read_table_set_EOS_params`, or the assigned `ghl_tabulated_free_beq_quantities` pointer (several run indirectly via initialization or cleanup).

Evidence (descriptive behavior): `GRHayL/EOS/Tabulated/NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T.c`; corroboration: `GRHayL/EOS/Tabulated/NRPyEOS_tabulated_get_index.c::NRPyEOS_tabulated_get_index_T`.

## Sources

- [GRHayL/include/ghl_nrpyeos_tabulated.h](../../../GRHayL/include/ghl_nrpyeos_tabulated.h) - NRPyEOS tabulated declarations
- [GRHayL/include/ghl_eos_functions.h](../../../GRHayL/include/ghl_eos_functions.h) - `ghl_tabulated_*` public pointers
- [GRHayL/EOS/Tabulated/NRPyEOS_initialize_tabulated_functions.c](../../../GRHayL/EOS/Tabulated/NRPyEOS_initialize_tabulated_functions.c) - pointer assignments
- [GRHayL/EOS/Tabulated/interpolators/make.code.defn](../../../GRHayL/EOS/Tabulated/interpolators/make.code.defn) - built interpolator files
- [GRHayL/EOS/Tabulated/NRPyEOS_tabulated_get_index.c](../../../GRHayL/EOS/Tabulated/NRPyEOS_tabulated_get_index.c) - `NRPyEOS_tabulated_get_index_T`
- [GRHayL/EOS/Tabulated/NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T.c](../../../GRHayL/EOS/Tabulated/NRPyEOS_tabulated_compute_Ye_of_rho_beq_constant_T.c) - beta-equilibrium `rho`-map caches and cleanup
- [GRHayL/EOS/Tabulated/NRPyEOS_tabulate_enthalpy.c](../../../GRHayL/EOS/Tabulated/NRPyEOS_tabulate_enthalpy.c) - checked and legacy enthalpy tabulation
- [Unit_Tests/unit_test_tabulated_eos.c](../../../Unit_Tests/unit_test_tabulated_eos.c) - direct tabulated wrapper tests

## See Also

- Depends on: [Tabulated Interpolation And Bounds](tabulated-interpolation-and-bounds.md)
- Depends on: [Tabulated Table Layout](tabulated-table-layout.md)
- Validated by: [EOS Tests And Reference Data](tests-and-reference-data.md)
