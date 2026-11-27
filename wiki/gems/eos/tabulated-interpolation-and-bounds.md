# Tabulated Interpolation And Bounds

> Uniform logarithmic grid requirements, direct and auxiliary-inversion interpolation, temperature recovery, bounds-check versus clamp helpers, and `NRPyEOS_tabulated_compute_enthalpy_and_cs2`. · Status: confirmed
> Up: [Equation Of State (EOS)](index.md)

## Summary

GRHayL's tabulated EOS interpolators assume strictly increasing uniform `logrho`, `logtemp`, and `ye` axes, derive inverse spacings from the first interval, and perform eight-point trilinear interpolation. Auxiliary wrappers (`from_eps`, `from_P`, `from_S`, `from_h`) invert temperature via `NRPyEOS_findtemp_from_any` with bisection fallback before direct interpolation, and hard-error table-bound checks are separate from configured-bound clamps.

## Detail

### Uniform Grid Assumption And Caller Input Domain

`NRPyEOS_read_table_set_EOS_params.c` (via `NRPyEOS_stellarcollapse.c::check_axis` and `validate_biased_axis`) validates that `logrho`, `logtemp`, and `ye` are finite, strictly increasing, uniform within `1e-8` index units of their first spacing, and have at least two points (logarithmic physical `rho` and `T`, linear `Y_e`), taking inverse spacings strictly from the first cell—nonuniform grids are unsupported.

**Caller preconditions (`NRPyEOS_tabulated_helpers.h`):**
- Supply finite coordinates and auxiliary targets with positive `rho` and `T`: hard-error checks (`NRPyEOS_checkbounds`, `NRPyEOS_checkbounds_kt0_noTcheck`) use ordered comparisons without finiteness checks, so `NaN` can bypass them into cell indexing.
- Auxiliary inversions take `*T` as an in/out initial guess: supply a finite positive in-table seed before `log(*T)`. Pressure and enthalpy targets must be positive; energy targets require `eps + energy_shift > 0` (source rejects negative shifted energy but does not reject exact zero before `log`).
- Supply finite positive `eos->root_finding_precision`. Direct shared-helper callers (`NRPyEOS_from_rho_Ye_T_interpolate_n_quantities`, `NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities`) must pass valid keys and `0 <= n <= NRPyEOS_ntablekeys` (only the upper bound is checked).

Evidence (descriptive behavior): `GRHayL/EOS/Tabulated/interpolators/NRPyEOS_tabulated_helpers.h::NRPyEOS_checkbounds`; corroboration: `GRHayL/EOS/Tabulated/interpolators/NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities.c`.

### Temperature Recovery, Bounds Paths, And Enthalpy/Sound-Speed Mutation

- **Temperature recovery (`NRPyEOS_findtemp_from_any` / `NRPyEOS_bisection`):** Auxiliary wrappers check `(rho, Y_e)` bounds first, solve for `T` to `eos->root_finding_precision`, write `*T`, and then invoke direct `(rho, Y_e, T)` interpolation (applying log/energy-shift conversions for `P`, `h`, `eps`). If bisection is reached and the target is at or below the interpolated value at table minimum temperature, bisection returns `table_T_min = exp(eos->table_logT[0])` (not narrower configured `eos->T_min`) with success status—so success does not guarantee exact inversion at the cold boundary. If the final direct interpolation fails, `*T` stays mutated while other outputs stay unassigned.
- **Hard errors versus clamps:** `NRPyEOS_checkbounds*` check **physical table bounds** and return error codes without mutating inputs. `NRPyEOS_enforce_table_bounds_rho_Ye_{T,eps,S,P}` (`NRPyEOS_enforce_table_bounds.c`) clamp caller variables in place to **configured `ghl_eos_parameters` bounds** (which may be narrower).
- **Enthalpy & sound speed (`NRPyEOS_tabulated_compute_enthalpy_and_cs2.c`):** clamps `prims->rho`, `Y_e`, and `temperature` via `ghl_tabulated_enforce_bounds_rho_Ye_T`, interpolates `press`, `eps`, and `cs2`, and computes/writes enthalpy even if interpolation returns an error (leaving primitive mutations); callers must check the returned status.

Evidence (descriptive behavior): `GRHayL/EOS/Tabulated/NRPyEOS_enforce_table_bounds.c::NRPyEOS_enforce_table_bounds_rho_Ye_T`; corroboration: `GRHayL/EOS/Tabulated/NRPyEOS_tabulated_compute_enthalpy_and_cs2.c`.

## Sources

- [GRHayL/include/ghl_nrpyeos_tabulated.h](../../../GRHayL/include/ghl_nrpyeos_tabulated.h) - tabulated wrapper declarations
- [GRHayL/EOS/Tabulated/interpolators/NRPyEOS_tabulated_helpers.h](../../../GRHayL/EOS/Tabulated/interpolators/NRPyEOS_tabulated_helpers.h) - `NRPyEOS_checkbounds`, `NRPyEOS_findtemp_from_any`, `NRPyEOS_bisection`
- [GRHayL/EOS/Tabulated/interpolators/NRPyEOS_from_rho_Ye_T_interpolate_n_quantities.c](../../../GRHayL/EOS/Tabulated/interpolators/NRPyEOS_from_rho_Ye_T_interpolate_n_quantities.c) - direct `(rho,Y_e,T)` helper
- [GRHayL/EOS/Tabulated/interpolators/NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities.c](../../../GRHayL/EOS/Tabulated/interpolators/NRPyEOS_from_rho_Ye_aux_find_T_and_interpolate_n_quantities.c) - auxiliary temperature-recovery helper
- [GRHayL/EOS/Tabulated/NRPyEOS_enforce_table_bounds.c](../../../GRHayL/EOS/Tabulated/NRPyEOS_enforce_table_bounds.c) - in-place configured-bound clamps
- [GRHayL/EOS/Tabulated/NRPyEOS_tabulated_compute_enthalpy_and_cs2.c](../../../GRHayL/EOS/Tabulated/NRPyEOS_tabulated_compute_enthalpy_and_cs2.c) - tabulated `ghl_compute_h_and_cs2`
- [Unit_Tests/unit_test_tabulated_eos.c](../../../Unit_Tests/unit_test_tabulated_eos.c) - interpolator and bounds tests

## See Also

- Depends on: [Tabulated Table Layout](tabulated-table-layout.md)
- Implements: [Tabulated Interpolator Catalog](tabulated-interpolator-catalog.md)
- Validated by: [EOS Tests And Reference Data](tests-and-reference-data.md)
- See also: [CompOSE EOS Leakage Callback](../neutrinos/compose-eos-leakage-callback.md)
