# Con2Prim Limits And Conversions

> Pre-recovery conservative limits, ordinary primitive guesses, post-recovery primitive limits/`u0`, `utilde` limiting, contractions, and primitive-to-conservative conversion. · Status: confirmed
> Up: [Con2Prim](index.md)

## Summary

Con2Prim separates recovery-side helpers (`ghl_apply_conservative_limits`, `ghl_guess_primitives`, `ghl_enforce_primitive_limits_and_compute_u0`, `ghl_limit_utilde_and_compute_v`, `ghl_compute_SU_Bsq_Ssq_BdotS`) from forward primitive-to-conservative routines (`ghl_compute_conservs`, `ghl_compute_conservs_and_Tmunu`). Callers must match densitized versus undensitized conservative inputs and check error returns before consuming updated primitives.

## Detail

### Pre-Recovery Limits And Ordinary Primitive Guesses

- **`ghl_apply_conservative_limits` (`apply_conservative_limits.c`):** operates in place on densitized `tau` and `SD`, requiring only `prims->BU` from primitives, `eos->tau_atm`/`press_atm`, and `params->psi6threshold` (the high-`psi6` branch runs only when `B2 >= eos->press_atm*1e-32` and `sqrt_detgamma > psi6threshold`; at smaller `B2` the weak-field momentum limit applies at any `sqrt_detgamma`). Sets `diagnostics->tau_fix` on atmosphere-floor, magnetic-energy, or high-`psi6` `tau` adjustments and `Stilde_fix` on `SD` rescaling. For tabulated EOS, negative energy or `tau` is not automatically invalid; admissibility depends on local table bounds at `(rho, Y_e)`.
- **`ghl_guess_primitives` (`guess_primitives.c`):** expects **undensitized** `cons_undens` (and `prims->BU` for tabulated EOS); takes no diagnostics pointer:
  - *Simple/hybrid:* sets `rho` from undensitized density, `u0 = 1`, `vU = -betaU`, and cold `press`/`eps`; leaves `BU`, `entropy`, `Y_e`, and `temperature` untouched.
  - *Tabulated:* calls `ghl_compute_SU_Bsq_Ssq_BdotS`, seeds temperature at `T_max`, enforces `(rho, Y_e, eps)` bounds via `ghl_tabulated_enforce_bounds_rho_Ye_eps`, inverts `eps` for `(press, entropy, temperature)`, and runs `ghl_limit_utilde_and_compute_v`. Invalid/nonfinite algebraic inputs fall back to the atmosphere guess of `set_tabulated_atmosphere_guess` (shared with `ghl_c2p_nn_guess_primitives`), whose electron fraction is `cons_undens->Y_e / cons_undens->rho` (conserved `D Y_e` divided by `D`) clamped to `[Y_e_min, Y_e_max]` when `D` is finite and positive and the quotient is finite, and `Y_e_atm` otherwise; an EOS-inversion failure preserves bounded `(rho, Y_e, eps)`, keeps `T_max`, substitutes finite atmosphere `press`/`entropy`, and completes velocities unless a later quantity is nonfinite.

Evidence (descriptive behavior): `apply_conservative_limits.c::ghl_apply_conservative_limits`; corroboration: `guess_primitives.c::ghl_guess_primitives` and `tabulated_primitive_guess_helpers.c::set_tabulated_atmosphere_guess`.

### Post-Recovery Limits, `utilde` Limiter, And Forward Conversions

- **`ghl_enforce_primitive_limits_and_compute_u0` (`enforce_primitive_limits_and_compute_u0.c`):** takes initialized in/out OR accumulator `*speed_limited`:
  - *Hybrid:* clamps `rho`, floors `press >= P_cold`, applies a pressure ceiling only when bounded `rho < 100 * rho_atm` or `sqrt_detgamma > psi6threshold` (`1e5 * P_cold` in high-`psi6` branch, `100 * P_cold` otherwise; no ceiling at higher density outside high-`psi6`), recomputes `eps`, and recomputes `entropy` when `evolve_entropy` is true.
  - *Simple:* clamps `rho` and `press`, recomputes `eps` (and `entropy` if `evolve_entropy`).
  - *Tabulated:* ignores `params->evolve_temp`; clamps `(rho, Y_e, temperature)` and recomputes `press`, `eps` (and `entropy` if `evolve_entropy`), or returns `ghl_error_used_disabled_hdf5` under `GHL_DISABLE_HDF5`.
  - Runs `ghl_limit_v_and_compute_u0` only after thermodynamic updates succeed; errors do not roll back earlier field edits.
- **`ghl_limit_utilde_and_compute_v` (`limit_utilde_and_compute_v.c`):** unchecked internal helper taking valid ADM metric, finite `utU[3]`, and finite `max_Lorentz_factor >= 1`; rescales `utU` only when the factor is `< 1` (never increasing speed), writes `prims->u0` and `prims->vU`, and returns whether rescaling occurred (without certifying finiteness).
- **Contractions & forward conversions:** `ghl_compute_SU_Bsq_Ssq_BdotS` (`compute_SU_Bsq_Ssq_BdotS.c`, declared in `ghl.h`) computes raised `SU`, `Bsq`, `Ssq`, `BdotS` from undensitized conservatives and `BU`, rescaling only a local copy of `SD`. `ghl_compute_conservs` (`compute_conservs.c`) and `ghl_compute_conservs_and_Tmunu` (`compute_conservs_and_Tmunu.c`, writing the ten upper-triangle `T4` components inline without calling `ghl_compute_TDNmunu`) require valid `prims->eps` and `prims->u0` and emit densitized conservatives.

Evidence (descriptive behavior): `enforce_primitive_limits_and_compute_u0.c::ghl_enforce_primitive_limits_and_compute_u0`; corroboration: `limit_utilde_and_compute_v.c` and `compute_SU_Bsq_Ssq_BdotS.c`.

## Sources

- [GRHayL/include/ghl_con2prim.h](../../../GRHayL/include/ghl_con2prim.h) - limit and conversion declarations
- [GRHayL/Con2Prim/apply_conservative_limits.c](../../../GRHayL/Con2Prim/apply_conservative_limits.c) - conservative floor/momentum limits
- [GRHayL/Con2Prim/guess_primitives.c](../../../GRHayL/Con2Prim/guess_primitives.c) - ordinary simple/hybrid/tabulated guesses
- [GRHayL/Con2Prim/Tabulated/tabulated_primitive_guess_helpers.c](../../../GRHayL/Con2Prim/Tabulated/tabulated_primitive_guess_helpers.c) - shared tabulated guess completion and atmosphere fallback
- [GRHayL/Con2Prim/enforce_primitive_limits_and_compute_u0.c](../../../GRHayL/Con2Prim/enforce_primitive_limits_and_compute_u0.c) - post-recovery limits and `u0`
- [GRHayL/Con2Prim/limit_utilde_and_compute_v.c](../../../GRHayL/Con2Prim/limit_utilde_and_compute_v.c) - `utilde` speed limiting
- [GRHayL/Con2Prim/compute_SU_Bsq_Ssq_BdotS.c](../../../GRHayL/Con2Prim/compute_SU_Bsq_Ssq_BdotS.c) - shared magnetic/momentum contractions
- [GRHayL/Con2Prim/compute_conservs.c](../../../GRHayL/Con2Prim/compute_conservs.c) - forward conservative conversion
- [GRHayL/Con2Prim/compute_conservs_and_Tmunu.c](../../../GRHayL/Con2Prim/compute_conservs_and_Tmunu.c) - combined conservative + `T_{\mu\nu}` conversion

## See Also

- Depends on: [Velocity And `u0`](../../core/velocity-u0.md)
- Depends on: [Stress-Energy And `smallb`](../../core/stress-energy-smallb.md)
- See also: [Recovery Flow](recovery-flow.md)
- See also: [Neural-Network Primitive Guess](neural-network-primitive-guess.md)
