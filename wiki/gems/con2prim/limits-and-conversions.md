# Con2Prim Limits And Conversions

## Purpose

This page routes helper contracts for Con2Prim limiting, primitive guesses,
velocity limiting, and primitive-to-conservative conversion. Source, headers,
Doxygen, and tests remain the authority; this page avoids copying derivations.

Read with [wiki/gems/con2prim.md](../con2prim.md), [docs/raw/Con2Prim.dox](../../../docs/raw/Con2Prim.dox), and
[GRHayL/include/ghl_con2prim.h](../../../GRHayL/include/ghl_con2prim.h).

## Recovery-Side Helpers

These helpers prepare or finish conservative-to-primitive recovery. They are not
primitive-to-conservative conversion routines.

### `ghl_apply_conservative_limits`

Route: [GRHayL/Con2Prim/apply_conservative_limits.c](../../../GRHayL/Con2Prim/apply_conservative_limits.c).

Contract:
- Inputs: `ghl_parameters`, `ghl_eos_parameters`, ADM metric, primitives, conservatives, and diagnostics.
- Required primitive field before recovery: `BU`. The source comments state this is the only primitive expected before the Con2Prim solve.
- Required conservative fields: densitized `rho`, `tau`, and `SD`.
- Required EOS bounds/atmosphere fields: `tau_atm` and `press_atm`.
- For tabulated EOS recovery, a negative energy or `tau` is not by itself
  invalid. Admissibility depends on the table's local bounds at the state's
  density and composition. Solver-specific sign clipping must not be inferred
  from the generic atmosphere-limit route.
- Required parameter: `psi6threshold`, used with `metric_adm->sqrt_detgamma` to choose the high-`psi6` momentum/energy limiting branch.
- Writes conservative outputs in place: limited `tau` and possibly rescaled `SD`.
- Writes diagnostics: `tau_fix` when the atmosphere floor, magnetic-energy, or
  high-`psi6` correction changes `tau`, and `Stilde_fix` when `SD` is rescaled.
- Caller should initialize diagnostics first with `ghl_initialize_diagnostics`; the helper only sets flags true when a fix occurs.

Tests and fixtures:
- Direct test: [Unit_Tests/unit_test_apply_conservative_limits.c](../../../Unit_Tests/unit_test_apply_conservative_limits.c).
- Shared fixture generator: [Unit_Tests/data_gen/unit_test_data_con2prim_multi_method_hybrid.c](../../../Unit_Tests/data_gen/unit_test_data_con2prim_multi_method_hybrid.c).
- Legacy recovery comparison route: [Unit_Tests/unit_test_ET_Legacy_primitives.c](../../../Unit_Tests/unit_test_ET_Legacy_primitives.c).

### `ghl_guess_primitives`

Route: [GRHayL/Con2Prim/guess_primitives.c](../../../GRHayL/Con2Prim/guess_primitives.c).

For the tabulated neural-network fallback path, route to
[neural-network primitive guess](neural-network-primitive-guess.md). That page
owns `ghl_c2p_nn_guess_primitives`, `c2p_nn`, and HDF5 model routing; this
section only covers the ordinary `ghl_guess_primitives` helper.

Contract:
- Inputs: `ghl_parameters`, `ghl_eos_parameters`, ADM metric, undensitized conservatives, and primitives.
- Required conservative fields: undensitized `rho`, `tau`, `SD`, and `Y_e` where the tabulated path needs electron fraction.
- Required primitive field for tabulated guess: `BU`, because the Palenzuela-style estimate computes magnetic contractions before velocity construction.
- Simple/hybrid EOS path sets `rho` from undensitized conservative density, sets `u0 = 1`, sets `vU = -betaU`, and computes cold `press` and `eps`.
- That simple/hybrid path is a partial update: incoming `BU`, entropy, `Y_e`,
  and temperature remain unchanged.
- Tabulated EOS path computes metric/magnetic contractions through `ghl_compute_SU_Bsq_Ssq_BdotS`, uses `T_max` as the temperature guess, enforces table bounds on `rho`, `Y_e`, and `eps`, computes `press`, `entropy`, and `temperature` from `eps`, then calls `ghl_limit_utilde_and_compute_v`.
- Invalid or nonfinite algebraic inputs produce an atmosphere guess that retains
  finite conservative composition at positive density within EOS bounds. An EOS
  inversion error instead preserves the bounded `rho`, `Y_e`, and `eps` seed,
  keeps `T_max`, supplies finite atmosphere pressure and entropy, and continues
  velocity construction. A later nonfinite or singular derived quantity still
  replaces that partial seed with the same composition-preserving atmosphere guess.
- Required tabulated EOS bounds: table-backed `rho`, `Y_e`, and `eps` bounds exposed through `ghl_tabulated_enforce_bounds_rho_Ye_eps`; required temperature bound/guess: `T_max`.
- Required parameter for tabulated guess: `max_Lorentz_factor` through `ghl_limit_utilde_and_compute_v`.
- Diagnostics: no diagnostics pointer is accepted. The tabulated path invokes utilde limiting but does not store `speed_limited`.

Tests and fixtures:
- Fixture route: `Unit_Tests/data_gen/unit_test_data_con2prim_multi_method_hybrid.c`.
- Recovery test consuming these fixtures: [Unit_Tests/unit_test_con2prim_multi_method_hybrid.c](../../../Unit_Tests/unit_test_con2prim_multi_method_hybrid.c). [Unit_Tests/unit_test_hybrid_failure.c](../../../Unit_Tests/unit_test_hybrid_failure.c) uses hard-coded cases and reads no fixture.

### `ghl_enforce_primitive_limits_and_compute_u0`

Route: [GRHayL/Con2Prim/enforce_primitive_limits_and_compute_u0.c](../../../GRHayL/Con2Prim/enforce_primitive_limits_and_compute_u0.c).

Contract:
- Inputs: `ghl_parameters`, `ghl_eos_parameters`, ADM metric, primitives, and an initialized `speed_limited` in/out OR accumulator. Initialize the flag to false at the start of a logical group; subsequent calls preserve an accumulated true value.
- Hybrid EOS: clamps `rho` to EOS bounds, recomputes cold pressure/`eps`, floors `press` at cold pressure, and applies a pressure ceiling only when bounded `rho < 100*rho_atm` or `sqrt_detgamma > psi6threshold`. The ceiling is `1e5*P_cold` in the high-`psi6` branch and `100*P_cold` otherwise; higher-density states outside that branch have no pressure ceiling. It recomputes `eps`, and recomputes `entropy` when `evolve_entropy` is true.
- Simple EOS: clamps `rho` and `press` to EOS bounds, recomputes `eps`, and recomputes `entropy` when `evolve_entropy` is true.
- Tabulated EOS: when HDF5 support is enabled, enforces `rho`, `Y_e`, and `temperature` table bounds, then recomputes `press` and `eps`; it also recomputes `entropy` when `evolve_entropy` is true. When HDF5 is disabled, returns `ghl_error_used_disabled_hdf5`.
- The helper does not branch on `params->evolve_temp`; it uses `temperature` as the tabulated EOS thermodynamic input.
- After successful thermodynamic processing: calls `ghl_limit_v_and_compute_u0`, enforcing `max_Lorentz_factor`, computing `u0`, and accumulating `speed_limited`.
- Possible error routes include unknown EOS type, disabled HDF5 for tabulated EOS, propagated tabulated EOS errors, and `u0` singularity from the velocity limiter. EOS errors return before velocity limiting. Earlier primitive changes are not rolled back; failure does not establish a valid updated `u0`.

Tests and fixtures:
- Direct test: [Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c](../../../Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c).
- Shared fixture generator: `Unit_Tests/data_gen/unit_test_data_con2prim_multi_method_hybrid.c`.
- Error-code route for `u0` limiting: [Unit_Tests/unit_test_code_error.c](../../../Unit_Tests/unit_test_code_error.c).

### `ghl_limit_utilde_and_compute_v`

Route: [GRHayL/Con2Prim/limit_utilde_and_compute_v.c](../../../GRHayL/Con2Prim/limit_utilde_and_compute_v.c).

Contract:
- Internal helper declared in `GRHayL/include/ghl_con2prim.h`.
- Inputs: `ghl_parameters`, ADM metric, mutable spatial `utU[3]`, and primitives.
- Caller must supply a valid ADM metric/lapse/shift, finite utilde with usable norm/scaling arithmetic and a nonnegative metric norm, and finite `max_Lorentz_factor >= 1`. The helper performs no finite/domain validation and has no error return.
- Computes the utilde norm using `metric_adm->gammaDD`. The roundoff-slack
  trigger rescales `utU` only when the computed factor is below one, so the
  helper never increases the supplied speed for usable inputs.
- Writes primitive outputs: `u0` and `vU`.
- Returns `true` when rescaling occurred, otherwise `false`; false does not certify finite valid outputs.
- Used by the tabulated primitive guess path; unlike `ghl_enforce_primitive_limits_and_compute_u0`, it does not write a diagnostics struct by itself.

## Conversion Helpers

These helpers convert initialized primitives to conservative variables. They do
not recover primitives from conservatives and should be kept separate from C2P
solver flow.

### `ghl_compute_conservs`

Route: [GRHayL/Con2Prim/compute_conservs.c](../../../GRHayL/Con2Prim/compute_conservs.c).

Contract:
- Inputs: ADM metric, ADM auxiliary metric, initialized primitives, and conservative output.
- Required primitive fields include `rho`, `press`, `eps`, `u0`, `vU`, `BU`, `entropy`, and `Y_e`.
- `eps` and `u0` are explicit prerequisites in source comments. `u0` is normally produced by `ghl_enforce_primitive_limits_and_compute_u0` or `ghl_limit_v_and_compute_u0`.
- Uses `metric_aux->g4DD` and core magnetic helper `ghl_compute_smallb_and_b2`.
- Outputs densitized conservative groups: density `rho`, energy `tau`, momentum `SD`, entropy, and `Y_e`.
- Recommended after primitive limits so primitive and conservative states stay consistent.

Tests and fixtures:
- Direct test coverage is in [Unit_Tests/unit_test_compute_conservs_and_Tmunu.c](../../../Unit_Tests/unit_test_compute_conservs_and_Tmunu.c).
- Fixture generator route: `Unit_Tests/data_gen/unit_test_data_con2prim_multi_method_hybrid.c`.
- Repo-local `ET_Legacy` conservative fixture route: [Unit_Tests/data_gen/unit_test_data_ET_Legacy_conservs.c](../../../Unit_Tests/data_gen/unit_test_data_ET_Legacy_conservs.c) and [Unit_Tests/unit_test_ET_Legacy_conservs.c](../../../Unit_Tests/unit_test_ET_Legacy_conservs.c).

### `ghl_compute_conservs_and_Tmunu`

Route: [GRHayL/Con2Prim/compute_conservs_and_Tmunu.c](../../../GRHayL/Con2Prim/compute_conservs_and_Tmunu.c).

Contract:
- Inputs match `ghl_compute_conservs`, plus a `ghl_stress_energy` output.
- Required primitive fields include `eps` and `u0`; the test explicitly stores `u0` because the helper consumes it.
- Outputs the same conservative groups as `ghl_compute_conservs`.
- Also outputs the lower-index stress-energy tensor route `Tmunu->T4` for the ten stored components returned by `ghl_return_stress_energy`.
- Reuses intermediate four-velocity and magnetic quantities rather than calling the standalone core stress-energy helper.
- Use this when conservatives and `Tmunu` are needed from the same initialized primitive state.

Tests and fixtures:
- Direct test: `Unit_Tests/unit_test_compute_conservs_and_Tmunu.c`.
- Fixture generator route: `Unit_Tests/data_gen/unit_test_data_con2prim_multi_method_hybrid.c`.
- Repo-local `ET_Legacy` routes: `Unit_Tests/data_gen/unit_test_data_ET_Legacy_conservs.c`, [Unit_Tests/data_gen/unit_test_data_ET_Legacy_primitives.c](../../../Unit_Tests/data_gen/unit_test_data_ET_Legacy_primitives.c), `Unit_Tests/unit_test_ET_Legacy_conservs.c`, and `Unit_Tests/unit_test_ET_Legacy_primitives.c`.

### Internal Contraction Route

`ghl_compute_SU_Bsq_Ssq_BdotS` lives in
[GRHayL/Con2Prim/compute_SU_Bsq_Ssq_BdotS.c](../../../GRHayL/Con2Prim/compute_SU_Bsq_Ssq_BdotS.c). It is a shared contraction helper
for tabulated primitive guesses and Palenzuela-style quantities:
- Inputs: ADM metric, undensitized conservatives, and primitives with `BU`.
- Outputs: raised `SU`, `Bsq`, `Ssq`, and `BdotS`.
- It may locally rescale a copy of `SD` before computing contractions, but it does not update the caller's conservative struct.

The tabulated NN path reuses the same aux quantities (`q`, `r`, `s`, `t`) and
finishes through `ghl_tabulated_primitive_guess_from_x`; keep detailed NN
schema and retry policy on
[neural-network primitive guess](neural-network-primitive-guess.md).

## Direct Ground Truth

- Public declarations and diagnostics fields: `GRHayL/include/ghl_con2prim.h`.
- Con2Prim parameters: [GRHayL/include/ghl.h](../../../GRHayL/include/ghl.h) and [GRHayL/GRHayL_Core/initialize_params.c](../../../GRHayL/GRHayL_Core/initialize_params.c).
- C2P docs: `docs/raw/Con2Prim.dox`.
- Source equations and step details: the source files named in each section above.
