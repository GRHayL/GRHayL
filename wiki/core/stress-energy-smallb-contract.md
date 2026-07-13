# Stress-Energy/smallb Contract

## Purpose

Route Core magnetic 4-vector and stress-energy helper contracts. This page
maps prerequisites, written components, and test routes for
`ghl_compute_smallb_and_b2`, `ghl_compute_TDNmunu`, `ghl_compute_TUPmunu`,
`ghl_initialize_stress_energy`, and `ghl_return_stress_energy`; source,
headers, Doxygen, and tests remain authority
(`GRHayL/GRHayL_Core/compute_smallb_and_b2.c`,
`GRHayL/GRHayL_Core/compute_TDNmunu.c`,
`GRHayL/GRHayL_Core/compute_TUPmunu.c`, `GRHayL/include/ghl.h`,
`docs/raw/GRHayL_Core.dox`).

## Prerequisites

Stress-energy compute paths require initialized primitive and metric state:
- `prims->u0` must be valid before `ghl_compute_smallb_and_b2`,
  `ghl_compute_TDNmunu`, or `ghl_compute_TUPmunu`; Doxygen comments and source
  use it directly (`ghl_compute_smallb_and_b2` in
  `GRHayL/GRHayL_Core/compute_smallb_and_b2.c`,
  `ghl_compute_TDNmunu` in `GRHayL/GRHayL_Core/compute_TDNmunu.c`, and
  `ghl_compute_TUPmunu` in `GRHayL/GRHayL_Core/compute_TUPmunu.c`).
- `prims->eps` must be valid for `ghl_compute_TDNmunu` and
  `ghl_compute_TUPmunu`, because both compute enthalpy from `eps`, pressure, and
  density (`ghl_compute_TDNmunu` in
  `GRHayL/GRHayL_Core/compute_TDNmunu.c` and `ghl_compute_TUPmunu` in
  `GRHayL/GRHayL_Core/compute_TUPmunu.c`).
- Initialized ADM metric is required for lapse, shift, spatial metric, and
  magnetic normalization in `ghl_compute_smallb_and_b2` (`ghl_metric_quantities`
  in `GRHayL/include/ghl.h` and `ghl_compute_smallb_and_b2` in
  `GRHayL/GRHayL_Core/compute_smallb_and_b2.c`).
- Initialized ADM auxiliaries are required for 4-metric and inverse 4-metric
  slots used by stress-energy routines (`ghl_ADM_aux_quantities` and stress-
  energy declarations in `GRHayL/include/ghl.h`, plus the compute routines in
  `GRHayL/GRHayL_Core/compute_TDNmunu.c` and
  `GRHayL/GRHayL_Core/compute_TUPmunu.c`).
- The `uDN` argument to `ghl_compute_smallb_and_b2` must be the lowered
  4-velocity consistent with the supplied metric and primitives. The helper
  trusts this caller-provided representation and does not derive or validate
  it (`GRHayL/include/ghl.h` and
  `GRHayL/GRHayL_Core/compute_smallb_and_b2.c`).
- Stress-energy enthalpy evaluation divides by `prims->rho`; callers must
  provide finite thermodynamic values and nonzero density. No Core stress
  routine checks these conditions (`GRHayL/GRHayL_Core/compute_TDNmunu.c` and
  `GRHayL/GRHayL_Core/compute_TUPmunu.c`).

## `smallb` and `b2`

`ghl_compute_smallb_and_b2` returns the magnetic 4-vector `smallb[0..3]`
(`b^mu`) plus scalar `smallb2` (`b^2`) through output arguments declared in
`GRHayL/include/ghl.h` and written by `ghl_compute_smallb_and_b2` in
`GRHayL/GRHayL_Core/compute_smallb_and_b2.c`.

This helper uses GRHayL's magnetic-field rescaling convention: GRHayL stores
magnetic quantities rescaled by `sqrt(4*pi)`, documented in
`docs/raw/derivation.md` and reflected by `ONE_OVER_SQRT_4PI` in
`GRHayL/include/ghl.h`. The helper comments in
`GRHayL/GRHayL_Core/compute_smallb_and_b2.c` also call out the same
normalization relative to the Duez et al. expressions.

The helper returns `void`, overwrites all four `smallb` entries and `smallb2`,
and reports no error for null, non-finite, inconsistent, or singular inputs.

Equation details belong in `docs/raw/derivation.md` as source context, especially
magnetic rescaling, `b^mu`, and ideal-MHD stress-energy background. Do not copy
derivations into this KB page.

## `TDNmunu` vs `TUPmunu`

`ghl_compute_TDNmunu` computes lower-index stress-energy components using the
lowered 4-velocity, lowered `smallb`, and `metric_aux->g4DD` in
`GRHayL/GRHayL_Core/compute_TDNmunu.c`. Use it when callers need the
`T_{\mu\nu}` route documented by the source comments.

`ghl_compute_TUPmunu` computes upper-index stress-energy components using the
upper 4-velocity, `smallb`, and `metric_aux->g4UU` in
`GRHayL/GRHayL_Core/compute_TUPmunu.c`. Use it when callers need the
`T^{\mu\nu}` route documented by the source comments.

Both routines populate the ten documented/returned symmetric components in
`Tmunu->T4`: `[0][0]`, `[0][1]`, `[0][2]`, `[0][3]`, `[1][1]`, `[1][2]`,
`[1][3]`, `[2][2]`, `[2][3]`, and `[3][3]` in
`GRHayL/GRHayL_Core/compute_TDNmunu.c` and
`GRHayL/GRHayL_Core/compute_TUPmunu.c`. They do not independently write every
symmetric matrix slot, so callers must not infer fresh writes to entries such as
`[1][0]` from these compute routines alone.

Both compute routines return `void` and perform no validation or finite-result
check. On entry, the six lower-triangle slots not named above may be
uninitialized; on return, they retain their prior contents.
`ghl_initialize_stress_energy` can build a mirrored matrix from ten supplied
values, but it is not part of either compute routine. Calling it before a
compute does not mirror the newly computed upper triangle; callers needing a
fully symmetric computed matrix must explicitly mirror those results or
unpack/reinitialize after computation.

## Packing and Return Helpers

`ghl_initialize_stress_energy` copies ten component inputs into a
`ghl_stress_energy` and mirrors off-diagonal entries into both symmetric slots
(`ghl_initialize_stress_energy` in `GRHayL/include/ghl.h` and
`GRHayL/GRHayL_Core/initialize_stress_energy.c`). This helper can initialize
either raised or lowered stress-energy data; the struct itself is just `T4[4][4]`
(`ghl_stress_energy` in `GRHayL/include/ghl.h`).

`ghl_return_stress_energy` returns the same ten documented components from
`Tmunu->T4`; it does not return every matrix slot (`ghl_return_stress_energy`
in `GRHayL/include/ghl.h` and
`GRHayL/GRHayL_Core/return_stress_energy.c`).

## Tests

`Unit_Tests/unit_test_compute_conservs_and_Tmunu.c` stores `u0` as explicit
input, initializes metric and ADM auxiliaries, assigns `prims.u0`, runs the
combined conservative/stress-energy path, and also calls standalone
`ghl_compute_TDNmunu` for comparison against fixture data in
`Unit_Tests/unit_test_compute_conservs_and_Tmunu.c`. This is normal-path
coverage for returned/stored stress-energy components, not proof that every
matrix slot is independently written.

`Unit_Tests/pert_test_fail_stress_energy.c` checks the same ten named
components used by the initialize/return contract (`Ttt`, `Ttx`, `Tty`, `Ttz`,
`Txx`, `Txy`, `Txz`, `Tyy`, `Tyz`, `Tzz`) and reports failures component by
component (`Unit_Tests/pert_test_fail_stress_energy.c`).

No direct fixture assertion for `ghl_compute_smallb_and_b2` or
`ghl_compute_TUPmunu` is visible in the listed Core test routes. Their use in
other computation paths is indirect coverage, not a standalone check of all
outputs or invalid-input behavior.
