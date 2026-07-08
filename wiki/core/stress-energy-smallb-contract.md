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
  use it directly (`GRHayL/GRHayL_Core/compute_smallb_and_b2.c:44`,
  `GRHayL/GRHayL_Core/compute_smallb_and_b2.c:61`,
  `GRHayL/GRHayL_Core/compute_TDNmunu.c:22`,
  `GRHayL/GRHayL_Core/compute_TUPmunu.c:22`).
- `prims->eps` must be valid for `ghl_compute_TDNmunu` and
  `ghl_compute_TUPmunu`, because both compute enthalpy from `eps`, pressure, and
  density (`GRHayL/GRHayL_Core/compute_TDNmunu.c:22`,
  `GRHayL/GRHayL_Core/compute_TDNmunu.c:35`,
  `GRHayL/GRHayL_Core/compute_TUPmunu.c:22`,
  `GRHayL/GRHayL_Core/compute_TUPmunu.c:35`).
- Initialized ADM metric is required for lapse, shift, spatial metric, and
  magnetic normalization in `ghl_compute_smallb_and_b2`
  (`GRHayL/include/ghl.h:154`,
  `GRHayL/GRHayL_Core/compute_smallb_and_b2.c:61`,
  `GRHayL/GRHayL_Core/compute_smallb_and_b2.c:96`).
- Initialized ADM auxiliaries are required for 4-metric and inverse 4-metric
  slots used by stress-energy routines (`GRHayL/include/ghl.h:167`,
  `GRHayL/include/ghl.h:502`,
  `GRHayL/GRHayL_Core/compute_TDNmunu.c:46`,
  `GRHayL/GRHayL_Core/compute_TUPmunu.c:46`).

## `smallb` and `b2`

`ghl_compute_smallb_and_b2` returns the magnetic 4-vector `smallb[0..3]`
(`b^mu`) plus scalar `smallb2` (`b^2`) through output arguments declared in
`GRHayL/include/ghl.h:514` and written in
`GRHayL/GRHayL_Core/compute_smallb_and_b2.c:69`,
`GRHayL/GRHayL_Core/compute_smallb_and_b2.c:84`, and
`GRHayL/GRHayL_Core/compute_smallb_and_b2.c:100`.

This helper uses GRHayL's magnetic-field rescaling convention: GRHayL stores
magnetic quantities rescaled by `sqrt(4*pi)`, documented in
`docs/raw/derivation.md:79` and reflected by `ONE_OVER_SQRT_4PI` in
`GRHayL/include/ghl.h:13`. The helper comments also call out the same
normalization relative to the Duez et al. expressions
(`GRHayL/GRHayL_Core/compute_smallb_and_b2.c:8`,
`GRHayL/GRHayL_Core/compute_smallb_and_b2.c:63`,
`GRHayL/GRHayL_Core/compute_smallb_and_b2.c:71`).

Equation details belong in `docs/raw/derivation.md` as source context, especially
magnetic rescaling, `b^mu`, and ideal-MHD stress-energy background
(`docs/raw/derivation.md:79`, `docs/raw/derivation.md:398`,
`docs/raw/derivation.md:489`). Do not copy derivations into this KB page.

## `TDNmunu` vs `TUPmunu`

`ghl_compute_TDNmunu` computes lower-index stress-energy components using the
lowered 4-velocity, lowered `smallb`, and `metric_aux->g4DD`
(`GRHayL/GRHayL_Core/compute_TDNmunu.c:29`,
`GRHayL/GRHayL_Core/compute_TDNmunu.c:46`,
`GRHayL/GRHayL_Core/compute_TDNmunu.c:59`,
`GRHayL/GRHayL_Core/compute_TDNmunu.c:65`). Use it when callers need the
`T_{\mu\nu}` route documented by the source comments.

`ghl_compute_TUPmunu` computes upper-index stress-energy components using the
upper 4-velocity, `smallb`, and `metric_aux->g4UU`
(`GRHayL/GRHayL_Core/compute_TUPmunu.c:29`,
`GRHayL/GRHayL_Core/compute_TUPmunu.c:39`,
`GRHayL/GRHayL_Core/compute_TUPmunu.c:52`,
`GRHayL/GRHayL_Core/compute_TUPmunu.c:62`). Use it when callers need the
`T^{\mu\nu}` route documented by the source comments.

Both routines populate the ten documented/returned symmetric components in
`Tmunu->T4`: `[0][0]`, `[0][1]`, `[0][2]`, `[0][3]`, `[1][1]`, `[1][2]`,
`[1][3]`, `[2][2]`, `[2][3]`, and `[3][3]`
(`GRHayL/GRHayL_Core/compute_TDNmunu.c:65`,
`GRHayL/GRHayL_Core/compute_TUPmunu.c:62`). They do not independently write
every symmetric matrix slot, so callers must not infer fresh writes to entries
such as `[1][0]` from these compute routines alone.

## Packing and Return Helpers

`ghl_initialize_stress_energy` copies ten component inputs into a
`ghl_stress_energy` and mirrors off-diagonal entries into both symmetric slots
(`GRHayL/include/ghl.h:470`,
`GRHayL/GRHayL_Core/initialize_stress_energy.c:20`,
`GRHayL/GRHayL_Core/initialize_stress_energy.c:27`). This helper can initialize
either raised or lowered stress-energy data; the struct itself is just `T4[4][4]`
(`GRHayL/include/ghl.h:181`,
`GRHayL/GRHayL_Core/initialize_stress_energy.c:7`).

`ghl_return_stress_energy` returns the same ten documented components from
`Tmunu->T4`; it does not return every matrix slot
(`GRHayL/include/ghl.h:483`,
`GRHayL/GRHayL_Core/return_stress_energy.c:35`,
`GRHayL/GRHayL_Core/return_stress_energy.c:42`).

## Tests

`Unit_Tests/unit_test_compute_conservs_and_Tmunu.c` stores `u0` as explicit
input, initializes metric and ADM auxiliaries, assigns `prims.u0`, runs the
combined conservative/stress-energy path, and also calls standalone
`ghl_compute_TDNmunu` for comparison against fixture data
(`Unit_Tests/unit_test_compute_conservs_and_Tmunu.c:93`,
`Unit_Tests/unit_test_compute_conservs_and_Tmunu.c:217`,
`Unit_Tests/unit_test_compute_conservs_and_Tmunu.c:226`,
`Unit_Tests/unit_test_compute_conservs_and_Tmunu.c:229`,
`Unit_Tests/unit_test_compute_conservs_and_Tmunu.c:272`). This is normal-path
coverage for returned/stored stress-energy components, not proof that every
matrix slot is independently written.

`Unit_Tests/pert_test_fail_stress_energy.c` checks the same ten named
components used by the initialize/return contract (`Ttt`, `Ttx`, `Tty`, `Ttz`,
`Txx`, `Txy`, `Txz`, `Tyy`, `Tyz`, `Tzz`) and reports failures component by
component (`Unit_Tests/pert_test_fail_stress_energy.c:10`,
`Unit_Tests/pert_test_fail_stress_energy.c:16`,
`Unit_Tests/pert_test_fail_stress_energy.c:70`).
