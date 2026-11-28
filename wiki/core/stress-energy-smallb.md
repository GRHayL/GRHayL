# Stress-Energy And `smallb`

> Core comoving magnetic 4-vector (`smallb`, `b^2`), covariant/contravariant stress-energy (`T_{\mu\nu}`, `T^{\mu\nu}`), upper-triangle write rules, and packing helpers. · Status: confirmed
> Up: [GRHayL Core](index.md)

## Summary

Core provides `ghl_compute_smallb_and_b2`, `ghl_compute_TDNmunu`, `ghl_compute_TUPmunu`, `ghl_initialize_stress_energy`, and `ghl_return_stress_energy`. All return `void` without input validation. `ghl_compute_smallb_and_b2` assumes rescaled `BU` (`B / sqrt(4*pi)`) and valid `u0`; the two `Tmunu` compute routines also need nonzero `rho` and initialized ADM auxiliaries and write only the ten upper-triangle components of `T4[4][4]`; the packing helpers only copy components.

## Detail

### Prerequisites And `smallb` / `b2`

`ghl_compute_smallb_and_b2` needs valid `prims->u0`, `prims->vU`, rescaled `prims->BU`, lowered `uDN`, and the lapse, shift, and `gammaDD` in `ghl_metric_quantities`; it reads no `ghl_ADM_aux_quantities`, `rho`, or `eps`. `ghl_compute_TDNmunu` and `ghl_compute_TUPmunu` additionally need finite `prims->press` and `prims->eps`, nonzero `prims->rho` (enthalpy divides by `rho`), and `ghl_ADM_aux_quantities` (`g4DD` for `TDNmunu`; `g4DD` and `g4UU` for `TUPmunu`). `ghl_initialize_stress_energy` and `ghl_return_stress_energy` only copy components.

`ghl_compute_smallb_and_b2` (`compute_smallb_and_b2.c`) takes caller-supplied lowered four-velocity `uDN` (trusted without validation) and already-rescaled `prims->BU` (divided by `sqrt(4*pi)`, matching `ONE_OVER_SQRT_4PI` and `docs/raw/derivation.md`), overwriting all four entries of `smallb[0..3]` (`b^\mu`) and scalar `smallb2` (`b^2`).

Evidence (descriptive behavior): `compute_smallb_and_b2.c::ghl_compute_smallb_and_b2`; corroboration: `derivation.md`.

### `TDNmunu`, `TUPmunu`, Packing, And Coverage Gaps

- **`ghl_compute_TDNmunu` (`compute_TDNmunu.c`):** computes lower-index `T_{\mu\nu}` from lowered `uDN`, lowered `smallb`, and `metric_aux->g4DD`.
- **`ghl_compute_TUPmunu` (`compute_TUPmunu.c`):** computes upper-index `T^{\mu\nu}` from upper `uUP`, upper `smallb`, and `metric_aux->g4UU`.
- **Upper-triangle-only writes:** both compute routines populate only the ten entries `[0][0]`, `[0][1]`, `[0][2]`, `[0][3]`, `[1][1]`, `[1][2]`, `[1][3]`, `[2][2]`, `[2][3]`, `[3][3]` of `Tmunu->T4[4][4]`. The six lower-triangle entries retain prior/uninitialized contents; calling `ghl_initialize_stress_energy` beforehand does **not** mirror newly computed upper-triangle values afterward.
- **Packing/unpacking:** `ghl_initialize_stress_energy` (`initialize_stress_energy.c`) copies ten symmetric components into `T4[4][4]` and mirrors off-diagonals; `ghl_return_stress_energy` (`return_stress_energy.c`) reads the same ten upper-triangle components back.
- **Test coverage:** `Unit_Tests/unit_test_compute_conservs_and_Tmunu.c` (with `pert_test_fail_stress_energy.c`) checks combined conservative/stress-energy output and standalone `ghl_compute_TDNmunu`, exercising `ghl_compute_smallb_and_b2` indirectly. Neither `ghl_compute_TUPmunu` nor standalone `ghl_compute_smallb_and_b2` has a direct reference-data assertion.

Evidence (descriptive behavior): `compute_TDNmunu.c::ghl_compute_TDNmunu`; corroboration: `compute_TUPmunu.c::ghl_compute_TUPmunu` and `unit_test_compute_conservs_and_Tmunu.c`.

## Sources

- [GRHayL/include/ghl.h](../../GRHayL/include/ghl.h) - `ghl_stress_energy` and stress-energy declarations
- [GRHayL/GRHayL_Core/compute_smallb_and_b2.c](../../GRHayL/GRHayL_Core/compute_smallb_and_b2.c) - `ghl_compute_smallb_and_b2`
- [GRHayL/GRHayL_Core/compute_TDNmunu.c](../../GRHayL/GRHayL_Core/compute_TDNmunu.c) - covariant `T_{\mu\nu}` upper-triangle computation
- [GRHayL/GRHayL_Core/compute_TUPmunu.c](../../GRHayL/GRHayL_Core/compute_TUPmunu.c) - contravariant `T^{\mu\nu}` upper-triangle computation
- [GRHayL/GRHayL_Core/initialize_stress_energy.c](../../GRHayL/GRHayL_Core/initialize_stress_energy.c) - mirrored 10-component packer
- [GRHayL/GRHayL_Core/return_stress_energy.c](../../GRHayL/GRHayL_Core/return_stress_energy.c) - 10-component unpacker
- [docs/raw/derivation.md](../../docs/raw/derivation.md) - magnetic rescaling and GRMHD stress-energy equations
- [Unit_Tests/unit_test_compute_conservs_and_Tmunu.c](../../Unit_Tests/unit_test_compute_conservs_and_Tmunu.c) - `ghl_compute_TDNmunu` validation
- [Unit_Tests/pert_test_fail_stress_energy.c](../../Unit_Tests/pert_test_fail_stress_energy.c) - 10-component perturbation comparator

## See Also

- Depends on: [Metric And ADM Helpers](metric-adm.md)
- Depends on: [Variables And Conventions](../physics/variables-and-conventions.md)
- Validated by: [Core Tests And Reference Data](tests-and-reference-data.md)
- See also: [Evolution Equation Map](../physics/evolution-equation-map.md)
