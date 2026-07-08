# Atmosphere Prescription Contract

## Purpose

This page records the active Atmosphere prescription contract and the known
source drift around the unused radial-falloff routine. Source, Doxygen, tests,
and build lists remain the authority; this page routes readers to those files
without copying full APIs or implementation blocks.

Read with [wiki/gems/atmosphere.md](../atmosphere.md) for the module router.

## Active Prescription

`ghl_set_prims_to_constant_atm` is the only Atmosphere prescription that is both
built and documented.

- Public declaration:
  [GRHayL/include/ghl_atmosphere.h](../../../GRHayL/include/ghl_atmosphere.h).
- Built source:
  [GRHayL/Atmosphere/set_prims_to_constant_atm.c](../../../GRHayL/Atmosphere/set_prims_to_constant_atm.c).
- Build list:
  [GRHayL/Atmosphere/make.code.defn](../../../GRHayL/Atmosphere/make.code.defn).
- Doxygen source:
  [docs/raw/Atmosphere.dox](../../../docs/raw/Atmosphere.dox).

The constant prescription resets thermodynamic and velocity primitive fields
from `ghl_eos_parameters`:

- `rho = eos->rho_atm`
- `press = eos->press_atm`
- `eps = eos->eps_atm`
- `entropy = eos->entropy_atm`
- `Y_e = eos->Y_e_atm`
- `temperature = eos->T_atm`
- `vU[0..2] = 0`

The primitive struct also contains `u0` and `BU[0..2]` in
[GRHayL/include/ghl.h](../../../GRHayL/include/ghl.h). The constant Atmosphere
source does not assign either field. Callers that need a complete primitive
state must preserve magnetic fields separately and compute or refresh `u0`
through the primitive-limit path.

## Unsupported Radial Falloff

`ghl_set_prims_to_radial_falloff_atm` is a source-present drift candidate, not a
supported Atmosphere prescription.

- It is declared in
  [GRHayL/include/ghl_atmosphere.h](../../../GRHayL/include/ghl_atmosphere.h).
- Its source file exists at
  [GRHayL/Atmosphere/set_prims_to_radial_falloff_atm.c](../../../GRHayL/Atmosphere/set_prims_to_radial_falloff_atm.c).
- It is absent from
  [GRHayL/Atmosphere/make.code.defn](../../../GRHayL/Atmosphere/make.code.defn),
  so the Atmosphere gem build list does not compile it.
- It is absent from
  [docs/raw/Atmosphere.dox](../../../docs/raw/Atmosphere.dox), which says only
  one prescription is provided and links only the constant-density routine.

The radial-falloff source leaves the core atmosphere values unresolved: the
assignments for `rho`, `press`, and `eps` are commented or unset, while
`entropy`, `Y_e`, `temperature`, and `vU[0..2]` are assigned from the same
atmosphere-style EOS fields used by the constant prescription. Treat this file
as unfinished until declaration, build list, documentation, implementation, and
tests agree.

## EOS Field Coupling

Atmosphere values live in `ghl_eos_parameters`, declared in
[GRHayL/include/ghl.h](../../../GRHayL/include/ghl.h). EOS initialization fills
the atmosphere and floor fields that the active Atmosphere prescription later
copies into primitives.

- Simple EOS: `ghl_initialize_simple_eos` accepts `rho_atm` and `press_atm`,
  stores them on `ghl_eos_parameters`, then derives `eps_atm`,
  `entropy_atm`, and `tau_atm` from the one-piece ideal-fluid setup in
  [GRHayL/GRHayL_Core/initialize_eos.c](../../../GRHayL/GRHayL_Core/initialize_eos.c).
- Hybrid EOS: `ghl_initialize_hybrid_eos` accepts `rho_atm`, stores it, then
  computes `press_atm`, `eps_atm`, `entropy_atm`, and `tau_atm` through the
  hybrid cold-pressure/entropy helpers in
  [GRHayL/GRHayL_Core/initialize_eos.c](../../../GRHayL/GRHayL_Core/initialize_eos.c).
- Tabulated EOS: `ghl_initialize_tabulated_eos` accepts `rho_atm`, `Y_e_atm`,
  and `T_atm`, clamps min/max bounds against table bounds, stores the
  atmosphere values on `ghl_eos_parameters`, and computes `press_atm`,
  `eps_atm`, and `entropy_atm` from the table-backed `(rho,Y_e,T)` route in
  [GRHayL/GRHayL_Core/initialize_eos.c](../../../GRHayL/GRHayL_Core/initialize_eos.c).

For wrapper-level EOS setup, the public
`ghl_initialize_*_eos_functions_and_params` functions in
[GRHayL/include/ghl.h](../../../GRHayL/include/ghl.h) route into the same
initialization source.

## Tests And Coverage

Direct constant-atmosphere coverage lives in
[Unit_Tests/unit_test_grhayl_core_test_suite.c](../../../Unit_Tests/unit_test_grhayl_core_test_suite.c).
That test poisons primitive fields, calls `ghl_set_prims_to_constant_atm`, and
checks the simple and hybrid EOS cases.

Do not claim direct tabulated Atmosphere test coverage from that file. Its
tabulated setup branch is commented out, includes a TODO about adding
table-based default checks, and the Atmosphere loop currently stops before the
tabulated case.

[Unit_Tests/unit_test_tabulated_eos.c](../../../Unit_Tests/unit_test_tabulated_eos.c)
covers tabulated EOS initialization and table behavior, not a direct Atmosphere
reset call.

## Impact Links Only

Primitive-limit and recovery pages are downstream impact routes, not extra
Atmosphere prescription contracts:

- [wiki/gems/con2prim/limits-and-conversions.md](../con2prim/limits-and-conversions.md)
  documents `ghl_enforce_primitive_limits_and_compute_u0`, including the final
  `u0` computation path and EOS-bound enforcement after primitives already
  exist.
- [wiki/gems/con2prim/recovery-flow.md](../con2prim/recovery-flow.md) routes
  recovery order and shows where primitive limiting happens after recovery.

## Source Of Truth

- [AGENTS.md](../../../AGENTS.md)
- [plan_synth.md](../../../plan_synth.md)
- [wiki/index.md](../../index.md)
- [wiki/catalog.md](../../catalog.md)
- [wiki/gems/atmosphere.md](../atmosphere.md)
- [GRHayL/include/ghl_atmosphere.h](../../../GRHayL/include/ghl_atmosphere.h)
- [GRHayL/Atmosphere/set_prims_to_constant_atm.c](../../../GRHayL/Atmosphere/set_prims_to_constant_atm.c)
- [GRHayL/Atmosphere/set_prims_to_radial_falloff_atm.c](../../../GRHayL/Atmosphere/set_prims_to_radial_falloff_atm.c)
- [GRHayL/Atmosphere/make.code.defn](../../../GRHayL/Atmosphere/make.code.defn)
- [GRHayL/include/ghl.h](../../../GRHayL/include/ghl.h)
- [GRHayL/GRHayL_Core/initialize_eos.c](../../../GRHayL/GRHayL_Core/initialize_eos.c)
- [docs/raw/Atmosphere.dox](../../../docs/raw/Atmosphere.dox)
- [Unit_Tests/unit_test_grhayl_core_test_suite.c](../../../Unit_Tests/unit_test_grhayl_core_test_suite.c)
- [Unit_Tests/unit_test_tabulated_eos.c](../../../Unit_Tests/unit_test_tabulated_eos.c)
- [Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c](../../../Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c)
