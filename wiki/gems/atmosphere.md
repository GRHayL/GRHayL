# Atmosphere Gem

## Purpose

Atmosphere routines reset primitive variables to prescribed low-density floor states.

## Read First

- `docs/raw/Atmosphere.dox`
- `GRHayL/include/ghl_atmosphere.h`
- `GRHayL/include/ghl.h`
- `wiki/physics/variables-and-conventions.md`

## Public Headers

- `GRHayL/include/ghl_atmosphere.h`
- `GRHayL/include/ghl.h`

Key active public surface:
- `ghl_set_prims_to_constant_atm`

Declared/unbuilt drift candidate:
- `ghl_set_prims_to_radial_falloff_atm` is declared in
  `GRHayL/include/ghl_atmosphere.h` and has a source file, but
  `GRHayL/Atmosphere/make.code.defn` does not compile it and
  `docs/raw/Atmosphere.dox` says only one prescription is provided.

## Implementation Paths

- Active build path: `GRHayL/Atmosphere/set_prims_to_constant_atm.c`
- Drift candidate: `GRHayL/Atmosphere/set_prims_to_radial_falloff_atm.c`
- `GRHayL/Atmosphere/make.code.defn`

## Test Paths

- Direct constant-atmosphere coverage: `Unit_Tests/unit_test_grhayl_core_test_suite.c`
- Nearby indirect coverage: `Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c`, `Unit_Tests/unit_test_apply_conservative_limits.c`, `Unit_Tests/unit_test_hybrid_failure.c`

## Key Contracts

- Atmosphere values come from `ghl_eos_parameters` fields such as density, pressure, internal energy, entropy, temperature, and `Y_e` atmosphere settings.
- Reset primitives must remain compatible with `ghl_enforce_primitive_limits_and_compute_u0`.
- Magnetic fields and velocity handling should preserve the assumptions expected by downstream conservative conversion.

## Common Edit Routes

- Change atmosphere values: update EOS parameter initialization/bounds and primitive-limit tests.
- Add prescription: add implementation, declaration, build-list entry, Doxygen entry, and direct or indirect unit coverage.

## Drift Risks

- Atmosphere floors affect Con2Prim recovery, conservative limiting, and EOS bounds.
- Public header changes affect GRHayLib through `GRHayLib.h`.

## Do Not Duplicate

Keep prescription-level details in `docs/raw/Atmosphere.dox` and source. This page routes edits and tests.
