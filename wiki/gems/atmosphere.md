# Atmosphere Gem

## Purpose

Atmosphere routines reset primitive variables to prescribed low-density floor states.

## Read First

- [Atmosphere prescription contract](atmosphere/prescription-contract.md)
- `docs/raw/Atmosphere.dox`
- `GRHayL/include/ghl_atmosphere.h`
- `GRHayL/include/ghl.h`
- `wiki/physics/variables-and-conventions.md`

## Public Headers

- `GRHayL/include/ghl_atmosphere.h`
- `GRHayL/include/ghl.h`

Built public surface:
- `ghl_set_prims_to_constant_atm`; see
  [prescription contract](atmosphere/prescription-contract.md).

Declared/source-present but unbuilt surface:
- `ghl_set_prims_to_radial_falloff_atm` is declared in
  `GRHayL/include/ghl_atmosphere.h` and has a source file, but
  `GRHayL/Atmosphere/make.code.defn` does not compile it and
  `docs/raw/Atmosphere.dox` says only one prescription is provided. No repo
  call or test was found. This establishes current manifest absence, not
  intended removal, internal status, or future support. Route radial questions to
  [prescription contract](atmosphere/prescription-contract.md).

`GRHayL/include/make.code.defn` lists `ghl_atmosphere.h` for installation.
Installed declaration, built definition, documentation, calls, and tests are
separate evidence; only the constant routine currently closes those surfaces.

## Implementation Paths

- Active build path: `GRHayL/Atmosphere/set_prims_to_constant_atm.c`
- Unbuilt/incomplete source path:
  `GRHayL/Atmosphere/set_prims_to_radial_falloff_atm.c`
- `GRHayL/Atmosphere/make.code.defn`

## Test Paths

- Direct constant-atmosphere coverage: `Unit_Tests/unit_test_grhayl_core_test_suite.c`
- Legacy branch coverage: `Unit_Tests/unit_test_ET_Legacy_primitives.c`
- Nearby downstream coverage: `Unit_Tests/unit_test_enforce_primitive_limits_and_compute_u0.c`, `Unit_Tests/unit_test_apply_conservative_limits.c`, `Unit_Tests/unit_test_hybrid_failure.c`

## Key Contracts

- Constant reset copies density, pressure, internal energy, entropy,
  temperature, and `Y_e` from `ghl_eos_parameters`; it zeros all three
  transport velocities.
- It does not write magnetic fields or `u0`. Existing magnetic values remain;
  `u0` may be stale until recomputed.
- Simple/hybrid EOS initialization does not populate `Y_e_atm` or `T_atm`,
  although constant reset reads both. Callers using those EOS families must
  initialize the fields or accept that no valid value is established.

## Common Edit Routes

- Change atmosphere values: update EOS parameter initialization/bounds and primitive-limit tests.
- Add/activate prescription: reconcile implementation, declaration,
  build-list entry, Doxygen entry, call surface, and direct unit coverage.

## Drift Risks

- Atmosphere floors affect Con2Prim recovery, conservative limiting, and EOS bounds.
- Public header changes affect GRHayLib through `GRHayLib.h`.

## Do Not Duplicate

Keep prescription-level details in `docs/raw/Atmosphere.dox` and source. This page routes edits and tests.
