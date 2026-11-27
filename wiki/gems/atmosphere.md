# Atmosphere

> Active `ghl_set_prims_to_constant_atm` floor reset, EOS field coupling, test coverage, and the unbuilt radial-falloff routine. · Status: stale
> Up: [GRHayL Gems](index.md)

## Summary

`GRHayL/Atmosphere/` resets primitive variables to low-density floor values stored in `ghl_eos_parameters`. `ghl_set_prims_to_constant_atm` is the sole declared, manifest-built, documented, and tested prescription; `ghl_set_prims_to_radial_falloff_atm` is header-declared and present as incomplete C source, but omitted from `GRHayL/Atmosphere/make.code.defn`, Doxygen, and tests.

## Detail

### Constant Atmosphere Reset And EOS Field Coupling

`ghl_set_prims_to_constant_atm` (`GRHayL/Atmosphere/set_prims_to_constant_atm.c`, declared in installed `GRHayL/include/ghl_atmosphere.h`) writes:
- `prims->rho = eos->rho_atm`, `press = eos->press_atm`, `eps = eos->eps_atm`, `entropy = eos->entropy_atm`, `Y_e = eos->Y_e_atm`, `temperature = eos->T_atm`;
- `prims->vU[0..2] = 0`.

It leaves `prims->BU[0..2]` and `prims->u0` untouched, returns `void`, and performs no pointer, EOS-family, or finiteness validation. Callers needing a complete primitive state preserve `BU` separately and recompute `u0`.

- **Production caller:** `GRHayL/Con2Prim/Tabulated/tabulated_primitive_guess_helpers.c` invokes `ghl_set_prims_to_constant_atm` on atmosphere fallback, then sets `vU = -betaU` and `u0 = lapseinv` ([Limits And Conversions](con2prim/limits-and-conversions.md), [Neural-Network Primitive Guess](con2prim/neural-network-primitive-guess.md)).
- **EOS initialization coupling:** `GRHayL/GRHayL_Core/initialize_eos.c` populates `rho_atm`, `press_atm`, `eps_atm`, `entropy_atm`, and `tau_atm` for simple/hybrid EOS while storing deterministic zero placeholders in unused `Y_e_atm` and `T_atm` so constant reset always reads defined values; tabulated EOS derives `press_atm`, `eps_atm`, and `entropy_atm` from table-clamped `(rho_atm, Y_e_atm, T_atm)`.

Evidence (descriptive behavior): `GRHayL/Atmosphere/set_prims_to_constant_atm.c::ghl_set_prims_to_constant_atm`; corroboration: `GRHayL/Con2Prim/Tabulated/tabulated_primitive_guess_helpers.c` and `GRHayL/GRHayL_Core/initialize_eos.c`.

### Unbuilt Radial-Falloff Routine And Test Coverage

Claim status: stale; contradiction: CONTR-0009. See [CONTR-0009](../contradictions.md#contr-0009).

| Routine | Header & manifest status | Current behavior |
| --- | --- | --- |
| `ghl_set_prims_to_constant_atm` | Installed (`ghl_atmosphere.h`), built (`make.code.defn`), documented (`Atmosphere.dox`), called in production and tests. (Ev: `GRHayL/Atmosphere/make.code.defn::set_prims_to_constant_atm.c`) | Full thermodynamic reset and zero velocity; preserves `BU` and `u0`. |
| `ghl_set_prims_to_radial_falloff_atm` | Installed (`ghl_atmosphere.h`), source-present, **omitted** from `make.code.defn`, `Atmosphere.dox`, callers, and tests. (Ev: `GRHayL/Atmosphere/set_prims_to_radial_falloff_atm.c::ghl_set_prims_to_radial_falloff_atm`) | Incomplete: `r` is unused; `rho`, `press`, `eps` writes are commented out; only copies `entropy`/`Y_e`/`temperature` and zeros `vU`. |

**Test coverage:**
- `Unit_Tests/unit_test_grhayl_core_test_suite.c` directly verifies `rho`, `press`, `eps`, `entropy`, and `vU = 0` for simple and hybrid EOS; its tabulated branch is commented out (`TODO`), and its `Y_e`/`temperature` checks sit inside an unreached `eos_type == 2` block.
- `Unit_Tests/unit_test_c2p_nn_guess.c` indirectly asserts tabulated atmosphere fields (including `Y_e` and `temperature`), `vU = -betaU`, `u0 = lapseinv`, and unchanged `BU` through synthetic tabulated metadata and fake callbacks.
- `Unit_Tests/unit_test_ET_Legacy_primitives.c` calls constant reset conditionally when `cons.rho <= 0.0` (with a `TODO` for explicit reset assertions). No test or caller references `ghl_set_prims_to_radial_falloff_atm`.

Evidence (descriptive behavior): `GRHayL/include/ghl_atmosphere.h::ghl_set_prims_to_radial_falloff_atm`; corroboration: `Unit_Tests/unit_test_grhayl_core_test_suite.c`.

## Sources

- [GRHayL/include/ghl_atmosphere.h](../../GRHayL/include/ghl_atmosphere.h) - constant and radial-falloff declarations
- [GRHayL/Atmosphere/make.code.defn](../../GRHayL/Atmosphere/make.code.defn) - `SRCS` (`set_prims_to_constant_atm.c` only)
- [GRHayL/Atmosphere/set_prims_to_constant_atm.c](../../GRHayL/Atmosphere/set_prims_to_constant_atm.c) - constant atmosphere implementation
- [GRHayL/Atmosphere/set_prims_to_radial_falloff_atm.c](../../GRHayL/Atmosphere/set_prims_to_radial_falloff_atm.c) - incomplete unbuilt radial-falloff source
- [GRHayL/Con2Prim/Tabulated/tabulated_primitive_guess_helpers.c](../../GRHayL/Con2Prim/Tabulated/tabulated_primitive_guess_helpers.c) - tabulated primitive-guess caller
- [GRHayL/GRHayL_Core/initialize_eos.c](../../GRHayL/GRHayL_Core/initialize_eos.c) - EOS atmosphere field setup and zero placeholders
- [docs/raw/Atmosphere.dox](../../docs/raw/Atmosphere.dox) - constant atmosphere Doxygen documentation
- [Unit_Tests/unit_test_grhayl_core_test_suite.c](../../Unit_Tests/unit_test_grhayl_core_test_suite.c) - direct simple/hybrid constant-atmosphere checks
- [Unit_Tests/unit_test_c2p_nn_guess.c](../../Unit_Tests/unit_test_c2p_nn_guess.c) - indirect tabulated atmosphere fallback checks

## See Also

- Depends on: [EOS Dispatch](../core/eos-dispatch.md)
- See also: [Limits And Conversions](con2prim/limits-and-conversions.md)
- See also: [Public API Map](../architecture/public-api-map.md)
- See also: [Contradictions](../contradictions.md#contr-0009)
