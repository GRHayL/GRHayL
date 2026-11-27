# Characteristic Speeds

> Direction-split characteristic speed kernels (`ghl_calculate_characteristic_speed_dirn{0,1,2}`), sign convention, EOS callback mutation, and Induction coupling. · Status: confirmed
> Up: [Flux_Source](index.md)

## Summary

`GRHayL/Flux_Source/ghl_calculate_characteristic_speed_dirn{0,1,2}.c` compute non-negative left-going (`cmin`) and right-going (`cmax`) coordinate-speed magnitudes `dx^i/dt` at cell faces in both default and `--disable-hdf5` builds. Every direction provides a legacy `void` entry point (which aborts on EOS errors) and a `_checked` variant returning `ghl_error_codes_t`, while mutable primitive arguments may be clamped in place by production tabulated `ghl_compute_h_and_cs2`.

## Detail

### Public API, Sign Convention, And EOS Callback Mutation

`GRHayL/include/ghl_flux_source.h` declares legacy and `_checked` pairs for `dirn0`, `dirn1`, and `dirn2` (built by `GRHayL/Flux_Source/make.code.defn`, generated from `generate_flux_source.py`). Each takes reconstructed face states `prims_r` and `prims_l`, `eos`, `metric_face`, and output pointers `cmin_dirn*`, `cmax_dirn*`:
- **Sign convention:** `cmin` is a **non-negative** left-going speed magnitude (the negated minimum of `0` and the left-going roots), not a signed minimum eigenvalue. Branch-free expressions can leave tiny negative roundoff residues that HLLE kernels floor to `0` within tolerance. A negative characteristic discriminant is clipped to `0` without returning an error (`docs/raw/Flux_Source.dox`).
- **Preconditions:** Supply finite ADM geometry with positive lapse, finite physical primitives with rescaled `BU` and initialized positive `u0`, and usable EOS enthalpy/sound speed. `_checked` propagates EOS errors only; it does not validate input physics or output finiteness.
- **Mutable primitive arguments:** Both states call global `ghl_compute_h_and_cs2`. Production tabulated dispatch clamps `rho`, `Y_e`, and `temperature` to configured EOS bounds and recomputes `press` and `eps` in place ([Tabulated Interpolation And Bounds](../eos/tabulated-interpolation-and-bounds.md)). On error, `_checked` leaves `cmin`/`cmax` unchanged while retaining any primitive mutation from the first or failing callback—pass copies if face states must remain immutable.

Evidence (descriptive behavior): `GRHayL/include/ghl_flux_source.h::ghl_calculate_characteristic_speed_dirn0_checked`; corroboration: `GRHayL/Flux_Source/ghl_calculate_characteristic_speed_dirn0.c`.

### Induction Coupling And Test Coverage

`ghl_HLL_vars` in `GRHayL/include/ghl_induction.h` expects `c1_min`, `c1_max`, `c2_min`, `c2_max` under this direction/sign convention, but no repository test feeds Flux_Source characteristic-speed outputs into Induction HLL routines (Induction HLL generators fill random signed speeds).

- `unit_test_ET_Legacy_flux_source.c` calls all three legacy speed wrappers with a test-local `compute_h_and_cs2` callback.
- `unit_test_hybrid_flux.c` and `unit_test_tabulated_flux.c` replay precomputed speed arrays from generators, while `unit_test_hybrid_flux.c` additionally injects first- and second-callback failures into every `_checked` speed direction to verify unchanged outputs.
- `unit_test_tabulated_eos_compose.c` validates analytic magnetized speed bounds in all three directions and, in direction 0 only, table-bound clamping of `rho`, `Y_e`, and `temperature` with finite `press` and `eps`.

Evidence (descriptive behavior): `Unit_Tests/unit_test_hybrid_flux.c`; corroboration: `Unit_Tests/unit_test_tabulated_eos_compose.c`.

## Sources

- [GRHayL/include/ghl_flux_source.h](../../../GRHayL/include/ghl_flux_source.h) - legacy and `_checked` speed declarations
- [GRHayL/Flux_Source/make.code.defn](../../../GRHayL/Flux_Source/make.code.defn) - built `dirn{0,1,2}` sources
- [GRHayL/Flux_Source/ghl_calculate_characteristic_speed_dirn0.c](../../../GRHayL/Flux_Source/ghl_calculate_characteristic_speed_dirn0.c) - `dirn0` speed kernel
- [GRHayL/Flux_Source/generate_flux_source.py](../../../GRHayL/Flux_Source/generate_flux_source.py) - `generate_speeds()`
- [GRHayL/include/ghl_induction.h](../../../GRHayL/include/ghl_induction.h) - `ghl_HLL_vars` speed coupling
- [docs/raw/Flux_Source.dox](../../../docs/raw/Flux_Source.dox) - discriminant clipping and speed definitions
- [Unit_Tests/unit_test_hybrid_flux.c](../../../Unit_Tests/unit_test_hybrid_flux.c) - checked callback-error tests
- [Unit_Tests/unit_test_tabulated_eos_compose.c](../../../Unit_Tests/unit_test_tabulated_eos_compose.c) - production tabulated speed and mutation checks

## See Also

- Depends on: [Tabulated Interpolation And Bounds](../eos/tabulated-interpolation-and-bounds.md)
- Implements: [Generated NRPy Boundary](generated-nrpy-boundary.md)
- See also: [HLLE Flux Variant Matrix](hlle-flux-variant-matrix.md)
- See also: [HLL Flux](../induction/hll-flux.md)
