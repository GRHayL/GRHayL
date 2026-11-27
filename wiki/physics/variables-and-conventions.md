# Variables And Conventions

> GRMHD primitive, conserved, magnetic, densitized, metric, and EOS variable conventions mapped to `ghl.h` and `docs/raw/derivation.md`. · Status: confirmed
> Up: [Physics](index.md)

## Summary

GRHayL evolves Valencia-form GRMHD equations using transport velocities `vU[i] = u^i / u^0`, magnetic fields rescaled by `1/sqrt(4*pi)`, densitized hydrodynamic conserved variables weighted by `sqrt_detgamma`, and family-dependent EOS entropy and atmosphere floors. `docs/raw/derivation.md` supplies the mathematical derivation and `GRHayL/include/ghl.h` defines the authoritative C struct layouts.

## Detail

### Primitive And Conserved Variables

- **Primitives (`ghl_primitive_quantities`):** stores `rho`, `press`, `eps`, `u0`, `vU[3]`, `BU[3]`, `Y_e`, `temperature`, and `entropy`. Transport velocity `vU[i] = u^i / u^0` relates to Eulerian velocity `(vU[i] + betaU[i]) / lapse` and Lorentz factor `W = lapse * u0`; see [Velocity And `u0`](../core/velocity-u0.md). Packing via `ghl_initialize_primitives` copies supplied fields without computing `u0`, so callers must set a consistent `u0` before conversion or stress-energy calls ([Struct Pack And Unpack](../core/struct-pack-unpack.md)). Constant atmosphere resets copy all scalar composition fields; simple/hybrid EOS setup zeroes the unused `Y_e_atm` and `T_atm` placeholders.
- **Conservatives (`ghl_conservative_quantities`):** stores `rho`, `tau`, `SD[3]`, `Y_e`, and `entropy`. `ghl_compute_conservs` writes `rho = sqrt_detgamma * W * prims.rho`, `Y_e = cons.rho * prims.Y_e`, and `entropy = sqrt_detgamma * W * prims.entropy`—so conserved `Y_e` is mass-weighted rather than a primitive fraction, and `entropy` carries `sqrt_detgamma * W` rather than `cons.rho`. `ghl_undensitize_conservatives` divides only by `sqrt_detgamma`; the same struct holds densitized or undensitized values depending on the called routine. Magnetic fields are evolved via vector-potential induction rather than a `BU` member in `ghl_conservative_quantities`.
- **Conservative energy floor `tau_atm`:** simple and hybrid EOS initialization compute `rho_atm * eps_atm`, whereas tabulated EOS uses `rho_min * eps_min`.

Evidence (scientific/interface requirement): `GRHayL/include/ghl.h::ghl_primitive_quantities`; corroboration: `GRHayL/Con2Prim/compute_conservs.c::ghl_compute_conservs` and `GRHayL/GRHayL_Core/initialize_eos.c`.

### Magnetic Field Rescaling And Densitized Quantities

`docs/raw/derivation.md` defines GRHayL's magnetic field as the unscaled (Gaussian) field divided by `sqrt(4*pi)` (Lorentz-Heaviside normalization; `ONE_OVER_SQRT_4PI` in `ghl.h`). Callers must supply already-rescaled `BU`: `initialize_primitives.c` copies `BU` verbatim and `compute_smallb_and_b2.c` consumes it without rescaling. Do not reintroduce `4*pi` factors in stress-energy, Flux_Source, or Induction routines.

Densitized hydrodynamic variables carry `sqrt_detgamma` from `ghl_metric_quantities`. In Induction (`docs/raw/Induction.dox`), `tilde{Phi}` (`phitilde`) is the densitized scalar potential and `tilde{B}` (`Btilde`) is the densitized magnetic field. The `psi6` argument in Induction routines denotes physical `sqrt_detgamma` (equating it to a BSSN conformal factor `e^{6\phi}` assumes a unit-determinant conformal metric): `ghl_HLL_flux_with_B` densitizes staggered `B` internally using `psi6`, whereas `ghl_HLL_flux_with_Btilde` expects already-densitized input ([HLL Flux](../gems/induction/hll-flux.md)).

Evidence (scientific/interface requirement): `GRHayL/include/ghl.h::ONE_OVER_SQRT_4PI`; corroboration: `GRHayL/GRHayL_Core/compute_smallb_and_b2.c::ghl_compute_smallb_and_b2` and `GRHayL/Induction/HLL_flux_with_B.c::ghl_HLL_flux_with_B`.

### Metric, Extrinsic Curvature, And EOS Variables

- **Metric structs:** `ghl_metric_quantities` stores lapse, inverse lapse values, shift `betaU[3]`, spatial 3-metric `gammaDD`, inverse `gammaUU`, `detgamma`, and `sqrt_detgamma`; `ghl_ADM_aux_quantities` stores 4-metric `g4DD` and inverse `g4UU`; `ghl_extrinsic_curvature` stores `K[3][3]`. Inline raise/lower helpers in `ghl_metric_helpers.h` require matching array indexing, and geometric source terms expect caller-supplied metric derivatives. Callers must supply positive-definite spatial metrics; see [Metric And ADM Helpers](../core/metric-adm.md) for current nonpositive-determinant handling (`CONTR-0011`).
- **EOS variables & flags:** `ghl_eos_parameters` stores atmosphere, minimum, and maximum bounds for density, pressure, specific internal energy, entropy, temperature, and `Y_e`.
  - `evolve_entropy` controls selected entropy updates (primitive limiting, hybrid Font recovery), not general entropy computation or entropy-variant HLLE flux calls.
  - `evolve_temp` controls temperature inversion versus supplied-temperature evaluation in tabulated Palenzuela recovery; other solvers do not uniformly read it ([Solver Matrix](../gems/con2prim/solver-matrix.md)).
  - Primitive `entropy` is EOS-family dependent: hybrid helpers compute `P / rho^(Gamma - 1)` using the active polytropic index, whereas tabulated entropy is interpolated from the HDF5 table—never interchange them.
  - Tabulated EOS dispatch installs `ghl_tabulated_free_beq_quantities` alongside interpolation pointers; see [EOS Dispatch](../core/eos-dispatch.md).

Evidence (descriptive behavior): `GRHayL/include/ghl.h::ghl_eos_parameters`; corroboration: `GRHayL/GRHayL_Core/initialize_eos.c::ghl_initialize_eos_functions`.

## Sources

- [docs/raw/derivation.md](../../docs/raw/derivation.md) - GRMHD variables, `sqrt(4*pi)` magnetic rescaling, and densitized conservative definitions
- [docs/raw/Induction.dox](../../docs/raw/Induction.dox) - `tilde{Phi}` and `tilde{B}` definitions
- [GRHayL/include/ghl.h](../../GRHayL/include/ghl.h) - `ghl_primitive_quantities`, `ghl_conservative_quantities`, `ghl_metric_quantities`, `ghl_eos_parameters`, `ONE_OVER_SQRT_4PI`
- [GRHayL/include/ghl_metric_helpers.h](../../GRHayL/include/ghl_metric_helpers.h) - inline metric raise/lower helpers
- [GRHayL/GRHayL_Core/initialize_eos.c](../../GRHayL/GRHayL_Core/initialize_eos.c) - `tau_atm` setup and `ghl_initialize_eos_functions`
- [GRHayL/GRHayL_Core/compute_smallb_and_b2.c](../../GRHayL/GRHayL_Core/compute_smallb_and_b2.c) - rescaled `BU` consumption
- [GRHayL/Con2Prim/compute_conservs.c](../../GRHayL/Con2Prim/compute_conservs.c) - densitized `rho`, `Y_e`, `entropy` formulas
- [GRHayL/Con2Prim/undensitize_conservatives.c](../../GRHayL/Con2Prim/undensitize_conservatives.c) - division by `sqrt_detgamma`
- [GRHayL/Induction/HLL_flux_with_B.c](../../GRHayL/Induction/HLL_flux_with_B.c) - internal `psi6` densitization

## See Also

- Depends on: [Evolution Equation Map](evolution-equation-map.md)
- Implements: [Shared Parameters And Enums](../core/shared-parameters-and-enums.md)
- See also: [Velocity And `u0`](../core/velocity-u0.md)
- See also: [Metric And ADM Helpers](../core/metric-adm.md)
- See also: [EOS Dispatch](../core/eos-dispatch.md)
