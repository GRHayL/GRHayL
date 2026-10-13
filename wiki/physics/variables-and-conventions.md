# Variables and Conventions

This page orients agents to variable meaning and source locations. Use [docs/raw/derivation.md](../../docs/raw/derivation.md) for derivations and [GRHayL/include/ghl.h](../../GRHayL/include/ghl.h) for struct authority.

## Primitive Variables

Physics orientation:
- `docs/raw/derivation.md` defines primitive variables as rest-mass density, pressure, transport velocity, and magnetic field.
- GRHayL stores primitive point data in `ghl_primitive_quantities` in `GRHayL/include/ghl.h`.
- Stored primitive fields include `rho`, `press`, `eps`, `u0`, `vU[3]`, `BU[3]`, `Y_e`, `temperature`, and `entropy`.

Primary code paths:
- Initialization/unpacking: [GRHayL/GRHayL_Core/initialize_primitives.c](../../GRHayL/GRHayL_Core/initialize_primitives.c), [GRHayL/GRHayL_Core/return_primitives.c](../../GRHayL/GRHayL_Core/return_primitives.c)
- Primitive bounds and `u0`: [GRHayL/Con2Prim/enforce_primitive_limits_and_compute_u0.c](../../GRHayL/Con2Prim/enforce_primitive_limits_and_compute_u0.c)
- Recovery route: [Con2Prim recovery flow](../gems/con2prim/recovery-flow.md)
- Limit route: [Con2Prim limits and conversions](../gems/con2prim/limits-and-conversions.md)
- Atmosphere reset: [GRHayL/Atmosphere/](../../GRHayL/Atmosphere)
- Reconstruction inputs: [GRHayL/Reconstruction/](../../GRHayL/Reconstruction)

Contracts:
- `vU[i] = u^i/u^0` is the transport velocity used in evolution fluxes.
  The Eulerian velocity is `(vU[i] + betaU[i])/lapse`, and the Lorentz
  factor is `W = lapse*u0`; see the [Core velocity contract](../core/velocity-u0-contract.md).
- Primitive packing copies the supplied thermodynamic and vector fields but
  leaves `u0` unchanged. Initialize a consistent `u0` before conversion or
  stress-energy evaluation; see [packing](../core/struct-pack-unpack-contract.md).
- `BU` follows GRHayL's rescaled magnetic-field convention.
- Primitive bounds come from `ghl_eos_parameters` and `ghl_parameters`.
- Constant atmosphere copies all scalar composition fields. Simple/hybrid EOS
  initialization defines the unused `Y_e_atm` and `T_atm` placeholders as zero.

## Conservative Variables

Physics orientation:
- `docs/raw/derivation.md` defines conservative evolution variables including densitized density, energy, momentum, and magnetic field.
- GRHayL stores hydrodynamic conservative point data in `ghl_conservative_quantities` in `GRHayL/include/ghl.h`.
- Stored conservative fields include `rho`, `tau`, `SD[3]`, plus evolved scalar entries `Y_e` and `entropy`.

Primary code paths:
- Initialization/unpacking: [GRHayL/GRHayL_Core/initialize_conservatives.c](../../GRHayL/GRHayL_Core/initialize_conservatives.c), [GRHayL/GRHayL_Core/return_conservatives.c](../../GRHayL/GRHayL_Core/return_conservatives.c)
- Primitive-to-conservative: [GRHayL/Con2Prim/compute_conservs.c](../../GRHayL/Con2Prim/compute_conservs.c), [GRHayL/Con2Prim/compute_conservs_and_Tmunu.c](../../GRHayL/Con2Prim/compute_conservs_and_Tmunu.c)
- Conservative limits: [GRHayL/Con2Prim/apply_conservative_limits.c](../../GRHayL/Con2Prim/apply_conservative_limits.c)
- Conservative-to-primitive recovery: [GRHayL/Con2Prim/](../../GRHayL/Con2Prim)
- Recovery and conversion routes: [Con2Prim recovery flow](../gems/con2prim/recovery-flow.md),
  [Con2Prim limits and conversions](../gems/con2prim/limits-and-conversions.md)

Contracts:
- Densitized hydrodynamic fields depend on `sqrt_detgamma`.
- `compute_conservs` writes `rho = sqrt_detgamma*W*prims.rho`,
  `Y_e = cons.rho*prims.Y_e`, and
  `entropy = sqrt_detgamma*W*prims.entropy`. Conservative `Y_e` is therefore
  not a primitive electron fraction, and entropy uses a different density
  weighting. [Undensitization](../../GRHayL/Con2Prim/undensitize_conservatives.c)
  divides these fields by `sqrt_detgamma`; it does not recover their primitive
  meanings. The same conservative struct carries either representation, so
  follow each routine's contract; conservative packing copies values without
  converting their representation.
- Conservative momentum and energy include magnetic contributions through primitive `BU` and metric data.
- Magnetic evolution is handled through induction/vector-potential paths, not by a `BU` member in `ghl_conservative_quantities`.
- `tau_atm` is used as a conservative energy floor, but initialization is
  family-dependent: simple/hybrid use `rho_atm*eps_atm`; tabulated uses
  `rho_min*eps_min`, matching the family-specific floor semantics documented in
  [`ghl.h`](../../GRHayL/include/ghl.h) and implemented in
  [`initialize_eos.c`](../../GRHayL/GRHayL_Core/initialize_eos.c).

## Magnetic Field Rescaling

Physics orientation:
- [`docs/raw/derivation.md`](../../docs/raw/derivation.md) defines GRHayL's
  magnetic field as the unscaled field divided by `sqrt(4*pi)`. Primitive `BU`
  must already use that normalization;
  [`initialize_primitives.c`](../../GRHayL/GRHayL_Core/initialize_primitives.c)
  copies it and
  [`compute_smallb_and_b2.c`](../../GRHayL/GRHayL_Core/compute_smallb_and_b2.c)
  consumes it without another conversion.
- `GRHayL/include/ghl.h` defines `ONE_OVER_SQRT_4PI`, matching that convention.
- The derivation assumes rescaled magnetic quantities everywhere after the convention is introduced.

Primary code paths:
- Stress-energy helpers: [GRHayL/GRHayL_Core/compute_smallb_and_b2.c](../../GRHayL/GRHayL_Core/compute_smallb_and_b2.c), [GRHayL/GRHayL_Core/compute_TDNmunu.c](../../GRHayL/GRHayL_Core/compute_TDNmunu.c), [GRHayL/GRHayL_Core/compute_TUPmunu.c](../../GRHayL/GRHayL_Core/compute_TUPmunu.c)
- Flux/source terms: [GRHayL/Flux_Source/](../../GRHayL/Flux_Source); route characteristic speeds, HLLE
  variants, and source terms through [Flux_Source hub](../gems/flux-source.md)
- Induction: [GRHayL/Induction/](../../GRHayL/Induction); for `Btilde` and vector-potential HLL
  densitization, use [Induction HLL flux contract](../gems/induction/hll-flux-contract.md)

Contracts:
- Do not reintroduce `4*pi` factors into GRMHD evolution paths unless the convention itself changes.
- Induction and flux/source paths must use consistent magnetic normalization.

## Densitized Quantities

Physics orientation:
- `docs/raw/derivation.md` defines densitized conservative variables with `sqrt(gamma)` factors.
- `ghl_metric_quantities` stores `detgamma` and `sqrt_detgamma` in `GRHayL/include/ghl.h`.
- [docs/raw/Induction.dox](../../docs/raw/Induction.dox) defines `tilde{Phi}` as a densitized scalar potential and `tilde{B}` as the densitized magnetic field.

Primary code paths:
- Metric setup: [GRHayL/GRHayL_Core/initialize_metric.c](../../GRHayL/GRHayL_Core/initialize_metric.c), [GRHayL/GRHayL_Core/compute_ADM_auxiliaries.c](../../GRHayL/GRHayL_Core/compute_ADM_auxiliaries.c)
- Conservative conversion: `GRHayL/Con2Prim/compute_conservs.c`, [GRHayL/Con2Prim/undensitize_conservatives.c](../../GRHayL/Con2Prim/undensitize_conservatives.c)
- Con2Prim densitization route: [Con2Prim recovery flow](../gems/con2prim/recovery-flow.md)
- Induction flux variants: [GRHayL/Induction/HLL_flux_with_B.c](../../GRHayL/Induction/HLL_flux_with_B.c), [GRHayL/Induction/HLL_flux_with_Btilde.c](../../GRHayL/Induction/HLL_flux_with_Btilde.c); route `Btilde` and `psi6` details through [Induction HLL flux contract](../gems/induction/hll-flux-contract.md)
- Gauge RHS: [GRHayL/Induction/calculate_phitilde_rhs.c](../../GRHayL/Induction/calculate_phitilde_rhs.c); route `tildePhi`/`phitilde` details through [Induction gauge RHS contract](../gems/induction/gauge-rhs-contract.md)

Contracts:
- `psi6` arguments in these APIs mean the physical `sqrt_detgamma` factor.
  Identifying it with the sixth power of a BSSN conformal factor assumes a
  unit-determinant conformal metric; see
  [Induction interpolation](../gems/induction/interpolation-and-staggering-contract.md).
- `ghl_HLL_flux_with_B` takes undensitized staggered `B` and densitizes internally with `psi6`.
- `ghl_HLL_flux_with_Btilde` expects already densitized magnetic input.

## Metric and Helper Conventions

Physics orientation:
- `ghl_metric_quantities` stores lapse, inverse lapse values, shift, 3-metric, inverse 3-metric, determinant, and square root determinant.
- `ghl_ADM_aux_quantities` stores 4-metric and inverse 4-metric.
- `ghl_extrinsic_curvature` stores ADM extrinsic curvature.

Primary code paths:
- Public helpers: [GRHayL/include/ghl_metric_helpers.h](../../GRHayL/include/ghl_metric_helpers.h)
- Metric packing: `GRHayL/GRHayL_Core/initialize_metric.c`, [GRHayL/GRHayL_Core/enforce_detgij_and_initialize_ADM_metric.c](../../GRHayL/GRHayL_Core/enforce_detgij_and_initialize_ADM_metric.c)
- ADM auxiliaries: `GRHayL/GRHayL_Core/compute_ADM_auxiliaries.c`
- Source terms: [GRHayL/Flux_Source/ghl_calculate_source_terms.c](../../GRHayL/Flux_Source/ghl_calculate_source_terms.c); caller-owned
  metric derivative routing lives in
  [Flux_Source source-term contract](../gems/flux-source/source-terms-contract.md)

Contracts:
- Raise/lower helpers work with arrays matching GRHayL metric array layout.
- Source terms require metric derivatives supplied by callers.
- Metric helper changes can affect stress-energy, Con2Prim, flux/source, and induction interpolation.
- Physics callers require a positive-definite spatial metric. Current Core
  initializers do not return an error for zero/negative determinants, and the
  determinant-enforced path can preserve a negative sign after warning.

## EOS Variables

Physics orientation:
- `ghl_eos_parameters` in `GRHayL/include/ghl.h` stores atmosphere, minimum, and maximum bounds for density, pressure, internal energy, entropy, temperature, and `Y_e`.
- Hybrid EOS routines provide cold/thermal pressure and entropy-function support.
- Tabulated EOS routines use density, `Y_e`, and temperature or other table variables to interpolate pressure, internal energy, entropy, sound speed, and chemical potentials.

Primary code paths:
- Function pointers: [GRHayL/include/ghl_eos_functions.h](../../GRHayL/include/ghl_eos_functions.h)
- Hybrid EOS: [GRHayL/EOS/Hybrid/](../../GRHayL/EOS/Hybrid)
- Tabulated EOS: [GRHayL/EOS/Tabulated/](../../GRHayL/EOS/Tabulated); route built interpolation wrappers
  through [tabulated interpolator catalog](../gems/eos/tabulated-interpolator-catalog.md)
- Neutrino leakage EOS inputs: [GRHayL/Neutrinos/NRPyLeakage/](../../GRHayL/Neutrinos/NRPyLeakage); route chemical
  potentials, free-nucleon fractions, units, species ordering, and source-term
  signs through the [physics and EOS contract](../gems/neutrinos/physics-and-eos-contract.md)

Contracts:
- `evolve_entropy` controls selected entropy updates, including primitive
  limiting and hybrid Font recovery; it is not a global gate on entropy
  computations or entropy-variant flux calls.
- `evolve_temp` selects temperature inversion versus supplied-temperature
  evaluation in tabulated Palenzuela recovery. Other solvers do not uniformly
  consult this flag; see [solver matrix](../gems/con2prim/solver-matrix.md).
  Time evolution itself belongs to the host application.
- Primitive entropy is EOS-family dependent: the hybrid helper computes
  `P/rho^(Gamma-1)` using the active polytropic index, while tabulated entropy
  comes from the table. Do not interchange those conventions; see
  [hybrid EOS](../gems/eos/hybrid-piecewise-polytrope.md) and the
  [table adapter](../gems/eos/stellarcollapse-table-adapter.md).
- `Y_e`, temperature, and entropy are especially coupled between tabulated EOS, Con2Prim, flux entropy variants, and neutrino leakage.
- Function-pointer declaration/storage does not by itself establish
  initialization; tabulated dispatch setup assigns
  `ghl_tabulated_free_beq_quantities` to the concrete cleanup implementation.
