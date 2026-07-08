# Variables and Conventions

This page orients agents to variable meaning and source locations. Use `docs/raw/derivation.md` for derivations and `GRHayL/include/ghl.h` for struct authority.

## Primitive Variables

Physics orientation:
- `docs/raw/derivation.md` defines primitive variables as rest-mass density, pressure, transport velocity, and magnetic field.
- GRHayL stores primitive point data in `ghl_primitive_quantities` in `GRHayL/include/ghl.h`.
- Stored primitive fields include `rho`, `press`, `eps`, `u0`, `vU[3]`, `BU[3]`, `Y_e`, `temperature`, and `entropy`.

Primary code paths:
- Initialization/unpacking: `GRHayL/GRHayL_Core/initialize_primitives.c`, `GRHayL/GRHayL_Core/return_primitives.c`
- Primitive bounds and `u0`: `GRHayL/Con2Prim/enforce_primitive_limits_and_compute_u0.c`
- Recovery route: [Con2Prim recovery flow](../gems/con2prim/recovery-flow.md)
- Limit route: [Con2Prim limits and conversions](../gems/con2prim/limits-and-conversions.md)
- Atmosphere reset: `GRHayL/Atmosphere/`
- Reconstruction inputs: `GRHayL/Reconstruction/`

Contracts:
- `vU` is the transport velocity used in evolution fluxes.
- `BU` follows GRHayL's rescaled magnetic-field convention.
- Primitive bounds come from `ghl_eos_parameters` and `ghl_parameters`.

## Conservative Variables

Physics orientation:
- `docs/raw/derivation.md` defines conservative evolution variables including densitized density, energy, momentum, and magnetic field.
- GRHayL stores hydrodynamic conservative point data in `ghl_conservative_quantities` in `GRHayL/include/ghl.h`.
- Stored conservative fields include `rho`, `tau`, `SD[3]`, plus evolved scalar entries `Y_e` and `entropy`.

Primary code paths:
- Initialization/unpacking: `GRHayL/GRHayL_Core/initialize_conservatives.c`, `GRHayL/GRHayL_Core/return_conservatives.c`
- Primitive-to-conservative: `GRHayL/Con2Prim/compute_conservs.c`, `GRHayL/Con2Prim/compute_conservs_and_Tmunu.c`
- Conservative limits: `GRHayL/Con2Prim/apply_conservative_limits.c`
- Conservative-to-primitive recovery: `GRHayL/Con2Prim/`
- Recovery and conversion routes: [Con2Prim recovery flow](../gems/con2prim/recovery-flow.md),
  [Con2Prim limits and conversions](../gems/con2prim/limits-and-conversions.md)

Contracts:
- Densitized hydrodynamic fields depend on `sqrt_detgamma`.
- Conservative momentum and energy include magnetic contributions through primitive `BU` and metric data.
- Magnetic evolution is handled through induction/vector-potential paths, not by a `BU` member in `ghl_conservative_quantities`.

## Magnetic Field Rescaling

Physics orientation:
- `docs/raw/derivation.md` states that GRHayL uses a magnetic field rescaled by `sqrt(4*pi)`.
- `GRHayL/include/ghl.h` defines `ONE_OVER_SQRT_4PI`, matching that convention.
- The derivation assumes rescaled magnetic quantities everywhere after the convention is introduced.

Primary code paths:
- Stress-energy helpers: `GRHayL/GRHayL_Core/compute_smallb_and_b2.c`, `GRHayL/GRHayL_Core/compute_TDNmunu.c`, `GRHayL/GRHayL_Core/compute_TUPmunu.c`
- Flux/source terms: `GRHayL/Flux_Source/`; route characteristic speeds, HLLE
  variants, and source terms through [Flux_Source hub](../gems/flux-source.md)
- Induction: `GRHayL/Induction/`; for `Btilde` and vector-potential HLL
  densitization, use [Induction HLL flux contract](../gems/induction/hll-flux-contract.md)

Contracts:
- Do not reintroduce `4*pi` factors into GRMHD evolution paths unless the convention itself changes.
- Induction and flux/source paths must use consistent magnetic normalization.

## Densitized Quantities

Physics orientation:
- `docs/raw/derivation.md` defines densitized conservative variables with `sqrt(gamma)` factors.
- `ghl_metric_quantities` stores `detgamma` and `sqrt_detgamma` in `GRHayL/include/ghl.h`.
- `docs/raw/Induction.dox` defines `tilde{Phi}` as a densitized scalar potential and `tilde{B}` as the densitized magnetic field.

Primary code paths:
- Metric setup: `GRHayL/GRHayL_Core/initialize_metric.c`, `GRHayL/GRHayL_Core/compute_ADM_auxiliaries.c`
- Conservative conversion: `GRHayL/Con2Prim/compute_conservs.c`, `GRHayL/Con2Prim/undensitize_conservatives.c`
- Con2Prim densitization route: [Con2Prim recovery flow](../gems/con2prim/recovery-flow.md)
- Induction flux variants: `GRHayL/Induction/HLL_flux_with_B.c`, `GRHayL/Induction/HLL_flux_with_Btilde.c`; route `Btilde` and `psi6` details through [Induction HLL flux contract](../gems/induction/hll-flux-contract.md)
- Gauge RHS: `GRHayL/Induction/calculate_phitilde_rhs.c`; route `tildePhi`/`phitilde` details through [Induction gauge RHS contract](../gems/induction/gauge-rhs-contract.md)

Contracts:
- `psi6` and `sqrt_detgamma` are common names for the same densitization factor in these paths.
- `ghl_HLL_flux_with_B` takes undensitized staggered `B` and densitizes internally with `psi6`.
- `ghl_HLL_flux_with_Btilde` expects already densitized magnetic input.

## Metric and Helper Conventions

Physics orientation:
- `ghl_metric_quantities` stores lapse, inverse lapse values, shift, 3-metric, inverse 3-metric, determinant, and square root determinant.
- `ghl_ADM_aux_quantities` stores 4-metric and inverse 4-metric.
- `ghl_extrinsic_curvature` stores ADM extrinsic curvature.

Primary code paths:
- Public helpers: `GRHayL/include/ghl_metric_helpers.h`
- Metric packing: `GRHayL/GRHayL_Core/initialize_metric.c`, `GRHayL/GRHayL_Core/enforce_detgij_and_initialize_ADM_metric.c`
- ADM auxiliaries: `GRHayL/GRHayL_Core/compute_ADM_auxiliaries.c`
- Source terms: `GRHayL/Flux_Source/ghl_calculate_source_terms.c`; caller-owned
  metric derivative routing lives in
  [Flux_Source source-term contract](../gems/flux-source/source-terms-contract.md)

Contracts:
- Raise/lower helpers work with arrays matching GRHayL metric array layout.
- Source terms require metric derivatives supplied by callers.
- Metric helper changes can affect stress-energy, Con2Prim, flux/source, and induction interpolation.

## EOS Variables

Physics orientation:
- `ghl_eos_parameters` in `GRHayL/include/ghl.h` stores atmosphere, minimum, and maximum bounds for density, pressure, internal energy, entropy, temperature, and `Y_e`.
- Hybrid EOS routines provide cold/thermal pressure and entropy-function support.
- Tabulated EOS routines use density, `Y_e`, and temperature or other table variables to interpolate pressure, internal energy, entropy, sound speed, and chemical potentials.

Primary code paths:
- Function pointers: `GRHayL/include/ghl_eos_functions.h`
- Hybrid EOS: `GRHayL/EOS/Hybrid/`
- Tabulated EOS: `GRHayL/EOS/Tabulated/`; route built interpolation wrappers
  through [tabulated interpolator catalog](../gems/eos/tabulated-interpolator-catalog.md)
- Neutrino leakage EOS inputs: `GRHayL/Neutrinos/NRPyLeakage/`

Contracts:
- Entropy evolution is controlled by `ghl_parameters::evolve_entropy`.
- Temperature evolution is controlled by `ghl_parameters::evolve_temp`.
- `Y_e`, temperature, and entropy are especially coupled between tabulated EOS, Con2Prim, flux entropy variants, and neutrino leakage.
