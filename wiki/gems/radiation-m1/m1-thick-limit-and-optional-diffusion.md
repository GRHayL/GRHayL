# Thick limit and optional diffusion flux

The M1 library contains a public thick-limit scalar and an optional
energy-flux diffusion correction. These helpers are useful for experiments or
downstream methods that explicitly select them. They are not part of the
canonical neutrino four-point blended Rusanov route: that route uses metric
light-cone speeds, does not apply a speed cap, and rejects its diffusion-policy
flag when enabled.

## Thick-limit comoving energy

[`ghl_m1_compute_Jthick`](../../../GRHayL/Radiation/ghl_m1_Jthick.c) computes
the thick-limit comoving energy estimate from an Eulerian state and fluid
velocity:

```text
J_thick = 3 / (2 W^2 + 1)
          * ((2 W^2 - 1) E - 2 W^2 F_i V^i).
```

The helper returns a validity flag. A successful return means that the
calculation completed; it does not certify that the candidate is usable. The
candidate must be finite and strictly positive, and callers must inspect
`Jthick_is_valid`, before passing it to the diffusion correction. A nonpositive
or nonfinite candidate is reported as invalid rather than repaired into a
diffusion state.

## Optional diffusion inputs

[`ghl_m1_compute_diffusion_flux`](../../../GRHayL/Radiation/ghl_m1_diffusion_flux.c)
is an energy-flux helper. The neutrino wrapper
`ghl_m1_compute_neutrino_diffusion_flux` obtains `chi_tr` from
`rates_face->kappa_tr` and otherwise uses the same interface. The caller
supplies:

- an already computed densitized thin-side numerical energy flux (the public
  parameter is historically named `hll_flux_tildeE`);
- `E_star`, the interface energy used for the shift term;
- valid left and right `J_thick` values;
- a face gradient of `J_thick`;
- face `W` and `V^i`;
- transport opacity `chi_tr`, diffusion coefficient `D_face`, and proper
  normal length `delta_l`.

The diffusion kernel does not construct a face coefficient from cell rates. If
the host wants the available harmonic choice, it can use
[`ghl_m1_compute_harmonic_diffusion_coefficient`](../../../GRHayL/Radiation/ghl_m1_utils.c)
and validate `D_face` before the call. Likewise, if a downstream pair model needs a
partner-dependent inverse-pair opacity, that total face opacity must be
assembled by the host; the scalar `kappa_tr` field alone is the independent
absorption/scattering contribution in the pair contract.

## Optical-depth blend

The optional helper forms

```text
tau_face = chi_tr * delta_l.
```

It returns the original input energy flux unchanged when `tau_face` is below
`m1_params->zeta_min`, either `J_thick` input is invalid, or the supplied face
velocity fails the timelike consistency check. In these no-op cases the blend
factor is `a_face = 1` when that output is requested.

These no-op cases return `ghl_success` and publish the input flux and blend
factor. Every non-success return, including direction, metric, gradient, and
late arithmetic failures, leaves each non-NULL output argument unchanged.

When the inputs are valid, it uses

```text
a = tanh(1 / tau_face)
J_face = 0.5 * (J_thick_L + J_thick_R)
F_adv = (4/3) W^2 V^d J_face
F_diff = W D_face * (gamma^di grad_i J_face + V^d V^i grad_i J_face)
F_asym = F_adv - F_diff.
```

The asymptotic coordinate energy flux is

```text
f_asym_tilde = sqrt_detgamma
               * (alpha F_asym - beta^d max(E_star, E_floor)).
```

The published corrected energy flux is

```text
f_tilde = a * f_input_tilde + (1 - a) * f_asym_tilde.
```

Thus the input numerical flux dominates at small optical depth and the
thick-limit asymptotic flux dominates as `tau_face` grows. The helper validates
the spatial metric, finite gradients, positive length, and face velocity
normalization before publishing a corrected value; active-path arithmetic is
also completed and checked before either output is published.

## Scope of the correction

This helper corrects only the energy flux. It does not replace the five
component four-point transport operation, does not update number flux, and
does not provide a coupled E/F/N diffusion scheme. It also does not alter
source terms or the local implicit update.

The current four-point implementation explicitly rejects a request to combine
this correction with its own face blend. A host must therefore choose and
document an explicit downstream method boundary; silently applying both would
make the stored flux and its diagnostics impossible to interpret.

## Geometry and current policy

`delta_l` is a proper face-normal length, distinct from the four-point
`delta_x` coordinate spacing. See the
[geometry and wave-speed leaf](m1-geometry-and-wave-speed-policy.md). The
diffusion path checks that `W_face` agrees with the metric norm of `V_face` and
that the velocity remains subluminal. The current helper uses
`rates_face->kappa_tr` for `chi_tr`; it does not apply the canonical
four-point opacity suppression factor.

## Evidence and historical context

The focused implementation check is
[`unit_test_m1_diffusion_flux.c`](../../../Unit_Tests/unit_test_m1_diffusion_flux.c).
That test replays the six existing seeded radiation baseline/perturbed RNG
pairs from [`jthick_thcm1.fixture`](../../../Unit_Tests/data/jthick_thcm1.fixture).
Each endpoint is evaluated by GRHayL and compared directly with the frozen
THC_M1 `calc_Pthick` result projected to comoving `J`; the paired replay also
compares `G(x1)-G(x0)` with `T(x1)-T(x0)` using the existing fixture policy.
This is scalar thick-limit evidence, not evidence for full-grid thick-limit
convergence. The checked-in replay agrees for all 12 endpoint comparisons and
all 6 paired-response comparisons: six baseline and six perturbed states, using
the existing `pointwise_a1_a2_v1` policy. The optional
Fick helper is not cross-code compared because THC_M1 has no one-to-one public
operation for that caller-prepared corrected flux; its local boundary and
arithmetic checks remain the appropriate coverage. The shared source and API are
[`ghl_m1_diffusion_flux.c`](../../../GRHayL/Radiation/ghl_m1_diffusion_flux.c),
[`ghl_m1_Jthick.c`](../../../GRHayL/Radiation/ghl_m1_Jthick.c), and
[`ghl_m1.h`](../../../GRHayL/include/ghl_m1.h).

The methods whitepaper's photon-era speed-cap and diffusion discussion is
valuable background for why these helpers exist, but current neutrino callers
must follow the separate canonical-path restriction in the
[integration contract](../../../GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md).
