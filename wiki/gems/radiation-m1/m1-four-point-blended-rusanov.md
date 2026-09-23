# Four-point blended Rusanov transport

`ghl_m1_compute_neutrino_four_point_transport_flux` is the canonical
neutrino face operation. It combines a two-state symmetric Rusanov low flux
with a centered high flux using four stencil states. The implementation is
componentwise and applies the same policy to `N`, `E`, and each covariant
`F_i`; it does not construct an HLL intermediate state.

The current implementation is in
[`ghl_m1_four_point_blended_rusanov.c`](../../../GRHayL/Radiation/ghl_m1_four_point_blended_rusanov.c),
with the public contract and component ordering in
[`ghl_m1.h`](../../../GRHayL/include/ghl_m1.h).

## Inputs and high/low fluxes

The stencil is `{j-1, j, j+1, j+2}`. Physical fluxes `f_L` and `f_R` belong
to cells `j` and `j+1`, respectively. The scalar speed supplied for each
side is nonnegative; the operation uses

```text
a = max(speed_L, speed_R).
```

For each transport component `q`, the high and low candidates are

```text
f_high[q] = 0.5 * (f_L[q] + f_R[q])

f_low[q]  = 0.5 * (f_L[q] + f_R[q])
              - 0.5 * a * (U_R[q] - U_L[q]).
```

Here `U_L` and `U_R` are the middle stencil states `j` and `j+1`. The high
candidate is the centered physical-flux average; the low candidate adds
Rusanov dissipation using the state jump. The caller remains responsible for
computing physical E/F and number fluxes from repaired states and closures.

## Componentwise four-point limiter

For each component, form

```text
d_minus = U[j]   - U[j-1]
d_center = U[j+1] - U[j]
d_plus  = U[j+2] - U[j+1].
```

With `theta = m1_params->minmod_theta`, the limiter is:

- if all three differences have the same nonzero sign,

  ```text
  phi = min(1, theta*d_minus/d_center, theta*d_plus/d_center);
  ```

- if both outer differences have a nonzero sign opposite to `d_center`, set
  `sawtooth = true`;
- otherwise set `phi = 0` and leave `sawtooth = false`.

The implementation accepts `0 <= theta <= 2`. `theta = 0` yields `phi = 0`
for every stencil, selecting the fully dissipative low-order flux; it is
decided before the ratios are formed, so an overflowing `d_minus/d_center`
cannot turn `0 * inf` into a `NaN` that `fmin` would resolve to `phi = 1`.
A zero difference does not qualify as a same-sign or opposite-sign pattern,
so it takes the zero-limiter case. The decision is made independently for all five components; a sawtooth
in one component does not change another component's limiter.

## Opacity suppression and final blend

The face opacity and the coordinate spacing form the optical width

```text
tau = kappa_face * delta_x.
```

`delta_x` is coordinate spacing in the selected direction and must use units
consistent with the inverse-length opacity. The suppression factor is

```text
A = 1                                         if tau <= 1
    max(mindiss, min(1, 1/tau))                if tau > 1.
```

For the ordinary limited branch, the final component is

```text
f_num = f_high - A * (1 - phi) * (f_high - f_low).
```

For a sawtooth pattern, the dissipation factor is `1` rather than `A`:

```text
f_num = f_high - (1 - phi) * (f_high - f_low).
```

Thus smooth, optically thin data can retain the high flux. In the ordinary
limited branch, optically thick data has `A >= mindiss`, but the applied
Rusanov contribution is `A * (1 - phi)` and can therefore be smaller than
`mindiss` (or zero when `phi = 1`). The sawtooth branch uses the full
Rusanov difference factor. The diagnostic record reports each component's
`phi`, each sawtooth flag, the common opacity suppression, and the selected
face speed.

## Densitization and policy guards

The pointwise API validates a coherent positive face metric, computes the
undensitized blend, and multiplies every component by
`metric_face->sqrt_detgamma` exactly once. Its output is therefore a
densitized face flux. The volume-weighted companion performs the same
componentwise blend on host-prepared operands but performs no metric
inspection or additional densitization. See the
[finite-volume boundary](m1-finite-volume-and-face-flux.md).

The canonical operation rejects `diffusion_correction_enabled = true`. It
also does not apply an optical-depth wave-speed cap. Metric light-cone speed
inputs are supplied by the host; see the
[geometry and speed policy](m1-geometry-and-wave-speed-policy.md). Separate
diffusion is an explicitly requested helper and is documented in the
[optional diffusion leaf](m1-thick-limit-and-optional-diffusion.md).

The generic two-state Rusanov and HLL-related shared helpers are not a second
canonical neutrino route. In particular, do not replace this four-point
operation with a historical photon HLL recipe merely because the older
methods whitepaper describes one.

## Failure behavior

The implementation validates all stencil states, physical fluxes, speeds,
opacity, spacing, limiter controls, and intermediate products before
publishing output. The pointwise and prepared paths keep candidate values
local; an error leaves the output flux and optional diagnostics unchanged.
This matters to hosts that reuse an output buffer across faces.

Relevant focused evidence is in
[`unit_test_m1_thcm1_blended_rusanov.c`](../../../Unit_Tests/unit_test_m1_thcm1_blended_rusanov.c),
[`unit_test_m1_neutrino_rusanov_flux.c`](../../../Unit_Tests/unit_test_m1_neutrino_rusanov_flux.c),
and the stored transport fixtures under
[`Unit_Tests/data/m1_thcm1`](../../../Unit_Tests/data/m1_thcm1). The
[M1 test guide](../../../Unit_Tests/README.m1.md) explains the pointwise and
variable-volume evidence boundary.

## Source authority

Use the [M1 integration contract](../../../GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md),
the [shared build manifest](../../../GRHayL/Radiation/make.code.defn), and the
[neutrino build manifest](../../../GRHayL/Radiation/Neutrinos/make.code.defn)
to distinguish compiled current behavior from historical whitepaper plans.
