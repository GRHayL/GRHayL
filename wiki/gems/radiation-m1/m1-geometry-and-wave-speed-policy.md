# M1 geometry and wave-speed policy

The current neutrino transport speed policy is deliberately simple: use the
metric light cone at the face, without an optical-depth cap. Geometry is still
part of the numerical contract because the lapse, shift, spatial metric,
inverse metric, determinant, coordinate spacing, and proper normal spacing
must not be mixed accidentally.

## Coherent face geometry

Every face operation must use one coherent
[`ghl_metric_quantities`](../../../GRHayL/include/ghl_m1.h) record. The
current kernels validate that the spatial metric is symmetric positive
definite and expect `gammaDD` and `gammaUU` to be mutually consistent. A
pointwise face flux also requires a finite positive `sqrt_detgamma`.

The face record supplies:

- `alpha`, the lapse;
- `betaU[i]`, the contravariant shift;
- `gammaDD[i][j]`, the spatial metric;
- `gammaUU[i][j]`, its inverse; and
- `sqrt_detgamma`, the spatial-volume factor used for one pointwise flux
  densitization.

Do not average `gammaDD`, `gammaUU`, and `sqrt_detgamma` independently if that
breaks inverse consistency. The face metric is also the metric used by the
physical flux, speed, and optional diffusion calculations at that face.

## Raw light-cone speeds

For coordinate direction `d`,
[`ghl_m1_compute_raw_lightcone_speeds`](../../../GRHayL/Radiation/ghl_m1_wavespeeds.c)
returns the signed coordinate speeds

```text
s_minus_raw = -beta^d - alpha * sqrt(gamma^dd)
s_plus_raw  = -beta^d + alpha * sqrt(gamma^dd).
```

The public `ghl_m1_compute_wavespeeds` wrapper additionally clips these
endpoints to the HLL envelope with `min(0, s_minus_raw)` and
`max(0, s_plus_raw)`. That clipping helper is useful to shared compatibility
callers; it is not the canonical neutrino four-point interface.

The four-point API receives the adjacent nonnegative speed scalars accepted by
its contract and takes their maximum for the Rusanov dissipation. A host must
map its raw light-cone estimate to that scalar consistently, without applying
an optical-depth cap. The current four-point implementation rejects negative,
nonfinite speed inputs and records the selected maximum in its diagnostics.

There is no closure-dependent M1 eigenvalue solve in this transport path. The
full M1 closure is still needed to construct the physical E/F fluxes, but it
does not replace the metric light-cone speed policy.

## Coordinate spacing versus proper normal spacing

The four-point opacity suppression uses `delta_x`, the coordinate spacing in
the selected direction:

```text
tau_four_point = kappa_face * delta_x.
```

It is not the proper face-normal thickness. When a proper normal length is
needed, the current helper
[`ghl_m1_compute_face_normal_delta_l`](../../../GRHayL/Radiation/ghl_m1_utils.c)
uses the face geometry to provide the corresponding `delta_l`. In the
coordinate-direction convention used by the methods contract, the geometric
relationship is

```text
delta_l = delta_x / sqrt(gamma^dd).
```

Use `delta_x` for the four-point limiter's opacity factor and `delta_l` for
the optional thick-limit diffusion optical depth. Substituting one for the
other changes the numerical method on a non-unit or non-Euclidean metric.

## Densitization and index placement

Radiation states use covariant `F_i`, while the physical energy flux uses the
raised component `F^d = gamma^di F_i`. The momentum physical flux uses the
mixed pressure component `P^d_i`. The shift terms are applied before
densitization:

```text
f_E^d   = alpha F^d - beta^d E
f_F_i^d = alpha P^d_i - beta^d F_i.
```

The pointwise four-point operation then multiplies the completed numerical
flux by `sqrt_detgamma` once. A host using volume-prepared operands follows a
different contract and must not add this factor in the library; see the
[finite-volume leaf](m1-finite-volume-and-face-flux.md).

## What is intentionally not combined

The current canonical sequence is:

1. construct repaired states and closures;
2. construct physical E/F and number fluxes;
3. derive adjacent metric light-cone speed inputs;
4. call four-point blended Rusanov; and
5. apply the returned flux in the host divergence.

An optical-depth speed cap is not inserted between steps 3 and 4. A separate
diffusion correction is not automatically composed with the four-point
operation; requesting that policy through the four-point API is rejected. The
optional diffusion helper has its own proper-length, face-velocity, and
thick-limit contract, documented in the
[thick-limit leaf](m1-thick-limit-and-optional-diffusion.md).

The older photon methods whitepaper discusses capped speeds and HLL-oriented
choices. Those are historical context only where they conflict with the
current [M1 integration contract](../../../GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md).

## Implementation and evidence

The shared speed and physical-flux kernels are compiled by
[`Radiation/make.code.defn`](../../../GRHayL/Radiation/make.code.defn); the
neutrino transport wrapper is compiled by
[`Radiation/Neutrinos/make.code.defn`](../../../GRHayL/Radiation/Neutrinos/make.code.defn).
The public declarations and comments in
[`ghl_m1.h`](../../../GRHayL/include/ghl_m1.h) are the interface authority.

Focused local checks include
[`unit_test_m1_neutrino_rusanov_flux.c`](../../../Unit_Tests/unit_test_m1_neutrino_rusanov_flux.c),
[`unit_test_m1_thcm1_blended_rusanov.c`](../../../Unit_Tests/unit_test_m1_thcm1_blended_rusanov.c),
and the generic shared flux test
[`unit_test_rusanov_flux.c`](../../../Unit_Tests/unit_test_rusanov_flux.c).
They exercise local operations and fixture comparisons; they do not prove a
complete downstream mesh evolution or a continuum convergence result.
