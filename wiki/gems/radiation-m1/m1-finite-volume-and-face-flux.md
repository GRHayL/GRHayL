# M1 finite-volume and face-flux boundary

This leaf describes the numerical boundary between a downstream finite-volume
driver and the host-neutral Radiation M1 kernels. It is about data ownership
and flux preparation; the four-point blending arithmetic is described in the
[four-point Rusanov leaf](m1-four-point-blended-rusanov.md).

## What the library evolves

For one grey neutrino species, the undensitized local state is ordered
`{N, E, Fx, Fy, Fz}`:

| component | meaning | storage convention |
| --- | --- | --- |
| `N` | Eulerian radiation number density | scalar |
| `E` | Eulerian radiation energy density | scalar |
| `F_i` | covariant Eulerian energy-flux components | lower spatial index |

The public order is fixed by
`ghl_m1_neutrino_transport_component_t` in
[`ghl_m1.h`](../../../GRHayL/include/ghl_m1.h). A host may store conservative
variables with a volume or `sqrt(det(gamma))` factor, but every local transport
kernel documents whether its inputs and outputs are weighted. Do not infer
that convention from a neighboring API.

## Physical fluxes

The caller supplies undensitized physical fluxes. In coordinate direction
`d`, the energy and momentum fluxes are

```text
f_E^d   = alpha F^d - beta^d E
f_F_i^d = alpha P^d_i - beta^d F_i
```

The neutrino number current supplies the corresponding number flux `n^d`, so

```text
f_N^d = alpha n^d - beta^d N.
```

These are coordinate physical fluxes, not numerical face fluxes. The shared
physical-flux implementation is
[`ghl_m1_compute_physical_flux`](../../../GRHayL/Radiation/ghl_m1_rusanov_flux.c),
and number-flux construction is implemented in
[`ghl_m1_neutrino_number_flux.c`](../../../GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_number_flux.c).
The caller must provide a closure consistent with each cell state before it
constructs the E/F fluxes.

## Pointwise four-point call

The canonical pointwise operation is
[`ghl_m1_compute_neutrino_four_point_transport_flux`](../../../GRHayL/include/ghl_m1.h).
The host assembles and owns:

- four neighboring, undensitized states in stencil order
  `{j-1, j, j+1, j+2}`;
- the undensitized physical fluxes at `j` and `j+1`;
- the two adjacent nonnegative, uncapped metric light-cone speed scalars;
- face opacity, coordinate spacing, face metric, and configured limiter values.

The kernel applies the low Rusanov flux, componentwise four-point limiter,
opacity suppression, and high/low blend to all five components. It then
multiplies the final face result by the face
`sqrt_detgamma` exactly once. The output is a densitized face flux in the same
five-component order. The host, not the kernel, takes the flux divergence.

The operation is pointwise and has no grid loop, reconstruction, ghost-zone,
boundary, AMR, or time-integrator knowledge. These ownership rules are part of
the [M1 integration contract](../../../GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md).

## Variable-volume prepared-face call

When the finite-volume discretization gives the four stencil cells and the
face different volume factors, prepare the nonlinear operands before calling
[`ghl_m1_compute_neutrino_four_point_volume_weighted_transport_flux`](../../../GRHayL/include/ghl_m1.h):

```text
U_prepared[j]   = V_cell[j] * U[j]
F_prepared[L/R] = V_face    * F[L/R]
```

Use one consistent `V_face` for both physical flux operands. This API performs
the same limiter, opacity suppression, and Rusanov blend, but it does not read
a metric and does not apply another `sqrt_detgamma`. Its result is already in
the caller's prepared face units.

The host must check that volumes and all products are finite and positive
where required. It must not pass prepared operands to the pointwise API, and
must not multiply the returned prepared flux by a second metric or volume
factor. The active variable-volume fixture adapter follows exactly this rule;
see [`README.m1.md`](../../../Unit_Tests/README.m1.md) and
[`unit_test_m1_thcm1_blended_rusanov.c`](../../../Unit_Tests/unit_test_m1_thcm1_blended_rusanov.c).

When all cell volume factors equal the face factor, the prepared operands share
a common volume scale. The prepared result carries that face volume factor,
while the pointwise result carries `sqrt_detgamma`. Their returned fluxes are
directly comparable only when those factors match; otherwise convert them to
the same units first.
Keeping preparation in the host makes the mesh measure, storage convention,
and finite-volume update explicit instead of hiding them in a Radiation API.

## Stage flow owned by the host

A typical explicit transport stage is:

1. reconstruct admissible cell states and repair the E/F realizability cone;
2. compute closures, comoving moments, number currents, and physical fluxes;
3. compute metric light-cone speed inputs and face opacity;
4. assemble four-point or prepared-face operands and call the canonical face
   kernel;
5. apply the returned face fluxes in the host's finite-volume divergence;
6. add geometry/source contributions according to the host's stage policy; and
7. pass the transport-predicted state to the local neutrino source update.

The list is a boundary description, not a shipped driver. Reconstruction,
RK/IMEX staging, mesh traversal, face sharing, boundary conditions, matter
recovery, and rate refresh remain downstream responsibilities. The two active
build manifests show which kernels are actually compiled:
[`Radiation/make.code.defn`](../../../GRHayL/Radiation/make.code.defn) and
[`Radiation/Neutrinos/make.code.defn`](../../../GRHayL/Radiation/Neutrinos/make.code.defn).

## Do not combine incompatible conventions

The canonical neutrino face path is four-point blended Rusanov with metric
light-cone speeds. A generic two-state Rusanov helper exists for shared
callers, and legacy HLL helpers remain part of the shared surface, but neither
defines the canonical neutrino four-point operation. An optical-depth speed
cap or separate diffusion correction must not be silently inserted into this
call. The four-point API rejects its diffusion-policy flag when enabled.

The [current integration contract](../../../GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md)
and the source implementation
[`ghl_m1_four_point_blended_rusanov.c`](../../../GRHayL/Radiation/ghl_m1_four_point_blended_rusanov.c)
are authoritative when a historical methods document describes a different
HLL, cap, or densitization sequence.
