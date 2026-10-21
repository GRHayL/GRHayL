# M1 public API quick reference

This leaf indexes the public shared-M1 surface declared in
[`GRHayL/include/ghl_m1.h`](../../../GRHayL/include/ghl_m1.h): initialization
and solver/closure controls, and the shared helper kernels the neutrino path
and hosts may call. The installed header remains the complete declaration
surface; this page routes it. Neutrino source/transport
entry points are covered by the [neutrino M1 contract](neutrino-m1-contract.md) and
the [Radiation M1 hub](../radiation-m1.md); rate-provider operations are owned by the
[rate-provider contract](rate-provider-contract.md).

## Source ownership

`ghl_m1.h` is an installed public header listed in
[`GRHayL/include/make.code.defn`](../../../GRHayL/include/make.code.defn).
The shared-M1 implementations live in `GRHayL/Radiation/` and are compiled
through [`GRHayL/Radiation/make.code.defn`](../../../GRHayL/Radiation/make.code.defn);
the [API and build boundary](api-build-boundary.md) gives the full manifest
routing. The private headers `ghl_m1_closure_private.h`,
`ghl_m1_rusanov_private.h`, and `ghl_m1_utils.h` are build inputs, not
installed declarations, and Radiation-private `_validated` kernel variants
are not part of the installed surface.

## Parameter initialization and solver/closure controls

- `ghl_m1_initialize` and `ghl_m1_initialize_with_newton_tolerances`
  initialize the canonical M1 method bundle: metric light-cone transport
  speeds, the squared-ratio realizability rescale, and the primary full
  four-dimensional closure with its flagged pressure/covariance admissibility
  fallback.
- `ghl_m1_set_newton_tolerances` updates the mixed relative/absolute Newton
  tolerances on an initialized bundle; `ghl_m1_set_closure_solver_controls`
  sets the bracketed scalar Minerbo root-solver interval tolerance and
  iteration bound; `ghl_m1_set_closure_residual_tolerance` sets the maximum
  normalized Minerbo consistency residual a finite candidate may publish.

Initialization and the parameter setters validate their parameter values.
Callers supply coherent metric data and initialized parameter structs;
individual kernels perform the local checks described by their contracts,
without a universal per-call validation of either complete struct. The
canonical transport controls in the bundle
(`minmod_theta` and `mindiss`) are defined by the
[four-point blended Rusanov](m1-four-point-blended-rusanov.md) leaf.

## Shared public helper surface

- Realizability: `ghl_m1_apply_energy_floor` and
  `ghl_m1_realizability_repair` — see
  [realizability repair](m1-realizability-repair.md).
- Closure: `ghl_m1_compute_closure_with_primitives` is the public closure
  entry point; `ghl_m1_compute_closure_decomposition_diagnostic` is
  diagnostic-only — see
  [four-dimensional Minerbo closure](m1-four-dimensional-minerbo-closure.md).
- Moments and stress: `ghl_m1_compute_comoving_moments`,
  `ghl_m1_compute_stress_energy`, and `ghl_m1_compute_diagnostics` — see
  [comoving moments and source projection](m1-comoving-moments-and-source-projection.md)
  and [diagnostics and verification](m1-diagnostics-and-verification.md).
- Sources: `ghl_m1_compute_geometry_sources` and
  `ghl_m1_compute_matter_coupling_sources` — see
  [3+1 radiation equations](m1-3plus1-radiation-equations.md).
- Wave speeds and face geometry: `ghl_m1_compute_raw_lightcone_speeds`,
  `ghl_m1_clip_hll_speeds`, `ghl_m1_compute_wavespeeds`, and
  `ghl_m1_compute_face_normal_delta_l` — see
  [geometry and wave-speed policy](m1-geometry-and-wave-speed-policy.md).
- Fluxes: `ghl_m1_compute_physical_flux`, `ghl_m1_compute_rusanov_flux`, the
  undensitized component helper `ghl_calculate_Rusanov_flux`, and
  `ghl_m1_compute_number_rusanov_flux` — see
  [finite volume and face flux](m1-finite-volume-and-face-flux.md) and the
  [API and build boundary](api-build-boundary.md) for the Radiation-owned
  Rusanov boundary against Flux_Source HLLE.
- Four-point transport stages: `ghl_m1_compute_neutrino_four_point_transport_flux`,
  `ghl_m1_compute_neutrino_four_point_volume_weighted_transport_flux`, plus
  the componentwise stages `ghl_m1_compute_four_point_flux_limiter`,
  `ghl_m1_compute_four_point_opacity_suppression`, and
  `ghl_m1_compute_four_point_blended_flux` — see
  [four-point blended Rusanov](m1-four-point-blended-rusanov.md) and
  [host stages and volume-weighted integration](m1-host-stage-and-volume-weighted-integration.md).
- Optional thick-limit/diffusion helpers: `ghl_m1_compute_Jthick`,
  `ghl_m1_compute_harmonic_diffusion_coefficient`, and
  `ghl_m1_compute_diffusion_flux` — see
  [thick limit and optional diffusion](m1-thick-limit-and-optional-diffusion.md).
  These are not selected by the canonical four-point transport operation.
- Newton solver driver: `ghl_m1_newton_solve_4d` with the
  `ghl_m1_newton_callbacks` bundle and optional stage observer, the
  `ghl_m1_newton_weighted_merit` mixed merit, and
  `ghl_m1_newton_project_admissible` — see
  [implicit Newton and failure policy](m1-implicit-newton-and-failure-policy.md).

## Header-inline neutrino adapters

[ghl_m1.h](../../../GRHayL/include/ghl_m1.h) also supplies caller-compiled
adapters for the shared E/F kernels. `ghl_m1_neutrino_project_rad_state` copies
E/F from a non-null neutrino state; it does not copy or repair `N`. The
status-returning wrappers reject a null state before calling the corresponding
shared operation, whose remaining prerequisites and output contract apply.

| Neutrino adapter | Shared operation |
| --- | --- |
| `ghl_m1_compute_neutrino_closure` | `ghl_m1_compute_closure_with_primitives` |
| `ghl_m1_compute_neutrino_comoving_moments` | `ghl_m1_compute_comoving_moments` |
| `ghl_m1_compute_neutrino_stress_energy` | `ghl_m1_compute_stress_energy` |
| `ghl_m1_compute_neutrino_geometry_sources` | `ghl_m1_compute_geometry_sources` |
| `ghl_m1_compute_neutrino_Jthick` | `ghl_m1_compute_Jthick` |
| `ghl_m1_compute_neutrino_diagnostics` | `ghl_m1_compute_diagnostics` |

## Prerequisites and output/failure boundaries

- Follow each entry point's header contract for required pointers, directions,
  scalar inputs, state and closure validation. Hosts supply coherent metric
  quantities and an initialized parameter bundle; a supplied closure must
  correspond to the metric, primitives, and E/F state being consumed.
- Kernels that consume a supplied closure, including the neutrino
  `..._from_closure` number-current/source variants and shared moment, physical
  flux, and stress helpers, check finite pressure entries, trace, and symmetry
  through `ghl_m1_validate_closure_tensor` in
  [ghl_m1_utils.h](../../../GRHayL/Radiation/ghl_m1_utils.h). They do not repeat
  construction's pressure or Eulerian covariance positivity guards. Supply a
  closure produced by `ghl_m1_compute_closure_with_primitives`, or equivalently
  validated, for the consumed inputs. Other state and scalar checks remain
  operation-specific.
- Publication and diagnostic behavior are operation-specific. Pointwise
  arithmetic helpers retain their declared unchanged-output error contracts.
  The [Newton driver](../../../GRHayL/Radiation/ghl_m1_newton.c) can update
  diagnostics on failure while publishing its state only on success. The
  [neutrino source dispatcher](neutrino-m1-contract.md#failure-and-fallback-behavior)
  instead initializes its state output to the source base and its exchange
  packet to zero. A blanket unchanged-output rule does not cover those records.
- A completed `ghl_m1_compute_Jthick` call can still report an unusable
  thick-limit candidate through its validity flag; success means the
  arithmetic completed, not that the optional diffusion route is available.
- A finite full-four-dimensional closure candidate that fails pressure or Eulerian covariance
  positivity is replaced by the built-in Eulerian Minerbo admissibility
  fallback flagged in the closure record; the fallback is diagnostic-only
  metadata, not a user-selectable method.

## Test routes

The shared helper surface is exercised by
[`unit_test_m1_error_handling.c`](../../../Unit_Tests/unit_test_m1_error_handling.c)
(initialization, closure, repair, stress, diagnostics, geometry/matter
sources, invalid inputs, and unchanged-output rejection),
[`unit_test_m1_closure_fallback.c`](../../../Unit_Tests/unit_test_m1_closure_fallback.c)
(admissibility-fallback routes),
[`unit_test_m1_fd_jacobian.c`](../../../Unit_Tests/unit_test_m1_fd_jacobian.c)
(Newton driver, callback validation, projection, and merit behavior),
[`unit_test_m1_diffusion_flux.c`](../../../Unit_Tests/unit_test_m1_diffusion_flux.c)
(`Jthick` and diffusion helpers), and
[`unit_test_m1_neutrino_seeded_invariants.c`](../../../Unit_Tests/unit_test_m1_neutrino_seeded_invariants.c)
(pointwise closure/moments/stress/geometry/speeds; that owner holds the
exact-zero and seeded-state cases). The
[M1 test guide](../../../docs/raw/Radiation_unit_tests.md) and the
[tests and fixtures](tests-and-fixtures.md) leaf route the remaining
transport and fixture-replay owners.
