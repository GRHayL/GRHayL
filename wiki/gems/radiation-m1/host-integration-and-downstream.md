# Radiation M1 Host Boundary

`GRHayL/Radiation` is a host-neutral, pointwise library. It does not own mesh
loops, reconstruction, ghost zones, boundaries, AMR, schedules, stage
sequencing, matter variables, or `Con2Prim`.

## Host-owned transport

The host stores or converts its evolved variables, assembles the four-point
stencil, prepares the two adjacent physical fluxes and uncapped speeds, and
calls `ghl_m1_compute_neutrino_four_point_transport_flux` for every canonical
neutrino M1 face. The operation returns all five components in
`{N,E,Fx,Fy,Fz}` and densitizes the face flux exactly once. The host performs
the final divergence. The operation always uses the metric light-cone speeds,
the four-point blended Rusanov formula, and no separate diffusion correction.

The generic two-state Rusanov route remains available to other shared callers,
but it is not the canonical neutrino M1 transport route. No host policy choice
can replace the canonical face operation with an optical-depth speed cap or a
separate diffusion correction.

## Host-owned source stage

For each species, the host supplies frozen primitives, rates, the pre-transport
state, the transport-predicted state, `dt`, and the conserved baryon-density
normalization. `ghl_m1_solve_neutrino_source_update` returns a temporary state,
an exchange packet, and diagnostics. It never updates matter or tracks a stage.

`dt` is the coordinate-time stage timestep. The library applies the metric
lapse once, as `dt_alpha = alpha * dt`, for source-policy decisions, predictors,
explicit updates, and the implicit solve. The host must not pre-apply the lapse.
Before first use, the host must initialize its caller-owned neutrino diagnostics
and provider cache with their public initializers; one cache/diagnostics record
must be used per independent cell or thread unless access is synchronized.

The host must:

1. solve all three species before publishing any result;
2. reject provider or solve failures without a partial update;
3. aggregate the packets and compute one admissible limiter interval;
4. intersect the interval with EOS and recovery restrictions;
5. apply the same `theta` to every radiation, matter, and `Y_e` increment;
6. recover matter primitives on a temporary candidate; and
7. publish the coupled result only after all checks succeed.

`dN_rad_total` is distinct from `dL_rad_cc`. The canonical `Y_e` increment uses
the charged-current exchange; total radiation-number exchange remains a
separate diagnostic and is not folded into `Y_e`.

The local implicit source update retains the existing finite-difference
Jacobian and Newton solver, including its line search and fallback substepping.
Any source-regime handling is internal to that source-update contract and does
not select a different transport, repair, or closure method.

## Host implementation boundary

No host driver is included in this checkout. A downstream host is responsible
for the species loop, packet aggregation, admissibility interval, temporary
matter recovery, rollback, and one shared limiter scalar described above. Those
responsibilities remain outside the GRHayL build and dependency boundary.

The library documentation and focused checks therefore make no grid-level,
schedule-level, AMR, complete-evolution, or physical-validation claim.

## Local contract sources

- [`M1_INTEGRATION_CONTRACT.md`](../../../GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md)
- [`ghl_m1.h`](../../../GRHayL/include/ghl_m1.h)
