# Host stage and volume-weighted integration

Radiation M1 is a host-neutral collection of local and face operations. The
downstream finite-volume or method-of-lines driver owns the mesh, storage,
reconstruction, stage order, boundaries, and matter update. The library does
not ship a three-species grid driver.

The authoritative boundary is
[`M1_INTEGRATION_CONTRACT.md`](../../../GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md);
the public declarations are in
[`ghl_m1.h`](../../../GRHayL/include/ghl_m1.h).

## One transport/source stage

A host stage should make the following objects explicit:

| object | owner | requirement |
| --- | --- | --- |
| cell and face arrays, ghost zones, mesh measures | downstream host | no Radiation-owned mesh/grid state or traversal; process-wide closure counters and failure snapshots are diagnostics only |
| reconstructed four-state stencil | downstream host | states are undensitized and in `{N,E,Fx,Fy,Fz}` order |
| closure, moments, number current, physical flux | local Radiation calls and host assembly | inputs must share a coherent metric and repaired state |
| provider rate bundle | provider/host boundary | one validated frozen bundle per species and stage |
| face numerical flux | Radiation face kernel | pointwise result is densitized exactly once, or prepared result is already volume-weighted |
| finite-volume divergence and transport predictor | downstream host | produces the `state_transport` source base |
| local source update and exchange packet | Radiation pointwise source call | no EOS, `Con2Prim`, matter recovery, or rate refresh inside the call |
| multi-species limiter and publication | downstream host | one common scalar limiter before matter recovery/publication |

The practical sequence is:

1. obtain or refresh provider rates at the host's declared stage boundary;
2. undensitize cell data and reconstruct admissible states;
3. repair E/F states before closure and physical-flux evaluation;
4. compute cell closures, comoving moments, number currents, and physical
   E/F/N fluxes;
5. build coherent face geometry, speed inputs, and transport opacity;
6. call the canonical four-point face operation;
7. apply face fluxes and geometry terms through the host's finite-volume/RK
   update to form `state_transport`;
8. call the local source update for each species, using the transport state as
   its source base; and
9. aggregate exchange packets, apply the common admissibility interval, recover
   matter, and publish only a fully accepted stage.

This is an integration pattern, not an API that GRHayL executes for the host.
The [finite-volume face-flux leaf](m1-finite-volume-and-face-flux.md) covers
the face operands and densitization; the
[source-branch leaf](m1-source-update-branches-and-rollback.md) covers local
source publication.

## Pointwise versus prepared faces

Use
[`ghl_m1_compute_neutrino_four_point_transport_flux`](../../../GRHayL/include/ghl_m1.h)
when the face metric supplies the one pointwise `sqrt_detgamma` factor. Its
state and physical-flux operands are undensitized and its output is a
densitized face flux.

Use
[`ghl_m1_compute_neutrino_four_point_volume_weighted_transport_flux`](../../../GRHayL/include/ghl_m1.h)
when the host's finite-volume measure is part of the nonlinear stencil
operation. Prepare every stencil state and both physical fluxes as

```text
U_prepared[j]   = V_cell[j] * U[j]
F_prepared[L/R] = V_face    * F[L/R].
```

The prepared API performs no metric inspection and no additional
densitization. It is the host's responsibility to validate finite positive
volumes and products, use the same face volume on both sides, and insert the
returned prepared flux into its own divergence. Passing weighted operands to
the pointwise API is a double-weighting error.

The varying-volume fixture adapter demonstrates this preparation in
[`unit_test_m1_thcm1_blended_rusanov.c`](../../../Unit_Tests/unit_test_m1_thcm1_blended_rusanov.c);
the input layout and stored-reference boundary are described in
[`README.m1.md`](../../../Unit_Tests/README.m1.md).

## Species and pair scheduling

The grey radiation state has three species: `nue`, `anue`, and already-summed
`nux`. The provider returns validated frozen bundles in that species order.
When electron-flavor pair processes are active, call the paired
`nue/anue` source operation together. It keeps pair candidates and exchange
packets private until both species succeed. The heavy-flavor species uses the
single-species path.

After local success, aggregate the three radiation increments and all matter
and lepton packets. Intersect the radiation realizability/positivity interval
with matter energy-momentum, electron-fraction, EOS, and `Con2Prim` intervals.
Apply one scalar `theta` to every species and all coupled increments, then do
temporary matter recovery. A failure leaves all caller outputs unchanged. The
library does not own this multi-species limiter or publication step.

## Provider and cache boundary

The rate provider owns EOS/table lookup, unit conversion, weak-equilibrium
targets, channel selection, bounds/failure policy, and cache contents. The
host chooses refresh timing and passes the resulting frozen bundle into local
Radiation calls. Provider cache and neutrino diagnostics are caller-owned
records; do not share mutable records across concurrent calls without an
external synchronization policy.

The provider API is declared in
[`ghl_neutrino_rate_provider.h`](../../../GRHayL/include/ghl_neutrino_rate_provider.h)
and its implementation is
[`ghl_neutrino_rate_provider.c`](../../../GRHayL/Radiation/Neutrinos/ghl_neutrino_rate_provider.c).
Both shared and neutrino source manifests are listed in
[`Radiation/make.code.defn`](../../../GRHayL/Radiation/make.code.defn) and
[`Radiation/Neutrinos/make.code.defn`](../../../GRHayL/Radiation/Neutrinos/make.code.defn).

## What the host must not assume

The local calls do not provide AMR refluxing, boundary conditions, global CFL
selection, mesh connectivity, stage counters, matter evolution, or framework
integration. A successful pointwise result proves only that local inputs met
the operation's contract. The host must supply the missing global pieces and
must record which source, transport, volume, and publication conventions it
uses.

The implementation and ownership details are maintained in the
[current integration contract](../../../GRHayL/Radiation/M1_INTEGRATION_CONTRACT.md),
not in the historical photon implementation whitepaper's old sixteen-function
API list.
