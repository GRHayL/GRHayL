# Radiation M1 Traceability

This page records only repository-local implementation evidence. A
declaration is part of the production surface only when its implementation is
listed by an active make manifest. A test is evidence for the behavior it
actually exercises; it is not evidence of a complete host evolution or
physical validation.

## Production boundary

The shared M1 implementation is listed in
`GRHayL/Radiation/make.code.defn`. The neutrino implementation is listed in
`GRHayL/Radiation/Neutrinos/make.code.defn`. The installed declarations are in
`GRHayL/include/ghl_m1.h` and `GRHayL/include/ghl_radiation.h`.

The canonical M1 production surfaces are:

| Public surface | Implementation | Contract |
| --- | --- | --- |
| `ghl_m1_compute_closure_with_primitives` and its Minerbo helpers | `GRHayL/Radiation/ghl_m1_closure.c` and `GRHayL/Radiation/ghl_m1_utils.h` | Primary direct four-dimensional closure; every published tensor obeys the radiation trace identity. Finite non-PSD candidates and invalid exact-zero-flux tensors use a flagged Eulerian Minerbo admissibility fallback. |
| `ghl_m1_compute_neutrino_four_point_transport_flux` and its limiter, opacity, and blend helpers | `GRHayL/Radiation/ghl_m1_four_point_blended_rusanov.c` | Canonical four-point host-prepared face operation; metric light-cone speeds; `{N,E,Fx,Fy,Fz}`; exactly-once face densitization; no separate diffusion correction. |
| `ghl_m1_solve_neutrino_source_update` and `ghl_m1_try_neutrino_explicit_thin_update` | `GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_source_update.c` and the existing neutrino source/implicit files | Frozen-input, transactional source update; transport state is the source base; existing finite-difference/Newton implicit solver; separate total-number and charged-current exchange. |
| `ghl_m1_solve_neutrino_pair_source_update` | `GRHayL/Radiation/Neutrinos/ghl_m1_neutrino_pair_source.c` | Independent charged-current/scattering stage followed by a shared electron-pair number reaction and grey E/F solve; equal pair number increments, zero pair electron-fraction exchange, and atomic two-species publication. The approximation is specified in `PAIR_SOURCE_MODEL.md`. |
| `ghl_neutrino_rate_provider_initialize_nrpyleakage` and `ghl_neutrino_rate_provider_compute_cell` | `GRHayL/Radiation/Neutrinos/ghl_neutrino_rate_provider.c` and the private `ghl_m1_nrpyleakage_kernel.c` adapter | Production table-backed NRPyLeakage rates with separated electron pair number/energy emissivities, aggregate heavy-flavor rates, stable finite-tail evaluation, positive equilibrium targets before detailed-balance division, and final bundle validation. The adapter consumes the NRPyLeakage-owned source-private blocking, beta-moment, reaction-shift, and bremsstrahlung helpers; it does not copy those formulas or expose them as installed API. |
| Prepared-transport low-flux wrapper | Private `ghl_m1_compute_neutrino_four_point_rusanov_flux` in `GRHayL/Radiation/ghl_m1_four_point_blended_rusanov.c` | Componentwise Rusanov wrapper for the volume-weighted operands of the prepared four-point face operation; it remains private because its operand units differ from the undensitized public helper. It shares the unit-agnostic scalar arithmetic in `GRHayL/Flux_Source/ghl_calculate_Rusanov_flux.c` through the private `GRHayL/Flux_Source/ghl_rusanov_private.h` declaration. |

Other source files present under `GRHayL/Radiation/` are not production code
unless their manifest lists them. In particular, a parked or diagnostic file
does not become a public implementation merely because a declaration exists.

## Claim boundary

The current evidence supports the listed installed API operations. It does not
support claims about a particular grid framework, schedule, AMR
implementation, production downstream host, or full-evolution equivalence. It
also does not turn a table-backed or table-free client check into physical
validation; downstream consumers own any such campaign.

## NRPyLeakage/M1 helper boundary

The private M1 adapter is permitted to include the canonical NRPyLeakage
source-private helpers directly:

- [`NRPyLeakage_nucleon_blocking.h`](../Neutrinos/NRPyLeakage/NRPyLeakage_nucleon_blocking.h)
- [`NRPyLeakage_rate_helpers.h`](../Neutrinos/NRPyLeakage/NRPyLeakage_rate_helpers.h)
- [`M1 NRPyLeakage manifest`](Neutrinos/make.code.defn)
- [`M1 NRPyLeakage adapter`](Neutrinos/ghl_m1_nrpyleakage_kernel.c)

This is a tracked cross-gem source dependency, not a new public interface: the
helpers remain NRPyLeakage-owned, are not duplicated in Radiation, and are not
installed through `GRHayL/include`. The M1 rate bundle remains tau-free. In
particular, NRPyLeakage optical-depth suppression, leakage luminosity/source
assembly, and legacy leakage finite-output fallback policies remain outside
the M1 adapter; M1 retains its own validation and provider recovery contract.
