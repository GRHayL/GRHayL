# Neutrinos

> Router for NRPyLeakage neutrino opacities, source terms, luminosities, optical depths, nucleon blocking, CompOSE EOS adaptation, and reference tests.

| Page | Go here when... |
| --- | --- |
| [Neutrinos API And Data](api-and-data.md) | You need `ghl_radiation.h` structs, `NRPyLeakage_*` entry points, input/unit preconditions, error codes, or the HDF5/EOS call boundary. |
| [Physics Model And Sources](physics-model-and-sources.md) | You need grey leakage reactions, `nue`/`anue`/`nux` and `[0]`/`[1]` slot meanings, effective rates, Valencia source projections, or cgs suppression-ratio units. |
| [EOS Inputs And Chemical Potentials](eos-inputs-and-chemical-potentials.md) | You need the six-output tabulated EOS callback mapping, `mu_e`/`muhat` MeV conventions, `munu` vs. `mu_e-muhat`, or gradual-underflow requirements. |
| [Nucleon Blocking And EOS Conventions](nucleon-blocking-and-eos-conventions.md) | You need `NRPyLeakage_nucleon_blocking.h`, density-derived kinetic degeneracies, `expm1` overlap algebra, Kirchhoff reaction pairing, or BNS_NURATES/DD2 qualification results. |
| [Implementation Flow](implementation-flow.md) | You need the `GRHayL/Neutrinos/NRPyLeakage/` build list, minimal link sets, `EnsureFinite`/sanitizer paths, or per-file execution flows. |
| [Generator Provenance](generator-provenance.md) | You need the ancestral NRPy+ notebooks, `tmp_*` common-subexpression origin, hand-applied C corrections, or FDINT/ILEAS attribution. |
| [CompOSE EOS Conversion And Regularization](compose-eos-conversion-and-regularization.md) | You need offline StellarCollapse HDF5 conversion, the `sro-sly4-sna-141-regularized-v1` profile, table preflight checks, or raw vs. generated CompOSE decoding. |
| [CompOSE EOS Leakage Callback](compose-eos-leakage-callback.md) | You need runtime `(rho, Y_e, T)` to `(T, n_b, Y_q)` mapping, chemical-potential zero points, safe six-output callback rules, or CompOSE validation steps. |
| [Neutrinos Tests And Reference Data](tests-and-reference-data.md) | You need `nrpyleakage_main.h`, SLy4 table setup, binary comparison tolerances, table-free physics checks, or CI matrix commands. |

Back to [GRHayL Gems](../index.md).
