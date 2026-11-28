# Glossary

> Canonical spellings for recurring GRHayL code, physics, and documentation terms. · Status: confirmed

## Summary

Use these spellings when routing, citing, or writing KB leaves. Terms link to their defining pages where GRHayL defines them.

## Terms

| Term | Meaning |
| --- | --- |
| GRHayL | General Relativistic Hydrodynamic Library: C99 GRMHD library (`extern "C"` headers) of `GRHayL_Core` and physics gems; see [Public API Map](architecture/public-api-map.md). |
| chalice | Shared core module (`GRHayL/GRHayL_Core/`) defining structs, metrics, stress-energy helpers, EOS dispatch, errors, and I/O; see [Core](core/index.md). |
| gem | Modular physics subdirectory under `GRHayL/` (`Atmosphere`, `Con2Prim`, `EOS`, `Flux_Source`, `Induction`, `Neutrinos`, `Reconstruction`); see [Gems](gems/index.md). |
| ADM | Arnowitt-Deser-Misner 3+1 spatial metric $\gamma_{ij}$, lapse $\alpha$, shift $\beta^i$, and extrinsic curvature $K_{ij}$; see [Metric And ADM Helpers](core/metric-adm.md). |
| Valencia formulation | Flux-conservative GRMHD equations for $\rho_*$, $\tilde{\tau}$, $\tilde{S}_i$, and $\tilde{Y}_e = \rho_* Y_e$; see [Evolution Equation Map](physics/evolution-equation-map.md). |
| `BU` | Undensitized, $1/\sqrt{4\pi}$-rescaled $B^i$ in `ghl_primitive_quantities` ($\tilde{B}^i = \sqrt{\gamma} B^i$); see [Variables](physics/variables-and-conventions.md). |
| `Btilde` | Densitized magnetic field $\tilde{B}^i = \sqrt{\gamma} B^i$ used in Induction (`ghl_HLL_flux_with_Btilde`); see [Variables](physics/variables-and-conventions.md). |
| `smallb` | Fluid-frame magnetic four-vector $b^\mu$ and norm $b^2 = b_\mu b^\mu$ computed by Core helpers; see [Stress-Energy And Smallb](core/stress-energy-smallb.md). |
| `phitilde` | Densitized scalar electromagnetic gauge variable $\tilde{\Phi} = \sqrt{\gamma}\,\Phi$ in generalized Lorenz gauge; see [Gauge RHS](gems/induction/gauge-rhs.md). |
| `vU` | Transport three-velocity $v^i = u^i / u^0$ stored in `ghl_primitive_quantities`, distinct from Eulerian velocity; see [Velocity u0](core/velocity-u0.md). |
| `Con2Prim` | Conservative-to-primitive root-finding gem (Noble, Palenzuela, Newman, Font) for hybrid and tabulated EOS; see [Con2Prim](gems/con2prim/index.md). |
| Neural-network primitive guess | Optional tabulated Con2Prim retry guess (flag `enable_neural_net_c2p`, HDF5 group `grhayl_nn_c2p`); see [Neural-Network Guess](gems/con2prim/neural-network-primitive-guess.md). |
| `NRPyEOS` | Hybrid piecewise-polytrope and tabulated StellarCollapse EOS routines under `GRHayL/EOS/{Hybrid,Tabulated}/`; see [EOS](gems/eos/index.md). |
| Hybrid EOS | Piecewise-polytropic cold component plus thermal $\Gamma_{\rm th}$ ideal-gas EOS; see [Hybrid Polytrope](gems/eos/hybrid-piecewise-polytrope.md). |
| StellarCollapse | O'Connor-Ott 19-field HDF5 tabulated EOS format read by `NRPyEOS_stellarcollapse_read_table`; see [StellarCollapse Adapter](gems/eos/stellarcollapse-table-adapter.md). |
| CompOSE | External nuclear EOS database (`(T, n_b, Y_q)`) converted via `compose_to_grhayl.py`; see [CompOSE Conversion](gems/neutrinos/compose-eos-conversion-and-regularization.md). |
| `SRO` | Schneider-Roberts-Ott nuclear EOS family (e.g., CompOSE table 141 SLy4 SNA profile); see [StellarCollapse Adapter](gems/eos/stellarcollapse-table-adapter.md). |
| `NRPyLeakage` | Grey neutrino-leakage module under `GRHayL/Neutrinos/NRPyLeakage/` computing opacities, sources, luminosities, and depths; see [Neutrinos](gems/neutrinos/index.md). |
| `muhat` | Neutron-proton chemical-potential difference $\widehat{\mu} = \mu_n - \mu_p$ (MeV) for leakage; see [EOS Inputs](gems/neutrinos/eos-inputs-and-chemical-potentials.md). |
| `PLM` | Piecewise linear reconstruction (`minmod`, `MC`, `superbee`) in `GRHayL/Reconstruction/`; see [PLM Limiters](gems/reconstruction/plm-limiters.md). |
| `PPM` | Piecewise parabolic reconstruction with Colella-Woodward steepening and flattening in `GRHayL/Reconstruction/`; see [PPM Flow](gems/reconstruction/ppm-flow.md). |
| `WENO-z` | Fifth-order weighted essentially non-oscillatory reconstruction in `GRHayL/Reconstruction/`; see [WENO-z Method](gems/reconstruction/wenoz-method.md). |
| `HLLE` | Harten-Lax-van Leer-Einfeldt 1D GRMHD flux solver in `GRHayL/Flux_Source/` (distinct from Induction 2D HLL); see [HLLE Matrix](gems/flux-source/hlle-flux-variant-matrix.md). |
| `GHL_DISABLE_HDF5` | Preprocessor guard set by `./configure --disable-hdf5` excluding HDF5 table I/O and table tests; see [Build And Configure](architecture/build-and-configure.md). |
| `make.code.defn` | Recursive Cactus-compatible source manifest (`SRCS`, `SUBDIRS`, `#! INCS`) parsed by `configure`; see [Build And Configure](architecture/build-and-configure.md). |
| `GRHayLib` | Cactus / Einstein Toolkit thorn under `implementations/GRHayLib/` compiling GRHayL via symlinks; see [Implementations](implementations/index.md). |
| `ET_Legacy` | Legacy Einstein Toolkit / `IllinoisGRMHD` comparison tests and `GRHayL/TestData` binaries; see [ET Legacy](tests/et-legacy-comparison.md). |
| `GRHayLET` | External Einstein Toolkit evolution thorn suite calling `GRHayLib` parameters and schedules. |
| `BNS_NURATES` | External spectral neutrino-reaction library for blocking qualification; see [Nucleon Blocking](gems/neutrinos/nucleon-blocking-and-eos-conventions.md). |
| `FDINT` | External scalar Fermi-Dirac minimax fit library ported under BSD-3 into `NRPyLeakage_nucleon_blocking.h`; see [Generator Provenance](gems/neutrinos/generator-provenance.md). |

## See Also

- See also: [KB Rules](SCHEMA.md)
- See also: [Workflows](workflows.md)
- See also: [Public API Map](architecture/public-api-map.md)
