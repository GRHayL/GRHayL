# Equation Of State (EOS)

> Router for simple, hybrid piecewise-polytrope, and tabulated StellarCollapse EOS initialization, table layout, interpolation, bounds, and tests.

| Page | Go here when... |
| --- | --- |
| [EOS Initialization And Dispatch](initialization-and-dispatch.md) | You need simple/hybrid/tabulated setup, parameter validation rules, global function-pointer dispatch, or `GHL_DISABLE_HDF5` behavior. |
| [Hybrid Piecewise-Polytrope EOS](hybrid-piecewise-polytrope.md) | You need piecewise-polytrope continuity constants (`K_ppoly`, `eps_integ_const`, `p_ppoly`), cold/thermal helpers, or breakpoint bounds. |
| [Tabulated Table Layout](tabulated-table-layout.md) | You need HDF5 table loading flow, `NRPyEOS_keys`, loader acceptance checks, energy shift, sound-speed adjustment, or `NRPyEOS_free_memory`. |
| [StellarCollapse Table Adapter](stellarcollapse-table-adapter.md) | You need `ghl_eos_table_stellarcollapse`, HDF5 helpers, producer chemical-potential references (legacy SLy4 vs. SRO-141), or CompOSE regularization interface. |
| [Tabulated Interpolation And Bounds](tabulated-interpolation-and-bounds.md) | You need uniform-log-grid requirements, `(rho,Y_e,T)` and auxiliary temperature-recovery wrappers, bounds checks vs. clamps, or `ghl_compute_h_and_cs2`. |
| [Tabulated Interpolator Catalog](tabulated-interpolator-catalog.md) | You need built interpolator wrapper groups by input tuple, registry audit, beta-equilibrium `rho`-map caches, or derived enthalpy/sound-speed routines. |
| [EOS Tests And Reference Data](tests-and-reference-data.md) | You need `unit_test_piecewise_polytrope`, `unit_test_tabulated_eos`, `unit_test_tabulated_eos_compose`, CompOSE Python tests, sample tables, or coverage gaps. |

Back to [GRHayL Gems](../index.md).
