# Global Catalog

> Keyword and scope index for normal leaves across `wiki/`. · Status: confirmed

Search this table with `rg -n -i -F '<term>' wiki/catalog.md` when router titles do not resolve a topic; never read the full table.

| Page | Scope | Query terms |
| --- | --- | --- |
| [Build And Configure](architecture/build-and-configure.md) | Configure flags, make.code.defn, HDF5 filter, style. | configure, make.code.defn, GHL_DISABLE_HDF5, parser, clang-format |
| [CI Workflows](architecture/ci-workflows.md) | GitHub Actions matrix, actions, coverage, codecov. | github-actions, compile_GRHayL, codecov.yml, run_tests.sh, paths-ignore |
| [Public API Map](architecture/public-api-map.md) | Public headers, install lists, conditional symbols. | ghl.h, ghl_con2prim.h, ghl_eos_functions.h, ghl_radiation.h, IHDS |
| [Generated Boundaries](architecture/generated-boundaries.md) | Doxygen outputs, Flux_Source, NRPyLeakage, Makefile. | Doxyfile, DoxygenLayout.xml, generate_flux_source.py, generate_makefile.sh |
| [Change Impact](architecture/change-impact.md) | Subsystem dependency map for code, header, doc edits. | change impact, header edit, EOS change, GRHayLib, TestData |
| [Variables And Conventions](physics/variables-and-conventions.md) | GRMHD primitives, conservatives, Btilde, 3+1 signs. | derivation.md, prims, cons, BU, Btilde, vU, ADM, Valencia |
| [Evolution Equation Map](physics/evolution-equation-map.md) | Flux-conservative RHSs, induction, leakage closure. | flux, source terms, induction, Lorenz gauge, phitilde, leakage |
| [Shared Parameters And Enums](core/shared-parameters-and-enums.md) | Core manifest, structs, enums, non-Core ghl.h audit. | ghl_parameters, ghl_eos_parameters, ghl_error_codes_t, GRHayL_Core |
| [Struct Pack And Unpack](core/struct-pack-unpack.md) | Pack/unpack helpers for primitives, conservatives, B. | ghl_initialize_primitives, ghl_return_primitives, ghl_return_conservatives |
| [Metric And ADM Helpers](core/metric-adm.md) | 3-metric, inverse, detgamma, ADM extrinsic curvature. | ghl_initialize_metric, ghl_ADM_aux_quantities, detgamma, CONTR-0011 |
| [Velocity u0](core/velocity-u0.md) | Four-velocity u0, Lorentz factor, speed limiting. | ghl_compute_u0, ghl_limit_v_and_compute_u0, W_max, vU |
| [Stress-Energy And Smallb](core/stress-energy-smallb.md) | Fluid-frame b^mu, b^2, T^{\mu\nu}, T_{\mu\nu}. | ghl_compute_smallb_and_b2, ghl_compute_Tmunu, smallb, b2 |
| [EOS Dispatch](core/eos-dispatch.md) | Core EOS initialization and function-pointer dispatch. | ghl_initialize_eos_functions, ghl_eos_type, initialize_eos.c |
| [Errors, IO, And Debug Utilities](core/errors-io-debug-utilities.md) | Error codes, abort macros, debug prints, info output. | ghl_error, ghl_abort_if_error, ghl_info, ghl_read_error_codes |
| [Core Tests And Reference Data](core/tests-and-reference-data.md) | Unit tests covering Core metric, u0, Tmunu, errors. | unit_test_metric, unit_test_u0, unit_test_Tmunu, unit_test_code_error |
| [Atmosphere](gems/atmosphere.md) | Density/pressure/velocity atmosphere reset and falloff. | ghl_set_prims_to_constant_atm, rho_atm, radial_falloff, CONTR-0009 |
| [Con2Prim Solver Matrix](gems/con2prim/solver-matrix.md) | Solver methods, EOS compatibility, backup chains. | Noble2D, Noble1D, Palenzuela1D, Newman1D, Font1D, CerdaDuran |
| [Con2Prim Recovery Flow](gems/con2prim/recovery-flow.md) | Primary/backup dispatch, diagnostics, entropy fallback. | ghl_con2prim_multi_method, ghl_con2prim_diagnostics, backup_routine |
| [Con2Prim Limits And Conversions](gems/con2prim/limits-and-conversions.md) | Prim2Con, undensitize, S/tau bounds, atmosphere check. | ghl_compute_conservs, ghl_undensitize_conservs, ghl_enforce_primitive_limits |
| [Neural-Network Primitive Guess](gems/con2prim/neural-network-primitive-guess.md) | Optional C NN initial guess ABI and fallback path. | GHL_NN_C2P_API_VERSION, grhayl_nn_c2p, neural network |
| [Con2Prim Internal Numerics](gems/con2prim/internal-numerics-and-root-finding.md) | Newton-Raphson, Brent, Noble/Palenzuela/Newman math. | con2prim_CerdaDuran3D.cc, Noble, Palenzuela, Newman, Font1D |
| [Con2Prim Tests And Reference Data](gems/con2prim/tests-and-reference-data.md) | Con2Prim C unit tests, generators, Newman history. | unit_test_con2prim, con2prim_newman_history.h, con2prim-testdata-ref |
| [EOS Initialization And Dispatch](gems/eos/initialization-and-dispatch.md) | Ideal, Hybrid, Simple, Tabulated EOS init & dispatch. | ghl_initialize_hybrid_eos, ghl_initialize_tabulated_eos, eos_type |
| [Hybrid Piecewise Polytrope](gems/eos/hybrid-piecewise-polytrope.md) | Piecewise-polytropic cold + thermal Gamma_th EOS. | neos, rho_ppoly, Gamma_ppoly, K_ppoly, Gamma_th, P_min, P_max |
| [Tabulated Table Layout](gems/eos/tabulated-table-layout.md) | Tabulated EOS structs, HDF5 guards, memory cleanup. | ghl_tabulated_free_memory, GHL_DISABLE_HDF5, energy_shift, CONTR-0010 |
| [StellarCollapse Table Adapter](gems/eos/stellarcollapse-table-adapter.md) | 19-field HDF5 reader, cgs-to-code, energy conventions. | NRPyEOS_stellarcollapse_read_table, CGS_TO_CODE_DENSITY, SLy4, SRO-141 |
| [Tabulated Interpolation And Bounds](gems/eos/tabulated-interpolation-and-bounds.md) | Trilinear log-grid lookup, T-inversion, domain errors. | NRPyEOS_get_interp_indices, T_from_eps, T_from_P, T_from_S, bounds |
| [Tabulated Interpolator Catalog](gems/eos/tabulated-interpolator-catalog.md) | Complete map of NRPyEOS tabulated interpolator calls. | NRPyEOS_P_and_eps_from_rho_Ye_T, muhat, cs2, entropy, interpolators |
| [EOS Tests And Reference Data](gems/eos/tests-and-reference-data.md) | Hybrid, Simple, Tabulated, and CompOSE EOS unit tests. | unit_test_hybrid_eos, unit_test_tabulated_eos, unit_test_tabulated_eos_compose |
| [Characteristic Speeds](gems/flux-source/characteristic-speeds.md) | Fast magnetosonic speeds cmin/cmax and dispersion. | ghl_compute_characteristic_speeds, cmin, cmax, fast magnetosonic |
| [HLLE Flux Variant Matrix](gems/flux-source/hlle-flux-variant-matrix.md) | HLLE flux entry points for Ideal/Hybrid/Tabulated/S. | ghl_calculate_HLLE_fluxes, entropy flux, Ye flux, dirn |
| [Source Terms](gems/flux-source/source-terms.md) | Valencia GRMHD geometric source terms for S_i and tau. | ghl_calculate_source_terms, S_tilde_rhs, tau_tilde_rhs, metric derivs |
| [Generated NRPy Flux/Source Boundary](gems/flux-source/generated-nrpy-boundary.md) | Flux_Source Python NRPy 2 generator vs checked-in C. | generate_flux_source.py, calculate_HLLE_fluxes, NRPy 2 |
| [Flux_Source Tests And Reference Data](gems/flux-source/tests-and-reference-data.md) | HLLE flux, speed, and source-term unit tests. | unit_test_HLLE_flux, unit_test_characteristic_speeds, test_compute_h_and_cs2 |
| [Induction HLL Flux](gems/induction/hll-flux.md) | Vector-potential HLL electric-field fluxes and Btilde. | ghl_calculate_HLL_fluxes, ghl_HLL_vars, A_i, Btilde |
| [Induction Interpolation And Staggering](gems/induction/interpolation-and-staggering.md) | Staggered A_i, edge-centered E_i, interpolation stencils. | stagger, A_x, A_y, A_z, electric field, edge interpolation |
| [Induction Gauge RHS](gems/induction/gauge-rhs.md) | Generalized Lorenz gauge phitilde RHS and damping. | ghl_calculate_phitilde_rhs, phitilde, xi_lorenz, Lorenz gauge |
| [Induction Tests And Reference Data](gems/induction/tests-and-reference-data.md) | Induction HLL, interpolation, and gauge unit tests. | unit_test_induction, HLL_flux, phitilde_rhs, induction tests |
| [Induction Verification Workflows](gems/induction/verification-workflows.md) | Focused build and verification steps for Induction. | Induction verification, make tests, disable-hdf5 |
| [Neutrinos API And Data](gems/neutrinos/api-and-data.md) | Radiation structs, NRPyLeakage calls, errors, HDF5. | ghl_radiation.h, ghl_nrpyleakage.h, robust_isfinite, NRPYLEAKAGE_FD_OR_RETURN |
| [Neutrinos Physics Model And Sources](gems/neutrinos/physics-model-and-sources.md) | Grey leakage rates, [0]/[1] slots, Valencia closure. | R_source, Q_source, nue, anue, nux, 6.0/NRPyLeakage_c_light |
| [Neutrinos EOS Inputs And Chemical Potentials](gems/neutrinos/eos-inputs-and-chemical-potentials.md) | Six-output EOS callback, mu_e/muhat, munu, underflow. | muhat, mu_e, munu, NRPyLeakage_Q_npmass, gradual underflow |
| [Nucleon Blocking And EOS Conventions](gems/neutrinos/nucleon-blocking-and-eos-conventions.md) | Density-derived blocking, expm1 overlaps, DD2 tests. | NRPyLeakage_nucleon_blocking.h, BNS_NURATES, DD2, ILEAS, FDINT |
| [Neutrinos Implementation Flow](gems/neutrinos/implementation-flow.md) | NRPyLeakage build list, link sets, per-file C flows. | NRPyLeakage_rate_helpers.h, EnsureFinite, PathOfLeastResistance |
| [Neutrinos Generator Provenance](gems/neutrinos/generator-provenance.md) | Ancestral NRPy+ notebooks, tmp_* CSE, C corrections. | Tabulated_EOS_IllinoisGRMHD, tmp_0, NRPyLeakage_ZL, Fukushima |
| [CompOSE EOS Conversion And Regularization](gems/neutrinos/compose-eos-conversion-and-regularization.md) | Offline SRO-141 converter, preflight, CompOSE decode. | compose_to_grhayl.py, sro-sly4-sna-141-regularized-v1, eoscompose.h5, Q3 |
| [CompOSE EOS Leakage Callback](gems/neutrinos/compose-eos-leakage-callback.md) | Runtime (rho,Ye,T)->(T,nb,Yq) map, 6 outputs, checks. | m_ref, Y_q, index_thermo, index_yi, EOSmaker, mu_n, mu_p |
| [Neutrinos Tests And Reference Data](gems/neutrinos/tests-and-reference-data.md) | nrpyleakage_main.h, SLy4 replays, table-free physics. | unit_test_nrpyleakage_physics, optically_thin_gas, constant_density_sphere |
| [Face And Stencil Conventions](gems/reconstruction/face-and-stencil-conventions.md) | Left-face vs centered two-face stencils, read-only dox. | Reconstruction.dox, r_face, l_face, stencil, reconstruction |
| [PLM Limiters](gems/reconstruction/plm-limiters.md) | Minmod, MC, and superbee piecewise-linear limiters. | ghl_minmod, ghl_mc, ghl_superbee, PLM |
| [PPM Flow](gems/reconstruction/ppm-flow.md) | Colella-Woodward PPM steepening, flattening, pressure. | ghl_ppm_reconstruction, steepen, flatten, tilde_f, PPM |
| [WENO-z Method](gems/reconstruction/wenoz-method.md) | Fifth-order WENO-z smoothness indicators and weights. | ghl_wenoz_reconstruction, WENO-z, smoothness indicator |
| [Reconstruction Tests And Reference Data](gems/reconstruction/tests-and-reference-data.md) | PLM/PPM/WENO-z unit tests and ET_Legacy binaries. | unit_test_reconstruction, ppm_reconstruction, CONTR-0012 |
| [GRHayLib Cactus Build](implementations/grhaylib-cactus-build.md) | GRHayLib thorn symlinks, CCL build, CarpetX note. | implementations/GRHayLib, interface.ccl, make.code.defn, CarpetX |
| [GRHayLib Runtime Parameters](implementations/grhaylib-runtime-parameters.md) | GRHayLib param.ccl, schedule.ccl, global EOS/params. | param.ccl, schedule.ccl, GRHayLib_Initialize, GRHayLET |
| [GRHayLib Verification](implementations/grhaylib-verification.md) | Static parity limits vs Cactus build/runtime proof. | GRHayLib verification, Cactus build, CONTR-0005 |
| [Unit Test Catalog](tests/unit-test-catalog.md) | Index of C unit tests and Python CompOSE tests. | unit_test_, Unit_Tests/compose, CONTR-0003 |
| [Data Generators And Helpers](tests/data-generators-and-helpers.md) | Data generators, shared test helpers, sample table. | Unit_Tests/data_gen, ghl_test_helpers, pert_test_fail, randomize |
| [Runner And Build Outputs](tests/runner-and-build-outputs.md) | .github/run_tests.sh, TestData fetches, EXIT cleanup. | run_tests.sh, TestData, EXIT trap, CONTR-0003 |
| [Reference Data And Test Drivers](tests/reference-data-and-test-drivers.md) | Perturbation tolerances, binary replay, driver rules. | ghl_pert_test_fail, TestData, CONTR-0003, CONTR-0004 |
| [HDF5 Sample Tables](tests/hdf5-sample-tables.md) | Tracked simple_table.h5 and downloaded SLy4/LS220. | simple_table.h5, SLy4, LS220, HDF5 table |
| [ET Legacy Comparison](tests/et-legacy-comparison.md) | ET_Legacy regression drivers and pinned TestData ref. | et-legacy-testdata-ref, ET_Legacy, CONTR-0003, CONTR-0004, CONTR-0012 |
| [Expected Failure And Error Keys](tests/expected-failure-and-error-keys.md) | unit_test_code_error.c keys 0..90 and error paths. | unit_test_code_error, error keys, ghl_error_codes_t, CONTR-0010 |
| [Unit Test Coverage And Gap Matrix](tests/unit-test-coverage-and-gap-matrix.md) | Module-by-module test coverage and unverified gaps. | coverage matrix, test gaps, untested paths, no-HDF5 |
